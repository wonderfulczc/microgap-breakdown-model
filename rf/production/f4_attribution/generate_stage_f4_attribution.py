#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import resource
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.rf.jefimenko.axisymmetric import rotate_axisymmetric_series  # noqa: E402
from streamer_rf.rf.jefimenko.constants import C0  # noqa: E402
from streamer_rf.rf.jefimenko.diagnostics import current_moment_radiation_approx, near_far_audit  # noqa: E402
from streamer_rf.rf.jefimenko.observer import Observer  # noqa: E402
from streamer_rf.rf.jefimenko.solver import evaluate_waveform  # noqa: E402
from streamer_rf.rf.source.adapters import load_petsc_stage4_field_source_csv  # noqa: E402
from streamer_rf.rf.source.schema import SourceSeries  # noqa: E402
from streamer_rf.rf.spectral.bands import PROJECT_BANDS, STANDARD_BANDS, Band, classify_band, integrate_band  # noqa: E402
from streamer_rf.rf.spectral.cwt import band_energy_vs_time, morlet_cwt  # noqa: E402
from streamer_rf.rf.spectral.fft import compute_one_sided_spectrum  # noqa: E402
from streamer_rf.rf.spectral.fluence import spectral_fluence_from_Erad, vector_power_esd  # noqa: E402
from streamer_rf.rf.spectral.sensitivity import compare_rf_mesh_levels  # noqa: E402
from streamer_rf.rf.spectral.trust import build_trust_report  # noqa: E402


SOURCE_ROOT = ROOT / "results/stage_f4/source_gate_recovery"
OUT_ROOT = ROOT / "rf/production/f4_attribution"
J_PROVENANCE = "CONTINUITY_CONSISTENT_FINITE_VOLUME_FLUX"
SOLVER_VERSION = "eb1dd95+stage-f4-attribution"
OBSERVER_X_M = 0.20
N_PHI = 4
MAX_SOURCE_SNAPSHOTS = 32
MAX_OBSERVER_SAMPLES = 32


@dataclass(frozen=True)
class DatasetSpec:
    case_id: str
    source_dir: Path
    role: str
    voltage_state: str
    geometry_id: str
    field_value_Vpm: float
    mesh_label: str
    mesh_trust_override_Hz: float | None = None


def first_time(path: Path) -> float:
    with path.open(newline="") as handle:
        return float(next(csv.DictReader(handle))["time_s"])


def field_files(path: Path) -> list[Path]:
    files = [p for p in path.glob("fields_*.csv") if p.name != "fields_final.csv"]
    return sorted(files, key=first_time)


def choose_source_files(files: list[Path], *, max_count: int = MAX_SOURCE_SNAPSHOTS) -> list[Path]:
    if len(files) <= max_count:
        return files
    idx = np.unique(np.round(np.linspace(0, len(files) - 1, max_count)).astype(int))
    return [files[int(i)] for i in idx]


def read_metrics(source_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    scalar = pd.read_csv(source_dir / "scalar_history.csv")
    cm = pd.read_csv(source_dir / "current_moment.csv")
    coll = pd.read_csv(source_dir / "collision_metrics.csv")
    return scalar, cm, coll


def source_validity(scalar: pd.DataFrame) -> dict[str, float | bool]:
    continuity_abs = (scalar["plasma_charge_derivative_A"] - scalar["outer_boundary_current_A"]).abs()
    scale = np.maximum.reduce(
        [
            scalar["plasma_charge_derivative_A"].abs().to_numpy(),
            scalar["outer_boundary_current_A"].abs().to_numpy(),
            np.full(len(scalar), 1e-30),
        ]
    )
    rel = continuity_abs.to_numpy() / scale
    rel_skip = rel[1:] if len(rel) > 1 else rel
    no_nan = bool(np.all(np.isfinite(scalar.select_dtypes(include=[float, int]).to_numpy())))
    return {
        "no_nan": no_nan,
        "conservation_residual": float(scalar["conservation_residual"].max()),
        "continuity_abs_A_max": float(continuity_abs.max()),
        "continuity_abs_A_median": float(continuity_abs.median()),
        "continuity_rel_skip_first_max": float(np.max(rel_skip)) if len(rel_skip) else math.nan,
        "continuity_valid": bool(
            no_nan
            and float(scalar["conservation_residual"].max()) < 1e-2
            and float(continuity_abs.max()) < 1.0
            and float(np.max(rel_skip)) < 1e-4
        ),
    }


def load_axisymmetric_series(spec: DatasetSpec) -> SourceSeries:
    records = [
        load_petsc_stage4_field_source_csv(
            path,
            case_id=spec.case_id,
            solver_version=SOLVER_VERSION,
            geometry_id=spec.geometry_id,
            voltage_state=spec.voltage_state,
            photoionization="on",
        )
        for path in choose_source_files(field_files(spec.source_dir))
    ]
    return SourceSeries(tuple(records))


def observer_for_series(series: SourceSeries) -> Observer:
    z = series.records[0].columns["z_center"]
    return Observer("F4_far_field_reference", OBSERVER_X_M, 0.0, float(np.mean(z)))


def valid_observer_times(series3d: SourceSeries, observer: Observer, *, dt_hint_s: float) -> np.ndarray:
    rec = series3d.records[0]
    xyz = np.column_stack((rec.columns["x_center"], rec.columns["y_center"], rec.columns["z_center"]))
    distances = np.linalg.norm(xyz - observer.position[None, :], axis=1)
    tmin = series3d.times[1] + float(distances.max()) / C0
    tmax = series3d.times[-2] + float(distances.min()) / C0
    if tmax <= tmin:
        raise RuntimeError("no valid retarded-time observer support")
    n = max(5, int(np.floor((tmax - tmin) / dt_hint_s)) + 1)
    n = min(n, MAX_OBSERVER_SAMPLES)
    return tmin + np.arange(n, dtype=float) * ((tmax - tmin) / (n - 1))


def current_moment_comparison(cm: pd.DataFrame, observer: Observer, waveform: pd.DataFrame) -> dict[str, float]:
    t = cm["time"].to_numpy(dtype=float)
    Mz = cm["I_CM_electron"].to_numpy(dtype=float)
    dMdt = np.gradient(Mz, t)
    Rvec = observer.position - np.array([0.0, 0.0, observer.z_m])
    approx_z = np.array([current_moment_radiation_approx(np.array([0.0, 0.0, dm]), Rvec)[2] for dm in dMdt])
    delay = float(np.linalg.norm(Rvec) / C0)
    q_time = waveform["time_s"].to_numpy(dtype=float) - delay
    approx_interp = np.interp(q_time, t, approx_z, left=np.nan, right=np.nan)
    full = waveform["Ez_dJ"].to_numpy(dtype=float)
    valid = np.isfinite(approx_interp)
    if np.count_nonzero(valid) < 3:
        return {"correlation": math.nan, "peak_time_difference_s": math.nan, "normalized_waveform_error": math.nan}
    a = approx_interp[valid]
    b = full[valid]
    corr = float(np.corrcoef(a, b)[0, 1]) if np.std(a) > 0 and np.std(b) > 0 else math.nan
    ta = q_time[valid][int(np.argmax(np.abs(a)))] + delay
    tb = waveform["time_s"].to_numpy(dtype=float)[valid][int(np.argmax(np.abs(b)))]
    err = float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300))
    return {"correlation": corr, "peak_time_difference_s": float(tb - ta), "normalized_waveform_error": err}


def p_stage_labels(scalar: pd.DataFrame, coll: pd.DataFrame, cm: pd.DataFrame) -> pd.DataFrame:
    t0 = float(scalar["time_s"].iloc[0])
    t1 = float(scalar["time_s"].iloc[-1])
    Ne = scalar["total_electrons"].to_numpy(dtype=float)
    ne = scalar["ne_max_m_3"].to_numpy(dtype=float)
    times = scalar["time_s"].to_numpy(dtype=float)
    cm_times = cm["time"].to_numpy(dtype=float)
    dM = np.abs(np.gradient(cm["I_CM_electron"].to_numpy(dtype=float), cm_times))
    dM_interp = np.interp(times, cm_times, dM)
    onset_mask = (Ne >= 1.5 * Ne[0]) | (ne >= 1.5 * ne[0])
    inception_start = float(times[np.argmax(onset_mask)]) if np.any(onset_mask) else float(times[len(times) // 3])
    head = np.interp(times, coll["time"].to_numpy(dtype=float), coll["left_inner_z"].to_numpy(dtype=float))
    head0 = float(head[0])
    prop_mask = np.abs(head - head0) >= 4.0e-5
    prop_start = float(times[np.argmax(prop_mask)]) if np.any(prop_mask) else float(times[2 * len(times) // 3])
    if prop_start <= inception_start:
        prop_start = float(times[min(len(times) - 1, np.searchsorted(times, inception_start) + max(2, len(times) // 8))])
    rows = [
        {
            "dataset": "F4-P",
            "stage": "AVALANCHE",
            "start_time_s": t0,
            "end_time_s": inception_start,
            "evidence": "Ne or ne_max before 1.5x growth threshold",
        },
        {
            "dataset": "F4-P",
            "stage": "INCEPTION",
            "start_time_s": inception_start,
            "end_time_s": prop_start,
            "evidence": "Ne/ne_max growth established; head displacement below 40 um",
        },
        {
            "dataset": "F4-P",
            "stage": "PROPAGATION",
            "start_time_s": prop_start,
            "end_time_s": t1,
            "evidence": f"head displacement={abs(float(head[-1])-head0):.3e} m; peak |dM/dt|={float(dM_interp.max()):.3e} A m s^-1",
        },
    ]
    return pd.DataFrame([r for r in rows if r["end_time_s"] > r["start_time_s"]])


def c_stage_labels(scalar: pd.DataFrame, coll: pd.DataFrame, event_time_s: float) -> pd.DataFrame:
    t0 = float(scalar["time_s"].iloc[0])
    t1 = float(scalar["time_s"].iloc[-1])
    if event_time_s <= 0.0:
        d = coll["d_head"].to_numpy(dtype=float)
        event_time_s = float(coll["time"].iloc[int(np.nanargmin(d))])
    encounter_start = max(t0, event_time_s - 2.0e-10)
    collision_start = max(t0, event_time_s - 5.0e-11)
    post_start = min(t1, event_time_s + 5.0e-11)
    rows = [
        {"dataset": "F4-C", "stage": "PRE_INTERACTION", "start_time_s": t0, "end_time_s": encounter_start, "evidence": "before event-time minus 0.2 ns"},
        {"dataset": "F4-C", "stage": "ENCOUNTER", "start_time_s": encounter_start, "end_time_s": collision_start, "evidence": "heads approaching event-time window"},
        {"dataset": "F4-C", "stage": "COLLISION", "start_time_s": collision_start, "end_time_s": post_start, "evidence": f"stage4_run event_time_s={event_time_s:.6e}"},
        {"dataset": "F4-C", "stage": "POST_COLLISION", "start_time_s": post_start, "end_time_s": t1, "evidence": "early post-event response"},
    ]
    return pd.DataFrame([r for r in rows if r["end_time_s"] > r["start_time_s"]])


def termination_event_time(source_dir: Path, scalar: pd.DataFrame) -> float:
    term = source_dir / "termination.csv"
    if term.exists():
        try:
            value = float(pd.read_csv(term)["event_time"].iloc[0])
            if value > 0.0:
                return value
        except Exception:
            pass
    coll = pd.read_csv(source_dir / "collision_metrics.csv")
    return float(coll["time"].iloc[int(np.argmin(coll["d_head"].to_numpy(dtype=float)))])


def trusted_band_status(trust_low: float, trust_high: float) -> dict[str, str]:
    return {band.name: classify_band(band, trust_low, trust_high) for band in (*STANDARD_BANDS, *PROJECT_BANDS)}


def trusted_band_energy(
    frequency: np.ndarray,
    density: np.ndarray,
    band: Band,
    status: str,
) -> float:
    if status != "TRUSTED":
        return math.nan
    return integrate_band(frequency, density, band)


def analyze_dataset(spec: DatasetSpec, stage_labels: pd.DataFrame, *, mesh_trust_Hz: float | None) -> dict[str, object]:
    started = time.perf_counter()
    scalar, cm, coll = read_metrics(spec.source_dir)
    validity = source_validity(scalar)
    source2d = load_axisymmetric_series(spec)
    source3d = rotate_axisymmetric_series(source2d, n_phi=N_PHI)
    obs = observer_for_series(source2d)
    dt_hint = float(np.median(np.diff(source2d.times)))
    obs_times = valid_observer_times(source3d, obs, dt_hint_s=dt_hint)
    waveform_samples = evaluate_waveform(source3d, [obs], obs_times, source_manifest_id=spec.case_id, chunk_size=200_000)
    waveform = pd.DataFrame([sample.to_row() for sample in waveform_samples])
    waveform.to_csv(OUT_ROOT / f"{spec.case_id}_jefimenko_waveform.csv", index=False)
    Erad = waveform[["Ex_dJ", "Ey_dJ", "Ez_dJ"]].to_numpy(dtype=float)
    fluence = spectral_fluence_from_Erad(waveform["time_s"].to_numpy(dtype=float), Erad, window="hann", far_field_valid=True)
    spectrum = fluence["spectrum"]
    power_esd = vector_power_esd(spectrum)
    duration = float(waveform["time_s"].iloc[-1] - waveform["time_s"].iloc[0])
    if mesh_trust_Hz is None:
        mesh_trust_Hz = float(spectrum.nyquist_Hz)
    trust = build_trust_report(
        source_valid=bool(validity["continuity_valid"]),
        current_provenance=J_PROVENANCE,
        continuity_status="PASS" if validity["continuity_valid"] else "FAIL",
        remap_status="CONSERVATIVE",
        dt_s=float(spectrum.dt_s),
        duration_s=duration,
        derivative_trust_frequency_Hz=0.65 * float(spectrum.nyquist_Hz),
        interpolation_trust_frequency_Hz=0.70 * float(spectrum.nyquist_Hz),
        mesh_trust_frequency_Hz=float(mesh_trust_Hz),
        sampling_trust_frequency_Hz=0.90 * float(spectrum.nyquist_Hz),
        far_field_valid=True,
    )
    bands = trusted_band_status(trust.trusted_frequency_low_Hz, trust.trusted_frequency_high_Hz)
    trust_payload = trust.to_dict()
    trust_payload["band_status"] = bands
    (OUT_ROOT / f"{spec.case_id}_rf_trust_report.json").write_text(json.dumps(trust_payload, indent=2, sort_keys=True) + "\n")

    freq_low = max(trust.trusted_frequency_low_Hz, 1.0 / max(duration, 1e-30))
    freq_high = min(trust.trusted_frequency_high_Hz, spectrum.nyquist_Hz)
    if freq_high > freq_low:
        cwt_freq = np.geomspace(freq_low, freq_high, 48)
        cwt = morlet_cwt(waveform["time_s"].to_numpy(dtype=float), Erad[:, 2], cwt_freq, n_cycles=5.0)
        cwt_bands = {
            band.name: (band.f_low_Hz, band.f_high_Hz)
            for band in (*STANDARD_BANDS, *PROJECT_BANDS)
            if classify_band(band, trust.trusted_frequency_low_Hz, trust.trusted_frequency_high_Hz) != "UNTRUSTED"
        }
        cwt_energy_t = band_energy_vs_time(cwt, cwt_bands) if cwt_bands else {}
    else:
        cwt = None
        cwt_energy_t = {}

    rows: list[dict[str, object]] = []
    for label in stage_labels.to_dict("records"):
        mask_time = (waveform["time_s"] >= label["start_time_s"] + OBSERVER_X_M / C0) & (
            waveform["time_s"] <= label["end_time_s"] + OBSERVER_X_M / C0
        )
        source_mask = (scalar["time_s"] >= label["start_time_s"]) & (scalar["time_s"] <= label["end_time_s"])
        cm_mask = (cm["time"] >= label["start_time_s"]) & (cm["time"] <= label["end_time_s"])
        if source_mask.any():
            E_range = f"{float(scalar.loc[source_mask, 'E_max_V_m'].min()):.6e}:{float(scalar.loc[source_mask, 'E_max_V_m'].max()):.6e}"
            Ne_range = f"{float(scalar.loc[source_mask, 'total_electrons'].min()):.6e}:{float(scalar.loc[source_mask, 'total_electrons'].max()):.6e}"
        else:
            E_range = ""
            Ne_range = ""
        head = coll["left_inner_z"].to_numpy(dtype=float)
        head_v = np.gradient(head, coll["time"].to_numpy(dtype=float))
        head_mask = (coll["time"] >= label["start_time_s"]) & (coll["time"] <= label["end_time_s"])
        head_range = (
            f"{float(np.nanmin(head_v[head_mask])):.6e}:{float(np.nanmax(head_v[head_mask])):.6e}"
            if np.any(head_mask)
            else ""
        )
        dM = np.gradient(cm["I_CM_electron"].to_numpy(dtype=float), cm["time"].to_numpy(dtype=float))
        peak_dMdt = float(np.nanmax(np.abs(dM[cm_mask]))) if cm_mask.any() else math.nan
        band_energy_cols: dict[str, float] = {}
        for band in (*STANDARD_BANDS, *PROJECT_BANDS):
            status = bands[band.name]
            if status != "UNTRUSTED" and cwt is not None and band.name in cwt_energy_t and np.any(mask_time):
                band_energy_cols[f"{band.name}_energy"] = float(
                    np.trapezoid(cwt_energy_t[band.name][mask_time.to_numpy()], waveform.loc[mask_time, "time_s"].to_numpy())
                )
            else:
                band_energy_cols[f"{band.name}_energy"] = math.nan
        trusted_mask = (spectrum.frequency_Hz >= trust.trusted_frequency_low_Hz) & (
            spectrum.frequency_Hz <= trust.trusted_frequency_high_Hz
        )
        dominant = float(spectrum.frequency_Hz[trusted_mask][int(np.argmax(power_esd[trusted_mask]))]) if np.any(trusted_mask) else math.nan
        peak_cwt = math.nan
        if cwt is not None and np.any(mask_time):
            sub = np.where(cwt.coi_mask[:, mask_time.to_numpy()], cwt.power[:, mask_time.to_numpy()], -np.inf)
            if np.any(np.isfinite(sub)):
                peak_cwt = float(cwt.frequency_Hz[np.unravel_index(int(np.nanargmax(sub)), sub.shape)[0]])
        rows.append(
            {
                "dataset": label["dataset"],
                "stage": label["stage"],
                "start_time_s": label["start_time_s"],
                "end_time_s": label["end_time_s"],
                "Emax_range": E_range,
                "Ne_range": Ne_range,
                "head_velocity_range": head_range,
                "peak_dMdt": peak_dMdt,
                "trusted_low": trust.trusted_frequency_low_Hz,
                "trusted_high": trust.trusted_frequency_high_Hz,
                "VHF_status": bands["VHF"],
                "UHF_status": bands["UHF"],
                "SHF_status": bands["SHF"],
                "60_90MHz_status": bands["60-90 MHz"],
                "100_200MHz_status": bands["100-200 MHz"],
                "200_300MHz_status": bands["200-300 MHz"],
                "300_500MHz_status": bands["300-500 MHz"],
                "500_1000MHz_status": bands["500 MHz-1 GHz"],
                "1_3GHz_status": bands["1-3 GHz"],
                "3_10GHz_status": bands["3-10 GHz"],
                **band_energy_cols,
                "dominant_trusted_frequency": dominant,
                "peak_CWT_trusted_frequency": peak_cwt,
                "stage_frequency_attribution": "SUPPORTED"
                if trust.scientific_rf_valid
                and trust.trusted_frequency_high_Hz > trust.trusted_frequency_low_Hz
                and np.any(mask_time)
                else "NOT_RESOLVED",
                "evidence": label["evidence"],
            }
        )

    stage_summary = pd.DataFrame(rows)
    stage_summary.to_csv(OUT_ROOT / f"{spec.case_id}_stage_rf_summary.csv", index=False)
    stage_labels.to_csv(OUT_ROOT / f"{spec.case_id}_stage_labels.csv", index=False)
    waveform_relation = current_moment_comparison(cm, obs, waveform)
    audit = near_far_audit(source3d, obs.position, source_timescale_s=max(float(np.median(np.diff(cm["time"]))), 1e-30))
    compact = {
        "case_id": spec.case_id,
        "role": spec.role,
        "source_dir": str(spec.source_dir.relative_to(ROOT)),
        "J_provenance": J_PROVENANCE,
        "source_duration_s": float(scalar["time_s"].iloc[-1] - scalar["time_s"].iloc[0]),
        "source_output_dt_median_s": float(np.median(np.diff([first_time(p) for p in field_files(spec.source_dir)]))),
        "source_snapshots": len(field_files(spec.source_dir)),
        "observer_waveform_samples": len(waveform),
        "observer_waveform_duration_s": duration,
        "rf_trust_report": trust_payload,
        "remap_note": "native 2D axisymmetric source uses an unchanged structured partition; no AMR remap is required",
        "band_status": bands,
        "validity": validity,
        "fluence_parseval_relative_error": float(fluence["parseval_relative_error"]),
        "M_dMdt_vs_jefimenko": waveform_relation,
        "near_far_audit": audit,
        "runtime_s": time.perf_counter() - started,
    }
    (OUT_ROOT / f"{spec.case_id}_summary.json").write_text(json.dumps(compact, indent=2, sort_keys=True) + "\n")
    return compact


def mesh_trust_from_collision() -> dict[str, object]:
    fine_dir = SOURCE_ROOT / "F4-C-interaction-collision"
    coarse_dir = SOURCE_ROOT / "F4-C-interaction-collision-coarse40um"
    if not coarse_dir.exists() or not (coarse_dir / "current_moment.csv").exists():
        return {"status": "NOT_RUN", "mesh_trust_frequency_Hz": math.inf}
    fine = pd.read_csv(fine_dir / "current_moment.csv")
    coarse = pd.read_csv(coarse_dir / "current_moment.csv")
    tf0 = max(float(fine["time"].min()), float(coarse["time"].min()))
    tf1 = min(float(fine["time"].max()), float(coarse["time"].max()))
    dt = max(float(np.median(np.diff(fine["time"]))), float(np.median(np.diff(coarse["time"]))))
    n = max(8, int(np.floor((tf1 - tf0) / dt)) + 1)
    t = np.linspace(tf0, tf1, n)
    fine_m = np.interp(t, fine["time"], fine["I_CM_electron"])
    coarse_m = np.interp(t, coarse["time"], coarse["I_CM_electron"])
    fine_d = np.gradient(fine_m, t)
    coarse_d = np.gradient(coarse_m, t)
    result = compare_rf_mesh_levels(t, coarse_d, t, fine_d, tolerance=0.50, rolling_bins=3)
    return {
        "status": "COARSE40UM_VS_FINE20UM",
        "overlap_start_s": tf0,
        "overlap_end_s": tf1,
        "mesh_trust_frequency_Hz": float(result["mesh_trust_frequency_Hz"]),
        "tolerance": result["tolerance"],
        "coarse_grid_um": 40.0,
        "fine_grid_um": 20.0,
    }


def main() -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    p_spec = DatasetSpec(
        case_id="F4-P-stage4-left-isolated",
        source_dir=SOURCE_ROOT / "F4-P-streamer-propagation",
        role="single continuous S4-LEFT-ISOLATED avalanche/inception/propagation",
        voltage_state="4.8e6 V/m frozen background field",
        geometry_id="S4-LEFT-ISOLATED frozen 20 um axisymmetric propagation case",
        field_value_Vpm=4.8e6,
        mesh_label="20um",
    )
    c_spec = DatasetSpec(
        case_id="F4-C-stage5-highfield-collision",
        source_dir=SOURCE_ROOT / "F4-C-interaction-collision",
        role="S5-HIGHFIELD-COLLISION event-centered collision window",
        voltage_state="6.4e6 V/m frozen high-field background field",
        geometry_id="S5-HIGHFIELD-COLLISION frozen 20 um axisymmetric head-on collision case",
        field_value_Vpm=6.4e6,
        mesh_label="20um",
    )
    p_scalar, p_cm, p_coll = read_metrics(p_spec.source_dir)
    c_scalar, c_cm, c_coll = read_metrics(c_spec.source_dir)
    p_labels = p_stage_labels(p_scalar, p_coll, p_cm)
    c_event = termination_event_time(c_spec.source_dir, c_scalar)
    c_labels = c_stage_labels(c_scalar, c_coll, c_event)
    mesh = mesh_trust_from_collision()
    p_result = analyze_dataset(p_spec, p_labels, mesh_trust_Hz=None)
    c_result = analyze_dataset(c_spec, c_labels, mesh_trust_Hz=float(mesh["mesh_trust_frequency_Hz"]))
    combined = pd.concat(
        [
            pd.read_csv(OUT_ROOT / f"{p_spec.case_id}_stage_rf_summary.csv"),
            pd.read_csv(OUT_ROOT / f"{c_spec.case_id}_stage_rf_summary.csv"),
        ],
        ignore_index=True,
    )
    combined.to_csv(OUT_ROOT / "stage_rf_summary.csv", index=False)
    payload = {
        "STAGE_F4_ATTRIBUTION_STATUS": "PASS_CANDIDATE"
        if any(combined["stage_frequency_attribution"] == "SUPPORTED")
        else "BLOCKED",
        "CURRENT_SOURCE_PROVENANCE": J_PROVENANCE,
        "VHF_ATTRIBUTION": "NOT_RESOLVED",
        "UHF_ATTRIBUTION": "NOT_RESOLVED",
        "GHz_SHF_ATTRIBUTION": "SUPPORTED"
        if any((combined["stage_frequency_attribution"] == "SUPPORTED") & (combined["SHF_status"] == "PARTIALLY_TRUSTED"))
        else "NOT_RESOLVED",
        "THREE_D_SCIENTIFIC_POLARIZATION": "DEFERRED",
        "mesh_trust": mesh,
        "datasets": [p_result, c_result],
        "peak_rss_kb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (OUT_ROOT / "stage_f4_attribution_summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": payload["STAGE_F4_ATTRIBUTION_STATUS"], "mesh": mesh}, indent=2))


if __name__ == "__main__":
    main()

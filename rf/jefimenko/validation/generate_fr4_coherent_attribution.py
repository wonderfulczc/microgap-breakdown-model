from __future__ import annotations

import json
import math
import resource
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.rf.coherent_attribution import (  # noqa: E402
    coherent_attribution,
    cwt_linearity_error,
    mechanism_radiative_fields,
    uniform_resample,
)
from streamer_rf.rf.mechanisms import compute_mechanism_decomposition  # noqa: E402
from streamer_rf.rf.spectral.bands import Band, PROJECT_BANDS, STANDARD_BANDS  # noqa: E402
from streamer_rf.rf.spectral.cwt import morlet_cwt  # noqa: E402


FR3 = ROOT / "rf/jefimenko/validation/fr3_mechanisms"
F4 = ROOT / "rf/production/f4_attribution"
OUT = ROOT / "rf/jefimenko/validation/fr4_coherent_attribution"
OBSERVER_VECTOR_M = np.array([0.2, 0.0, 0.0], dtype=float)
N_CWT_FREQ = 48
N_CYCLES = 5.0

MECHANISMS = {
    "HEAD_CHARGE_EVOLUTION": (
        "dM_charge_evolution_x_Am_s",
        "dM_charge_evolution_y_Am_s",
        "dM_charge_evolution_z_Am_s",
    ),
    "HEAD_ACCELERATION": (
        "dM_head_acceleration_x_Am_s",
        "dM_head_acceleration_y_Am_s",
        "dM_head_acceleration_z_Am_s",
    ),
    "CURRENT_MOMENT_REDISTRIBUTION": (
        "dM_redistribution_x_Am_s",
        "dM_redistribution_y_Am_s",
        "dM_redistribution_z_Am_s",
    ),
}


def _finite_or_none(x: float) -> float | None:
    return float(x) if math.isfinite(float(x)) else None


def _vector_columns(df: pd.DataFrame, cols: tuple[str, str, str]) -> np.ndarray:
    return df.loc[:, list(cols)].to_numpy(dtype=float)


def _total_dM(df: pd.DataFrame) -> np.ndarray:
    return df[["dM_total_x_Am_s", "dM_total_y_Am_s", "dM_total_z_Am_s"]].to_numpy(dtype=float)


def _mechanism_terms(df: pd.DataFrame) -> dict[str, np.ndarray]:
    return {name: _vector_columns(df, cols) for name, cols in MECHANISMS.items()}


def _mechanism_fields(df: pd.DataFrame) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray]:
    terms = _mechanism_terms(df)
    fields = mechanism_radiative_fields(terms, OBSERVER_VECTOR_M)
    rec = sum(fields.values(), np.zeros_like(next(iter(fields.values()))))
    total = mechanism_radiative_fields({"total": _total_dM(df)}, OBSERVER_VECTOR_M)["total"]
    return fields, rec, total


def _l2(a: np.ndarray, b: np.ndarray, mask: np.ndarray) -> float:
    m = np.asarray(mask, dtype=bool)
    return float(np.linalg.norm((a - b)[m]) / max(np.linalg.norm(b[m]), 1e-300)) if np.any(m) else math.nan


def _max_abs(a: np.ndarray, b: np.ndarray, mask: np.ndarray) -> float:
    m = np.asarray(mask, dtype=bool)
    return float(np.max(np.abs((a - b)[m]))) if np.any(m) else math.nan


def _stage4_bands(trust: dict[str, object]) -> list[Band]:
    low = float(trust["trusted_frequency_low_Hz"])
    high = float(trust["trusted_frequency_high_Hz"])
    bands = [Band("TRUSTED_CONTINUOUS", low, high)]
    for band in (*STANDARD_BANDS, *PROJECT_BANDS):
        if band.name in {"VHF", "UHF", "1-3 GHz", "60-90 MHz", "100-200 MHz", "200-300 MHz", "300-500 MHz", "500 MHz-1 GHz"}:
            continue
        lo = max(low, band.f_low_Hz)
        hi = min(high, band.f_high_Hz)
        if hi > lo:
            bands.append(Band(band.name, lo, hi))
    return bands


def _compute_cwts(time_s: np.ndarray, fields_z: dict[str, np.ndarray], freqs: np.ndarray) -> tuple[dict[str, object], object, float]:
    c = {name: morlet_cwt(time_s, sig, freqs, n_cycles=N_CYCLES) for name, sig in fields_z.items()}
    rec_sum = sum((v.coefficients for v in c.values()), np.zeros_like(next(iter(c.values())).coefficients))
    rec_direct = morlet_cwt(time_s, sum(fields_z.values(), np.zeros_like(next(iter(fields_z.values())))), freqs, n_cycles=N_CYCLES)
    closure = cwt_linearity_error(rec_sum, rec_direct.coefficients, rec_direct.coi_mask)
    return c, rec_direct, closure


def _stage_band_mask(time_s: np.ndarray, freq: np.ndarray, cwt_mask: np.ndarray, stage: dict[str, object], band: Band) -> np.ndarray:
    tmask = (time_s >= float(stage["start_time_s"])) & (time_s <= float(stage["end_time_s"]))
    fmask = (freq >= band.f_low_Hz) & (freq <= band.f_high_Hz)
    return cwt_mask & fmask[:, None] & tmask[None, :]


def _eta_rows(
    *,
    case_id: str,
    labels: pd.DataFrame,
    trust: dict[str, object],
    cwt_terms: dict[str, object],
    cwt_rec_direct: object,
    cwt_closure: float,
    derivative_tag: str,
) -> list[dict[str, object]]:
    freq = cwt_rec_direct.frequency_Hz
    rec_sum = sum((v.coefficients for v in cwt_terms.values()), np.zeros_like(cwt_rec_direct.coefficients))
    rows: list[dict[str, object]] = []
    for stage in labels.to_dict("records"):
        for band in _stage4_bands(trust):
            mask = _stage_band_mask(cwt_rec_direct.time_s, freq, cwt_rec_direct.coi_mask, stage, band)
            result = coherent_attribution(
                {name: c.coefficients for name, c in cwt_terms.items()},
                rec_sum,
                mask,
            )
            eta = result.eta
            contrib = result.coherent_contributions
            roles = result.interference_role
            rows.append(
                {
                    "case_id": case_id,
                    "stage": stage["stage"],
                    "band_name": band.name,
                    "f_low_Hz": band.f_low_Hz,
                    "f_high_Hz": band.f_high_Hz,
                    "eta_head_charge_evolution": eta.get("HEAD_CHARGE_EVOLUTION", math.nan),
                    "eta_head_acceleration": eta.get("HEAD_ACCELERATION", math.nan),
                    "eta_current_moment_redistribution": eta.get("CURRENT_MOMENT_REDISTRIBUTION", math.nan),
                    "eta_sum": result.eta_sum,
                    "eta_sum_error": result.eta_sum_error,
                    "coherent_charge": contrib.get("HEAD_CHARGE_EVOLUTION", math.nan),
                    "coherent_acceleration": contrib.get("HEAD_ACCELERATION", math.nan),
                    "coherent_redistribution": contrib.get("CURRENT_MOMENT_REDISTRIBUTION", math.nan),
                    "P_reconstructed": result.P_reconstructed,
                    "charge_interference_role": roles.get("HEAD_CHARGE_EVOLUTION", ""),
                    "acceleration_interference_role": roles.get("HEAD_ACCELERATION", ""),
                    "redistribution_interference_role": roles.get("CURRENT_MOMENT_REDISTRIBUTION", ""),
                    "trusted_coefficient_count": result.trusted_count,
                    "trusted_fraction": result.trusted_fraction,
                    "cwt_linearity_closure": cwt_closure,
                    "derivative_configuration": derivative_tag,
                    "attribution_status": "ATTRIBUTION_INTERPRET_WITH_CAUTION"
                    if result.status == "OK"
                    else "ATTRIBUTION_NOT_RESOLVED",
                    "status_reason": result.status,
                }
            )
    return rows


def _load_stage4() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object], dict[str, object]]:
    case = "F4-P-stage4-left-isolated"
    return (
        pd.read_csv(FR3 / f"{case}_mechanisms.csv"),
        pd.read_csv(F4 / f"{case}_stage_labels.csv"),
        json.loads((F4 / f"{case}_rf_trust_report.json").read_text()),
        json.loads((F4 / f"{case}_summary.json").read_text()),
    )


def _alternative_stage4_decomposition(df: pd.DataFrame) -> pd.DataFrame:
    t = df["time_s"].to_numpy(dtype=float)
    q = df["q_head_C"].to_numpy(dtype=float)
    pos = np.column_stack((np.zeros(t.size), np.zeros(t.size), df["head_position_z_m"].to_numpy(dtype=float)))
    M = df[["M_total_x_Am", "M_total_y_Am", "M_total_z_Am"]].to_numpy(dtype=float)
    decomp = compute_mechanism_decomposition(
        t,
        M,
        q,
        pos,
        head_valid=df["head_valid"].to_numpy(dtype=bool).copy(),
        derivative_radius=2,
    )
    alt = df.copy()
    alt["derivative_valid"] = decomp["derivative_valid"]
    for name, cols in MECHANISMS.items():
        key = {
            "HEAD_CHARGE_EVOLUTION": "dM_charge_evolution_Am_s",
            "HEAD_ACCELERATION": "dM_head_acceleration_Am_s",
            "CURRENT_MOMENT_REDISTRIBUTION": "dM_redistribution_Am_s",
        }[name]
        alt.loc[:, list(cols)] = decomp[key]
    alt.loc[:, ["dM_total_x_Am_s", "dM_total_y_Am_s", "dM_total_z_Am_s"]] = decomp["dM_total_Am_s"]
    return alt


def _analyze_stage4_configuration(
    df: pd.DataFrame,
    labels: pd.DataFrame,
    trust: dict[str, object],
    *,
    derivative_tag: str,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    valid = df["derivative_valid"].to_numpy(dtype=bool) & df["trust_valid"].to_numpy(dtype=bool) & df["head_valid"].to_numpy(dtype=bool)
    fields, rec, total = _mechanism_fields(df)
    td_l2 = _l2(rec, total, valid)
    td_max = _max_abs(rec, total, valid)
    tu, fields_u = uniform_resample(df["time_s"].to_numpy(dtype=float), np.stack([fields[k] for k in MECHANISMS], axis=1), valid)
    fields_z = {name: fields_u[:, i, 2] for i, name in enumerate(MECHANISMS)}
    freqs = np.geomspace(float(trust["trusted_frequency_low_Hz"]), float(trust["trusted_frequency_high_Hz"]), N_CWT_FREQ)
    cwt_terms, cwt_rec_direct, cwt_closure = _compute_cwts(tu, fields_z, freqs)
    rows = _eta_rows(
        case_id="F4-P-stage4-left-isolated",
        labels=labels,
        trust=trust,
        cwt_terms=cwt_terms,
        cwt_rec_direct=cwt_rec_direct,
        cwt_closure=cwt_closure,
        derivative_tag=derivative_tag,
    )
    diag = {
        "derivative_configuration": derivative_tag,
        "valid_input_samples": int(np.count_nonzero(valid)),
        "uniform_cwt_samples": int(tu.size),
        "uniform_time_start_s": float(tu[0]),
        "uniform_time_end_s": float(tu[-1]),
        "time_domain_rec_vs_total_L2": _finite_or_none(td_l2),
        "time_domain_rec_vs_total_max_abs": _finite_or_none(td_max),
        "cwt_linearity_closure": _finite_or_none(cwt_closure),
    }
    return rows, diag


def _robustness(primary: pd.DataFrame, alternative: pd.DataFrame, fr3_summary: dict[str, object]) -> pd.DataFrame:
    rows = []
    keys = [
        "eta_head_charge_evolution",
        "eta_head_acceleration",
        "eta_current_moment_redistribution",
    ]
    merged = primary.merge(
        alternative,
        on=["case_id", "stage", "band_name"],
        suffixes=("_primary", "_alternative"),
    )
    for row in merged.to_dict("records"):
        for key in keys:
            p = float(row[f"{key}_primary"]) if pd.notna(row[f"{key}_primary"]) else math.nan
            a = float(row[f"{key}_alternative"]) if pd.notna(row[f"{key}_alternative"]) else math.nan
            sign_ok = math.isfinite(p) and math.isfinite(a) and (p == 0.0 or a == 0.0 or math.copysign(1.0, p) == math.copysign(1.0, a))
            rows.append(
                {
                    "stage": row["stage"],
                    "band_name": row["band_name"],
                    "mechanism": key.removeprefix("eta_").upper(),
                    "eta_primary": p,
                    "eta_alternative": a,
                    "absolute_change": abs(p - a) if math.isfinite(p) and math.isfinite(a) else math.nan,
                    "sign_consistent": sign_ok,
                    "derivative_robustness": "ATTRIBUTION_NOT_ROBUST_TO_DERIVATIVE" if not sign_ok else "DERIVATIVE_SIGN_CONSISTENT",
                    "segmentation_robustness": "RECORDED_F_R3_SENSITIVITY_ONLY",
                    "segmentation_note": "F-R3 retained summary sensitivity, not threshold-specific compact time series; no volumetric reread in F-R4",
                }
            )
    return pd.DataFrame(rows)


def _apply_robustness_to_attribution(primary: pd.DataFrame, robust: pd.DataFrame) -> pd.DataFrame:
    out = primary.copy()
    for mech_key, mech_name in (
        ("eta_head_charge_evolution", "HEAD_CHARGE_EVOLUTION"),
        ("eta_head_acceleration", "HEAD_ACCELERATION"),
        ("eta_current_moment_redistribution", "CURRENT_MOMENT_REDISTRIBUTION"),
    ):
        status_col = mech_key.replace("eta_", "attribution_status_")
        robust_col = mech_key.replace("eta_", "derivative_robustness_")
        seg_col = mech_key.replace("eta_", "segmentation_robustness_")
        out[status_col] = "ATTRIBUTION_INTERPRET_WITH_CAUTION"
        out[robust_col] = "DERIVATIVE_SIGN_CONSISTENT"
        out[seg_col] = "RECORDED_F_R3_SENSITIVITY_ONLY"
        for idx, row in out.iterrows():
            sub = robust[
                (robust["stage"] == row["stage"])
                & (robust["band_name"] == row["band_name"])
                & (robust["mechanism"] == mech_name)
            ]
            if sub.empty:
                out.loc[idx, status_col] = "ATTRIBUTION_NOT_RESOLVED"
                out.loc[idx, robust_col] = "MISSING_DERIVATIVE_ROBUSTNESS"
                continue
            rob = str(sub.iloc[0]["derivative_robustness"])
            out.loc[idx, robust_col] = rob
            out.loc[idx, seg_col] = str(sub.iloc[0]["segmentation_robustness"])
            if rob == "ATTRIBUTION_NOT_ROBUST_TO_DERIVATIVE":
                out.loc[idx, status_col] = "ATTRIBUTION_NOT_RESOLVED"
    status_cols = [c for c in out.columns if c.startswith("attribution_status_")]
    out["attribution_status"] = np.where(
        (out[status_cols] == "ATTRIBUTION_INTERPRET_WITH_CAUTION").all(axis=1),
        "ATTRIBUTION_INTERPRET_WITH_CAUTION",
        "ATTRIBUTION_NOT_RESOLVED",
    )
    return out


def _stage5_status(fr3_summary: dict[str, object]) -> dict[str, object]:
    stage5 = next(d for d in fr3_summary["datasets"] if d["case_id"] == "F4-C-stage5-highfield-collision")
    payload = {
        "case_id": "F4-C-stage5-highfield-collision",
        "MECHANISM_ATTRIBUTION": "NOT_RESOLVED",
        "reason": "F_R3_MECHANISM_DECOMPOSITION_NOT_RESOLVED",
        "FULL_MAXWELL_REFERENCE_PENDING": True,
        "formal_eta_available": False,
        "eta_fields": None,
        "F_R3_mechanism_status": stage5["mechanism_status"],
        "F_R3_closure": stage5["primary_metrics"]["normalized_rms_closure"],
        "trusted_frequency_low_Hz": stage5["trusted_frequency_low_Hz"],
        "trusted_frequency_high_Hz": stage5["trusted_frequency_high_Hz"],
        "current_moment_full_jefimenko_cross_check": stage5["current_moment_full_jefimenko_cross_check"],
    }
    (OUT / "fr4_stage5_status.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def main() -> None:
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    fr3_summary = json.loads((FR3 / "fr3_mechanism_summary.json").read_text())
    stage4_csv, labels, trust, f4_summary = _load_stage4()
    primary_rows, primary_diag = _analyze_stage4_configuration(stage4_csv, labels, trust, derivative_tag="primary_radius5_degree2")
    alt_csv = _alternative_stage4_decomposition(stage4_csv)
    alt_rows, alt_diag = _analyze_stage4_configuration(alt_csv, labels, trust, derivative_tag="alternative_radius2_degree2")
    primary = pd.DataFrame(primary_rows)
    alternative = pd.DataFrame(alt_rows)
    robust = _robustness(primary, alternative, fr3_summary)
    primary = _apply_robustness_to_attribution(primary, robust)
    primary.to_csv(OUT / "fr4_stage4_attribution.csv", index=False)
    robust.to_csv(OUT / "fr4_stage4_robustness.csv", index=False)
    stage5 = _stage5_status(fr3_summary)
    stage4_fr3 = next(d for d in fr3_summary["datasets"] if d["case_id"] == "F4-P-stage4-left-isolated")
    summary = {
        "F_R4_STATUS": "PASS_CANDIDATE",
        "scope": "coherent mechanism-stage-frequency attribution relative to reduced-order current-moment field",
        "F_R3_checkpoint_commit": "6f23e99",
        "F_R3_input_flags": {
            "Stage4": stage4_fr3["mechanism_status"],
            "Stage5": stage5["F_R3_mechanism_status"],
            "Stage5_FULL_MAXWELL_REFERENCE_PENDING": True,
        },
        "CWT_settings": {"function": "morlet_cwt", "n_cycles": N_CYCLES, "frequency_count": N_CWT_FREQ},
        "observer_vector_m": OBSERVER_VECTOR_M.tolist(),
        "Stage4": {
            "trust": {
                "trusted_frequency_low_Hz": trust["trusted_frequency_low_Hz"],
                "trusted_frequency_high_Hz": trust["trusted_frequency_high_Hz"],
                "band_status": trust["band_status"],
            },
            "primary_diagnostics": primary_diag,
            "alternative_diagnostics": alt_diag,
            "F_R3_time_domain_closure": stage4_fr3["primary_metrics"]["normalized_rms_closure"],
            "current_moment_full_jefimenko_cross_check": f4_summary["M_dMdt_vs_jefimenko"],
            "attribution_rows": int(len(primary)),
            "caution_rows": int((primary["attribution_status"] == "ATTRIBUTION_INTERPRET_WITH_CAUTION").sum()),
            "not_resolved_rows": int((primary["attribution_status"] == "ATTRIBUTION_NOT_RESOLVED").sum()),
            "derivative_robustness": {
                "not_robust_count": int((robust["derivative_robustness"] == "ATTRIBUTION_NOT_ROBUST_TO_DERIVATIVE").sum()),
                "total_count": int(len(robust)),
            },
            "formal_status": "ATTRIBUTION_NOT_RESOLVED"
            if (primary["attribution_status"] == "ATTRIBUTION_NOT_RESOLVED").any()
            else "ATTRIBUTION_INTERPRET_WITH_CAUTION",
        },
        "Stage5": stage5,
    }
    summary["runtime_s"] = time.perf_counter() - started
    summary["peak_rss_kb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    summary["new_output_bytes"] = sum(p.stat().st_size for p in OUT.glob("*") if p.is_file())
    (OUT / "fr4_attribution_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

from __future__ import annotations

from dataclasses import asdict, dataclass, field

from .bands import PROJECT_BANDS, classify_band


@dataclass(frozen=True)
class RFTrustReport:
    source_valid: bool
    current_provenance: str
    continuity_status: str
    remap_status: str
    dt_s: float
    duration_s: float
    df_Hz: float
    nyquist_Hz: float
    low_frequency_interpretation_bound_Hz: float
    derivative_trust_frequency_Hz: float
    interpolation_trust_frequency_Hz: float
    mesh_trust_frequency_Hz: float
    sampling_trust_frequency_Hz: float
    far_field_valid: bool
    trusted_frequency_low_Hz: float
    trusted_frequency_high_Hz: float
    limiting_factor: str
    scientific_rf_valid: bool
    band_status: dict[str, str] = field(default_factory=dict)

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def build_trust_report(
    *,
    source_valid: bool,
    current_provenance: str,
    continuity_status: str,
    remap_status: str,
    dt_s: float,
    duration_s: float,
    derivative_trust_frequency_Hz: float,
    interpolation_trust_frequency_Hz: float,
    mesh_trust_frequency_Hz: float,
    sampling_trust_frequency_Hz: float,
    far_field_valid: bool,
    minimum_cycles_for_interpretation: float = 3.0,
) -> RFTrustReport:
    df = 1.0 / duration_s
    nyquist = 0.5 / dt_s
    low = minimum_cycles_for_interpretation / duration_s
    candidates = {
        "derivative": derivative_trust_frequency_Hz,
        "interpolation": interpolation_trust_frequency_Hz,
        "mesh": mesh_trust_frequency_Hz,
        "sampling": sampling_trust_frequency_Hz,
        "nyquist": nyquist,
    }
    high = min(candidates.values())
    limiting = min(candidates, key=candidates.get)
    scientific = bool(source_valid and far_field_valid and remap_status == "CONSERVATIVE" and high > low)
    bands = {band.name: classify_band(band, low, high) for band in PROJECT_BANDS}
    return RFTrustReport(
        source_valid=source_valid,
        current_provenance=current_provenance,
        continuity_status=continuity_status,
        remap_status=remap_status,
        dt_s=dt_s,
        duration_s=duration_s,
        df_Hz=df,
        nyquist_Hz=nyquist,
        low_frequency_interpretation_bound_Hz=low,
        derivative_trust_frequency_Hz=derivative_trust_frequency_Hz,
        interpolation_trust_frequency_Hz=interpolation_trust_frequency_Hz,
        mesh_trust_frequency_Hz=mesh_trust_frequency_Hz,
        sampling_trust_frequency_Hz=sampling_trust_frequency_Hz,
        far_field_valid=far_field_valid,
        trusted_frequency_low_Hz=low,
        trusted_frequency_high_Hz=high,
        limiting_factor=limiting,
        scientific_rf_valid=scientific,
        band_status=bands,
    )


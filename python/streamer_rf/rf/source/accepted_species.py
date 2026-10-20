"""Accepted-interval all-mobile-species source; no displacement or terminal current."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from .schema import SourceMetadata, SourceRecord


def load_accepted_species_packet(manifest_path, species="total"):
    path = Path(manifest_path)
    m = json.loads(path.read_text())
    required = ("geometry_id", "mesh_id", "operator_id", "reference_plane_id", "time_layer",
                "species_model_id", "units", "source_hash", "source_definition", "files")
    if any(not m.get(k) for k in required):
        raise ValueError("SOURCE_IDENTITY_MISSING")
    if m.get("version") != "I1-ACCEPTED-SPECIES-FV-1.0" or m["time_layer"] != "ACCEPTED_INTERVAL_PRE_FLUX":
        raise ValueError("SOURCE_TIME_OR_VERSION_MISMATCH")
    if m["source_definition"] != "TOTAL_CONDUCTION_FV_ELECTRON_PLUS_POSITIVE_ION" or m.get("includes_displacement_current") is not False:
        raise ValueError("SOURCE_CURRENT_IDENTITY_MISMATCH")
    if m["units"] != {"J":"A/m2", "rho":"C/m3", "time":"s", "volume":"m3"}:
        raise ValueError("SOURCE_UNIT_MISMATCH")
    if not {'species_center_source.csv','species_resolved_current_faces.csv'} <= set(m['files']):
        raise ValueError('CONSERVATIVE_FACE_SOURCE_REQUIRED')
    for name, expected in m["files"].items():
        if Path(name).name != name or hashlib.sha256((path.parent/name).read_bytes()).hexdigest() != expected:
            raise ValueError("SOURCE_FILE_IDENTITY_MISMATCH")
    if not 0 <= m["pre_time_s"] < m["post_time_s"]:
        raise ValueError("INVALID_SOURCE_INTERVAL")
    with (path.parent/"species_center_source.csv").open() as f:
        data = list(csv.DictReader(f))
    def col(k): return np.asarray([float(r[k]) for r in data])
    if not np.all(col("time_pre_s") == m["pre_time_s"]) or not np.all(col("time_post_s") == m["post_time_s"]):
        raise ValueError("SOURCE_TIME_LAYER_MISMATCH")
    if species not in ("total", "electron", "positive_ion"):
        raise ValueError("UNSUPPORTED_SPECIES_CURRENT")
    er,ez,pr,pz = (col(k) for k in ("Je_r_A_m2","Je_z_A_m2","Jp_r_A_m2","Jp_z_A_m2"))
    jr,jz = (er+pr,ez+pz) if species=="total" else (er,ez) if species=="electron" else (pr,pz)
    meta = SourceMetadata(case_id=m["case_id"], solver="StageC_NATIVE_FV", solver_version=m["source_hash"],
        time_s=m["pre_time_s"], coordinate_system="axisymmetric_rz", pressure_Pa=m["pressure_Pa"],
        temperature_K=m["temperature_K"], geometry_id=m["geometry_id"], voltage_state="PRESCRIBED_500V_NOT_TENG",
        photoionization="OPAQUE_ABSORBING_METAL", source_definition=m["source_definition"]+":"+species,
        units=m["units"], extra=m)
    n=len(data);r=col("r_m");dr=col("dr_m");dz=col("dz_m")
    return SourceRecord(meta, dict(cell_id=col("cell_id").astype(int),level=np.zeros(n,dtype=int),
        x_center=r,y_center=np.zeros(n),z_center=col("z_m"),dx=dr,dy=2*np.pi*r,dz=dz,
        cell_volume=col("volume_m3"),rho=col("rho_C_m3"),Jx=jr,Jy=np.zeros(n),Jz=jz))

"""Backend-independent H1 contracts and conservative engineering cost estimates."""
import re
import numpy as np

C0=299792458.
REFERENCE_PLANE="EXTERNAL_CEXT_TO_GAP_CGAP_PARALLEL_GSP"


def validate_backend(record):
    if record.get("upstream_url")!="https://github.com/thliebig/openEMS-Project.git":
        raise ValueError("UNAPPROVED_BACKEND")
    for key in ("project","openEMS","CSXCAD","fparser"):
        if not re.fullmatch(r"[0-9a-f]{40}",record["commits"].get(key,"")):
            raise ValueError("INVALID_BACKEND_PIN")
    return record["installation_status"]=="WORKING_RUNTIME_VERIFIED"


def validate_geometry(g):
    if g["kind"] not in ("H1_REFERENCE_GEOMETRY","PRODUCTION_GEOMETRY"):
        raise ValueError("INVALID_GEOMETRY_KIND")
    ext=np.asarray(g["domain_extent_m"],float)
    if ext.shape!=(3,) or np.any(~np.isfinite(ext)|(ext<=0)):
        raise ValueError("INVALID_GEOMETRY_EXTENT")
    if g["kind"]=="PRODUCTION_GEOMETRY" and not g.get("stage_b_geometry_reference"):
        raise ValueError("PRODUCTION_ELECTRODE_GEOMETRY_PENDING_STAGE_B")
    if (g["reference_plane_id"]!=REFERENCE_PLANE or g["voltage_reference"]!="GAP_NODE_MINUS_GROUND"
            or g["positive_current_direction"]!="EXTERNAL_TO_GAP"):
        raise ValueError("PORT_CONVENTION_MISMATCH")
    strategy,partition=g["microgap_representation"],g["cgap_partition"]
    if (strategy,partition) not in [("PORT_EQUIVALENT_MICROGAP","UPSTREAM_LUMPED_CGAP"),
                                  ("EXPLICIT_MICROGAP_GEOMETRY","GEOMETRY_RESOLVED_CGAP")]:
        raise ValueError("CGAP_PARTITION_CONFLICT")
    if g.get("add_duplicate_lumped_cgap",False):
        raise ValueError("DUPLICATE_CGAP_FORBIDDEN")
    return True


def wavelength_cell(fmax_Hz, cells_per_wavelength=20, eps_r=1.):
    values=np.array([fmax_Hz,cells_per_wavelength,eps_r],float)
    if np.any(~np.isfinite(values)|(values<=0)) or cells_per_wavelength<10:
        raise ValueError("INVALID_MESH_FREQUENCY")
    return C0/fmax_Hz/np.sqrt(eps_r)/cells_per_wavelength


def graded_axis(length_m,edges_m,max_cell_m,min_cell_m,growth=1.4):
    """Preserve explicit edges; split large neighboring cells until graded.

    An irrelevant microscopic edge is rejected, not silently promoted into
    a global timestep constraint. This is a geometry policy, not gap resolution.
    """
    p=np.asarray([length_m,max_cell_m,min_cell_m,growth])
    edges=np.asarray(edges_m,float)
    if np.any(~np.isfinite(p)) or not 0<min_cell_m<=max_cell_m or length_m<=0 or not 1<growth<=2:
        raise ValueError("INVALID_MESH_POLICY")
    if edges.ndim!=1 or np.any(~np.isfinite(edges)|(edges<0)|(edges>length_m)):
        raise ValueError("INVALID_GEOMETRY_EDGE")
    points=np.unique(np.r_[0.,edges,length_m])
    for _ in range(100):
        widths=np.diff(points)
        if widths.min()<min_cell_m*(1-1e-12):
            raise ValueError("MICROSCOPIC_REFINEMENT_REJECTED")
        split=widths>max_cell_m*(1+1e-12)
        split[:-1]|=widths[:-1]>growth*widths[1:]*(1+1e-12)
        split[1:]|=widths[1:]>growth*widths[:-1]*(1+1e-12)
        if not np.any(split): return points
        points=np.sort(np.r_[points,((points[:-1]+points[1:])/2)[split]])
    raise ValueError("GRADING_DID_NOT_CONVERGE")


def cfl_dt(minimum_xyz_m):
    d=np.asarray(minimum_xyz_m,float)
    if d.shape!=(3,) or np.any(~np.isfinite(d)|(d<=0)):
        raise ValueError("INVALID_CELL_SIZE")
    return float(1/(C0*np.sqrt(np.sum(1/d**2))))


def mesh_report(axes_m, fmax_Hz, growth_limit=1.4, minimum_allowed_m=1e-4):
    if len(axes_m)!=3:
        raise ValueError("INVALID_MESH_AXES")
    widths=[np.diff(np.asarray(a,float)) for a in axes_m]
    if any(d.size<2 or np.any(~np.isfinite(d)|(d<=0)) for d in widths):
        raise ValueError("INVALID_MESH_AXES")
    low=[float(d.min()) for d in widths];high=[float(d.max()) for d in widths]
    growth=max(float(np.max(np.maximum(d[1:]/d[:-1],d[:-1]/d[1:]))) for d in widths)
    if min(low)<minimum_allowed_m*(1-1e-10):
        raise ValueError("MICROSCOPIC_REFINEMENT_REJECTED")
    if max(high)>wavelength_cell(fmax_Hz)*(1+1e-10) or growth>growth_limit*(1+1e-10):
        raise ValueError("MESH_POLICY_VIOLATION")
    return dict(cells=int(np.prod([len(d) for d in widths])),min_xyz_m=low,max_xyz_m=high,
        growth_max=growth,estimated_CFL_dt_s=cfl_dt(low),target_max_frequency_Hz=fmax_Hz)


def cost_estimate(count_xyz,minimum_xyz_m,windows_s,bytes_per_cell=256):
    n=np.asarray(count_xyz,float)
    if n.shape!=(3,) or np.any(~np.isfinite(n)|(n<=0)|(n!=np.floor(n))):
        raise ValueError("INVALID_CELL_COUNTS")
    if bytes_per_cell<=0 or any(not np.isfinite(t) or t<=0 for t in windows_s.values()):
        raise ValueError("INVALID_COST_INPUT")
    cells=int(np.prod(n)); dt=cfl_dt(minimum_xyz_m)
    steps={k:int(np.ceil(v/dt)) for k,v in windows_s.items()}
    return dict(cells=cells,estimated_CFL_dt_s=dt,estimated_memory_bytes=cells*bytes_per_cell,
        bytes_per_cell_assumption=bytes_per_cell,steps=steps,
        cell_updates={k:cells*v for k,v in steps.items()},
        status="ENGINEERING_ESTIMATE_NOT_MEASURED_RENNIGS_TIMESTEP")


def multiband_decision(short_window_s,long_window_s):
    if not np.isfinite(short_window_s+long_window_s) or not 0<short_window_s<=long_window_s:
        raise ValueError("INVALID_WINDOWS")
    return dict(window_cost_ratio=long_window_s/short_window_s,
        decision="MULTIBAND_FULLWAVE_RUNS_REQUIRED" if long_window_s/short_window_s>=10 else "REVIEW_BAND_COST",
        rule="ENGINEERING_COST_RATIO_10_NOT_A_PHYSICAL_THRESHOLD",final_bands="PENDING_H2_H3")

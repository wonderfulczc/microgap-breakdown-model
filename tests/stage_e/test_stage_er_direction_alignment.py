import importlib.util
import math
import sys
from pathlib import Path


MODULE = Path(__file__).resolve().parents[2] / "solver3d" / "afivo_reference" / "stage_e" / "er" / "stage_er_direction_alignment.py"
spec = importlib.util.spec_from_file_location("stage_er_direction_alignment", MODULE)
er = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = er
spec.loader.exec_module(er)


def assert_close(a, b, tol=1e-12):
    assert abs(a - b) <= tol


def test_parallel_field_and_trajectory_zero_delta():
    out = er.direction_alignment((0.0, 0.0, -2.0), (0.0, 0.0, -1.0))
    assert out["field_direction_valid"]
    assert out["trajectory_valid"]
    assert_close(out["delta_theta_deg"], 0.0)
    assert_close(out["theta_E_deg"], 0.0)
    assert_close(out["theta_head_deg"], 0.0)


def test_perpendicular_field_and_trajectory():
    out = er.direction_alignment((1.0, 0.0, 0.0), (0.0, 0.0, -1.0))
    assert_close(out["delta_theta_deg"], 90.0)


def test_antiparallel_field_and_trajectory():
    out = er.direction_alignment((0.0, 0.0, 1.0), (0.0, 0.0, -1.0))
    assert_close(out["delta_theta_deg"], 180.0)


def test_arbitrary_known_3d_vectors():
    a = (1.0, 2.0, 3.0)
    b = (-4.0, 5.0, -6.0)
    expected = math.degrees(math.acos(sum(x * y for x, y in zip(a, b)) /
                                     (math.sqrt(sum(x * x for x in a)) * math.sqrt(sum(y * y for y in b)))))
    assert_close(er.vector_angle_deg(a, b), expected)


def test_clamp_safety():
    assert er.clamp(1.0 + 1e-13) == 1.0
    assert er.clamp(-1.0 - 1e-13) == -1.0


def test_polarity_orientation():
    raw_against_propagation = (0.0, 0.0, 1.0)
    out = er.direction_alignment(raw_against_propagation, (0.0, 0.0, -1.0), polarity_sign=-1.0)
    assert_close(out["delta_theta_deg"], 0.0)


def test_zero_displacement_invalid():
    out = er.direction_alignment((0.0, 0.0, -1.0), (0.0, 0.0, 0.0))
    assert out["field_direction_valid"]
    assert not out["trajectory_valid"]
    assert out["status"] == "ZERO_HEAD_DISPLACEMENT"
    assert math.isnan(out["delta_theta_deg"])


def test_zero_field_invalid():
    out = er.direction_alignment((0.0, 0.0, 0.0), (0.0, 0.0, -1.0))
    assert not out["field_direction_valid"]
    assert out["status"] == "ZERO_FIELD_MAGNITUDE"


def test_nonfinite_input_invalid():
    out = er.direction_alignment((math.nan, 0.0, -1.0), (0.0, 0.0, -1.0))
    assert not out["field_direction_valid"]
    assert out["status"] == "NONFINITE_INPUT"


def test_trajectory_resolution_valid():
    valid, status, over_dx = er.trajectory_resolution_status(2.0e-6, 1.0e-6)
    assert valid
    assert status == "VALID"
    assert_close(over_dx, 2.0)


def test_trajectory_resolution_under_resolved():
    valid, status, over_dx = er.trajectory_resolution_status(0.25e-6, 1.0e-6)
    assert not valid
    assert status == "INSUFFICIENT_SPATIAL_DISPLACEMENT"
    assert_close(over_dx, 0.25)


def test_trajectory_resolution_bad_dx():
    valid, status, over_dx = er.trajectory_resolution_status(1.0e-6, 0.0)
    assert not valid
    assert status == "INVALID_SPATIAL_RESOLUTION"
    assert math.isnan(over_dx)

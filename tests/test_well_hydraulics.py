import numpy as np
import pytest

from src.hydraulics import (
    annulus_area, annulus_pressure_loss, dpdl_annulus, ecd_dynamic, mean_velocity,
)
from src.trajectory import minimum_curvature
from src.units import G
from src.well_hydraulics import ecd_profile

RHO, Q, PV, YP = 1200.0, 0.03, 0.025, 6.0


def _vertical(td=2000.0):
    return minimum_curvature([0, td / 2, td], [0, 0, 0], [0, 0, 0])


def _j_type():
    md = [0, 500, 600, 700, 800, 900, 1000, 1500, 2000]
    inc = [0, 0, 10, 20, 30, 40, 45, 45, 45]
    azi = [0, 45, 45, 45, 45, 45, 45, 45, 45]
    return minimum_curvature(md, inc, azi)


def test_vertical_uniform_ecd_is_rho_plus_gradient_over_g():
    iv = [(0.0, 2000.0, 0.25, 0.127)]
    p = ecd_profile(_vertical(), RHO, Q, iv, PV, YP)
    v = float(mean_velocity(Q, annulus_area(0.25, 0.127)))
    grad = dpdl_annulus(RHO, v, 0.25, 0.127, PV, YP)[0]
    assert np.allclose(p['ECD'].iloc[1:], RHO + grad / G)
    assert np.isnan(p['ECD'].iloc[0])


def test_matches_annulus_pressure_loss_at_td():
    iv = [(0.0, 1000.0, 0.25, 0.18), (1000.0, 2000.0, 0.25, 0.127)]
    p = ecd_profile(_vertical(), RHO, Q, iv, PV, YP)
    total, _ = annulus_pressure_loss(
        RHO, Q, [(1000.0, 0.25, 0.18), (1000.0, 0.25, 0.127)], PV, YP)
    assert np.isclose(p['APL'].iloc[-1], total)
    assert np.isclose(p['ECD'].iloc[-1], ecd_dynamic(RHO, total, 2000.0))


def test_deviated_uses_md_for_friction_and_tvd_for_hydrostatic():
    iv = [(0.0, 2000.0, 0.25, 0.127)]
    dev, vert = _j_type(), _vertical(2000.0)
    p_dev = ecd_profile(dev, RHO, Q, iv, PV, YP)
    p_ver = ecd_profile(vert, RHO, Q, iv, PV, YP)
    assert np.isclose(p_dev['APL'].iloc[-1], p_ver['APL'].iloc[-1])    # MD sama
    assert p_dev['TVD'].iloc[-1] < p_ver['TVD'].iloc[-1]               # TVD lebih kecil
    assert p_dev['ECD'].iloc[-1] > p_ver['ECD'].iloc[-1]               # ECD lebih besar
    assert np.isclose(p_dev['ECD'].iloc[-1],
                      RHO + p_dev['APL'].iloc[-1] / (G * p_dev['TVD'].iloc[-1]))


def test_apl_increases_with_depth():
    iv = [(0.0, 800.0, 0.22, 0.127), (800.0, 2000.0, 0.2159, 0.127)]
    p = ecd_profile(_j_type(), RHO, Q, iv, PV, YP)
    assert np.all(np.diff(p['APL']) >= 0)


def test_md_eval_interpolates():
    iv = [(0.0, 2000.0, 0.25, 0.127)]
    p = ecd_profile(_j_type(), RHO, Q, iv, PV, YP, md_eval=[750.0, 1250.0])
    assert len(p) == 2 and p['TVD'].is_monotonic_increasing


def test_invalid_inputs():
    s = _vertical()
    with pytest.raises(ValueError):
        ecd_profile(s, RHO, Q, [(0.0, 1000.0, 0.25, 0.127)], PV, YP)     # tidak sampai TD
    with pytest.raises(ValueError):
        ecd_profile(s, RHO, Q, [(100.0, 2000.0, 0.25, 0.127)], PV, YP)  # tidak mulai dari 0
    with pytest.raises(ValueError):
        ecd_profile(s, RHO, 0.0, [(0.0, 2000.0, 0.25, 0.127)], PV, YP)  # q = 0
    with pytest.raises(ValueError):
        ecd_profile(s, RHO, Q, [(0.0, 2000.0, 0.1, 0.127)], PV, YP)     # lubang < pipa
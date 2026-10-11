import numpy as np
import pytest

from src.hydraulics import (
    bingham_dpdl_annulus, bingham_dpdl_pipe, dpdl_annulus, dpdl_pipe,
    flow_rate_window, turbulent_dpdl_annulus, turbulent_dpdl_pipe,
)

LUMPUR = [(1200.0, 0.025, 6.0), (1600.0, 0.02, 20.0), (1100.0, 0.01, 2.0)]


@pytest.mark.parametrize('rho, pv, yp', LUMPUR)
def test_dpdl_pipe_naik_terhadap_kecepatan(rho, pv, yp):
    v = np.linspace(0.1, 8.0, 800)
    g = [dpdl_pipe(rho, x, 0.1086, pv, yp)[0] for x in v]
    assert np.all(np.diff(g) >= 0)


@pytest.mark.parametrize('rho, pv, yp', LUMPUR)
def test_dpdl_annulus_naik_terhadap_kecepatan(rho, pv, yp):
    v = np.linspace(0.1, 8.0, 800)
    g = [dpdl_annulus(rho, x, 0.25, 0.127, pv, yp)[0] for x in v]
    assert np.all(np.diff(g) >= 0)


def test_dpdl_batas_rezim_pipa():
    rho, d, pv, yp = 1200.0, 0.1086, 0.025, 6.0
    g, rezim = dpdl_pipe(rho, 0.2, d, pv, yp)           # kecepatan rendah
    assert rezim == 'laminar'
    assert np.isclose(g, bingham_dpdl_pipe(0.2, d, pv, yp))
    g, rezim = dpdl_pipe(rho, 6.0, d, pv, yp)           # kecepatan tinggi
    assert rezim == 'turbulent'
    assert np.isclose(g, turbulent_dpdl_pipe(rho, 6.0, d, pv))


def test_dpdl_batas_rezim_annulus():
    rho, pv, yp = 1200.0, 0.025, 6.0
    g, rezim = dpdl_annulus(rho, 0.2, 0.25, 0.127, pv, yp)
    assert rezim == 'laminar'
    assert np.isclose(g, bingham_dpdl_annulus(0.2, 0.25, 0.127, pv, yp))
    g, rezim = dpdl_annulus(rho, 6.0, 0.25, 0.127, pv, yp)
    assert rezim == 'turbulent'
    assert np.isclose(g, turbulent_dpdl_annulus(rho, 6.0, 0.25, 0.127, pv))


def test_validasi_reologi():
    with pytest.raises(ValueError):
        dpdl_pipe(1200.0, 1.0, 0.1, 0.0, 6.0)                 # pv = 0
    with pytest.raises(ValueError):
        dpdl_annulus(1200.0, 1.0, 0.25, 0.127, 0.025, -1.0)   # yp negatif


def test_flow_window_contiguous():
    d13 = 13 / 32 * 0.0254
    res = flow_rate_window(
        1200.0, [(1000, 0.1086), (200, 0.0714)],
        [(1000, 0.25, 0.127), (200, 0.25, 0.18)], [d13] * 3, 0.025, 6.0,
        tvd=1200.0, ecd_max=1300.0, p_max=20.0e6, d_cut=0.006, rho_cut=2500.0,
        ft_min=0.5, p_surface=1.0e5)
    assert 'contiguous' in res
    if res['window'] is not None:
        assert res['contiguous'] is True
import numpy as np
import pytest

from src.hydraulics import (
    annulus_area, bingham_dpdl_annulus, bingham_dpdl_pipe, ecd,
    hydrostatic_pressure, mean_velocity, pipe_area,
)
from src.units import (
    FT_TO_M, G, PSI_TO_PA, cp_to_pas, ft_to_m, in_to_m, kgm3_to_ppg,
    lbf100ft2_to_pa, pa_to_psi, ppg_to_kgm3, psi_to_pa,
)


def test_hydrostatic_analytic():
    # 1200 kg/m3 pada 3000 m: 1200 x 9.80665 x 3000 = 35.30 MPa
    assert np.isclose(hydrostatic_pressure(1200.0, 3000.0), 1200.0 * G * 3000.0)
    assert np.isclose(hydrostatic_pressure(1200.0, 3000.0), 35.304e6, rtol=1e-4)


def test_hydrostatic_field_rule():
    # Rumus lapangan: 0.052 x ppg x ft = psi (0.052 adalah pembulatan)
    p = hydrostatic_pressure(ppg_to_kgm3(10.0), ft_to_m(10000.0))
    assert np.isclose(pa_to_psi(p), 0.052 * 10.0 * 10000.0, rtol=2e-3)


def test_ecd_field_example():
    # MW 10 ppg, dP annulus 200 psi, TVD 10000 ft: ECD ~ 10 + 200/(0.052 x 10000) = 10.385 ppg
    e = ecd(ppg_to_kgm3(10.0), psi_to_pa(200.0), ft_to_m(10000.0))
    assert np.isclose(kgm3_to_ppg(e), 10.385, atol=0.01)


def test_ecd_no_friction_equals_mw():
    assert np.isclose(ecd(1200.0, 0.0, 2500.0), 1200.0)


def test_ecd_rejects_zero_tvd():
    with pytest.raises(ValueError):
        ecd(1200.0, 1.0e5, 0.0)

def test_pipe_area_and_velocity():
    a = pipe_area(0.1)
    assert np.isclose(a, np.pi / 4 * 0.01)
    assert np.isclose(mean_velocity(0.01, a), 1.2732, atol=1e-3)


def test_annulus_area():
    assert np.isclose(annulus_area(0.2, 0.12), np.pi / 4 * (0.04 - 0.0144))
    with pytest.raises(ValueError):
        annulus_area(0.1, 0.12)               # lubang lebih kecil dari pipa


def test_bingham_pipe_hand_calc():
    # D=0.1 m, v=1 m/s, PV=0.02 Pa.s, YP=5 Pa
    # viskos 32*0.02*1/0.01 = 64 ; luluh 16*5/(3*0.1) = 266.667 ; total 330.667 Pa/m
    assert np.isclose(bingham_dpdl_pipe(1.0, 0.1, 0.02, 5.0), 330.6667, atol=1e-3)


def test_bingham_annulus_hand_calc():
    # lubang 0.2 m, pipa 0.12 m (celah 0.08), v=0.5 m/s
    # viskos 48*0.02*0.5/0.0064 = 75 ; luluh 6*5/0.08 = 375 ; total 450 Pa/m
    assert np.isclose(bingham_dpdl_annulus(0.5, 0.2, 0.12, 0.02, 5.0), 450.0)


def test_bingham_pipe_matches_field_formula():
    # PV 20 cp, YP 10 lbf/100ft2, D 4 in, v 10 ft/s
    pv, yp = cp_to_pas(20.0), lbf100ft2_to_pa(10.0)
    d, v = in_to_m(4.0), ft_to_m(10.0)
    si_psi_ft = bingham_dpdl_pipe(v, d, pv, yp) * FT_TO_M / PSI_TO_PA
    field = 20.0 * 10.0 / (1500 * 4.0 ** 2) + 10.0 / (225 * 4.0)   # Bourgoyne, psi/ft
    assert np.isclose(si_psi_ft, field, rtol=5e-3)


def test_bingham_invalid_diameter():
    with pytest.raises(ValueError):
        bingham_dpdl_pipe(1.0, 0.0, 0.02, 5.0)
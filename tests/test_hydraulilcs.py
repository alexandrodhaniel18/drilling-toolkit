import numpy as np
import pytest

from src.hydraulics import ecd, hydrostatic_pressure
from src.units import (
    G, ft_to_m, kgm3_to_ppg, pa_to_psi, ppg_to_kgm3, psi_to_pa,
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
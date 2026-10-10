import numpy as np
import pytest

from src.mse import (
    bit_area, mse, mse_axial, mse_breakdown, mse_efficiency, mse_rotary,
)
from src.units import (
    PSI_TO_PA, fthr_to_ms, in_to_m, lbf_to_n, ftlbf_to_nm, mhr_to_ms,
    rpm_to_rps,
)


def test_unit_converters():
    assert np.isclose(lbf_to_n(1.0), 4.4482216, rtol=1e-7)
    assert np.isclose(ftlbf_to_nm(1.0), 1.3558179, rtol=1e-7)
    assert np.isclose(fthr_to_ms(3600.0), 0.3048)
    assert np.isclose(mhr_to_ms(36.0), 0.01)
    assert np.isclose(rpm_to_rps(120.0), 2.0)


def test_bit_area_hand_calc():
    assert np.isclose(bit_area(0.2), np.pi / 4 * 0.04)


def test_mse_axial_hand_calc():
    # WOB 50 kN, bit 0.2 m: 50000 / 0.0314159 = 1.5915 MPa
    assert np.isclose(mse_axial(50000.0, 0.2), 1.59155e6, rtol=1e-5)


def test_mse_rotary_hand_calc():
    # T 5000 N.m, 2 rev/s, ROP 0.005 m/s, bit 0.2 m: 2*pi*5000*2/(0.0314159*0.005) = 400 MPa
    assert np.isclose(mse_rotary(5000.0, 2.0, 0.005, 0.2), 4.0e8, rtol=1e-9)


def test_mse_total_and_breakdown():
    b = mse_breakdown(50000.0, 5000.0, 2.0, 0.005, 0.2)
    assert np.isclose(b['total'], 401.5915e6, rtol=1e-6)
    assert np.isclose(b['total'], b['axial'] + b['rotary'])
    assert np.isclose(b['rotary_share'], 0.99604, atol=1e-5)
    assert np.isclose(mse(50000.0, 5000.0, 2.0, 0.005, 0.2), b['total'])


def test_mse_matches_field_formula_teale():
    # WOB 20000 lbf, T 5000 ft.lbf, 100 rpm, ROP 30 ft/jam, bit 8.5 in
    # Rumus lapangan (psi): WOB/A + 120*pi*N*T/(A*ROP), A dalam in2
    a_in2 = np.pi / 4 * 8.5 ** 2
    field_psi = 20000.0 / a_in2 + 120 * np.pi * 100 * 5000 / (a_in2 * 30.0)
    si_psi = mse(lbf_to_n(20000.0), ftlbf_to_nm(5000.0), rpm_to_rps(100.0),
                 fthr_to_ms(30.0), in_to_m(8.5)) / PSI_TO_PA
    assert np.isclose(si_psi, field_psi, rtol=1e-6)
    assert np.isclose(field_psi, 111079.1, atol=0.1)


def test_mse_behaviour():
    base = mse_rotary(5000.0, 2.0, 0.005, 0.2)
    assert np.isclose(mse_rotary(5000.0, 2.0, 0.010, 0.2), base / 2)    # ROP 2x -> rotari 1/2
    assert np.isclose(mse(0.0, 5000.0, 2.0, 0.005, 0.2), base)           # tanpa WOB: murni rotari
    assert mse(50000.0, 5000.0, 2.0, 0.004, 0.2) > mse(50000.0, 5000.0, 2.0, 0.005, 0.2)


def test_mse_efficiency():
    assert np.isclose(mse_efficiency(2.0e8, 2.0e8), 1.0)
    assert np.isclose(mse_efficiency(401.5915e6, 100.0e6), 0.24901, atol=1e-5)


def test_mse_array_input():
    r = mse(np.array([50000.0, 60000.0]), np.array([5000.0, 6000.0]),
            np.array([2.0, 2.0]), np.array([0.005, 0.006]), 0.2)
    assert r.shape == (2,)
    assert np.allclose(r, [401.5915e6, 401.9099e6], rtol=1e-6)


def test_mse_invalid_inputs():
    with pytest.raises(ValueError):
        mse(50000.0, 5000.0, 2.0, 0.0, 0.2)         # ROP nol
    with pytest.raises(ValueError):
        mse(50000.0, 5000.0, 2.0, 0.005, 0.0)       # diameter nol
    with pytest.raises(ValueError):
        mse(-1.0, 5000.0, 2.0, 0.005, 0.2)          # WOB negatif
    with pytest.raises(ValueError):
        mse(50000.0, -5.0, 2.0, 0.005, 0.2)         # torsi negatif
    with pytest.raises(ValueError):
        mse_efficiency(2.0e8, 0.0)
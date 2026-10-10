import numpy as np
import pytest

from src.hydraulics import (
    annulus_area, annulus_pressure_loss, bingham_dpdl_annulus,
    bingham_dpdl_pipe, bit_pressure_drop, dpdl_annulus, dpdl_pipe, ecd,
    ecd_dynamic, hydraulic_power, hydrostatic_pressure, mean_velocity,
    pipe_area, pipe_pressure_loss, pump_pressure, reynolds_annulus,
    reynolds_pipe, total_flow_area, turbulent_dpdl_annulus,
    turbulent_dpdl_pipe,
)
from src.units import (
    FT_TO_M, G, PSI_TO_PA, cp_to_pas, ft_to_m, gpm_to_m3s, in_to_m,
    kgm3_to_ppg, lbf100ft2_to_pa, nozzle32_to_m, pa_to_psi, ppg_to_kgm3,
    psi_to_pa,
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

def test_reynolds_pipe_hand_calc():
    # rho 1200, v 1, D 0.1, PV 0.02, YP 5
    # mu_e = 0.02 + 5*0.1/6 = 0.103333 ; Re = 1200*1*0.1/0.103333 = 1161.29
    assert np.isclose(reynolds_pipe(1200, 1.0, 0.1, 0.02, 5.0), 1161.29, atol=0.01)


def test_reynolds_annulus_hand_calc():
    # rho 1200, v 0.5, lubang 0.2, pipa 0.12 (celah 0.08)
    # mu_e = 0.02 + 5*0.08/(8*0.5) = 0.12 ; Re = 1200*0.5*0.8165*0.08/0.12 = 326.6
    assert np.isclose(reynolds_annulus(1200, 0.5, 0.2, 0.12, 0.02, 5.0), 326.60, atol=0.01)


def test_turbulent_pipe_hand_calc():
    # rho 1200, v 3, D 0.1, PV 0.02: Re = 18000, f = 0.006829, dP/dL = 1475.07 Pa/m
    assert np.isclose(turbulent_dpdl_pipe(1200, 3.0, 0.1, 0.02), 1475.07, atol=0.01)


def test_turbulent_annulus_hand_calc():
    # rho 1200, v 2, lubang 0.2, pipa 0.12, PV 0.02 -> 1235.51 Pa/m
    assert np.isclose(turbulent_dpdl_annulus(1200, 2.0, 0.2, 0.12, 0.02), 1235.51, atol=0.01)


def test_reynolds_pipe_matches_field_formula():
    # 10 ppg, v 10 ft/s, D 4 in, PV 20 cp, YP 10 lbf/100ft2
    rho = ppg_to_kgm3(10.0)
    pv, yp = cp_to_pas(20.0), lbf100ft2_to_pa(10.0)
    re_si = reynolds_pipe(rho, ft_to_m(10.0), in_to_m(4.0), pv, yp)
    re_field = 928.0 * 10.0 * 10.0 * 4.0 / (20.0 + 6.66 * 10.0 * 4.0 / 10.0)
    assert np.isclose(re_si, re_field, rtol=5e-3)


def test_turbulent_pipe_matches_field_formula():
    # 10 ppg, v 15 ft/s, D 4 in, PV 20 cp. Field: f*rho*v^2/(25.8*d) psi/ft
    f = 0.0791 / (928.0 * 10.0 * 15.0 * 4.0 / 20.0) ** 0.25
    field = f * 10.0 * 15.0 ** 2 / (25.8 * 4.0)
    si = turbulent_dpdl_pipe(ppg_to_kgm3(10.0), ft_to_m(15.0), in_to_m(4.0),
                             cp_to_pas(20.0)) * FT_TO_M / PSI_TO_PA
    assert np.isclose(si, field, rtol=5e-3)


def test_regime_selection():
    g1, r1 = dpdl_pipe(1200, 1.0, 0.1, 0.02, 5.0)      # Re 1161 -> laminar
    assert r1 == 'laminar'
    assert np.isclose(g1, bingham_dpdl_pipe(1.0, 0.1, 0.02, 5.0))
    g2, r2 = dpdl_pipe(1200, 3.0, 0.1, 0.02, 5.0)      # Re 7535 -> turbulen
    assert r2 == 'turbulent'
    assert np.isclose(g2, turbulent_dpdl_pipe(1200, 3.0, 0.1, 0.02))
    _, r3 = dpdl_annulus(1200, 0.5, 0.2, 0.12, 0.02, 5.0)   # Re 327 -> laminar
    _, r4 = dpdl_annulus(1200, 2.0, 0.2, 0.12, 0.02, 5.0)   # Re 3484 -> turbulen
    assert (r3, r4) == ('laminar', 'turbulent')


def test_velocity_must_be_positive():
    with pytest.raises(ValueError):
        reynolds_pipe(1200, 0.0, 0.1, 0.02, 5.0)

def test_apl_single_section_hand_calc():
    # rho 1200, q 0.02, lubang 0.2, pipa 0.12 (celah 0.08), L 1000, PV 0.02, YP 5
    # A = 0.020106 -> v = 0.99472 ; dP/dL = 149.2 + 375 = 524.208 Pa/m
    total, det = annulus_pressure_loss(1200, 0.02, [(1000, 0.2, 0.12)], 0.02, 5.0)
    assert np.isclose(total, 524207.76, rtol=1e-6)
    assert det[0]['regime'] == 'laminar'
    assert np.isclose(det[0]['dpdl'] * 1000, total)


def test_apl_two_sections_sum():
    secs = [(1000, 0.25, 0.127), (200, 0.25, 0.18)]
    total, det = annulus_pressure_loss(1200, 0.03, secs, 0.025, 6.0)
    assert np.isclose(det[0]['dp'], 358019.44, rtol=1e-6)
    assert np.isclose(det[1]['dp'], 165012.71, rtol=1e-6)
    assert np.isclose(total, det[0]['dp'] + det[1]['dp'])
    assert [d['regime'] for d in det] == ['laminar', 'laminar']


def test_apl_turbulent_sections():
    secs = [(1000, 0.25, 0.127), (200, 0.25, 0.18)]
    total, det = annulus_pressure_loss(1200, 0.08, secs, 0.025, 6.0)
    assert [d['regime'] for d in det] == ['turbulent', 'turbulent']
    assert np.isclose(total, 1674127.1, rtol=1e-6)


def test_ecd_dynamic_hand_calc():
    # rho 1200, APL 2 MPa, TVD 2000 m -> 1200 + 2e6/(G*2000)
    assert np.isclose(ecd_dynamic(1200, 2.0e6, 2000), 1200 + 2.0e6 / (G * 2000))
    assert np.isclose(ecd_dynamic(1200, 2.0e6, 2000), 1301.97, atol=0.01)
    assert ecd_dynamic(1200, 0.0, 2000) == 1200       # tanpa gesekan, ECD = rho


def test_ecd_dynamic_from_apl_in_ppg():
    secs = [(1000, 0.25, 0.127), (200, 0.25, 0.18)]
    total, _ = annulus_pressure_loss(1200, 0.03, secs, 0.025, 6.0)
    e = ecd_dynamic(1200, total, 1200)
    assert np.isclose(e, 1244.45, atol=0.01)
    assert e > 1200                                    # ECD selalu >= densitas statis
    assert np.isclose(kgm3_to_ppg(e), 1244.45 / 119.8264, rtol=1e-3)


def test_apl_and_ecd_invalid_inputs():
    with pytest.raises(ValueError):
        annulus_pressure_loss(1200, 0.02, [], 0.02, 5.0)           # tanpa segmen
    with pytest.raises(ValueError):
        annulus_pressure_loss(1200, 0.02, [(1000, 0.1, 0.12)], 0.02, 5.0)
    with pytest.raises(ValueError):
        ecd_dynamic(1200, 1.0e5, 0.0)                              # TVD nol

def test_pipe_loss_laminar_hand_calc():
    # rho 1200, q 0.005, ID 0.1086, L 1000, PV 0.025, YP 6
    # v = 0.53978 ; dP/dL = 36.62 + 294.65 = 331.27 Pa/m
    total, det = pipe_pressure_loss(1200, 0.005, [(1000, 0.1086)], 0.025, 6.0)
    assert np.isclose(total, 331273.61, rtol=1e-6)
    assert det[0]['regime'] == 'laminar'


def test_pipe_loss_two_sections_turbulent():
    secs = [(1000, 0.1086), (200, 0.0714)]       # drill pipe + drill collar
    total, det = pipe_pressure_loss(1200, 0.03, secs, 0.025, 6.0)
    assert [d['regime'] for d in det] == ['turbulent', 'turbulent']
    assert np.isclose(det[0]['dp'], 1608572.54, rtol=1e-6)
    assert np.isclose(det[1]['dp'], 2358282.40, rtol=1e-6)
    assert np.isclose(total, 3966854.93, rtol=1e-6)


def test_bit_pressure_drop_hand_calc():
    # 1200*0.03^2 / (2*0.95^2*(4e-4)^2) = 1.08 / 2.888e-7 = 3.7396 MPa
    assert np.isclose(bit_pressure_drop(1200, 0.03, 4.0e-4), 3739612.19, rtol=1e-6)


def test_bit_drop_matches_field_formula():
    # 10 ppg, 500 gpm, 3 nozzle 13/32. Lapangan: 8.311e-5*rho*Q^2/(Cd^2*At^2) psi
    at_in2 = 3 * np.pi / 4 * (13 / 32) ** 2
    field = 8.311e-5 * 10.0 * 500.0 ** 2 / (0.95 ** 2 * at_in2 ** 2)
    tfa = total_flow_area([nozzle32_to_m(13)] * 3)
    si = bit_pressure_drop(ppg_to_kgm3(10.0), gpm_to_m3s(500.0), tfa) / PSI_TO_PA
    assert np.isclose(si, field, rtol=5e-3)


def test_nozzle_and_tfa():
    assert np.isclose(nozzle32_to_m(13), 0.0103188, atol=1e-6)
    assert np.isclose(total_flow_area([nozzle32_to_m(13)] * 3), 2.5088e-4, rtol=1e-4)


def test_pump_pressure_is_sum_of_parts():
    pipe = [(1000, 0.1086), (200, 0.0714)]
    ann = [(1000, 0.25, 0.127), (200, 0.25, 0.18)]
    noz = [nozzle32_to_m(13)] * 3
    r = pump_pressure(1200, 0.03, pipe, ann, noz, 0.025, 6.0, p_surface=1.0e5)
    assert np.isclose(r['pipe'], 3966854.93, rtol=1e-6)
    assert np.isclose(r['annulus'], 523032.15, rtol=1e-6)
    assert np.isclose(r['bit'], 9506396.29, rtol=1e-6)
    assert np.isclose(r['total'], 14096283.37, rtol=1e-6)
    assert np.isclose(r['total'], r['surface'] + r['pipe'] + r['bit'] + r['annulus'])
    assert np.isclose(r['bit'] / r['total'], 0.6744, atol=1e-4)


def test_hydraulic_power():
    assert np.isclose(hydraulic_power(14096283.37, 0.03), 422888.50, rtol=1e-6)
    assert np.isclose(hydraulic_power(9506396.29, 0.03), 285191.89, rtol=1e-6)


def test_bit_and_pipe_invalid_inputs():
    with pytest.raises(ValueError):
        total_flow_area([0.0103, 0.0])               # nozzle nol
    with pytest.raises(ValueError):
        bit_pressure_drop(1200, 0.03, 4.0e-4, cd=0.0)
    with pytest.raises(ValueError):
        pipe_pressure_loss(1200, 0.03, [], 0.025, 6.0)
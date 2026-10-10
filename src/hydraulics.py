"""Hidrolika lumpur bor: tekanan hidrostatik dan ECD (satuan SI).

Referensi: Bourgoyne dkk. (1986), Applied Drilling Engineering.
"""
import numpy as np

from src.units import G


def hydrostatic_pressure(mw, tvd):
    """Tekanan hidrostatik (Pa) = MW (kg/m3) x g x TVD (m).

    Bergantung pada TVD (kedalaman vertikal), bukan MD.
    """
    return np.asarray(mw, dtype=float) * G * np.asarray(tvd, dtype=float)


def ecd(mw, dp_annulus, tvd):
    """Equivalent circulating density (kg/m3).

    ECD = MW + dp_annulus / (g x TVD)
    dp_annulus: kehilangan tekanan gesekan di annulus (Pa) dari dasar sampai permukaan.
    """
    tvd = np.asarray(tvd, dtype=float)
    if np.any(tvd <= 0):
        raise ValueError('TVD harus > 0')
    return np.asarray(mw, dtype=float) + np.asarray(dp_annulus, dtype=float) / (G * tvd)

def _require_positive(x, name):
    """Ubah ke array float dan tolak nilai <= 0 (fungsi bantu, awalan _)."""
    x = np.asarray(x, dtype=float)
    if np.any(x <= 0):
        raise ValueError(f'{name} harus > 0')
    return x


def pipe_area(d):
    """Luas penampang dalam pipa (m2) dari diameter dalam d (m)."""
    d = _require_positive(d, 'diameter')
    return np.pi / 4 * d ** 2


def annulus_area(d_hole, d_pipe):
    """Luas penampang annulus (m2): lubang (atau ID casing) dikurangi OD pipa."""
    d_hole = _require_positive(d_hole, 'd_hole')
    d_pipe = _require_positive(d_pipe, 'd_pipe')
    if np.any(d_hole <= d_pipe):
        raise ValueError('d_hole harus > d_pipe')
    return np.pi / 4 * (d_hole ** 2 - d_pipe ** 2)


def mean_velocity(q, area):
    """Kecepatan rata-rata (m/s) = laju alir q (m3/s) / luas penampang (m2)."""
    return np.asarray(q, dtype=float) / _require_positive(area, 'area')


def bingham_dpdl_pipe(v, d, pv, yp):
    """Gradien tekanan (Pa/m) di dalam pipa, Bingham Plastic, aliran LAMINAR.

    v: kecepatan (m/s); d: diameter dalam (m); pv: PV (Pa.s); yp: YP (Pa).
    Pendekatan Bourgoyne dkk. (1986): 32*pv*v/d^2 + 16*yp/(3*d).
    """
    d = _require_positive(d, 'diameter')
    v = np.asarray(v, dtype=float)
    return 32.0 * pv * v / d ** 2 + 16.0 * yp / (3.0 * d)


def bingham_dpdl_annulus(v, d_hole, d_pipe, pv, yp):
    """Gradien tekanan (Pa/m) di annulus, Bingham Plastic, aliran LAMINAR.

    Pendekatan celah sempit: 48*pv*v/gap^2 + 6*yp/gap, gap = d_hole - d_pipe.
    """
    d_hole = _require_positive(d_hole, 'd_hole')
    d_pipe = _require_positive(d_pipe, 'd_pipe')
    gap = d_hole - d_pipe
    if np.any(gap <= 0):
        raise ValueError('d_hole harus > d_pipe')
    v = np.asarray(v, dtype=float)
    return 48.0 * pv * v / gap ** 2 + 6.0 * yp / gap

RE_CRIT = 2100.0                 # batas laminar/turbulen (penyederhanaan)
K_ANNULUS = np.sqrt(2.0 / 3.0)   # 0,8165: diameter ekuivalen annulus


def reynolds_pipe(rho, v, d, pv, yp):
    """Bilangan Reynolds Bingham di dalam pipa (tanpa satuan).

    rho: kg/m3; v: m/s; d: m; pv: Pa.s; yp: Pa.
    Viskositas efektif: mu_e = pv + yp*d/(6*v).
    """
    rho = _require_positive(rho, 'rho')
    v = _require_positive(v, 'v')
    d = _require_positive(d, 'diameter')
    mu_e = pv + yp * d / (6.0 * v)
    return rho * v * d / mu_e


def reynolds_annulus(rho, v, d_hole, d_pipe, pv, yp):
    """Bilangan Reynolds Bingham di annulus (tanpa satuan).

    Viskositas efektif: mu_e = pv + yp*gap/(8*v); diameter ekuivalen K*gap.
    """
    rho = _require_positive(rho, 'rho')
    v = _require_positive(v, 'v')
    d_hole = _require_positive(d_hole, 'd_hole')
    d_pipe = _require_positive(d_pipe, 'd_pipe')
    gap = d_hole - d_pipe
    if np.any(gap <= 0):
        raise ValueError('d_hole harus > d_pipe')
    mu_e = pv + yp * gap / (8.0 * v)
    return rho * v * K_ANNULUS * gap / mu_e


def _fanning_blasius(re):
    """Faktor gesekan Fanning, korelasi Blasius: 0,0791 / Re^0,25."""
    return 0.0791 / np.asarray(re, dtype=float) ** 0.25


def turbulent_dpdl_pipe(rho, v, d, pv):
    """Gradien tekanan (Pa/m) di pipa, aliran TURBULEN (Blasius, Re dari pv)."""
    rho = _require_positive(rho, 'rho')
    v = _require_positive(v, 'v')
    d = _require_positive(d, 'diameter')
    re = rho * v * d / pv
    return 2.0 * _fanning_blasius(re) * rho * v ** 2 / d


def turbulent_dpdl_annulus(rho, v, d_hole, d_pipe, pv):
    """Gradien tekanan (Pa/m) di annulus, aliran TURBULEN (Blasius)."""
    rho = _require_positive(rho, 'rho')
    v = _require_positive(v, 'v')
    d_hole = _require_positive(d_hole, 'd_hole')
    d_pipe = _require_positive(d_pipe, 'd_pipe')
    gap = d_hole - d_pipe
    if np.any(gap <= 0):
        raise ValueError('d_hole harus > d_pipe')
    d_eq = K_ANNULUS * gap
    re = rho * v * d_eq / pv
    return 2.0 * _fanning_blasius(re) * rho * v ** 2 / d_eq


def dpdl_pipe(rho, v, d, pv, yp):
    """Gradien tekanan di pipa; rezim dipilih otomatis. Hanya untuk SKALAR.

    Mengembalikan (dpdl dalam Pa/m, 'laminar' atau 'turbulent').
    """
    if reynolds_pipe(rho, v, d, pv, yp) < RE_CRIT:
        return float(bingham_dpdl_pipe(v, d, pv, yp)), 'laminar'
    return float(turbulent_dpdl_pipe(rho, v, d, pv)), 'turbulent'


def dpdl_annulus(rho, v, d_hole, d_pipe, pv, yp):
    """Gradien tekanan di annulus; rezim dipilih otomatis. Hanya untuk SKALAR."""
    if reynolds_annulus(rho, v, d_hole, d_pipe, pv, yp) < RE_CRIT:
        return float(bingham_dpdl_annulus(v, d_hole, d_pipe, pv, yp)), 'laminar'
    return float(turbulent_dpdl_annulus(rho, v, d_hole, d_pipe, pv)), 'turbulent'

def annulus_pressure_loss(rho, q, sections, pv, yp):
    """Kehilangan tekanan annulus total (Pa) dari beberapa segmen seri.

    rho: kg/m3; q: laju alir (m3/s); pv: Pa.s; yp: Pa.
    sections: daftar (panjang_MD_m, d_hole_m, d_pipe_m), satu tuple per segmen.
    Mengembalikan (total_Pa, rincian), rincian = daftar dict per segmen:
    length, v, re, regime, dpdl, dp.
    """
    if len(sections) == 0:
        raise ValueError('sections tidak boleh kosong')
    total = 0.0
    details = []
    for length, d_hole, d_pipe in sections:
        length = float(_require_positive(length, 'length'))
        v = float(mean_velocity(q, annulus_area(d_hole, d_pipe)))
        dpdl, regime = dpdl_annulus(rho, v, d_hole, d_pipe, pv, yp)
        dp = dpdl * length
        total += dp
        details.append({
            'length': length,
            'v': v,
            're': float(reynolds_annulus(rho, v, d_hole, d_pipe, pv, yp)),
            'regime': regime,
            'dpdl': dpdl,
            'dp': dp,
        })
    return total, details


def ecd_dynamic(rho, apl, tvd):
    """ECD dinamis (kg/m3) = rho + APL/(g*TVD). apl dalam Pa, tvd dalam m."""
    tvd = float(_require_positive(tvd, 'tvd'))
    return rho + apl / (G * tvd)

CD_BIT = 0.95                    # koefisien discharge nozzle (default API)


def pipe_pressure_loss(rho, q, sections, pv, yp):
    """Kehilangan tekanan di dalam drill string (Pa), segmen seri.

    sections: daftar (panjang_m, d_dalam_m). Mengembalikan (total_Pa, rincian).
    """
    if len(sections) == 0:
        raise ValueError('sections tidak boleh kosong')
    total = 0.0
    details = []
    for length, d in sections:
        length = float(_require_positive(length, 'length'))
        v = float(mean_velocity(q, pipe_area(d)))
        dpdl, regime = dpdl_pipe(rho, v, d, pv, yp)
        dp = dpdl * length
        total += dp
        details.append({
            'length': length,
            'v': v,
            're': float(reynolds_pipe(rho, v, d, pv, yp)),
            'regime': regime,
            'dpdl': dpdl,
            'dp': dp,
        })
    return total, details


def total_flow_area(nozzle_diameters):
    """TFA (m2) dari daftar diameter nozzle (m)."""
    d = _require_positive(nozzle_diameters, 'nozzle')
    return float(np.sum(np.pi / 4 * d ** 2))


def bit_pressure_drop(rho, q, tfa, cd=CD_BIT):
    """Penurunan tekanan di nozzle bit (Pa) = rho*q^2 / (2*cd^2*tfa^2)."""
    rho = _require_positive(rho, 'rho')
    q = _require_positive(q, 'q')
    tfa = _require_positive(tfa, 'tfa')
    cd = _require_positive(cd, 'cd')
    return rho * q ** 2 / (2.0 * cd ** 2 * tfa ** 2)


def hydraulic_power(p, q):
    """Daya hidrolik (W) = tekanan (Pa) x laju alir (m3/s)."""
    return float(p) * float(_require_positive(q, 'q'))


def pump_pressure(rho, q, pipe_sections, annulus_sections, nozzles,
                  pv, yp, p_surface=0.0, cd=CD_BIT):
    """Tekanan pompa total (Pa) dan rinciannya, dalam satu dict."""
    p_pipe, pipe_det = pipe_pressure_loss(rho, q, pipe_sections, pv, yp)
    p_ann, ann_det = annulus_pressure_loss(rho, q, annulus_sections, pv, yp)
    tfa = total_flow_area(nozzles)
    p_bit = float(bit_pressure_drop(rho, q, tfa, cd))
    return {
        'surface': p_surface, 'pipe': p_pipe, 'bit': p_bit, 'annulus': p_ann,
        'total': p_surface + p_pipe + p_bit + p_ann, 'tfa': tfa,
        'pipe_details': pipe_det, 'annulus_details': ann_det,
    }

CD_CUTTING = 1.5                 # koefisien hambat serpihan (rezim Newton)


def slip_velocity(d_cut, rho_cut, rho, cd=CD_CUTTING):
    """Kecepatan jatuh serpihan (m/s), rezim Newton, sumur vertikal.

    d_cut: diameter serpihan (m); rho_cut, rho: densitas serpihan, lumpur (kg/m3).
    """
    d_cut = _require_positive(d_cut, 'd_cut')
    rho = _require_positive(rho, 'rho')
    cd = _require_positive(cd, 'cd')
    rho_cut = np.asarray(rho_cut, dtype=float)
    if np.any(rho_cut <= rho):
        raise ValueError('rho_cut harus > rho')
    return np.sqrt(4.0 * G * d_cut * (rho_cut - rho) / (3.0 * cd * rho))


def transport_ratio(v_ann, v_slip):
    """Rasio transport FT = 1 - v_slip/v_ann. Negatif berarti serpihan jatuh."""
    v_ann = _require_positive(v_ann, 'v_ann')
    return 1.0 - np.asarray(v_slip, dtype=float) / v_ann


def min_flow_rate_hole_cleaning(rho, annulus_sections, d_cut, rho_cut,
                                ft_min=0.5, cd=CD_CUTTING):
    """Laju alir minimum (m3/s) agar FT >= ft_min di SEMUA segmen annulus."""
    if not 0.0 < ft_min < 1.0:
        raise ValueError('ft_min harus di antara 0 dan 1')
    if len(annulus_sections) == 0:
        raise ValueError('sections tidak boleh kosong')
    v_slip = float(slip_velocity(d_cut, rho_cut, rho, cd))
    v_min = v_slip / (1.0 - ft_min)
    q_min = 0.0
    for length, d_hole, d_pipe in annulus_sections:
        q_min = max(q_min, v_min * float(annulus_area(d_hole, d_pipe)))
    return q_min


def flow_rate_window(rho, pipe_sections, annulus_sections, nozzles, pv, yp,
                     tvd, ecd_max, p_max, d_cut, rho_cut, ft_min=0.5,
                     p_surface=0.0, q_grid=None):
    """Pindai laju alir dan tandai mana yang memenuhi tiga batas sekaligus.

    ecd_max (kg/m3), p_max (Pa). Mengembalikan dict: q, ecd, p_pump, ft,
    ok_clean, ok_ecd, ok_pump, ok, dan window = (q_min, q_max) atau None.
    """
    if not 0.0 < ft_min < 1.0:
        raise ValueError('ft_min harus di antara 0 dan 1')
    if q_grid is None:
        q_grid = np.linspace(0.005, 0.100, 191)
    q_grid = _require_positive(q_grid, 'q_grid')
    v_slip = float(slip_velocity(d_cut, rho_cut, rho))
    ecd_arr, p_arr, ft_arr = [], [], []
    for q in q_grid:
        r = pump_pressure(rho, q, pipe_sections, annulus_sections, nozzles,
                          pv, yp, p_surface)
        ecd_arr.append(ecd_dynamic(rho, r['annulus'], tvd))
        p_arr.append(r['total'])
        ft_arr.append(min(float(transport_ratio(d['v'], v_slip))
                          for d in r['annulus_details']))
    ecd_arr, p_arr, ft_arr = np.array(ecd_arr), np.array(p_arr), np.array(ft_arr)
    ok_clean = ft_arr >= ft_min
    ok_ecd = ecd_arr <= ecd_max
    ok_pump = p_arr <= p_max
    ok = ok_clean & ok_ecd & ok_pump
    window = (float(q_grid[ok].min()), float(q_grid[ok].max())) if ok.any() else None
    return {'q': q_grid, 'ecd': ecd_arr, 'p_pump': p_arr, 'ft': ft_arr,
            'ok_clean': ok_clean, 'ok_ecd': ok_ecd, 'ok_pump': ok_pump,
            'ok': ok, 'window': window}
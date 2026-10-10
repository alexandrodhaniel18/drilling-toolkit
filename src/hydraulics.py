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
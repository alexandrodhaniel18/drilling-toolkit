"""Hidrolika sepanjang sumur: menggabungkan trajectory (MD, TVD) dengan hidrolika.

Gesekan annulus bergantung pada panjang lintasan (MD), sedangkan hidrostatik
bergantung pada kedalaman vertikal (TVD). Modul ini memakai keduanya sekaligus.
Satuan SI: m, m3/s, Pa, kg/m3, Pa.s.
"""
import numpy as np
import pandas as pd

from src.hydraulics import annulus_area, dpdl_annulus, mean_velocity
from src.units import G


def _check_intervals(intervals):
    """Pastikan interval berurutan dari MD 0 tanpa celah; kembalikan MD paling dalam."""
    if len(intervals) == 0:
        raise ValueError('intervals tidak boleh kosong')
    prev_bottom = 0.0
    for top, bottom, _d_hole, _d_pipe in intervals:
        if not np.isclose(top, prev_bottom):
            raise ValueError('intervals harus berurutan, mulai dari MD 0, tanpa celah')
        if bottom <= top:
            raise ValueError('md_bottom harus > md_top')
        prev_bottom = bottom
    return prev_bottom


def ecd_profile(survey, rho, q, intervals, pv, yp, md_eval=None):
    """ECD pada setiap kedalaman (stasiun survei atau md_eval).

    survey: DataFrame berkolom MD dan TVD (hasil minimum_curvature).
    rho: densitas lumpur (kg/m3); q: laju alir (m3/s); pv (Pa.s); yp (Pa).
    intervals: daftar (md_atas, md_bawah, d_lubang, d_pipa) dalam meter, dari
        permukaan ke bawah, tanpa celah. d_lubang = diameter lubang atau ID casing.
    md_eval: MD yang dievaluasi; bawaannya MD stasiun survei.

    APL pada kedalaman MD_i = jumlah (gradien tekanan x panjang interval di atas
    MD_i). ECD_i = rho + APL_i / (g x TVD_i). TVD antar stasiun diinterpolasi linear.
    Mengembalikan DataFrame: MD, TVD, APL (Pa), ECD (kg/m3); ECD = NaN di TVD = 0.
    """
    if rho <= 0:
        raise ValueError('rho harus > 0')
    if q <= 0:
        raise ValueError('q harus > 0')
    md_max = _check_intervals(intervals)

    md_s = survey['MD'].to_numpy(dtype=float)
    tvd_s = survey['TVD'].to_numpy(dtype=float)
    md = md_s.copy() if md_eval is None else np.asarray(md_eval, dtype=float)
    if np.any(md < 0) or np.any(md > md_s[-1] + 1e-9):
        raise ValueError('md_eval di luar rentang survei')
    if np.any(md > md_max + 1e-9):
        raise ValueError('intervals tidak mencakup MD sedalam md_eval')
    tvd = np.interp(md, md_s, tvd_s)

    apl = np.zeros_like(md)
    for top, bottom, d_hole, d_pipe in intervals:
        v = float(mean_velocity(q, annulus_area(d_hole, d_pipe)))
        gradient, _regime = dpdl_annulus(rho, v, d_hole, d_pipe, pv, yp)
        apl += gradient * np.clip(md - top, 0.0, bottom - top)

    ecd = np.full_like(md, np.nan)
    ok = tvd > 0
    ecd[ok] = rho + apl[ok] / (G * tvd[ok])
    return pd.DataFrame({'MD': md, 'TVD': tvd, 'APL': apl, 'ECD': ecd})
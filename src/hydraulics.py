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
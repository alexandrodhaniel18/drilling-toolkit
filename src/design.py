"""Pembuatan survei sintetis untuk profil J-type (build and hold).

Konvensi sama dengan trajectory.py: meter dan derajat. Build dengan DLS
konstan dan azimuth konstan, sehingga inklinasi naik linear terhadap MD.
Referensi rumus: Bourgoyne dkk. (1986), Applied Drilling Engineering.
"""
import numpy as np
import pandas as pd


def design_j_well(kop, dls, inc_final, azi, td, step=30.0, dls_length=30.0):
    """Buat tabel survei (MD, INC, AZI) untuk sumur J-type.

    kop: kedalaman mulai membelok (m); dls: build rate (deg per dls_length);
    inc_final: inklinasi hold (deg); azi: azimuth (deg); td: MD akhir (m);
    step: jarak antar stasiun (m).
    """
    if dls <= 0:
        raise ValueError('dls harus > 0')
    if not (0 < inc_final <= 90):
        raise ValueError('inc_final harus di antara 0 dan 90 deg')

    eob = kop + inc_final / dls * dls_length      # MD akhir build
    if not (0 <= kop < eob < td):
        raise ValueError('syarat: 0 <= kop < eob < td')

    md = np.unique(np.concatenate((np.arange(0.0, td, step), [kop, eob, td])))
    inc = np.clip(dls * (md - kop) / dls_length, 0.0, inc_final)
    azi_col = np.full_like(md, float(azi))

    return pd.DataFrame({'MD': md, 'INC': inc, 'AZI': azi_col})
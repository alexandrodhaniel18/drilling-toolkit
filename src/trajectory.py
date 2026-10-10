"""Perhitungan trajectory sumur dengan metode minimum curvature.

Konvensi: inklinasi dari vertikal (deg), azimuth dari utara searah jarum
jam (deg), TVD positif ke bawah. Satuan panjang: meter.
Referensi: Sawaryn & Thorogood (2005), SPE 84246.
"""
import numpy as np
import pandas as pd


def minimum_curvature(md, inc_deg, azi_deg, dls_length=30.0):
    """Hitung TVD, N, E, dan DLS dari data survei.

    Stasiun pertama dianggap di permukaan (TVD = N = E = 0).
    dls_length: panjang acuan DLS dalam satuan MD (30 -> deg/30 m).
    """
    md = np.asarray(md, dtype=float)
    inc = np.radians(np.asarray(inc_deg, dtype=float))
    azi = np.radians(np.asarray(azi_deg, dtype=float))

    if not (md.shape == inc.shape == azi.shape):
        raise ValueError('md, inc, dan azi harus sama panjang')
    if md.size < 2:
        raise ValueError('minimal 2 stasiun survei')
    d_md = np.diff(md)
    if np.any(d_md <= 0):
        raise ValueError('MD harus naik tegas')

    i1, i2 = inc[:-1], inc[1:]
    a1, a2 = azi[:-1], azi[1:]

    cos_beta = np.cos(i2 - i1) - np.sin(i1) * np.sin(i2) * (1 - np.cos(a2 - a1))
    beta = np.arccos(np.clip(cos_beta, -1.0, 1.0))

    rf = np.ones_like(beta)
    curved = beta > 1e-8
    rf[curved] = 2.0 / beta[curved] * np.tan(beta[curved] / 2.0)

    d_n = d_md / 2 * (np.sin(i1) * np.cos(a1) + np.sin(i2) * np.cos(a2)) * rf
    d_e = d_md / 2 * (np.sin(i1) * np.sin(a1) + np.sin(i2) * np.sin(a2)) * rf
    d_tvd = d_md / 2 * (np.cos(i1) + np.cos(i2)) * rf
    dls = np.degrees(beta) / d_md * dls_length

    zero = np.array([0.0])
    return pd.DataFrame({
        'MD': md,
        'INC': np.degrees(inc),
        'AZI': np.degrees(azi),
        'TVD': np.concatenate((zero, np.cumsum(d_tvd))),
        'N': np.concatenate((zero, np.cumsum(d_n))),
        'E': np.concatenate((zero, np.cumsum(d_e))),
        'DLS': np.concatenate((zero, dls)),
    })

def vertical_section(df, section_azi_deg):
    """Vertical section: proyeksi (N, E) ke bidang berazimuth section_azi_deg."""
    th = np.radians(section_azi_deg)
    return df['N'] * np.cos(th) + df['E'] * np.sin(th)


def closure(df):
    """Jarak dan azimuth penutup (deg) di stasiun terakhir."""
    n, e = df['N'].iloc[-1], df['E'].iloc[-1]
    return float(np.hypot(n, e)), float(np.degrees(np.arctan2(e, n)) % 360)


def read_survey(path, dls_length=30.0):
    """Baca CSV berkolom MD, INC, AZI lalu hitung trajectory."""
    raw = pd.read_csv(path)
    missing = {'MD', 'INC', 'AZI'} - set(raw.columns)
    if missing:
        raise ValueError(f'kolom hilang: {missing}')
    return minimum_curvature(raw['MD'], raw['INC'], raw['AZI'], dls_length)

if __name__ == '__main__':
    md = [0, 500, 600, 700, 800, 900, 1000, 1500, 2000]
    inc = [0, 0, 10, 20, 30, 40, 45, 45, 45]
    azi = [0, 45, 45, 45, 45, 45, 45, 45, 45]
    print(minimum_curvature(md, inc, azi).round(2))
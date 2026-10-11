import csv

import numpy as np


def _require_positive(x, name):
    """Ubah ke array float dan tolak nilai <= 0."""
    x = np.asarray(x, dtype=float)
    if np.any(x <= 0):
        raise ValueError(f'{name} harus > 0')
    return x


def _require_nonnegative(x, name):
    """Ubah ke array float dan tolak nilai < 0 (nol diperbolehkan)."""
    x = np.asarray(x, dtype=float)
    if np.any(x < 0):
        raise ValueError(f'{name} tidak boleh negatif')
    return x


def bit_area(d_bit):
    """Luas dasar lubang (m2) dari diameter bit d_bit (m)."""
    d_bit = _require_positive(d_bit, 'd_bit')
    return np.pi / 4.0 * d_bit ** 2


def mse_axial(wob, d_bit):
    """Suku aksial MSE (Pa) = WOB / A. wob dalam N."""
    return _require_nonnegative(wob, 'wob') / bit_area(d_bit)


def mse_rotary(torque, rps, rop, d_bit):
    """Suku rotari MSE (Pa) = 2*pi*N*T / (A*ROP).

    torque: N.m; rps: putaran/detik; rop: m/s (harus > 0); d_bit: m.
    """
    torque = _require_nonnegative(torque, 'torque')
    rps = _require_nonnegative(rps, 'rps')
    rop = _require_positive(rop, 'rop')
    return 2.0 * np.pi * torque * rps / (bit_area(d_bit) * rop)


def mse(wob, torque, rps, rop, d_bit):
    """MSE total (Pa), rumus Teale (1965). Semua masukan SI; array didukung."""
    return mse_axial(wob, d_bit) + mse_rotary(torque, rps, rop, d_bit)


def mse_breakdown(wob, torque, rps, rop, d_bit):
    """Rincian MSE: dict berisi axial, rotary, total (Pa), dan rotary_share."""
    axial = mse_axial(wob, d_bit)
    rotary = mse_rotary(torque, rps, rop, d_bit)
    total = axial + rotary
    return {'axial': axial, 'rotary': rotary, 'total': total,
            'rotary_share': rotary / total}


def mse_efficiency(mse_value, ccs):
    """Efisiensi mekanis = CCS / MSE (tanpa satuan). Keduanya dalam Pa."""
    mse_value = _require_positive(mse_value, 'mse')
    ccs = _require_positive(ccs, 'ccs')
    return ccs / mse_value

# faktor ke SI: depth -> m, wob -> N, torque -> N.m, rpm -> putaran/detik, rop -> m/s
UNIT_FACTORS = {
    'depth': {'m': 1.0, 'ft': 0.3048},
    'wob': {'n': 1.0, 'kn': 1.0e3, 'lbf': 4.4482216152605,
            'klbf': 4448.2216152605},
    'torque': {'nm': 1.0, 'knm': 1.0e3, 'ftlbf': 1.3558179483314,
               'kftlbf': 1355.8179483314},
    'rpm': {'rpm': 1.0 / 60.0, 'rps': 1.0},
    'rop': {'m/s': 1.0, 'm/h': 1.0 / 3600.0, 'ft/h': 0.3048 / 3600.0},
}
OUTPUT_KEYS = {'depth': 'depth', 'wob': 'wob', 'torque': 'torque',
               'rpm': 'rps', 'rop': 'rop'}


def read_drilling_csv(path, columns, units):
    """Baca CSV data pemboran dan ubah ke SI.

    columns: {'depth': nama_kolom, 'wob': ..., 'torque': ..., 'rpm': ..., 'rop': ...}
    units:   {'depth': 'm', 'wob': 'kn', 'torque': 'knm', 'rpm': 'rpm', 'rop': 'm/h'}
    Satuan WAJIB disebut (tidak ada nilai bawaan) supaya tidak salah tebak.
    Sel kosong atau bukan angka menjadi NaN. Mengembalikan dict array SI dengan
    kunci depth (m), wob (N), torque (N.m), rps (putaran/s), rop (m/s).
    """
    for key in OUTPUT_KEYS:
        if key not in columns:
            raise ValueError(f'columns tidak memuat kunci {key!r}')
        if key not in units:
            raise ValueError(f'units tidak memuat kunci {key!r}')
        if units[key].lower() not in UNIT_FACTORS[key]:
            raise ValueError(
                f'satuan {units[key]!r} untuk {key} tidak dikenal; '
                f'pilihan: {sorted(UNIT_FACTORS[key])}')
    values = {key: [] for key in OUTPUT_KEYS}
    with open(path, newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        missing = [columns[k] for k in OUTPUT_KEYS
                   if columns[k] not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f'kolom tidak ada di file: {missing}')
        for row in reader:
            for key in OUTPUT_KEYS:
                try:
                    values[key].append(float(row[columns[key]]))
                except (TypeError, ValueError):
                    values[key].append(np.nan)
    data = {}
    for key, out in OUTPUT_KEYS.items():
        factor = UNIT_FACTORS[key][units[key].lower()]
        data[out] = np.array(values[key], dtype=float) * factor
    return data


def filter_on_bottom(data, min_wob=0.0, min_rps=0.0, min_rop=0.0, min_torque=0.0):
    """Buang baris yang bukan pemboran aktif (bit di atas dasar, sambungan, NaN).

    Baris dipertahankan jika semua nilainya terhingga dan wob, rps, rop, torque
    masing-masing LEBIH BESAR dari batas minimumnya.
    Mengembalikan (data_baru, mask), mask boolean sepanjang data awal.
    """
    mask = np.ones(len(data['depth']), dtype=bool)
    limits = {'wob': min_wob, 'rps': min_rps, 'rop': min_rop, 'torque': min_torque}
    for key in ('depth', 'wob', 'torque', 'rps', 'rop'):
        mask &= np.isfinite(data[key])
    for key, lim in limits.items():
        with np.errstate(invalid='ignore'):
            mask &= data[key] > lim
    return {k: v[mask] for k, v in data.items()}, mask


def mse_series(data, d_bit):
    """MSE (Pa) untuk setiap baris data (dict hasil read_drilling_csv)."""
    return mse(data['wob'], data['torque'], data['rps'], data['rop'], d_bit)


def rolling_median(x, window):
    """Median bergulir (window ganjil); tepi diperlakukan dengan nilai tepi."""
    x = np.asarray(x, dtype=float)
    if window < 1 or window % 2 == 0:
        raise ValueError('window harus bilangan ganjil >= 1')
    if x.ndim != 1:
        raise ValueError('x harus 1 dimensi')
    if window == 1:
        return x.copy()
    half = window // 2
    padded = np.pad(x, half, mode='edge')
    windows = np.lib.stride_tricks.sliding_window_view(padded, window)
    return np.median(windows, axis=1)
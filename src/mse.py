import csv

import numpy as np


def _require_positive(x, name):
    x = np.asarray(x, dtype=float)
    if np.any(x <= 0):
        raise ValueError(f'{name} harus > 0')
    return x


def _require_nonnegative(x, name):
    x = np.asarray(x, dtype=float)
    if np.any(x < 0):
        raise ValueError(f'{name} tidak boleh negatif')
    return x


def bit_area(d_bit):
    d_bit = _require_positive(d_bit, 'd_bit')
    return np.pi / 4.0 * d_bit ** 2


def mse_axial(wob, d_bit):
    return _require_nonnegative(wob, 'wob') / bit_area(d_bit)


def mse_rotary(torque, rps, rop, d_bit):
    torque = _require_nonnegative(torque, 'torque')
    rps = _require_nonnegative(rps, 'rps')
    rop = _require_positive(rop, 'rop')
    return 2.0 * np.pi * torque * rps / (bit_area(d_bit) * rop)


def mse(wob, torque, rps, rop, d_bit):
    return mse_axial(wob, d_bit) + mse_rotary(torque, rps, rop, d_bit)


def mse_breakdown(wob, torque, rps, rop, d_bit):
    axial = mse_axial(wob, d_bit)
    rotary = mse_rotary(torque, rps, rop, d_bit)
    total = axial + rotary
    return {'axial': axial, 'rotary': rotary, 'total': total, 'rotary_share': rotary / total}


def mse_efficiency(mse_value, ccs):
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


def baseline_trend(mse_values, window=101):
    """Tren dasar MSE = median bergulir lebar (Pa), untuk pembanding.

    window (ganjil) harus lebih dari 2x panjang anomali terpanjang yang ingin
    dideteksi (dalam jumlah titik), kalau tidak anomali ikut menarik tren.
    """
    return rolling_median(mse_values, window)


def find_high_mse_intervals(depth, mse_pa, trend, ratio=1.5, min_length=10.0,
                            rotary_share=None):
    """Cari interval kedalaman dengan MSE >= ratio x tren secara berurutan.

    depth (m) harus naik. Titik berurutan yang lolos syarat digabung menjadi
    satu interval; interval dibuang jika rentangnya (akhir - awal) < min_length (m).
    rotary_share (opsional, array) menambah kolom median porsi rotari.
    Mengembalikan daftar dict: start_depth, end_depth, length, points,
    mse_median, trend_median, ratio_median, rotary_share_median (atau None).
    """
    depth = np.asarray(depth, dtype=float)
    mse_pa = np.asarray(mse_pa, dtype=float)
    trend = np.asarray(trend, dtype=float)
    if not (depth.shape == mse_pa.shape == trend.shape) or depth.ndim != 1:
        raise ValueError('depth, mse_pa, dan trend harus 1 dimensi sama panjang')
    if rotary_share is not None:
        rotary_share = np.asarray(rotary_share, dtype=float)
        if rotary_share.shape != depth.shape:
            raise ValueError('rotary_share harus sama panjang dengan depth')
    if np.any(np.diff(depth) <= 0):
        raise ValueError('depth harus naik (strictly increasing)')
    if ratio <= 1.0:
        raise ValueError('ratio harus > 1')
    if min_length < 0:
        raise ValueError('min_length tidak boleh negatif')
    _require_positive(trend, 'trend')

    rel = mse_pa / trend
    flagged = rel >= ratio
    intervals = []
    i, n = 0, len(depth)
    while i < n:
        if not flagged[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and flagged[j + 1]:
            j += 1
        span = depth[j] - depth[i]
        if span >= min_length:
            sl = slice(i, j + 1)
            intervals.append({
                'start_depth': float(depth[i]),
                'end_depth': float(depth[j]),
                'length': float(span),
                'points': int(j - i + 1),
                'mse_median': float(np.median(mse_pa[sl])),
                'trend_median': float(np.median(trend[sl])),
                'ratio_median': float(np.median(rel[sl])),
                'rotary_share_median': (None if rotary_share is None
                                        else float(np.median(rotary_share[sl]))),
            })
        i = j + 1
    return intervals


def write_intervals_csv(intervals, path):
    """Tulis tabel ringkasan interval ke CSV (MSE dalam MPa). Mengembalikan path."""
    header = ['start_depth_m', 'end_depth_m', 'length_m', 'points',
              'mse_median_MPa', 'trend_median_MPa', 'ratio_median',
              'rotary_share_median']
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(header)
        for iv in intervals:
            share = iv['rotary_share_median']
            w.writerow([
                f"{iv['start_depth']:.1f}", f"{iv['end_depth']:.1f}",
                f"{iv['length']:.1f}", iv['points'],
                f"{iv['mse_median'] / 1e6:.2f}", f"{iv['trend_median'] / 1e6:.2f}",
                f"{iv['ratio_median']:.2f}",
                '' if share is None else f'{share:.3f}',
            ])
    return path
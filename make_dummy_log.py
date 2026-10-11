import csv
import os

import numpy as np


def make_dummy_log(path, n=400, seed=42):
    """Tulis log pemboran sintetis (BUKAN data lapangan) ke file CSV.

    Kolom: depth_m, wob_kN, torque_kNm, rpm, rop_mh.
    Ada formasi yang makin keras terhadap kedalaman, satu lapisan keras
    (stringer) di 1200-1230 m, beberapa baris sambungan (bit di atas dasar),
    dan dua sel kosong.
    """
    rng = np.random.default_rng(seed)
    depth = 1000.0 + np.arange(n, dtype=float)
    hard = np.clip(1.0 + (depth - 1000.0) / 800.0, 1.0, None)
    stringer = (depth >= 1200.0) & (depth <= 1230.0)
    hard = np.where(stringer, hard * 2.2, hard)

    wob = 60.0 + rng.normal(0.0, 2.0, n)                    # kN
    rpm = 110.0 + rng.normal(0.0, 3.0, n)
    rop = 36.0 / hard * (1.0 + rng.normal(0.0, 0.05, n))    # m/jam
    torque = 4.0 * hard ** 0.6 * (1.0 + rng.normal(0.0, 0.04, n))  # kN.m

    off = np.zeros(n, dtype=bool)
    off[::100] = True                                       # sambungan
    wob[off], rop[off], torque[off] = 0.0, 0.0, 1.0

    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['depth_m', 'wob_kN', 'torque_kNm', 'rpm', 'rop_mh'])
        for i in range(n):
            row = [depth[i], wob[i], torque[i], rpm[i], rop[i]]
            row = [f'{v:.3f}' for v in row]
            if i == 57:
                row[2] = ''                                 # sel kosong
            if i == 213:
                row[4] = ''
            w.writerow(row)
    return path


if __name__ == '__main__':
    print(make_dummy_log(os.path.join('data', 'dummy_drilling_log.csv')))

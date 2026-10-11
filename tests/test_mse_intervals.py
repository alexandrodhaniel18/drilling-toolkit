import csv

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pytest

from make_dummy_log import make_dummy_log
from src.mse import (
    baseline_trend, filter_on_bottom, find_high_mse_intervals, mse_rotary,
    mse_series, read_drilling_csv, write_intervals_csv,
)
from src.plots import plot_mse_intervals

COLS = {'depth': 'depth_m', 'wob': 'wob_kN', 'torque': 'torque_kNm',
        'rpm': 'rpm', 'rop': 'rop_mh'}
UNITS = {'depth': 'm', 'wob': 'kn', 'torque': 'knm', 'rpm': 'rpm', 'rop': 'm/h'}


def _synthetic():
    depth = np.arange(0.0, 100.0)            # 0..99 m, langkah 1 m
    m = np.full(100, 100.0e6)                # MSE dasar 100 MPa
    m[20:35] = 200.0e6                       # interval 20..34 m, rasio 2.0
    m[60:64] = 300.0e6                       # lonjakan pendek (rentang 3 m)
    m[80:96] = 140.0e6                       # rasio 1.4, di bawah ambang 1.5
    trend = np.full(100, 100.0e6)
    return depth, m, trend


def test_find_intervals_hand_calc():
    depth, m, trend = _synthetic()
    iv = find_high_mse_intervals(depth, m, trend, ratio=1.5, min_length=10.0)
    assert len(iv) == 1                      # lonjakan pendek dan rasio 1.4 tidak lolos
    a = iv[0]
    assert a['start_depth'] == 20.0 and a['end_depth'] == 34.0
    assert a['length'] == 14.0 and a['points'] == 15
    assert np.isclose(a['mse_median'], 200.0e6)
    assert np.isclose(a['trend_median'], 100.0e6)
    assert np.isclose(a['ratio_median'], 2.0)
    assert a['rotary_share_median'] is None


def test_min_length_boundary_and_ratio_threshold():
    depth, m, trend = _synthetic()
    assert len(find_high_mse_intervals(depth, m, trend, 1.5, min_length=14.0)) == 1   # sama dengan: lolos
    assert len(find_high_mse_intervals(depth, m, trend, 1.5, min_length=14.5)) == 0
    # ambang diturunkan ke 1.4: interval 80..95 (rentang 15 m) ikut lolos (>= ambang)
    iv = find_high_mse_intervals(depth, m, trend, 1.4, 10.0)
    assert [(a['start_depth'], a['end_depth']) for a in iv] == [(20.0, 34.0), (80.0, 95.0)]


def test_rotary_share_median():
    depth, m, trend = _synthetic()
    share = np.full(100, 0.9)
    share[20:35] = 0.95
    iv = find_high_mse_intervals(depth, m, trend, 1.5, 10.0, rotary_share=share)
    assert np.isclose(iv[0]['rotary_share_median'], 0.95)


def test_baseline_trend_ignores_short_anomaly():
    _, m, _ = _synthetic()
    flat = np.full(100, 100.0e6)
    flat[20:35] = 200.0e6                    # 15 titik < setengah jendela 51
    assert np.allclose(baseline_trend(flat, 51), 100.0e6)
    # jendela terlalu kecil (11) ikut tertarik anomali 15 titik
    assert baseline_trend(flat, 11)[27] == 200.0e6


def test_write_intervals_csv(tmp_path):
    depth, m, trend = _synthetic()
    iv = find_high_mse_intervals(depth, m, trend, 1.5, 10.0)
    p = write_intervals_csv(iv, str(tmp_path / 'iv.csv'))
    rows = list(csv.reader(open(p, encoding='utf-8')))
    assert rows[0][0] == 'start_depth_m' and rows[0][4] == 'mse_median_MPa'
    assert rows[1] == ['20.0', '34.0', '14.0', '15', '200.00', '100.00', '2.00', '']
    assert len(rows) == 2


def test_invalid_inputs():
    depth, m, trend = _synthetic()
    with pytest.raises(ValueError):
        find_high_mse_intervals(depth[::-1], m, trend)           # depth turun
    with pytest.raises(ValueError):
        find_high_mse_intervals(depth, m[:50], trend)            # beda panjang
    with pytest.raises(ValueError):
        find_high_mse_intervals(depth, m, trend, ratio=1.0)
    with pytest.raises(ValueError):
        find_high_mse_intervals(depth, m, trend, min_length=-1.0)
    with pytest.raises(ValueError):
        find_high_mse_intervals(depth, m, np.zeros(100))         # tren nol
    with pytest.raises(ValueError):
        baseline_trend(m, 10)                                    # jendela genap


def test_dummy_log_detects_stringer(tmp_path):
    p = make_dummy_log(str(tmp_path / 'd.csv'), n=400, seed=7)
    d, _ = filter_on_bottom(read_drilling_csv(p, COLS, UNITS))
    m = mse_series(d, 0.2159)
    rs = mse_rotary(d['torque'], d['rps'], d['rop'], 0.2159) / m
    iv = find_high_mse_intervals(d['depth'], m, baseline_trend(m, 101), 1.5, 10.0, rs)
    assert len(iv) == 1
    assert 1199 <= iv[0]['start_depth'] <= 1203
    assert 1228 <= iv[0]['end_depth'] <= 1231
    assert iv[0]['ratio_median'] > 3.0
    assert iv[0]['rotary_share_median'] > 0.99


def test_plot_mse_intervals():
    depth, m, trend = _synthetic()
    iv = find_high_mse_intervals(depth, m, trend, 1.4, 10.0)
    fig = plot_mse_intervals(depth, m, trend, iv, window=5)
    ax = fig.axes[0]
    assert len(ax.patches) == 2                              # satu arsiran per interval
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert 'Tren dasar' in labels and 'Interval MSE tinggi' in labels
    plt.close(fig)
    with pytest.raises(ValueError):
        plot_mse_intervals(depth, m, trend, iv, unit='bar')
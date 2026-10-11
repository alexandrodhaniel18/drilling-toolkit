import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pytest

from make_dummy_log import make_dummy_log
from src.mse import (
    filter_on_bottom, mse_series, read_drilling_csv, rolling_median,
)
from src.plots import plot_mse_depth

COLS = {'depth': 'depth_m', 'wob': 'wob_kN', 'torque': 'torque_kNm',
        'rpm': 'rpm', 'rop': 'rop_mh'}
UNITS = {'depth': 'm', 'wob': 'kn', 'torque': 'knm', 'rpm': 'rpm', 'rop': 'm/h'}


def _write(tmp_path, text):
    p = tmp_path / 'log.csv'
    p.write_text(text, encoding='utf-8')
    return str(p)


def test_read_converts_units(tmp_path):
    p = _write(tmp_path,
               'depth_m,wob_kN,torque_kNm,rpm,rop_mh\n'
               '1000,50,5,120,18\n')
    d = read_drilling_csv(p, COLS, UNITS)
    assert np.isclose(d['depth'][0], 1000.0)
    assert np.isclose(d['wob'][0], 50000.0)       # kN -> N
    assert np.isclose(d['torque'][0], 5000.0)     # kN.m -> N.m
    assert np.isclose(d['rps'][0], 2.0)           # 120 rpm -> 2 rev/s
    assert np.isclose(d['rop'][0], 0.005)         # 18 m/jam -> m/s


def test_read_field_units(tmp_path):
    p = _write(tmp_path, 'd,w,t,n,r\n3280.84,10,2,60,100\n')
    cols = {'depth': 'd', 'wob': 'w', 'torque': 't', 'rpm': 'n', 'rop': 'r'}
    units = {'depth': 'ft', 'wob': 'klbf', 'torque': 'kftlbf', 'rpm': 'rpm', 'rop': 'ft/h'}
    d = read_drilling_csv(p, cols, units)
    assert np.isclose(d['depth'][0], 1000.0, atol=0.01)
    assert np.isclose(d['wob'][0], 44482.216, rtol=1e-6)
    assert np.isclose(d['torque'][0], 2711.636, rtol=1e-6)
    assert np.isclose(d['rop'][0], 100 * 0.3048 / 3600)


def test_read_blank_cells_become_nan(tmp_path):
    p = _write(tmp_path,
               'depth_m,wob_kN,torque_kNm,rpm,rop_mh\n'
               '1000,50,,120,18\n'
               '1001,50,5,120,abc\n')
    d = read_drilling_csv(p, COLS, UNITS)
    assert np.isnan(d['torque'][0])
    assert np.isnan(d['rop'][1])


def test_read_invalid_arguments(tmp_path):
    p = _write(tmp_path, 'depth_m,wob_kN,torque_kNm,rpm,rop_mh\n1000,50,5,120,18\n')
    with pytest.raises(ValueError):
        read_drilling_csv(p, COLS, {**UNITS, 'rop': 'km/h'})      # satuan tidak dikenal
    with pytest.raises(ValueError):
        read_drilling_csv(p, {**COLS, 'rop': 'tidak_ada'}, UNITS)  # kolom tidak ada
    with pytest.raises(ValueError):
        read_drilling_csv(p, COLS, {k: v for k, v in UNITS.items() if k != 'wob'})


def test_filter_on_bottom(tmp_path):
    p = _write(tmp_path,
               'depth_m,wob_kN,torque_kNm,rpm,rop_mh\n'
               '1000,50,5,120,18\n'      # baik
               '1001,0,1,120,0\n'        # sambungan: wob 0, rop 0
               '1002,50,,120,18\n'       # torsi kosong
               '1003,50,5,120,18\n')     # baik
    d = read_drilling_csv(p, COLS, UNITS)
    f, mask = filter_on_bottom(d)
    assert mask.tolist() == [True, False, False, True]
    assert np.allclose(f['depth'], [1000.0, 1003.0])


def test_mse_series_hand_calc(tmp_path):
    p = _write(tmp_path,
               'depth_m,wob_kN,torque_kNm,rpm,rop_mh\n'
               '1000,50,5,120,18\n')
    d = read_drilling_csv(p, COLS, UNITS)
    assert np.isclose(mse_series(d, 0.2)[0], 401.5915e6, rtol=1e-6)


def test_rolling_median():
    x = np.array([1.0, 100.0, 3.0, 4.0, 5.0])
    assert np.allclose(rolling_median(x, 3), [1.0, 3.0, 4.0, 4.0, 5.0])
    assert np.allclose(rolling_median(x, 1), x)
    with pytest.raises(ValueError):
        rolling_median(x, 4)


def test_dummy_log_pipeline(tmp_path):
    p = make_dummy_log(str(tmp_path / 'dummy.csv'), n=400, seed=1)
    d = read_drilling_csv(p, COLS, UNITS)
    assert len(d['depth']) == 400
    f, mask = filter_on_bottom(d)
    assert 380 < len(f['depth']) < 400                     # sambungan + 2 sel kosong dibuang
    assert np.all(np.isfinite(f['rop'])) and np.all(f['rop'] > 0)
    m = mse_series(f, 0.2159)
    assert np.all(m > 0)
    inside = (f['depth'] >= 1200) & (f['depth'] <= 1230)
    assert np.median(m[inside]) > 1.5 * np.median(m[~inside])   # lapisan keras terdeteksi


def test_plot_mse_depth(tmp_path):
    depth = np.arange(1000.0, 1010.0)
    m = np.linspace(1.0e8, 2.0e8, 10)
    fig = plot_mse_depth(depth, m, window=3)
    ax = fig.axes[0]
    assert ax.yaxis_inverted()
    assert np.allclose(ax.lines[0].get_xdata(), m * 1e-6)
    assert np.allclose(ax.lines[0].get_ydata(), depth)
    assert len(ax.lines) == 2                              # mentah + median bergulir
    plt.close(fig)
    fig = plot_mse_depth(depth, m, unit='psi')
    assert 'psi' in fig.axes[0].get_xlabel()
    assert len(fig.axes[0].lines) == 1
    plt.close(fig)
    with pytest.raises(ValueError):
        plot_mse_depth(depth, m, unit='bar')
    with pytest.raises(ValueError):
        plot_mse_depth(depth, m[:5])

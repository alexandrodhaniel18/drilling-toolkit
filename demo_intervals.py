import os

import matplotlib.pyplot as plt

from make_dummy_log import make_dummy_log
from src.mse import (
    baseline_trend, filter_on_bottom, find_high_mse_intervals, mse_rotary,
    mse_series, read_drilling_csv, write_intervals_csv,
)
from src.plots import plot_mse_intervals
from src.units import in_to_m

PATH = os.path.join('data', 'dummy_drilling_log.csv')
if not os.path.exists(PATH):
    make_dummy_log(PATH)

columns = {'depth': 'depth_m', 'wob': 'wob_kN', 'torque': 'torque_kNm',
           'rpm': 'rpm', 'rop': 'rop_mh'}
units = {'depth': 'm', 'wob': 'kn', 'torque': 'knm', 'rpm': 'rpm', 'rop': 'm/h'}

D_BIT = in_to_m(8.5)
WINDOW = 101        # harus > 2x panjang anomali terpanjang (dalam titik)
RATIO = 1.5         # interval ditandai jika MSE >= 1.5 x tren
MIN_LENGTH = 10.0   # rentang kedalaman minimum interval (m)

data, mask = filter_on_bottom(read_drilling_csv(PATH, columns, units))
m = mse_series(data, D_BIT)
trend = baseline_trend(m, WINDOW)
share = mse_rotary(data['torque'], data['rps'], data['rop'], D_BIT) / m
intervals = find_high_mse_intervals(data['depth'], m, trend, RATIO, MIN_LENGTH, share)

os.makedirs('outputs', exist_ok=True)
write_intervals_csv(intervals, os.path.join('outputs', 'mse_intervals.csv'))
for iv in intervals:
    print(f"{iv['start_depth']:.0f}-{iv['end_depth']:.0f} m: MSE median "
          f"{iv['mse_median'] / 1e6:.0f} MPa, {iv['ratio_median']:.2f}x tren")

fig = plot_mse_intervals(data['depth'], m, trend, intervals, window=11,
                         title='Interval MSE tinggi (data dummy)')
fig.savefig(os.path.join('outputs', 'mse_intervals.png'), dpi=200, bbox_inches='tight')
plt.show()
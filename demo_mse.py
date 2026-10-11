import os

import matplotlib.pyplot as plt

from make_dummy_log import make_dummy_log
from src.mse import filter_on_bottom, mse_series, read_drilling_csv
from src.plots import plot_mse_depth
from src.units import in_to_m

PATH = os.path.join('data', 'dummy_drilling_log.csv')
if not os.path.exists(PATH):
    make_dummy_log(PATH)

columns = {'depth': 'depth_m', 'wob': 'wob_kN', 'torque': 'torque_kNm',
           'rpm': 'rpm', 'rop': 'rop_mh'}
units = {'depth': 'm', 'wob': 'kn', 'torque': 'knm', 'rpm': 'rpm', 'rop': 'm/h'}

raw = read_drilling_csv(PATH, columns, units)
data, mask = filter_on_bottom(raw)
print(f'{mask.sum()} dari {len(mask)} baris dipakai')

m = mse_series(data, in_to_m(8.5))
fig = plot_mse_depth(data['depth'], m, window=11,
                     title='MSE terhadap kedalaman (data dummy)')
os.makedirs('outputs', exist_ok=True)
fig.savefig(os.path.join('outputs', 'mse_depth.png'), dpi=200, bbox_inches='tight')
plt.show()

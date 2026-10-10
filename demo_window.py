from src.hydraulics import flow_rate_window
from src.plots import plot_flow_window
from src.units import nozzle32_to_m

pipe = [(1000, 0.1086), (200, 0.0714)]          # (panjang m, ID m)
ann = [(1000, 0.25, 0.127), (200, 0.25, 0.18)]  # (panjang m, d lubang m, d pipa m)
nozzles = [nozzle32_to_m(13)] * 3

res = flow_rate_window(1200, pipe, ann, nozzles, 0.025, 6.0,
                       tvd=1200, ecd_max=1300.0, p_max=20.0e6,
                       d_cut=0.006, rho_cut=2500, ft_min=0.5, p_surface=1.0e5)
print('Jendela layak (m3/s):', res['window'])

fig = plot_flow_window(res, 1200, 1300.0, 20.0e6, unit='gpm',
                       title='Jendela laju alir (data dummy)')
fig.savefig('outputs/ecd_window.png', dpi=200, bbox_inches='tight')

import matplotlib.pyplot as plt
plt.show()

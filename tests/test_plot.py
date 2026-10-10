import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pytest

from src.hydraulics import flow_rate_window
from src.plots import plot_flow_window
from src.units import GPM_TO_M3S, nozzle32_to_m


def _result(p_max=20.0e6):
    pipe = [(1000, 0.1086), (200, 0.0714)]
    ann = [(1000, 0.25, 0.127), (200, 0.25, 0.18)]
    return flow_rate_window(1200, pipe, ann, [nozzle32_to_m(13)] * 3, 0.025, 6.0,
                            tvd=1200, ecd_max=1300.0, p_max=p_max, d_cut=0.006,
                            rho_cut=2500, ft_min=0.5, p_surface=1.0e5)


def test_plot_structure_and_data():
    res = _result()
    fig = plot_flow_window(res, 1200, 1300.0, 20.0e6)
    assert len(fig.axes) == 3
    assert np.allclose(fig.axes[0].lines[0].get_xdata(), res['q'] * 1000.0)
    assert np.allclose(fig.axes[0].lines[0].get_ydata(), res['ecd'])
    assert np.allclose(fig.axes[1].lines[0].get_ydata(), res['p_pump'] / 1e6)
    assert np.allclose(fig.axes[2].lines[0].get_ydata(), res['ft'])
    assert len(fig.axes[0].patches) == 1          # jendela layak diarsir
    plt.close(fig)


def test_plot_gpm_axis():
    res = _result()
    fig = plot_flow_window(res, 1200, 1300.0, 20.0e6, unit='gpm')
    assert 'gpm' in fig.axes[2].get_xlabel()
    assert np.allclose(fig.axes[0].lines[0].get_xdata(), res['q'] / GPM_TO_M3S)
    plt.close(fig)


def test_plot_without_window_has_no_shading():
    res = _result(p_max=1.0e6)
    assert res['window'] is None
    fig = plot_flow_window(res, 1200, 1300.0, 1.0e6)
    assert len(fig.axes[0].patches) == 0
    plt.close(fig)


def test_plot_invalid_unit():
    with pytest.raises(ValueError):
        plot_flow_window(_result(), 1200, 1300.0, 20.0e6, unit='bpm')
import matplotlib.pyplot as plt
import numpy as np

from src.units import GPM_TO_M3S
from src.trajectory import vertical_section


def plot_well(df, section_azi_deg=45.0, save_path=None):
    """Plan view (kiri) dan section view (kanan)."""
    vs = vertical_section(df, section_azi_deg)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

    ax1.plot(df['E'], df['N'], '-o')
    ax1.set(title='Plan view', xlabel='Easting (m)', ylabel='Northing (m)')
    ax1.set_aspect('equal')
    ax1.grid(True)

    ax2.plot(vs, df['TVD'], '-o')
    ax2.invert_yaxis()
    ax2.set(title='Section view', xlabel='Vertical section (m)', ylabel='TVD (m)')
    ax2.grid(True)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150)


def plot_well_3d(df, save_path=None):
    """Lintasan sumur dalam 3D."""
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(projection='3d')
    ax.plot(df['E'], df['N'], df['TVD'], '-o')
    ax.invert_zaxis()
    ax.set(xlabel='Easting (m)', ylabel='Northing (m)', zlabel='TVD (m)')
    if save_path:
        fig.savefig(save_path, dpi=150)

def plot_flow_window(res, rho, ecd_max, p_max, ft_min=0.5, unit='lps', title=None):
    """Plot ECD, tekanan pompa, dan rasio transport terhadap laju alir.

    res: dict hasil flow_rate_window. rho: densitas statis (kg/m3);
    ecd_max (kg/m3); p_max (Pa); ft_min: batas rasio transport.
    unit: 'lps' (L/s) atau 'gpm'. Jendela layak diarsir hijau. Mengembalikan fig.
    """
    q = np.asarray(res['q'], dtype=float)
    if unit == 'lps':
        conv, xlabel = 1000.0, 'Laju alir (L/s)'
    elif unit == 'gpm':
        conv, xlabel = 1.0 / GPM_TO_M3S, 'Laju alir (gpm)'
    else:
        raise ValueError("unit harus 'lps' atau 'gpm'")
    x = q * conv

    fig, (ax_e, ax_p, ax_f) = plt.subplots(3, 1, figsize=(8, 9), sharex=True)

    ax_e.plot(x, res['ecd'], color='tab:blue', label='ECD dinamis')
    ax_e.axhline(ecd_max, color='tab:red', ls='--', label='ECD maksimum')
    ax_e.axhline(rho, color='gray', ls=':', label='Densitas statis')
    ax_e.set_ylabel('ECD (kg/m3)')

    ax_p.plot(x, np.asarray(res['p_pump']) / 1e6, color='tab:blue',
              label='Tekanan pompa')
    ax_p.axhline(p_max / 1e6, color='tab:red', ls='--', label='Tekanan pompa maksimum')
    ax_p.set_ylabel('Tekanan pompa (MPa)')

    ax_f.plot(x, res['ft'], color='tab:blue', label='Rasio transport (segmen terburuk)')
    ax_f.axhline(ft_min, color='tab:red', ls='--', label='FT minimum')
    ax_f.set_ylabel('Rasio transport FT')
    ax_f.set_xlabel(xlabel)

    if res['window'] is not None:
        lo, hi = res['window']
        for i, ax in enumerate((ax_e, ax_p, ax_f)):
            ax.axvspan(lo * conv, hi * conv, color='tab:green', alpha=0.15,
                       label='Jendela layak' if i == 0 else None)
    for ax in (ax_e, ax_p, ax_f):
        ax.grid(alpha=0.3)
        ax.legend(loc='best', fontsize=8)
    if title:
        ax_e.set_title(title)
    fig.tight_layout()
    return fig

def plot_mse_depth(depth, mse_pa, window=1, unit='mpa', title=None):
    """Plot MSE terhadap kedalaman (sumbu kedalaman ke bawah).

    depth: m; mse_pa: Pa. window: ganjil; >1 menambah garis median bergulir.
    unit: 'mpa' atau 'psi'. Mengembalikan fig.
    """
    from src.mse import rolling_median
    depth = np.asarray(depth, dtype=float)
    mse_pa = np.asarray(mse_pa, dtype=float)
    if depth.shape != mse_pa.shape:
        raise ValueError('depth dan mse_pa harus sama panjang')
    if unit == 'mpa':
        conv, xlabel = 1.0e-6, 'MSE (MPa)'
    elif unit == 'psi':
        conv, xlabel = 1.0 / 6894.757293168, 'MSE (psi)'
    else:
        raise ValueError("unit harus 'mpa' atau 'psi'")
    x = mse_pa * conv

    fig, ax = plt.subplots(figsize=(5, 8))
    ax.plot(x, depth, color='tab:blue', alpha=0.35, lw=1, label='MSE')
    if window > 1:
        ax.plot(rolling_median(x, window), depth, color='tab:blue', lw=2,
                label=f'Median bergulir ({window} titik)')
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Kedalaman (m)')
    ax.grid(alpha=0.3)
    ax.legend(loc='best', fontsize=8)
    if title:
        ax.set_title(title)
    fig.tight_layout()
    return fig
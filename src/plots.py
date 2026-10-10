import matplotlib.pyplot as plt

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
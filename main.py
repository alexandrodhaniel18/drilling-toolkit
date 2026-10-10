from pathlib import Path

import matplotlib.pyplot as plt

from src.plots import plot_well, plot_well_3d
from src.trajectory import closure, read_survey


def main():
    Path('outputs').mkdir(exist_ok=True)
    df = read_survey('data/survey_example.csv')
    print(df.round(2))
    print('Closure (m, deg):', closure(df))

    df.round(3).to_csv('outputs/trajectory_result.csv', index=False)
    plot_well(df, section_azi_deg=45.0, save_path='outputs/well_2d.png')
    plot_well_3d(df, save_path='outputs/well_3d.png')
    plt.show()


if __name__ == '__main__':
    main()
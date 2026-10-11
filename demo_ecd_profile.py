from pathlib import Path

import matplotlib.pyplot as plt

from src.design import design_j_well
from src.trajectory import minimum_curvature
from src.units import kgm3_to_ppg
from src.well_hydraulics import ecd_profile

RHO, Q, PV, YP = 1200.0, 0.03, 0.025, 6.0       # data dummy (Q = 0,03 m3/s ~ 475 gpm)


def main():
    Path('outputs').mkdir(exist_ok=True)
    design = design_j_well(kop=500, dls=3, inc_final=45, azi=45, td=2000)
    survey = minimum_curvature(design['MD'], design['INC'], design['AZI'])

    intervals = [(0.0, 800.0, 0.22, 0.127),       # casing (ID 0,22 m), data dummy
                 (800.0, 2000.0, 0.2159, 0.127)]  # open hole 8,5 in, data dummy

    prof = ecd_profile(survey, RHO, Q, intervals, PV, YP)
    prof['ECD_ppg'] = kgm3_to_ppg(prof['ECD'])
    ringkas = prof.round(2)
    print(ringkas.iloc[::10].to_string(index=False))
    print(ringkas.tail(1).to_string(index=False, header=False))   # baris TD

    fig, ax = plt.subplots(figsize=(5, 6))
    ax.plot(prof['ECD'], prof['TVD'], '-o', label='ECD')
    ax.axvline(RHO, linestyle='--', color='gray', label='MW statis')
    ax.invert_yaxis()
    ax.set(xlabel='ECD (kg/m3)', ylabel='TVD (m)', title='Profil ECD (data dummy)')
    ax.grid(True)
    ax.legend()
    fig.tight_layout()
    fig.savefig('outputs/ecd_profile.png', dpi=150)
    plt.show()


if __name__ == '__main__':
    main()
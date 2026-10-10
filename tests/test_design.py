import numpy as np
import pytest

from src.design import design_j_well
from src.trajectory import minimum_curvature


def test_j_well_matches_analytic():
    kop, dls, inc_f, azi, td = 500.0, 3.0, 45.0, 45.0, 2000.0
    sv = design_j_well(kop, dls, inc_f, azi, td)
    df = minimum_curvature(sv['MD'], sv['INC'], sv['AZI'])

    R = 180 * 30 / (np.pi * dls)
    I = np.radians(inc_f)
    eob = kop + inc_f / dls * 30

    row = df[np.isclose(df['MD'], eob)].iloc[0]          # baris di EOB
    dep_eob = np.hypot(row['N'], row['E'])
    assert np.isclose(row['TVD'], kop + R * np.sin(I))
    assert np.isclose(dep_eob, R * (1 - np.cos(I)))

    hold = td - eob
    end = df.iloc[-1]
    assert np.isclose(end['TVD'], row['TVD'] + hold * np.cos(I))
    assert np.isclose(np.hypot(end['N'], end['E']), dep_eob + hold * np.sin(I))


def test_invalid_design():
    with pytest.raises(ValueError):
        design_j_well(500, 3, 45, 45, td=900)    # EOB di 950 m melewati TD
    with pytest.raises(ValueError):
        design_j_well(500, 0, 45, 45, 2000)      # dls harus > 0
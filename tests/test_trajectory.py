import numpy as np
import pytest

from src.trajectory import closure, minimum_curvature, read_survey, vertical_section

def test_vertical_well():
    df = minimum_curvature([0, 500, 1000], [0, 0, 0], [0, 0, 0])
    assert np.allclose(df['TVD'], df['MD'])
    assert np.allclose(df['N'], 0) and np.allclose(df['E'], 0)


def test_straight_hold():
    # INC 30 deg, AZI 90 deg konstan: ke timur, TVD = MD cos30, E = MD sin30
    df = minimum_curvature([0, 1000], [30, 30], [90, 90])
    assert np.isclose(df['TVD'].iloc[-1], 1000 * np.cos(np.radians(30)))
    assert np.isclose(df['E'].iloc[-1], 1000 * np.sin(np.radians(30)))
    assert np.isclose(df['N'].iloc[-1], 0, atol=1e-9)


def test_quarter_circle_arc():
    R = 1000.0
    df = minimum_curvature([0, np.pi * R / 2], [0, 90], [0, 0])
    assert np.isclose(df['TVD'].iloc[-1], R)
    assert np.isclose(df['N'].iloc[-1], R)


def test_invariants():
    md = [0, 500, 600, 700, 800, 900, 1000, 1500, 2000]
    inc = [0, 0, 10, 20, 30, 40, 45, 45, 45]
    azi = [0, 45, 45, 45, 45, 45, 45, 45, 45]
    df = minimum_curvature(md, inc, azi)
    assert (df['MD'] >= df['TVD'] - 1e-9).all()           # MD tidak pernah < TVD
    dist = np.sqrt(df['N']**2 + df['E']**2 + df['TVD']**2)
    assert (dist <= df['MD'] + 1e-9).all()                # tali busur <= panjang busur
    assert (df['DLS'] >= 0).all()

def _j_type():
    md = [0, 500, 600, 700, 800, 900, 1000, 1500, 2000]
    inc = [0, 0, 10, 20, 30, 40, 45, 45, 45]
    azi = [0, 45, 45, 45, 45, 45, 45, 45, 45]
    return minimum_curvature(md, inc, azi)


def test_closure_j_type():
    dist, azi = closure(_j_type())
    assert np.isclose(azi, 45.0)                         # semua azimuth 45 deg
    assert np.isclose(dist, 642.54 * np.sqrt(2), atol=0.05)


def test_vertical_section_directions():
    df = _j_type()
    along = vertical_section(df, 45.0)                   # searah sumur
    across = vertical_section(df, 135.0)                 # tegak lurus sumur
    assert np.isclose(along.iloc[-1], closure(df)[0])
    assert np.allclose(across, 0, atol=1e-9)


def test_read_survey_csv(tmp_path):
    f = tmp_path / 'survey.csv'
    f.write_text('MD,INC,AZI\n0,0,0\n500,0,45\n600,10,45\n')
    df = read_survey(f)
    assert np.isclose(df['TVD'].iloc[-1], 599.49, atol=0.01)


def test_read_survey_missing_column(tmp_path):
    f = tmp_path / 'bad.csv'
    f.write_text('MD,INC\n0,0\n100,5\n')
    with pytest.raises(ValueError):
        read_survey(f)
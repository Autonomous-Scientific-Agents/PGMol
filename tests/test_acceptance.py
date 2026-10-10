from pathlib import Path

import numpy as np
import pytest

from pgmol import find_point_group, point_group, read_xyz

FIXTURES = sorted((Path(__file__).parents[1] / "core_xyz").glob("*.xyz"))


@pytest.mark.parametrize("fixture", FIXTURES, ids=lambda path: path.stem)
def test_every_public_fixture(fixture: Path) -> None:
    expected = fixture.stem.rsplit("_", 1)[1]
    assert find_point_group(fixture) == expected


@pytest.mark.parametrize("fixture", FIXTURES)
def test_invariant_to_rigid_motion_and_order(fixture: Path) -> None:
    symbols, coordinates = read_xyz(fixture)
    expected = fixture.stem.rsplit("_", 1)[1]
    rng = np.random.default_rng(1729)
    matrix, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    matrix[:, 0] *= -1  # include an orientation-reversing rigid frame change
    order = rng.permutation(len(symbols))
    moved = coordinates @ matrix + np.array([8.2, -4.1, 2.7])
    assert point_group([symbols[i] for i in order], moved[order]) == expected

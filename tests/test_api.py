import numpy as np
import pytest

from pgmol import point_group


def test_coordinate_input_is_not_mutated() -> None:
    coordinates = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    original = coordinates.copy()
    assert point_group(["H", "Cl"], coordinates) == "Cinfv"
    np.testing.assert_array_equal(coordinates, original)


def test_species_must_match() -> None:
    coordinates = [[0, 0, -1], [0, 0, 1]]
    assert point_group(["H", "H"], coordinates) == "Dinfh"
    assert point_group(["H", "F"], coordinates) == "Cinfv"


def test_tolerance_controls_near_linearity() -> None:
    coordinates = [[0, 0, -1], [0.06, 0, 0], [0, 0, 1]]
    assert point_group(["H", "O", "H"], coordinates, tolerance=0.1) == "Dinfh"
    assert point_group(["H", "O", "H"], coordinates, tolerance=0.01) == "C2v"


def test_pure_s6_group_is_not_reduced_to_inversion() -> None:
    coordinates: list[list[float]] = []
    symbols: list[str] = []
    for radius, height, twist, symbol in (
        (1.0, 1.0, 0.0, "H"),
        (1.5, 0.5, np.radians(15), "F"),
    ):
        for sign, offset in ((1, 0.0), (-1, np.pi)):
            for index in range(3):
                angle = index * 2 * np.pi / 3 + offset + twist
                coordinates.append(
                    [radius * np.cos(angle), radius * np.sin(angle), sign * height]
                )
                symbols.append(symbol)
    assert point_group(symbols, coordinates, tolerance=0.001) == "S6"


@pytest.mark.parametrize("tolerance", [0, -1, float("nan"), float("inf")])
def test_bad_tolerance(tolerance: float) -> None:
    with pytest.raises(ValueError, match="tolerance"):
        point_group(["He"], [[0, 0, 0]], tolerance)


@pytest.mark.parametrize(
    ("symbols", "coordinates", "message"),
    [
        ([], [], "at least one"),
        (["H"], [0, 0, 0], "shape"),
        (["H"], [[0, 0, np.inf]], "finite"),
        ([""], [[0, 0, 0]], "symbols"),
    ],
)
def test_bad_coordinate_input(symbols, coordinates, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        point_group(symbols, coordinates)

import numpy as np

from pgmol.symmetry import (
    _discover_operations,
    _matching_permutation,
    _perfect_matching,
)


def test_bipartite_matching_finds_non_greedy_assignment() -> None:
    assert _perfect_matching([[0, 1], [0]]) == [1, 0]
    assert _perfect_matching([[0], [0]]) is None


def test_matching_permutation_handles_coincident_atoms() -> None:
    coordinates = np.zeros((2, 3))
    permutation = _matching_permutation(
        ["H", "H"], coordinates, np.eye(3), tolerance=0.1
    )
    assert permutation is not None
    assert sorted(permutation) == [0, 1]


def test_operation_discovery_preserves_exact_identity() -> None:
    # Close same-species points let a small nonidentity rotation induce the
    # identity permutation within tolerance. The exact identity must win.
    angle = 0.01
    coordinates = np.array(
        [[np.cos(value), np.sin(value), 0.2] for value in (0.0, angle, 2 * angle)]
        + [[0.0, 0.0, -0.6]]
    )
    centered = coordinates - np.mean(coordinates, axis=0)
    operations = _discover_operations(["H"] * 4, centered, tolerance=0.1)
    assert any(np.array_equal(operation, np.eye(3)) for operation in operations)

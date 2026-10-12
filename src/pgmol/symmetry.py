"""Geometric symmetry operation discovery and Schoenflies classification."""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np
import numpy.typing as npt

Array = npt.NDArray[np.float64]
_IDENTITY = np.eye(3)


def _perfect_matching(adjacency: list[list[int]]) -> list[int] | None:
    """Return one row-to-column perfect matching, if one exists."""
    assigned: dict[int, int] = {}

    def augment(row: int, seen: set[int]) -> bool:
        for column in adjacency[row]:
            if column in seen:
                continue
            seen.add(column)
            if column not in assigned or augment(assigned[column], seen):
                assigned[column] = row
                return True
        return False

    if not all(augment(row, set()) for row in range(len(adjacency))):
        return None
    result = [0] * len(adjacency)
    for column, row in assigned.items():
        result[row] = column
    return result


def _matching_permutation(
    symbols: Sequence[str], coordinates: Array, operation: Array, tolerance: float
) -> tuple[int, ...] | None:
    """Find a species-preserving atom permutation induced by an operation."""
    transformed = coordinates @ operation.T
    limit2 = tolerance * tolerance
    permutation = [-1] * len(symbols)
    for species in sorted(set(symbols)):
        indices = np.asarray(
            [i for i, name in enumerate(symbols) if name == species], dtype=int
        )
        differences = transformed[indices, None, :] - coordinates[indices, :]
        distances2 = np.einsum("ijk,ijk->ij", differences, differences)

        # Symmetry-equivalent atoms are normally farther apart than tolerance.
        # This vectorized path avoids graph matching in that common case.
        nearest = np.argmin(distances2, axis=1)
        if np.all(distances2[np.arange(len(indices)), nearest] <= limit2) and len(
            np.unique(nearest)
        ) == len(indices):
            local_match = nearest.tolist()
        else:
            adjacency = [np.flatnonzero(row <= limit2).tolist() for row in distances2]
            if any(not choices for choices in adjacency):
                return None
            local_match = _perfect_matching(adjacency)
            if local_match is None:
                return None
        for row, column in enumerate(local_match):
            permutation[int(indices[row])] = int(indices[column])
    return tuple(permutation)


def _has_matching(
    symbols: Sequence[str], coordinates: Array, operation: Array, tolerance: float
) -> bool:
    """Test for a species-preserving perfect matching after an operation."""
    return _matching_permutation(symbols, coordinates, operation, tolerance) is not None


def _frame(first: Array, second: Array) -> Array:
    one = first / np.linalg.norm(first)
    residual = second - np.dot(second, one) * one
    two = residual / np.linalg.norm(residual)
    three = np.cross(one, two)
    return np.column_stack((one, two, three))


def _discover_operations(
    symbols: Sequence[str], coordinates: Array, tolerance: float
) -> list[Array]:
    """Construct all operations from the possible images of a non-collinear pair."""
    best: tuple[int, int] | None = None
    best_cross = 0.0
    for i, first in enumerate(coordinates):
        for j, second in enumerate(coordinates):
            size = float(np.linalg.norm(np.cross(first, second)))
            if size > best_cross:
                best = (i, j)
                best_cross = size
    if best is None or best_cross <= 1e-10:
        return [_IDENTITY]

    first_index, second_index = best
    source = _frame(coordinates[first_index], coordinates[second_index])
    operations: dict[tuple[int, ...], tuple[float, Array]] = {}
    radial_slack = 2.0 * tolerance
    for i, first in enumerate(coordinates):
        if symbols[i] != symbols[first_index]:
            continue
        if (
            abs(np.linalg.norm(first) - np.linalg.norm(coordinates[first_index]))
            > radial_slack
        ):
            continue
        for j, second in enumerate(coordinates):
            if symbols[j] != symbols[second_index]:
                continue
            if (
                abs(np.linalg.norm(second) - np.linalg.norm(coordinates[second_index]))
                > radial_slack
            ):
                continue
            if np.linalg.norm(np.cross(first, second)) <= 1e-10:
                continue
            target = _frame(first, second)
            for handedness in (1.0, -1.0):
                operation = target @ np.diag((1.0, 1.0, handedness)) @ source.T
                permutation = _matching_permutation(
                    symbols, coordinates, operation, tolerance
                )
                if permutation is not None:
                    key = (*permutation, 1 if handedness > 0.0 else -1)
                    transformed = coordinates @ operation.T
                    displacement = float(
                        np.max(
                            np.linalg.norm(
                                transformed - coordinates[np.asarray(permutation)],
                                axis=1,
                            )
                        )
                    )
                    previous = operations.get(key)
                    if previous is None or displacement < previous[0]:
                        operations[key] = (displacement, operation)
    identity = (*range(len(symbols)), 1)
    operations[identity] = (0.0, _IDENTITY)
    return [operation for _, operation in operations.values()]


def _rotation_order(operation: Array) -> int:
    cosine = float(np.clip((np.trace(operation) - 1.0) / 2.0, -1.0, 1.0))
    angle = math.acos(cosine)
    if angle < 1e-7:
        return 1
    estimate = max(2, int(round(2.0 * math.pi / angle)))
    for order in range(max(2, estimate - 1), estimate + 2):
        if np.allclose(np.linalg.matrix_power(operation, order), _IDENTITY, atol=1e-5):
            return order
    return estimate


def _axis(operation: Array, eigenvalue: float = 1.0) -> Array:
    values, vectors = np.linalg.eig(operation)
    index = int(np.argmin(np.abs(values - eigenvalue)))
    vector = np.real(vectors[:, index])
    return vector / np.linalg.norm(vector)


def _parallel(first: Array, second: Array) -> bool:
    return abs(float(np.dot(first, second))) > 1.0 - 1e-5


def _classify_operations(operations: Sequence[Array]) -> str:
    proper = [op for op in operations if np.linalg.det(op) > 0.0]
    improper = [op for op in operations if np.linalg.det(op) < 0.0]
    rotations = [(op, _rotation_order(op)) for op in proper]
    maximum = max(order for _, order in rotations)

    # The proper rotational subgroups uniquely identify polyhedral families.
    if len(proper) >= 60 and maximum == 5:
        return "Ih" if len(improper) else "I"
    if len(proper) >= 24 and maximum == 4:
        return "Oh" if len(improper) else "O"
    if len(proper) >= 12 and maximum == 3:
        if not improper:
            return "T"
        inversion = any(np.allclose(op, -_IDENTITY, atol=1e-5) for op in improper)
        return "Th" if inversion else "Td"

    principal_operation, n = max(rotations, key=lambda item: item[1])
    principal = _axis(principal_operation)
    perpendicular_twos = [
        op
        for op, order in rotations
        if order == 2 and not _parallel(_axis(op), principal)
    ]
    reflections = [op for op in improper if abs(np.trace(op) - 1.0) < 1e-5]
    inversion = any(np.allclose(op, -_IDENTITY, atol=1e-5) for op in improper)
    horizontal = any(_parallel(_axis(op, -1.0), principal) for op in reflections)
    vertical = any(not _parallel(_axis(op, -1.0), principal) for op in reflections)

    if n == 1:
        if inversion:
            return "Ci"
        if reflections:
            return "Cs"
        return "C1"
    if perpendicular_twos:
        if horizontal:
            return f"D{n}h"
        if improper:
            return f"D{n}d"
        return f"D{n}"
    if vertical:
        return f"C{n}v"
    if horizontal:
        return f"C{n}h"
    if improper:
        # A pure S_2n group has proper subgroup C_n. Some such groups also
        # contain inversion (S2), so the smallest improper order is not the
        # canonical group label.
        return f"S{2 * n}"
    return f"C{n}"


def classify(symbols: Sequence[str], coordinates: Array, tolerance: float) -> str:
    """Center coordinates, discover symmetry operations, and classify them."""
    if len(symbols) == 1:
        return "Kh"
    centered = coordinates - np.mean(coordinates, axis=0)
    _, singular, right = np.linalg.svd(centered, full_matrices=False)
    axis = right[0]
    perpendicular = centered - np.outer(centered @ axis, axis)
    if float(np.max(np.linalg.norm(perpendicular, axis=1))) <= tolerance:
        inversion = _has_matching(symbols, centered, -_IDENTITY, tolerance)
        return "Dinfh" if inversion else "Cinfv"
    return _classify_operations(_discover_operations(symbols, centered, tolerance))

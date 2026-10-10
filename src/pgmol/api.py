"""Small public entry points for point-group determination."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np
import numpy.typing as npt

from .symmetry import classify
from .xyz import read_xyz


def point_group(
    symbols: Sequence[str],
    coordinates: npt.ArrayLike,
    tolerance: float = 0.1,
) -> str:
    """Determine a Schoenflies point group from Cartesian coordinates.

    ``tolerance`` is the maximum displacement in ångströms allowed between
    every transformed atom and a one-to-one matching atom of the same element.
    Coordinates are not modified. The result is invariant to atom order and
    rigid translation, rotation, or reflection.
    """
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be a positive finite number")
    names = tuple(symbols)
    xyz = np.asarray(coordinates, dtype=float)
    if not names:
        raise ValueError("at least one atom is required")
    if xyz.shape != (len(names), 3):
        raise ValueError("coordinates must have shape (number of symbols, 3)")
    if not np.all(np.isfinite(xyz)):
        raise ValueError("coordinates must be finite")
    if any(not isinstance(name, str) or not name for name in names):
        raise ValueError("symbols must be nonempty strings")
    return classify(names, xyz.copy(), float(tolerance))


def find_point_group(path: str | Path, tolerance: float = 0.1) -> str:
    """Read an XYZ file and return its canonical Schoenflies point-group label.

    ``tolerance`` is a maximum species-aware atom matching distance in
    ångströms; see :func:`point_group` for full semantics.
    """
    symbols, coordinates = read_xyz(path)
    return point_group(symbols, coordinates, tolerance)

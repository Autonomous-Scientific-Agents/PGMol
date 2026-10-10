"""Strict parsing of ordinary XYZ coordinate files."""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import numpy.typing as npt

_ELEMENTS = frozenset(
    [
        "H",
        "He",
        "Li",
        "Be",
        "B",
        "C",
        "N",
        "O",
        "F",
        "Ne",
        "Na",
        "Mg",
        "Al",
        "Si",
        "P",
        "S",
        "Cl",
        "Ar",
        "K",
        "Ca",
        "Sc",
        "Ti",
        "V",
        "Cr",
        "Mn",
        "Fe",
        "Co",
        "Ni",
        "Cu",
        "Zn",
        "Ga",
        "Ge",
        "As",
        "Se",
        "Br",
        "Kr",
        "Rb",
        "Sr",
        "Y",
        "Zr",
        "Nb",
        "Mo",
        "Tc",
        "Ru",
        "Rh",
        "Pd",
        "Ag",
        "Cd",
        "In",
        "Sn",
        "Sb",
        "Te",
        "I",
        "Xe",
        "Cs",
        "Ba",
        "La",
        "Ce",
        "Pr",
        "Nd",
        "Pm",
        "Sm",
        "Eu",
        "Gd",
        "Tb",
        "Dy",
        "Ho",
        "Er",
        "Tm",
        "Yb",
        "Lu",
        "Hf",
        "Ta",
        "W",
        "Re",
        "Os",
        "Ir",
        "Pt",
        "Au",
        "Hg",
        "Tl",
        "Pb",
        "Bi",
        "Po",
        "At",
        "Rn",
        "Fr",
        "Ra",
        "Ac",
        "Th",
        "Pa",
        "U",
        "Np",
        "Pu",
        "Am",
        "Cm",
        "Bk",
        "Cf",
        "Es",
        "Fm",
        "Md",
        "No",
        "Lr",
        "Rf",
        "Db",
        "Sg",
        "Bh",
        "Hs",
        "Mt",
        "Ds",
        "Rg",
        "Cn",
        "Nh",
        "Fl",
        "Mc",
        "Lv",
        "Ts",
        "Og",
    ]
)
_ELEMENT = re.compile(r"^[A-Z][a-z]?$")


class XYZFormatError(ValueError):
    """Raised when an XYZ input is malformed or inconsistent."""


def read_xyz(path: str | Path) -> tuple[tuple[str, ...], npt.NDArray[np.float64]]:
    """Read a standard XYZ file, returning element symbols and ångström coordinates.

    The first line is the atom count and the second is an arbitrary comment. Exactly
    that many four-field atom records must follow. Blank trailing lines are accepted;
    extended-XYZ properties and extra atom-line fields are intentionally rejected.
    """
    source = Path(path)
    try:
        text = source.read_text(encoding="utf-8")
    except OSError as exc:
        raise XYZFormatError(f"cannot read XYZ file {source}: {exc}") from exc
    lines = text.splitlines()
    if len(lines) < 2:
        raise XYZFormatError("XYZ input must contain an atom count and comment line")
    try:
        count = int(lines[0].strip())
    except ValueError as exc:
        raise XYZFormatError("first line must be an integer atom count") from exc
    if count < 1:
        raise XYZFormatError("atom count must be positive")
    records = lines[2:]
    while records and not records[-1].strip():
        records.pop()
    if len(records) != count:
        raise XYZFormatError(f"expected {count} atom records, found {len(records)}")

    symbols: list[str] = []
    coordinates: list[list[float]] = []
    for number, line in enumerate(records, start=3):
        fields = line.split()
        if len(fields) != 4:
            raise XYZFormatError(f"line {number}: expected symbol and 3 coordinates")
        symbol = fields[0]
        if not _ELEMENT.fullmatch(symbol) or symbol not in _ELEMENTS:
            raise XYZFormatError(f"line {number}: invalid element symbol {symbol!r}")
        try:
            xyz = [float(value) for value in fields[1:]]
        except ValueError as exc:
            raise XYZFormatError(f"line {number}: coordinates must be numbers") from exc
        if not all(math.isfinite(value) for value in xyz):
            raise XYZFormatError(f"line {number}: coordinates must be finite")
        symbols.append(symbol)
        coordinates.append(xyz)
    return tuple(symbols), np.asarray(coordinates, dtype=float)

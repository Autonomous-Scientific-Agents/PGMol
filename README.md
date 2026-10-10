# PGMol

PGMol determines a molecule's Schoenflies point group from element identities
and Cartesian coordinates. It matches geometric symmetry operations directly;
filenames, formulas, and bonding assumptions are not used.

## Install

PGMol requires Python 3.10+ and NumPy.

```bash
python -m pip install .
# Development tools:
python -m pip install -e '.[dev]'
```

## Use

```python
from pgmol import find_point_group, point_group

print(find_point_group("water.xyz"))
print(point_group(["O", "H", "H"], [[0, 0, 0], [0.757, 0, 0.587], [-0.757, 0, 0.587]]))
```

```bash
pgmol water.xyz
pgmol --tolerance 0.05 water.xyz
# Also supported:
python -m pgmol water.xyz
```

Both calls return a canonical ASCII label such as `C2v`, `D3h`, `Td`, `Oh`,
`Ih`, `Cinfv`, `Dinfh`, or `Kh`. Finite groups are determined from detected
proper rotations, reflections, inversion, and improper operations.

## Tolerance

The default tolerance is 0.1 ångström. A candidate operation is accepted only
when **every** transformed atom has a one-to-one match of the same element no
farther than this distance. The same maximum-distance rule detects linearity.
This deterministic, absolute rule deliberately does not round coordinates or
scale tolerance with molecular size.

## XYZ format

A standard XYZ file contains an integer atom count, one required comment line,
and exactly that many records, each with an element symbol and three finite
Cartesian coordinates in ångströms:

```xyz
3
water
O  0.000  0.000  0.000
H  0.757  0.000  0.587
H -0.757  0.000  0.587
```

Blank lines at the end are accepted. Extra fields, extended-XYZ property
columns, blank atom records, fractional coordinates, and multiple concatenated
frames are not supported and produce a useful error.

## Testing

```bash
ruff check .
ruff format --check .
pytest --cov=pgmol --cov-branch --cov-report=term-missing
```

The test suite checks every `core_xyz` fixture through the public API as well
as rigid motions, reflected coordinate frames, atom permutations, malformed
input, tolerance behavior, and CLI paths.

## Limitations and performance

Atoms are treated as unlabeled points within each element; isotopes and bonding
are not represented. The exhaustive operation search uses possible images of
a non-collinear atom pair and species-aware bipartite matching. Its practical
cost is roughly cubic to quartic in the largest same-element population, so it
is intended for small and medium molecules rather than macromolecules. Very
loose tolerances can merge genuinely distinct geometries, as any
coordinate-only classifier must.

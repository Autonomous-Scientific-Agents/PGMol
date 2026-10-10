"""Command-line interface for PGMol."""

from __future__ import annotations

import argparse

from .api import find_point_group


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="determine a molecule's point group")
    parser.add_argument("xyz", help="path to an XYZ file")
    parser.add_argument(
        "-t",
        "--tolerance",
        type=float,
        default=0.1,
        help="maximum atom matching displacement in angstrom (default: 0.1)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        result = find_point_group(args.xyz, args.tolerance)
    except (OSError, ValueError) as exc:
        parser.exit(2, f"pgmol: error: {exc}\n")
    print(result)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())

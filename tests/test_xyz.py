from pathlib import Path

import pytest

from pgmol import XYZFormatError, read_xyz


def write(tmp_path: Path, text: str) -> Path:
    path = tmp_path / "input.xyz"
    path.write_text(text, encoding="utf-8")
    return path


def test_read_xyz_and_trailing_blanks(tmp_path: Path) -> None:
    symbols, xyz = read_xyz(write(tmp_path, "1\ncomment\nHe 1 2 3\n\n"))
    assert symbols == ("He",)
    assert xyz.tolist() == [[1.0, 2.0, 3.0]]


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("", "count and comment"),
        ("x\ncomment\n", "integer atom count"),
        ("0\ncomment\n", "positive"),
        ("2\ncomment\nH 0 0 0\n", "expected 2"),
        ("1\ncomment\nH 0 0\n", "expected symbol"),
        ("1\ncomment\nh 0 0 0\n", "invalid element"),
        ("1\ncomment\nXx 0 0 0\n", "invalid element"),
        ("1\ncomment\nH x 0 0\n", "must be numbers"),
        ("1\ncomment\nH nan 0 0\n", "must be finite"),
        ("1\ncomment\nH 0 0 0 extra\n", "expected symbol"),
    ],
)
def test_rejects_malformed_xyz(tmp_path: Path, text: str, message: str) -> None:
    with pytest.raises(XYZFormatError, match=message):
        read_xyz(write(tmp_path, text))


def test_missing_file_is_format_error(tmp_path: Path) -> None:
    with pytest.raises(XYZFormatError, match="cannot read"):
        read_xyz(tmp_path / "missing.xyz")

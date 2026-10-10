from pathlib import Path

import pytest

from pgmol.cli import main

ROOT = Path(__file__).parents[1]


def test_cli_success(capsys: pytest.CaptureFixture[str]) -> None:
    assert main([str(ROOT / "core_xyz/H2O_C2v.xyz")]) == 0
    assert capsys.readouterr().out == "C2v\n"


def test_cli_error(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    bad = tmp_path / "bad.xyz"
    bad.write_text("bad", encoding="utf-8")
    with pytest.raises(SystemExit) as error:
        main([str(bad)])
    assert error.value.code != 0
    assert "error" in capsys.readouterr().err

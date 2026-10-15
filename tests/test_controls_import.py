"""controls.py must import cleanly both ways it is used:
as `controls` (python/ on sys.path, the Research notebooks) and as
`python.controls` (repo root on sys.path, evoca_test.ipynb).
Regression: the eg/dyn activity probes did `from python.evoca_py import ...`
inside run_with_controls, which raised ModuleNotFoundError from Research/."""
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
pytest.importorskip("ipywidgets")


def _check(path_entry, modname, cwd):
    code = (f"import sys; sys.path.insert(0, {str(path_entry)!r}); "
            f"import importlib; m = importlib.import_module({modname!r}); "
            "assert m.EGENE_KEY_COUNT > 0 and m.DYN_ACT_COUNT > 0")
    r = subprocess.run([sys.executable, "-c", code], cwd=cwd,
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_import_as_toplevel_from_elsewhere(tmp_path):
    _check(ROOT / "python", "controls", tmp_path)


def test_import_as_package_from_root():
    _check(ROOT, "python.controls", ROOT)

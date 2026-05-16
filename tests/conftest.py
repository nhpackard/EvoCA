"""Shared pytest fixtures for the EvoCA test suite.

Run from repo root:  python3 -m pytest -q
The C library is rebuilt once per session if missing or stale.
"""
import os
import subprocess
import sys

import pytest

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(_ROOT, 'python'))


@pytest.fixture(scope='session', autouse=True)
def _build_lib():
    src = os.path.join(_ROOT, 'C', 'evoca.c')
    lib = os.path.join(_ROOT, 'C', 'libevoca.dylib')
    if (not os.path.exists(lib)
            or os.path.getmtime(lib) < os.path.getmtime(src)):
        subprocess.run(
            ['gcc', '-O2', '-Wall', '-fPIC', '-dynamiclib',
             '-o', lib, src], cwd=_ROOT, check=True)
    return lib


@pytest.fixture
def make_sim():
    """Factory returning a fresh, freed-on-teardown EvoCA instance."""
    from evoca_py import EvoCA
    created = []

    def _factory(N, **kw):
        s = EvoCA()
        s.init(N, **kw)
        created.append(s)
        return s

    yield _factory
    for s in created:
        try:
            s.free()
        except Exception:
            pass

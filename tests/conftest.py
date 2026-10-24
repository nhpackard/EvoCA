"""Shared pytest fixtures for the EvoCA test suite.

Run from repo root:  python3 -m pytest -q
The C library is rebuilt once per session if missing or stale.
"""
import os
import sys

import pytest

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(_ROOT, 'python'))


@pytest.fixture(scope='session', autouse=True)
def _build_lib():
    import evoca_build
    return evoca_build.build()


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

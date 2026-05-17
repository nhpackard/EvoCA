"""D2 regression: tax_lut and tax_ring must round-trip through .evoca
recipes (both metaparams_init and metaparams_final), so the genelife
ring-ladder A/B configs reproduce from a saved recipe.

Guards the exact bug D2 fixed: these knobs were absent from
params()/metaparams_final, so a saved recipe silently lost them.
"""
import os

import pytest

from evoca_py import EvoCA, import_run


@pytest.mark.parametrize("recipe_kind", ["init", "final"])
def test_tax_lut_ring_roundtrip(recipe_kind):
    tax_lut, tax_ring = 7e-4, 3e-3            # both non-default
    s = EvoCA()
    s.init(32, tax_lut=tax_lut, tax_ring=tax_ring)
    s.state(lut='gol', egenome='uniform')      # populate _state_params
    path = s.export_recipe(f"_pytest_d2_{recipe_kind}")
    s.free()
    try:
        assert os.path.exists(path)
        sim, _kw = import_run(path, recipe=recipe_kind)
        try:
            assert sim.tax_lut == pytest.approx(tax_lut), \
                f"tax_lut lost in metaparams_{recipe_kind}"
            assert sim.tax_ring == pytest.approx(tax_ring), \
                f"tax_ring lost in metaparams_{recipe_kind}"
            # And they reach the C side via init(**mp).
            assert sim._lib.evoca_get_tax_lut() == pytest.approx(tax_lut)
            assert sim._lib.evoca_get_tax_ring() == pytest.approx(tax_ring)
        finally:
            sim.free()
    finally:
        if os.path.exists(path):
            os.remove(path)


def test_params_includes_tax_lut_ring():
    s = EvoCA()
    s.init(16, tax_lut=1e-3, tax_ring=2e-3)
    try:
        p = s.params()
        assert p.get('tax_lut') == pytest.approx(1e-3)
        assert p.get('tax_ring') == pytest.approx(2e-3)
    finally:
        s.free()

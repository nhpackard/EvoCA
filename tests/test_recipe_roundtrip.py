"""D2/D4 regression: the full metaparameter set must round-trip
through .evoca recipes (both metaparams_init and metaparams_final),
so saved configs reproduce. D2 added tax_lut/tax_ring (genelife
ring-ladder A/B); D4 added mu_egenome/p_dup_egene/tax_per_egene
(egene-driven §7 cognition experiments). Guards the exact bug: these
knobs were absent from params()/metaparams_final, so a saved recipe
silently lost them.
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


# D4: the rest of the param set (egene-driven recipes, e.g. the §7
# cognition experiments) must also round-trip through metaparams_final.
_D4 = dict(mu_egenome=4e-3, p_dup_egene=0.5, tax_per_egene=6e-4)


@pytest.mark.parametrize("recipe_kind", ["init", "final"])
def test_full_param_set_roundtrip(recipe_kind):
    s = EvoCA()
    s.init(32, **_D4)
    s.state(lut='gol', egenome='uniform')
    path = s.export_recipe(f"_pytest_d4_{recipe_kind}")
    s.free()
    try:
        sim, _kw = import_run(path, recipe=recipe_kind)
        try:
            for k, v in _D4.items():
                assert getattr(sim, k) == pytest.approx(v), \
                    f"{k} lost in metaparams_{recipe_kind}"
                assert s_params_get(sim, k) == pytest.approx(v)
        finally:
            sim.free()
    finally:
        if os.path.exists(path):
            os.remove(path)


def s_params_get(sim, k):
    return sim.params().get(k)

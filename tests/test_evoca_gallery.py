"""Tests for python/evoca_gallery.py (Gallery / Entry / GalleryRun).

Each capability is entered through the door a notebook uses:
  - lookup: G[i], G['label'], G['cfg246'] (unique suffix), G[[...]], G[a:b], split, where
  - load / show: the overrides and the scan-row load-time corrections take effect
    on the sim that comes back (show: on the sim handed to run_with_controls)
  - run: parallel == serial; a cached result is returned without re-running;
    the cache key follows overrides and source content, not label or title
  - save / open: round trip; an unchanged gallery is not rewritten
"""
import json
import os
import sys

import numpy as np
import pytest

_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(_ROOT, 'python'))

import evoca_gallery as eg
from evoca_gallery import Entry, Gallery

PUSH = '2026-04-27_push_edges'
R102 = '2026-04-27_top_16_cfg102_score0.73.evoca'          # pre-A9 recipe: no food_cap key
R236 = '2026-04-27_top_1_cfg236_score0.90.evoca'           # m_scale 1.2
FILL = '2026-04-26_fill_dead_region_with_evolving_patches.evoca'


@pytest.fixture
def cache(tmp_path, monkeypatch):
    monkeypatch.setattr(eg, 'CACHE_DIR', str(tmp_path / 'cache'))
    return tmp_path / 'cache'


def small():
    return Gallery('t', [Entry(recipe=R102), Entry(scan=PUSH, cfg=246),
                         Entry(recipe=FILL, label='fill', title='mosaic'),
                         Entry(scan='2026-04-27_initial', cfg=49, meta={'layer': 2})])


# ── construction and lookup ──────────────────────────────────────────

def test_entry_needs_exactly_one_source():
    with pytest.raises(ValueError):
        Entry(label='x')
    with pytest.raises(ValueError):
        Entry(recipe=R102, scan=PUSH, cfg=1)
    with pytest.raises(ValueError):
        Entry(scan=PUSH)


def test_default_labels_match_between_recipe_and_scan_row():
    # a recipe exported from a scan row gets the same label as the row itself
    assert Entry(recipe=R102).label == 'push_edges.cfg102'
    assert Entry(scan=PUSH, cfg=102).label == 'push_edges.cfg102'


def test_lookup_by_index_label_and_suffix():
    G = small()
    assert G[1].label == 'push_edges.cfg246'
    assert G['push_edges.cfg246'] is G[1]
    assert G['cfg246'] is G[1]
    assert G['fill'] is G[2]
    assert G[-1].label == 'initial.cfg49'
    with pytest.raises(KeyError):
        G['cfg24']                    # not a '.'-separated suffix
    with pytest.raises(IndexError):
        G[4]


def test_ambiguous_suffix_raises():
    G = Gallery('t', [Entry(scan=PUSH, cfg=49), Entry(scan='2026-04-27_initial', cfg=49)])
    with pytest.raises(KeyError, match='matches'):
        G['cfg49']


def test_duplicate_labels_rejected():
    with pytest.raises(ValueError, match='duplicate'):
        Gallery('t', [Entry(recipe=R102), Entry(scan=PUSH, cfg=102)])


def test_sub_galleries():
    G = small()
    assert G[1:3].labels() == ['push_edges.cfg246', 'fill']
    assert G[['fill', 0]].labels() == ['fill', 'push_edges.cfg102']
    assert G.where(layer=2).labels() == ['initial.cfg49']
    sp = Gallery('p', [Entry(scan=PUSH, cfg=c, meta={'layer': l}) for c, l in ((1, 2), (2, 1), (3, 2))]).split('layer')
    assert list(sp) == [1, 2] and sp[2].labels() == ['push_edges.cfg1', 'push_edges.cfg3']


# ── load / show: the public door applies overrides and corrections ───

def test_load_scan_entry_applies_scan_row_corrections():
    sim, _ = Entry(scan=PUSH, cfg=246).load(N=16)
    try:
        assert sim.food_cap == 1                          # pre-A9 row loads the standard rule
        assert sim.mu_egene == pytest.approx(0.01)        # pre-rename column read as mu_egene
        assert sim.mu_egenome == 0.0
    finally:
        sim.free()


def test_load_recipe_entry_overrides():
    e = Entry(recipe=R236)
    assert e.effective_params()['food_cap'] == 1          # as written (RXD-D4 undecided)
    e0 = Entry(recipe=R236, overrides={'food_cap': 0})
    assert e0.effective_params()['food_cap'] == 0
    sim, _ = e0.load(N=16, food_cap=1, tax=0.5)           # call-time overrides win
    try:
        assert (sim.food_cap, sim.tax) == (1, 0.5)
    finally:
        sim.free()


def test_show_forwards_to_run_with_controls(monkeypatch):
    import controls
    got = {}
    def fake(sim, **kw):
        got.update(kw, N=sim.N, food_cap=sim.food_cap); sim.free()
    monkeypatch.setattr(controls, 'run_with_controls', fake)
    small().show('cfg102', N=32, food_cap=0)
    assert got['N'] == 32 and got['food_cap'] == 0 and got['paused'] is False and got['cell_px'] == 2


def test_table_shows_effective_food_cap():
    t = Gallery('t', [Entry(recipe=R236), Entry(scan=PUSH, cfg=236, label='row')]).table()
    assert list(t.food_cap) == [1, 1]                     # both doors load the standard rule (RXD-D4)


# ── run: parallel, cached ────────────────────────────────────────────

RUN = dict(N=32, n_steps=60, every=20, snaps=(60,), cms=(0, 1))


def test_parallel_equals_serial(cache):
    G = small()
    a = G.run(**RUN, workers=2, cache=False, progress=False)
    b = G.run(**RUN, workers=1, cache=False, progress=False)
    for i in range(len(G)):
        assert a[i].ts.equals(b[i].ts)
        assert np.array_equal(a[i].shots[60][1], b[i].shots[60][1])


def test_cache_hit_does_not_rerun(cache, monkeypatch):
    G = small()[:2]
    first = G.run(**RUN, workers=1, progress=False)
    assert len(os.listdir(cache)) == 2
    monkeypatch.setattr(eg, '_run_job', lambda job: pytest.fail('re-ran a cached entry'))
    again = G.run(**RUN, workers=1, progress=False)
    assert all(first[i].ts.equals(again[i].ts) for i in range(2))
    # a new label or title is the same run
    G2 = Gallery('u', [Entry(recipe=R102, label='other', title='x'), Entry(scan=PUSH, cfg=246)])
    G2.run(**RUN, workers=1, progress=False)


def test_cache_key_follows_overrides_and_spec(cache):
    e = Entry(recipe=R102)
    k = lambda ent, **s: eg._cache_key(ent, {**RUN, 'overrides': {}, **s})
    assert k(e) == k(Entry(recipe=R102, label='z', title='t'))
    assert k(e) != k(Entry(recipe=R102, overrides={'food_cap': 0}))
    assert k(e) != k(e, N=64)
    assert k(Entry(scan=PUSH, cfg=246)) != k(Entry(scan=PUSH, cfg=102))


def test_run_overrides_take_effect(cache):
    G = Gallery('t', [Entry(recipe=R236)])
    a = G.run(**RUN, metrics='alive', workers=1, progress=False, overrides={'food_cap': 0})
    b = G.run(**RUN, metrics='alive', workers=1, progress=False, overrides={'food_cap': 0, 'tax': 0.2})
    c = Gallery('t', [Entry(recipe=R236, overrides={'food_cap': 0})]).run(**RUN, metrics='alive', workers=1, progress=False)
    assert a[0].ts.equals(c[0].ts)
    assert not a[0].ts.equals(b[0].ts)                    # tax 0.2 kills the lattice
    assert len(os.listdir(cache)) == 2                    # a and c share one cache file


# ── save / open ──────────────────────────────────────────────────────

def test_save_open_roundtrip_and_no_churn(tmp_path, monkeypatch):
    monkeypatch.setattr(eg, 'galleries_dir', lambda ws: str(tmp_path / ws))
    G = small(); G.workstream = 'Reaction-diffusion'; G.description = 'd'
    p = G.save()
    saved = json.load(open(p))['saved']
    H = Gallery.open(p)
    assert H.labels() == G.labels() and H['fill'].title == 'mosaic' and H[3].meta == {'layer': 2}
    assert H.to_dict() == G.to_dict()
    G.save()
    assert json.load(open(p))['saved'] == saved           # unchanged: not rewritten
    G.entries[2].title = 'changed'; G.save()
    assert Gallery.open(p)['fill'].title == 'changed'

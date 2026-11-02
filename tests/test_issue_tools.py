"""tools/link_issue_refs.py and tools/build_board.py: issue-ID linking and the
aggregate ISSUES copies. Constructed fixtures, so each rule is exercised."""
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'tools'))
sys.path.insert(0, os.path.join(HERE, '..', 'python'))
import link_issue_refs as L  # noqa: E402
import build_board as B  # noqa: E402

DEFS = {'RES-N1': 'Research/Resource-driving/ISSUES.md',
        'CHEM-N1': 'Research/Chem/ISSUES.md',
        'CHEM-N10': 'Research/Chem/ISSUES.md'}


def test_bare_id_becomes_relative_link():
    out, n = L.link_text('See RES-N1 for the campaign.', 'Docs/x.md', DEFS)
    assert out == 'See [RES-N1](../Research/Resource-driving/ISSUES.md) for the campaign.'
    assert n == 1


def test_link_is_relative_to_the_referring_file():
    out, _ = L.link_text('RES-N1', 'Research/Chem/plans/p.md', DEFS)
    assert out == '[RES-N1](../../Resource-driving/ISSUES.md)'


def test_same_file_references_are_left_alone():
    out, n = L.link_text('- **CHEM-N1**: see also CHEM-N10.', 'Research/Chem/ISSUES.md', DEFS)
    assert n == 0 and '[' not in out


def test_protected_spans_untouched():
    text = ('`RES-N1` and [RES-N1](elsewhere.md) and\n```\nRES-N1\n```\n'
            'https://example.org/RES-N1 end')
    out, n = L.link_text(text, 'Docs/x.md', DEFS)
    assert out == text and n == 0


def test_id_boundaries():
    # CHEM-N10 must not be read as CHEM-N1 followed by "0"; undefined IDs stay bare
    out, _ = L.link_text('CHEM-N10, CHEM-N1, XYZ-Q9, SUB-CHEM-N1', 'Docs/x.md', DEFS)
    assert out == ('[CHEM-N10](../Research/Chem/ISSUES.md), [CHEM-N1](../Research/Chem/ISSUES.md), '
                   'XYZ-Q9, SUB-CHEM-N1')


def test_bold_reference_keeps_bold():
    out, _ = L.link_text('**RES-N1** first', 'Docs/x.md', DEFS)
    assert out == '**[RES-N1](../Research/Resource-driving/ISSUES.md)** first'


def test_idempotent():
    once, _ = L.link_text('RES-N1 and CHEM-N1', 'Docs/x.md', DEFS)
    twice, n = L.link_text(once, 'Docs/x.md', DEFS)
    assert twice == once and n == 0


def _tree(tmp_path, monkeypatch, files):
    for rel, text in files.items():
        p = tmp_path / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)
    monkeypatch.setattr(L, 'ROOT', str(tmp_path))
    monkeypatch.setattr(L.evoca_paths, 'WORKSTREAMS', ['Chem', 'Resource-driving'])
    monkeypatch.setattr(B.evoca_paths, 'RESEARCH', str(tmp_path / 'Research'))


def test_definitions_forms(tmp_path, monkeypatch):
    _tree(tmp_path, monkeypatch, {
        'Research/Chem/ISSUES.md': ('# Chem\n- **CHEM-D5**: x\n- **CHEM-N7 done** (date): y\n'
                                    '- **2026-10-23, CHEM-D4 decided: z**\n- mentions **RES-N1** in passing\n'),
        'Research/Resource-driving/ISSUES.md': '# R\n- **RES-N1**: the campaign\n'})
    d = L.definitions()
    assert d == {'CHEM-D5': 'Research/Chem/ISSUES.md', 'CHEM-N7': 'Research/Chem/ISSUES.md',
                 'CHEM-D4': 'Research/Chem/ISSUES.md', 'RES-N1': 'Research/Resource-driving/ISSUES.md'}


def test_duplicate_definition_is_refused(tmp_path, monkeypatch):
    _tree(tmp_path, monkeypatch, {
        'Research/Chem/ISSUES.md': '# C\n- **RES-N1**: wrong home\n',
        'Research/Resource-driving/ISSUES.md': '# R\n- **RES-N1**: right home\n'})
    with pytest.raises(SystemExit):
        L.definitions()


def test_aggregate_rebases_links_and_links_own_ids(tmp_path, monkeypatch):
    _tree(tmp_path, monkeypatch, {
        'Research/Chem/ISSUES.md': ('# Chem\n## Plans\n- [`plans/p.md`](plans/p.md): a plan\n'
                                    '## Open\n- **CHEM-N1**: x, see CHEM-N10 and RES-N1\n- **CHEM-N10**: y\n'),
        'Research/Resource-driving/ISSUES.md': '# R\n- **RES-N1**: z\n'})
    text = B.render_issues(B.render())
    assert '[`plans/p.md`](Chem/plans/p.md)' in text           # rebased into Research/
    assert '[CHEM-N10](Chem/ISSUES.md)' in text                # own-file ref linked in the copy
    assert '[RES-N1](Resource-driving/ISSUES.md)' in text
    assert '\n## Chem\n' in text and '\n### Open\n' in text      # workstream is a top-level section


def test_id_glued_to_a_word_is_not_linked():
    out, n = L.link_text('RES-N1_v2 and RES-N1-old', 'Docs/x.md', DEFS)
    assert out == 'RES-N1_v2 and RES-N1-old' and n == 0

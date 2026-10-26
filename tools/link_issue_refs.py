"""Turn bare issue-ID references (RES-N1, CHEM-D5, ...) in markdown files into
relative links to the Research/<Workstream>/ISSUES.md that defines them.

    python3 tools/link_issue_refs.py          # rewrite in place, print a summary
    python3 tools/link_issue_refs.py --check  # exit 1 if any file would change

An ID is *defined* by a list item in a workstream's ISSUES.md whose opening
bold span names it (`- **ID**: ...`, `- **ID done** ...`). Left alone: the definition sites, references inside the defining
file itself (a link to the same file goes nowhere new), IDs already inside a
link, inline code or fenced code, and IDs no ISSUES.md defines. Historical
records and generated files are skipped (FROZEN, GENERATED); run
tools/build_board.py afterwards to refresh the generated copies.
Idempotent: a second run changes nothing.
"""
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'python'))
import evoca_paths  # noqa: E402

ROOT = evoca_paths.REPO_ROOT
ID = r'[A-Z]{2,5}-[A-Z]\d+[a-z]?'
ID_RE = re.compile(r'(?<![\w/#-])(' + ID + r')(?![\w-])')
FROZEN = ('Docs/devlog/', 'Docs/Archive/', 'Research/research_board_archival.md')
GENERATED = ('Research/ISSUES.md',)
# spans to leave untouched: fenced code, inline code, existing links [text](target), bare URLs
PROTECT = re.compile(r'```.*?```|`[^`\n]*`|\[[^\]\n]*\]\([^)\n]*\)|https?://\S+', re.S)


def definitions():
    """ID -> repo-relative path of the ISSUES.md that defines it."""
    defs = {}
    for ws in evoca_paths.WORKSTREAMS:
        rel = f'Research/{ws}/ISSUES.md'
        p = os.path.join(ROOT, rel)
        if not os.path.exists(p):
            continue
        # a definition is a list item whose opening bold span names the ID
        # ("- **CHEM-D5**: ...", "- **CHEM-N7 done** ...", "- **2026-10-23, RXD-D4 decided: ...**")
        for m in re.finditer(r'^\s*- \*\*([^*\n]+)\*\*', open(p, encoding='utf-8').read(), re.M):
            for ident in ID_RE.findall(m.group(1)):
                if defs.get(ident, rel) != rel:
                    sys.exit(f'{ident} is defined in both {defs[ident]} and {rel}')
                defs[ident] = rel
    return defs


def markdown_files():
    out = subprocess.run(['git', '-C', ROOT, 'ls-files', '-co', '--exclude-standard', '*.md'],
                         capture_output=True, text=True, check=True).stdout
    for rel in out.splitlines():
        if rel.startswith(FROZEN + GENERATED) or '.ipynb_checkpoints' in rel:
            continue
        if os.path.exists(os.path.join(ROOT, rel)):
            yield rel


def link_text(text, rel, defs):
    """Return (new_text, n_links) for the markdown file at repo path rel."""
    here = os.path.dirname(rel)
    n = 0

    def link_ids(chunk):
        nonlocal n

        def one(m):
            nonlocal n
            ident = m.group(1)
            target = defs.get(ident)
            if target is None or target == rel:
                return ident
            n += 1
            return f'[{ident}]({os.path.relpath(target, here)})'
        return ID_RE.sub(one, chunk)

    out, last = [], 0
    for m in PROTECT.finditer(text):
        out.append(link_ids(text[last:m.start()]))
        out.append(m.group(0))
        last = m.end()
    out.append(link_ids(text[last:]))
    return ''.join(out), n


def main():
    defs = definitions()
    check = '--check' in sys.argv
    total, files = 0, 0
    for rel in markdown_files():
        path = os.path.join(ROOT, rel)
        text = open(path, encoding='utf-8').read()
        new, n = link_text(text, rel, defs)
        if n:
            files += 1
            total += n
            print(f'{n:4d}  {rel}')
            if not check:
                open(path, 'w', encoding='utf-8').write(new)
    print(f'{len(defs)} defined IDs; {total} references linked in {files} files'
          + (' (check only)' if check else ''))
    if check and total:
        sys.exit(1)


if __name__ == '__main__':
    main()

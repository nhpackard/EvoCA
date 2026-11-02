"""Build the one shared library (EvoCA + ChemCA + shared core).

The single Python statement of the build: tests/conftest.py and the
notebooks call build(), and run as a script it builds (if stale; --force
always):

    python3 python/evoca_build.py [--force]

The Makefile and CLAUDE.md carry the same command for the shell:

    gcc -O2 -Wall -fPIC -dynamiclib -o C/libevoca.dylib C/evoca.c C/core.c C/chemca.c
"""
import os
import platform
import subprocess
import sys

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

SOURCES = ('C/evoca.c', 'C/core.c', 'C/chemca.c')
HEADERS = ('C/evoca.h', 'C/core.h', 'C/chemca.h')
CFLAGS = ('-O2', '-Wall', '-fPIC')


class BuildError(subprocess.CalledProcessError):
    """CalledProcessError whose message shows the compiler's stderr."""

    def __str__(self):
        return f"{super().__str__()}\n{self.stderr}"


def lib_path(root=ROOT):
    ext = 'dylib' if platform.system() == 'Darwin' else 'so'
    return os.path.join(root, 'C', f'libevoca.{ext}')


def command(root=ROOT):
    """The gcc argument list (paths relative to root)."""
    shared = '-dynamiclib' if platform.system() == 'Darwin' else '-shared'
    return ['gcc', *CFLAGS, shared, '-o', os.path.relpath(lib_path(root), root), *SOURCES]


def stale(root=ROOT):
    """True if the library is missing or older than any source or header."""
    lib = lib_path(root)
    if not os.path.exists(lib):
        return True
    t = os.path.getmtime(lib)
    return any(os.path.getmtime(os.path.join(root, f)) > t
               for f in SOURCES + HEADERS if os.path.exists(os.path.join(root, f)))


def build(root=ROOT, force=False):
    """Rebuild the library if stale (or force). Returns its path; on failure
    raises CalledProcessError whose message includes gcc's stderr."""
    if force or stale(root):
        r = subprocess.run(command(root), cwd=root, capture_output=True, text=True)
        if r.returncode != 0:
            raise BuildError(r.returncode, r.args, r.stdout, r.stderr)
        if r.stderr.strip():
            print(r.stderr)
    return lib_path(root)


if __name__ == '__main__':
    print(build(force='--force' in sys.argv[1:]))

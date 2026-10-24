"""run_sim's step_hook: entered through run_sim itself. The hook sees every
step (t = 1 .. n_steps) on the live sim, its summary reaches the result, it
may not overwrite a run_sim key, and a passive hook leaves the run
bit-identical."""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'python'))

from evoca_explore import run_sim           # noqa: E402

P = dict(food_inc=0.016, m_scale=3.5, gdiff=0.08, mu_lut=0.03, mu_egene=0.003, tax=0.035,
         n_chem=2, wdiff=0.01, waste_sink=0.019, waste_disp=0.5, waste_tox=0.4)


class Counter:
    def __init__(self, clash=False):
        self.ts, self.started, self.alive, self.clash = [], 0, [], clash

    def start(self, sim):
        self.started += 1

    def step(self, sim, t):
        self.ts.append(t)
        self.alive.append(float(sim.get_alive().mean()))

    def summary(self):
        return {'extinct': 'x'} if self.clash else {'hook_steps': len(self.ts), 'hook_alive_last': self.alive[-1]}


def test_hook_sees_every_step_and_reaches_result():
    h = Counter()
    out = run_sim(P, n_steps=50, sample_every=25, N=32, seed=0, shadow=False, step_hook=h)
    assert h.started == 1 and h.ts == list(range(1, 51))
    assert out['hook_steps'] == 50
    assert out['hook_alive_last'] == pytest.approx(out['alive_density_final'])


def test_hook_may_not_overwrite_run_sim_keys():
    with pytest.raises(ValueError, match='extinct'):
        run_sim(P, n_steps=10, sample_every=5, N=32, seed=0, shadow=False, step_hook=Counter(clash=True))


def test_passive_hook_leaves_run_identical():
    a = run_sim(P, n_steps=60, sample_every=20, N=32, seed=1, shadow=True)
    b = run_sim(P, n_steps=60, sample_every=20, N=32, seed=1, shadow=True, step_hook=Counter())
    b = {k: v for k, v in b.items() if not k.startswith('hook_')}
    assert a.keys() == b.keys()
    for k in a:
        assert a[k] == b[k] or (isinstance(a[k], float) and np.isnan(a[k]) and np.isnan(b[k])), k

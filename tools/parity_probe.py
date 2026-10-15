"""Bit-for-bit parity probe for merge gates.

Usage:
  python3 tools/parity_probe.py run  <repo_root> <out.json> [food_cap]  # build root first (make)
  python3 tools/parity_probe.py cmp  <base.json> <new.json>

`run` records, for a viable pure-evo config (N=128, mutation + egene on):
run_sim scalars for 2 seeds, and SHA1 hashes of v/alive/f/F/egenome/LUT
hashes every 100 steps of a direct 1200-step trajectory. `cmp` reports
every differing common scalar, added/lost keys, and trajectory identity.
An additive branch must give BAD 0 and IDENTICAL; new keys are allowed.
Check the fixture is not extinct before trusting a pass (an extinct
fixture passes trivially).

food_cap (optional, branch `food-cap` onward): pass 0 to run the legacy
uncapped mouthful rule, which must be IDENTICAL to a pre-fix baseline.
The fixture has m_scale=1.5 > 1, so the default (capped, food_cap=1) path
is NOT expected to match a pre-fix baseline. Omit it on pre-fix trees.
"""
import sys, os, json, hashlib


def run(root, out, food_cap=None):
    root, out = os.path.abspath(root), os.path.abspath(out)
    sys.path.insert(0, os.path.join(root, 'python'))
    os.chdir(root)
    import numpy as np
    from evoca_explore import run_sim
    from evoca_py import EvoCA

    P = dict(food_inc=0.01, m_scale=1.5, gdiff=0.08, mu_lut=0.01,
             mu_egene=0.003, mu_egenome=0.03, p_dup_egene=1.0,
             tax=0.03, tax_per_egene=0.0, restricted_mu=False)
    if food_cap is not None:
        P['food_cap'] = int(food_cap)
    res = {}
    for seed in (1, 2):
        m = run_sim(P, n_steps=2000, sample_every=100, N=128, seed=seed, shadow=True)
        res[f'run_sim_{seed}'] = {
            k: (float(v) if isinstance(v, (int, float, np.floating, np.integer, bool, np.bool_))
                else repr(v)) for k, v in m.items()}

    def h(a):
        return hashlib.sha1(np.ascontiguousarray(a).tobytes()).hexdigest()[:16]

    N = 128
    sim = EvoCA()
    sim.init(N, **{k: P[k] for k in ('food_inc', 'm_scale', 'gdiff', 'mu_lut', 'mu_egenome', 'tax',
                                     'food_cap') if k in P})
    sim.set_seed(7)
    sim.state(lut='gol', egenome='uniform', egenome_value=0b000011)
    sim.set_v(np.random.default_rng(7).integers(0, 2, (N, N), dtype=np.uint8))
    sim.set_f_all(0.5); sim.set_F_all(1.0); sim.set_alive_halfplane(0)
    traj = []
    for t in range(1200):
        sim.step()
        if t % 100 == 99:
            traj.append([h(sim.get_v()), h(sim.get_alive()), h(sim.get_f()), h(sim.get_F()),
                         h(sim.get_egenome()), h(sim.get_lut_hashes())])
    res['traj'] = traj
    json.dump(res, open(out, 'w'), indent=1, sort_keys=True)
    print('extinct:', res['run_sim_1']['extinct'], res['run_sim_2']['extinct'],
          '| final alive density:', float(sim.get_alive().mean()))


def cmp(a_path, b_path):
    a, b = json.load(open(a_path)), json.load(open(b_path))
    bad = 0
    for k in ('run_sim_1', 'run_sim_2'):
        common = set(a[k]) & set(b[k])
        for c in sorted(common):
            if a[k][c] != b[k][c]:
                bad += 1
                print('DIFF', k, c, a[k][c], b[k][c])
        print(k, 'common', len(common), 'new keys', sorted(set(b[k]) - set(a[k])),
              'lost', sorted(set(a[k]) - set(b[k])))
    print('traj', 'IDENTICAL' if a['traj'] == b['traj'] else 'DIFFERENT')
    print('BAD', bad)
    return bad == 0 and a['traj'] == b['traj']


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        run(sys.argv[2], sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else None)
    else:
        sys.exit(0 if cmp(sys.argv[2], sys.argv[3]) else 1)

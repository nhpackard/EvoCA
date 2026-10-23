import json, sys, os, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from multiprocessing import Pool
ROOT = os.path.dirname(os.path.abspath(__file__))
while not os.path.exists(os.path.join(ROOT, 'C', 'evoca.h')):
    ROOT = os.path.dirname(ROOT)
sys.path.insert(0, os.path.join(ROOT, 'python'))
from evoca_py import EvoCA
from evoca_explore import init_kwargs
X = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(X, 'libevoca_2way.dylib')   # built from twoway.patch, see notes.md   # build: see notes.md
P0 = json.load(open(os.path.join(ROOT, 'Research', 'Chem', 'Runs', '2026-10-23_chem_fw_spiralpoint.evoca')))['metaparams_init']

def build(p, N, seed, lib=LIB):
    """build_sim's initialization, with the library path chosen."""
    rng = np.random.default_rng(seed)
    s = EvoCA(lib_path=lib); s.init(N, **init_kwargs(p)); s.set_seed(seed)
    s.state(lut='gol', egenome='uniform', egenome_value=0b000011)
    s.set_v(rng.integers(0, 2, (N, N), dtype=np.uint8))
    s.set_f_all(0.5); s.set_F_all(1.0); s.set_alive_halfplane(0)
    return s

def img6(s, N):
    c = s.colorize(np.zeros(N * N, np.int32), 6).view(np.uint32).reshape(N, N)
    return np.dstack([(c >> 16) & 255, (c >> 8) & 255, c & 255]).astype(np.uint8)

def mass_check():
    s = build(dict(P0, waste_sink=0.015, wdiff=0.01), 64, 0)
    m = lambda: s.get_W().astype(np.float64).sum() + s.get_w().astype(np.float64).sum()
    m0 = m(); f0 = s.get_chem_flux(1); worst = 0
    for _ in range(500):
        s.step(); f = s.get_chem_flux(1)
        worst = max(worst, abs((m() - m0) - ((f[0] - f0[0]) - (f[1] - f0[1]))))
    s.free(); return worst

def run(arg):
    ws, wd = arg
    s = build(dict(P0, waste_sink=ws, wdiff=wd), 256, 0)
    al = []; td = dd = 0
    for t in range(3000):
        s.step()
        if t >= 1500:
            td += s.get_tox_deaths_last(); dd += s._lib.evoca_get_deaths_last()
        if t % 10 == 0:
            a = s.get_alive().mean()
            if t >= 1500: al.append(a)
            if a == 0: break
    im = img6(s, 256); Wm = float(s.get_W().mean()); s.free()
    al = np.array(al) if al else np.zeros(1)
    return ws, wd, im, al.mean(), al.std(), td / max(dd, 1), Wm

if __name__ == '__main__':
    print("waste mass identity, worst |error| over 500 steps:", mass_check())
    WS = [0.010, 0.013, 0.015, 0.017, 0.019]; WD = [0.01, 0.08]
    with Pool(8) as pool:
        R = pool.map(run, [(w, d) for d in WD for w in WS])
    fig, axs = plt.subplots(len(WD), len(WS), figsize=(4 * len(WS), 4.4 * len(WD)))
    for k, (ws, wd, im, am, sd, fr, Wm) in enumerate(R):
        a = axs.flat[k]; a.imshow(im); a.axis('off')
        a.set_title(f"two-way: sink {ws} wdiff {wd}\nalive {am:.2f}±{sd:.3f} toxic {100*fr:.0f}% W {Wm:.2f}", fontsize=9)
        print(f"sink {ws} wdiff {wd}: alive {am:.3f} std {sd:.3f} toxic {100*fr:.0f}%  W {Wm:.2f}")
    plt.tight_layout(); plt.savefig(os.path.join(X, 'grid2.png'), dpi=50)

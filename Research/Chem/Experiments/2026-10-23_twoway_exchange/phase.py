import sys, os, importlib.util, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run2 import build, P0, X, ROOT
spec = importlib.util.spec_from_file_location('n5', os.path.join(ROOT, 'Research', 'Chem', 'Scans', '2026-10-23_N5_global_chem', 'scan.py'))
n5 = importlib.util.module_from_spec(spec); spec.loader.exec_module(n5)
CASES = [('one-way', 0, 0.015, 0.01), ('two-way', 1, 0.017, 0.01), ('two-way', 1, 0.015, 0.08)]
N = 512; B = 8   # finer blocks than the scan, to resolve pinwheels

def run(case):
    name, ex, ws, wd = case
    s = build(dict(P0, waste_sink=ws, wdiff=wd, waste_exchange=ex), N, 0)
    for _ in range(3000): s.step()
    out = []; snaps = []; prev = None
    for t in range(3000, 4501):
        if t % 250 in (0, 250 - 60):
            pos, q, n, th = n5.phase_singularities(s.get_F(), s.get_W(), b=B)
            if t % 250 == 0:
                if prev is not None:
                    pp, pq, pth = prev
                    out.append((len(pos), n5.persistent(pos, q, pp, pq, n),
                                n5.persistent(pos, q, (pp + n // 2) % n, pq, n)))
                if len(snaps) < 3: snaps.append((t, th, pos, q))
            else:
                prev = (pos, q, th)
        s.step()
    s.free()
    return case, np.array(out), snaps

if __name__ == '__main__':
    with Pool(3) as pool:
        R = pool.map(run, CASES)
    fig, axs = plt.subplots(len(R), 3, figsize=(15, 5.2 * len(R)))
    for r, ((name, _, ws, wd), o, snaps) in enumerate(R):
        ns, pers, null = o.mean(0)
        print(f"{name} sink {ws} wdiff {wd}: singularities {ns:.0f}, persist(60 steps) {pers:.0f}, chance {null:.0f}, ratio {pers/max(null,1):.2f}")
        for k, (t, th, pos, q) in enumerate(snaps):
            a = axs[r, k]; a.imshow(th, cmap='hsv'); a.axis('off')
            a.scatter(pos[:, 1], pos[:, 0], c=['k' if v > 0 else 'w' for v in q], s=10)
            a.set_title(f'{name} {ws}/{wd} t={t}: phase, {len(pos)} defects', fontsize=10)
    plt.tight_layout(); plt.savefig(os.path.join(X, 'phase.png'), dpi=45)

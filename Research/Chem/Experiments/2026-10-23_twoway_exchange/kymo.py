import sys, os, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run2 import build, img6, P0, X, ROOT
CASES = [('one-way', 0, 0.015, 0.01), ('two-way', 1, 0.017, 0.01),
         ('two-way', 1, 0.015, 0.08), ('two-way', 1, 0.019, 0.01)]
N = 512

def run(case):
    name, ex, ws, wd = case
    s = build(dict(P0, waste_sink=ws, wdiff=wd, waste_exchange=ex), N, 0)
    for _ in range(3000): s.step()
    rows = []; frames = []
    for t in range(600):
        s.step()
        im = img6(s, N)
        rows.append(im[N // 2])            # the middle row, every step
        if t % 40 == 0 and len(frames) < 6: frames.append(im[:256, :256])
    s.free()
    return case, np.array(rows), frames

if __name__ == '__main__':
    with Pool(4) as pool:
        R = pool.map(run, CASES)
    fig, axs = plt.subplots(len(R), 7, figsize=(28, 4.2 * len(R)),
                            gridspec_kw=dict(width_ratios=[2.2] + [1] * 6))
    for r, ((name, _, ws, wd), kymo, frames) in enumerate(R):
        axs[r, 0].imshow(kymo, aspect='auto'); axs[r, 0].set_ylabel('t - 3000')
        axs[r, 0].set_title(f'{name} sink {ws} wdiff {wd}: middle row (x) vs time', fontsize=11)
        for k, f in enumerate(frames):
            axs[r, k + 1].imshow(f); axs[r, k + 1].axis('off'); axs[r, k + 1].set_title(f't={3000 + 40 * k}')
    plt.tight_layout(); plt.savefig(os.path.join(X, 'kymo.png'), dpi=45)

import sys, json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import os
ROOT = os.path.dirname(os.path.abspath(__file__))
while not os.path.exists(os.path.join(ROOT, 'C', 'evoca.h')): ROOT = os.path.dirname(ROOT)
sys.path.insert(0, os.path.join(ROOT, 'python'))
import evoca_waves as ew
from evoca_explore import build_sim
p = json.load(open(os.path.join(ROOT, 'Research', 'Chem', 'Runs', '2026-10-23_chem_fw_spiralpoint.evoca')))['metaparams_init']
p.update(waste_exchange=1, waste_sink=0.017, wdiff=0.01, waste_tox=1.0)
N, B = 1024, 16
s = build_sim(p, N=N, seed=0)
for _ in range(2000): s.step()
Fs = Ws = 0
for _ in range(1000):
    s.step(); Fs = Fs + ew.block_means(s.get_F(), B); Ws = Ws + ew.block_means(s.get_W(), B)
cen = (Fs / 1000, Ws / 1000)
for _ in range(1000): s.step()
snaps = []
for k in range(6):
    th, amp = ew.phase_field(s.get_F(), s.get_W(), B, center=cen)
    pos, sg = ew.cores(th, amp)
    c = s.colorize(np.zeros(N * N, np.int32), 6).view(np.uint32).reshape(N, N)
    snaps.append((4000 + 30 * k, th, pos, sg, np.dstack([(c >> 16) & 255, (c >> 8) & 255, c & 255]).astype(np.uint8)))
    for _ in range(30): s.step()
s.free()
fig, axs = plt.subplots(2, 6, figsize=(36, 12))
for k, (t, th, pos, sg, img) in enumerate(snaps):
    axs[0, k].imshow(th, cmap='hsv'); axs[0, k].scatter(pos[:, 1], pos[:, 0], c=['k' if v > 0 else 'w' for v in sg], s=14)
    axs[0, k].set_title(f't={t}: phase (16x16 blocks), {len(pos)} defects'); axs[0, k].axis('off')
    axs[1, k].imshow(img); axs[1, k].axis('off'); axs[1, k].set_title(f't={t} colormode 6')
plt.tight_layout(); plt.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'phase_N1024.png'), dpi=30)

"""Side-by-side movie: EvoCA spiral recipe (left) and a ChemCA port (right),
N = 512, both rendered the same way: red = alive, green = food F (0..1),
blue = waste W (0..W_MAX).

    python3 s4_movie.py <chem case> [seed] [warmup] [steps] [every]

Writes figs/s4_<case>_s<seed>.mp4 (ffmpeg).
"""
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s3_spirals as S  # noqa: E402  (same constructors and cases)

W_MAX = 0.5


def render(F, W, A):
    img = np.empty(F.shape + (3,), np.uint8)
    img[..., 0] = np.where(A, 255, 0)
    img[..., 1] = (np.clip(F, 0, 1) * 255).astype(np.uint8)
    img[..., 2] = (np.clip(W / W_MAX, 0, 1) * 255).astype(np.uint8)
    return img


def main(case, seed=0, warmup=2000, steps=2000, every=4):
    left, lf = S.make('evoca', seed)
    lsim = left
    # ChemCA and EvoCA live in separate C worlds, so both can run at once
    right, rf = S.make(case, seed)
    N = S.N
    out = os.path.join(HERE, 'figs', f's4_{case}_s{seed}.mp4')
    ff = subprocess.Popen(['ffmpeg', '-loglevel', 'error', '-y', '-f', 'rawvideo', '-pix_fmt', 'rgb24',
                           '-s', f'{2 * N + 8}x{N}', '-r', '30', '-i', '-', '-vcodec', 'libx264',
                           '-pix_fmt', 'yuv420p', '-crf', '20', out], stdin=subprocess.PIPE)
    sep = np.full((N, 8, 3), 255, np.uint8)
    for t in range(1, warmup + steps + 1):
        lsim.step()
        right.step()
        if t > warmup and t % every == 0:
            frame = np.concatenate([render(*lf()), sep, render(*rf())], axis=1)
            ff.stdin.write(frame.tobytes())
    ff.stdin.close()
    ff.wait()
    lsim.free()
    right.free()
    print(out)


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], *(int(x) for x in a[1:]))

"""evoca_waves on constructed fields: a planted rotating spiral pair is found
as two defects of opposite sign that persist while the phase turns; defects
of independent noise frames match only at the chance level; R separates a
synchronized bulk oscillation from a spiral, given the time-mean center."""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'python'))

import evoca_waves as ew                    # noqa: E402

N, B = 256, 8
n = N // B
Y, X = np.mgrid[0:N, 0:N].astype(float)


def _pair(phase0, c1=(80.3, 70.7), c2=(176.3, 186.7), k=0.05):
    """A +1/-1 spiral pair (periodic-consistent: the two windings cancel),
    rotating with phase0; F, W = cycle position on the (F, W) circle."""
    th = (np.arctan2(Y - c1[0], X - c1[1]) - np.arctan2(Y - c2[0], X - c2[1])
          + k * np.hypot(Y - c1[0], X - c1[1]) + phase0)
    return 0.5 + 0.4 * np.cos(th), 0.5 + 0.4 * np.sin(th)


def _smooth_noise(rng, s=4):
    k = np.exp(-((np.minimum(X, N - X)) ** 2 + (np.minimum(Y, N - Y)) ** 2) / (2 * s * s))
    a = np.fft.ifft2(np.fft.fft2(rng.standard_normal((N, N))) * np.fft.fft2(k)).real
    return a / a.std()


def test_spiral_pair_cores_found_and_persist():
    th0, amp0 = ew.phase_field(*_pair(0.0), B)
    p0, s0 = ew.cores(th0, amp0)
    assert sorted(s0.tolist()) == [-1, 1]
    near = lambda p, c: np.max(np.abs(p - np.array(c) // B)) <= 1
    assert any(near(p, (80, 70)) for p in p0) and any(near(p, (176, 186)) for p in p0)
    th1, amp1 = ew.phase_field(*_pair(1.0), B)            # the phase turned; cores stay
    p1, s1 = ew.cores(th1, amp1)
    assert ew.matched(p1, s1, p0, s0, n) == 2
    assert ew.matched(p1, s1, (p0 + n // 2) % n, s0, n) == 0


def test_noise_defects_match_only_at_chance():
    """Independent noise frames (defects come and go): later-frame defects
    match earlier ones about as often as the half-lattice-shifted chance level."""
    rng = np.random.default_rng(0)
    frames = [ew.phase_field(0.5 + 0.3 * _smooth_noise(rng), 0.5 + 0.3 * _smooth_noise(rng), B)
              for _ in range(6)]
    obs = ch = 0
    for (ta, aa), (tb, ab) in zip(frames[:-1], frames[1:]):
        pa, sa = ew.cores(ta, aa); pb, sb = ew.cores(tb, ab)
        obs += ew.matched(pb, sb, pa, sa, n); ch += ew.matched(pb, sb, (pa + n // 2) % n, sa, n)
    assert obs > 20 and 0.7 < obs / ch < 1.3


def test_order_parameter_bulk_vs_spiral():
    rng = np.random.default_rng(1)
    phase = 1.3                                            # every block at one phase
    F = 0.5 + 0.4 * np.cos(phase) + 0.01 * rng.standard_normal((N, N))
    W = 0.5 + 0.4 * np.sin(phase) + 0.01 * rng.standard_normal((N, N))
    th, amp = ew.phase_field(F, W, B, center=(0.5, 0.5))   # the time-mean state
    assert ew.order_parameter(th, amp) > 0.9
    th, amp = ew.phase_field(F, W, B)                      # spatial-mean center: blind to it
    assert np.median(amp) < 0.05
    th, amp = ew.phase_field(*_pair(0.0), B, center=(0.5, 0.5))
    assert ew.order_parameter(th, amp) < 0.3


def test_period_of_first_peak():
    t = np.arange(600)
    per, peak = ew.period_of(np.sin(2 * np.pi * t / 12.5), 10)
    assert abs(per - 125) <= 10 and peak > 0.8

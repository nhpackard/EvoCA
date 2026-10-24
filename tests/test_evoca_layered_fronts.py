"""Layered-front statistics (evoca_waves.ring_spectrum, kw_spectrum) on
constructed fields: a rotating spiral and a target source give an isotropic,
prominent, moving ring at the planted wavelength; a plane wave gives a moving
but anisotropic one; a frozen mosaic an isotropic ring that does not move;
white noise no prominent ring. The space-time spectrum recovers the planted
wavelength and phase speed and puts a frozen pattern at w = 0."""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'python'))

import evoca_waves as ew                    # noqa: E402

N, LAM, T = 256, 32, 16                     # ring at k = N / LAM = 8
K = N // LAM
Y, X = np.mgrid[0:N, 0:N].astype(float)
DX = np.minimum(np.abs(X - N / 2), N - np.abs(X - N / 2))
DY = np.minimum(np.abs(Y - N / 2), N - np.abs(Y - N / 2))
R = np.hypot(DX, DY)


def _frames(f):
    return np.array([f(2 * np.pi * t / T) for t in range(T)])


def spiral():
    th = np.arctan2(Y - N / 2 - 0.3, X - N / 2 - 0.7)
    return _frames(lambda p: np.cos(th - 2 * np.pi * R / LAM + p))


def target():
    return _frames(lambda p: np.cos(2 * np.pi * R / LAM - p))


def plane():
    return _frames(lambda p: np.cos(2 * np.pi * X / LAM - p))


def mosaic(seed=0):
    rng = np.random.default_rng(seed)
    kr, _ = ew._radial_index(N)
    A = np.fft.ifft2(np.fft.fft2(rng.standard_normal((N, N))) * (np.abs(kr - K) < 1)).real
    return np.repeat(A[None], T, axis=0)


def noise(seed=1):
    return np.random.default_rng(seed).standard_normal((T, N, N))


def test_spiral_and_target_are_moving_isotropic_rings():
    for fr in (spiral(), target()):
        r = ew.ring_spectrum(fr, b=4)
        assert abs(r['k_peak'] - K) <= 1
        assert r['prominence'] > 10 and r['isotropy'] > 0.6 and r['moving'] > 0.9


def test_plane_wave_is_anisotropic():
    r = ew.ring_spectrum(plane(), b=4)
    assert r['k_peak'] == K and r['moving'] > 0.9 and r['isotropy'] < 0.2


def test_frozen_mosaic_ring_does_not_move():
    r = ew.ring_spectrum(mosaic(), b=4)
    assert abs(r['k_peak'] - K) <= 1 and r['prominence'] > 10 and r['isotropy'] > 0.6
    assert r['moving'] < 0.05


def test_reddened_noise_is_not_a_ring():
    """Smooth (low-pass) noise has its power at the largest scales: a high
    prominence at k = 1, but no ring. Seen on the CHEM-N1 frames, where the
    prominence alone ranked large-scale variance."""
    rng = np.random.default_rng(2)
    kr, _ = ew._radial_index(N)
    A = np.fft.ifft2(np.fft.fft2(rng.standard_normal((T, N, N))) / (1 + (kr / 2) ** 2)).real
    r = ew.ring_spectrum(A, b=4)
    assert not (r['contrast'] > 2)              # nan (k_peak = 1) or no rise


def test_rings_rise_above_the_low_k_dip():
    for fr in (spiral(), target(), mosaic()):
        assert ew.ring_spectrum(fr, b=4)['contrast'] > 10


def test_white_noise_has_no_prominent_ring():
    r = ew.ring_spectrum(noise(), b=4)
    assert r['prominence'] < 3


def test_single_frame_has_no_moving_fraction():
    r = ew.ring_spectrum(target()[0], b=4)
    assert abs(r['k_peak'] - K) <= 1 and np.isnan(r['moving'])


def test_kw_spectrum_recovers_wavelength_and_speed():
    period, steps = 16, 64                   # speed = LAM / period = 2 cells/step
    x = np.arange(N, dtype=float)
    ky = np.array([np.cos(2 * np.pi * (x / LAM - t / period)) for t in range(steps)])
    r = ew.kw_spectrum(ky)
    assert r['wavelength'] == LAM and abs(r['speed'] - 2.0) < 1e-9
    assert r['moving'] > 0.99 and r['peak_share'] > 0.9
    frozen = np.repeat(np.cos(2 * np.pi * x / LAM)[None], steps, axis=0)
    assert ew.kw_spectrum(frozen)['moving'] < 1e-9


def test_kw_spectrum_ignores_bulk_pulsing():
    """A frozen pattern on a uniform pulse: the pulse is k = 0 (bulk), not a
    moving front, so the moving fraction stays near 0."""
    steps, x = 64, np.arange(N, dtype=float)
    ky = np.array([np.cos(2 * np.pi * x / LAM) + 3 * np.cos(2 * np.pi * t / 16) for t in range(steps)])
    assert ew.kw_spectrum(ky)['moving'] < 1e-9



def _sources_kymo(c=2.0, period=110, steps=600, n=512, sources=(40, 170, 300, 420), seed=3):
    """Irregular excitation waves (sites fire with one period, fronts
    irregular in space): each site's phase is set by its nearest
    source, phi = t / period - |x - s| / (c * period), so fronts leave the
    sources at speed c (V shapes) and annihilate where neighbors meet (Lambda
    shapes). Sharp fronts (pulses), plus noise."""
    x = np.arange(n, dtype=float)
    d = np.min([np.minimum(np.abs(x - s), n - np.abs(x - s)) for s in sources], axis=0)
    t = np.arange(steps, dtype=float)[:, None]
    phi = (t / period - d[None, :] / (c * period)) % 1.0
    rng = np.random.default_rng(seed)
    return (phi < 0.2).astype(float) + 0.3 * rng.standard_normal((steps, n))


def test_site_period_of_irregular_fronts():
    per, peak = ew.site_period(_sources_kymo())
    assert abs(per - 110) <= 2 and peak > 0.3


def test_site_period_noise_and_frozen_have_none():
    rng = np.random.default_rng(4)
    assert np.isnan(ew.site_period(rng.standard_normal((600, 512)))[0])
    frozen = np.repeat(rng.standard_normal(512)[None], 600, axis=0)
    assert np.isnan(ew.site_period(frozen)[0])

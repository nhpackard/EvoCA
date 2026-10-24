"""Wave and spiral diagnostics for an oscillating medium (Chem CHEM-N1;
shared with Reaction-diffusion RXD-N6).

The cycle phase of a coarse block is theta = atan2(W - W0, F - F0), from its
mean food F and waste W minus a center (F0, W0), in absolute units (not
standardized). In the food -> colony -> waste -> food cycle each block goes
around the (F, W) plane once per period, so theta is its position in the
cycle. The center should be each block's TIME-averaged state (pass it as
`center`). The default, the frame's spatial mean, makes a synchronized bulk
oscillation invisible (every block sits at the mean), so it sees only
spatial heterogeneity.

- `phase_field`: theta and the cycle amplitude on the b x b coarse grid.
- `defects`: plaquettes where theta winds by +-2 pi (phase singularities).
  Noise makes many, in close +- pairs.
- `loop_winding`: the winding of theta around a square loop of half-size L.
  Measured on smooth random fields, it does NOT separate noise defects from
  cores once defects are dense: a large loop encloses many, and the net is
  +-1 about half the time. Kept for inspection, not used as a filter.
- `cores`: defects in surroundings with cycle amplitude >= amp_min.
- `matched`: defects of a later frame with a same-sign defect of an earlier
  frame nearby. A spiral core persists for many periods while the phase
  turns around it; noise defects come and go. The discriminator is the
  survival of matches at lags of several periods against a chance level
  (the earlier frame shifted by half the lattice).
- `order_parameter`: Kuramoto R = |mean exp(i theta)| over amplitude-masked
  blocks. R near 1 is bulk (synchronized) pulsing; waves and spirals give
  low R with large coherent phase domains.
- `period_of`: dominant period of a time series from its autocorrelation.

Layered fronts (N, 2026-10-24): rather than locate moving cores, detect the
spiral *arms*. Away from its core a spiral is a train of evenly spaced
fronts, so it puts a ring in the spatial power spectrum at k = N / lambda.
- `ring_spectrum`: from a stack of frames, the radial spectrum, its peak
  (wavelength), the ring's prominence over the low-k band, a crude count of
  coherent layers k / FWHM, the ring's isotropy (a ring for spirals and
  targets, two spots for a plane wave), and its moving fraction (ring power
  in the fluctuation about the time-mean field, against the total: about 1
  for moving fronts, about 0 for a frozen mosaic).
- `kw_spectrum`: from a space-time plot (rows = steps), the power at
  frequency w != 0 against all power at k != 0 (moving against stationary
  structure), and the wavelength and phase speed of the strongest moving
  mode.
Neither separates a spiral from a target source: both are layered moving
fronts. That needs a core or a handedness, which these deliberately skip.
"""
import numpy as np


def _wrap(d):
    return (d + np.pi) % (2 * np.pi) - np.pi


def block_means(A, b):
    """Mean of A over each b x b block: the (N // b)^2 coarse grid."""
    n = A.shape[0] // b
    return A[:n * b, :n * b].reshape(n, b, n, b).mean(axis=(1, 3))


def phase_field(F, W, b, center=None):
    """(theta, amplitude) on the (N // b)^2 coarse grid. center = (F0, W0),
    scalars or coarse-grid arrays (e.g. time-averaged block_means); None
    uses the frame's spatial means (see the module note)."""
    Fc, Wc = block_means(F, b), block_means(W, b)
    F0, W0 = (Fc.mean(), Wc.mean()) if center is None else center
    x, y = Fc - F0, Wc - W0
    return np.arctan2(y, x), np.hypot(x, y)


def defects(theta):
    """Integer winding (+1, -1, 0) of every plaquette (i, j)-(i+1, j+1), periodic."""
    t01 = np.roll(theta, -1, axis=1)
    t11 = np.roll(t01, -1, axis=0)
    t10 = np.roll(theta, -1, axis=0)
    wind = _wrap(t01 - theta) + _wrap(t11 - t01) + _wrap(t10 - t11) + _wrap(theta - t10)
    return np.rint(wind / (2 * np.pi)).astype(int)


def loop_winding(theta, i, j, L):
    """Winding of theta around the square loop of half-size L enclosing
    plaquette (i, j), i.e. the corners (i - L + 1 .. i + L) x (j - L + 1 .. j + L)."""
    n = theta.shape[0]
    lo_i, hi_i, lo_j, hi_j = i - L + 1, i + L, j - L + 1, j + L
    path = ([(lo_i, c) for c in range(lo_j, hi_j)] + [(r, hi_j) for r in range(lo_i, hi_i)] +
            [(hi_i, c) for c in range(hi_j, lo_j, -1)] + [(r, lo_j) for r in range(hi_i, lo_i, -1)])
    th = np.array([theta[r % n, c % n] for r, c in path])
    return int(np.rint(_wrap(np.roll(th, -1) - th).sum() / (2 * np.pi)))


def cores(theta, amp, amp_min=0.1):
    """Positions (k x 2) and signs of the defects whose 2 x 2 corner blocks
    have mean cycle amplitude >= amp_min."""
    q = defects(theta)
    a = (amp + np.roll(amp, -1, 0) + np.roll(amp, -1, 1) + np.roll(np.roll(amp, -1, 0), -1, 1)) / 4
    q[a < amp_min] = 0
    return np.argwhere(q != 0), q[q != 0]


def matched(pos_a, sgn_a, pos_b, sgn_b, n, tol=2):
    """Number of entries of a with a same-sign entry of b within tol coarse
    cells (periodic Chebyshev distance)."""
    if len(pos_a) == 0 or len(pos_b) == 0:
        return 0
    d = np.abs(pos_a[:, None, :] - pos_b[None, :, :])
    d = np.minimum(d, n - d).max(axis=2)
    return int(((d <= tol) & (sgn_a[:, None] == sgn_b[None, :])).any(axis=1).sum())


def order_parameter(theta, amp, amp_min=0.1):
    """Kuramoto R over blocks with amplitude >= amp_min (nan if none)."""
    m = amp >= amp_min
    return float(np.abs(np.exp(1j * theta[m]).mean())) if m.any() else float('nan')


def period_of(series, dt):
    """(period in steps, autocorrelation peak): the first local maximum of
    the autocorrelation after its first zero crossing; period nan if that
    peak is below 0.2."""
    x = np.asarray(series, float) - np.mean(series)
    if len(x) < 4 or x.std() == 0:
        return float('nan'), 0.0
    acf = np.correlate(x, x, 'full')[len(x) - 1:]
    acf /= acf[0]
    half = len(acf) // 2
    neg = np.where(acf[:half] < 0)[0]
    if len(neg) == 0:
        return float('nan'), 0.0
    for k in range(neg[0] + 1, half - 1):
        if acf[k] >= acf[k - 1] and acf[k] >= acf[k + 1]:
            peak = float(acf[k])
            return (float(k * dt) if peak >= 0.2 else float('nan')), peak
    return float('nan'), 0.0


def _radial_index(n):
    k = np.fft.fftfreq(n) * n
    kx, ky = np.meshgrid(k, k)
    return np.hypot(kx, ky), np.arctan2(ky, kx)


def _radial_mean(P, kr, kmax):
    idx = np.rint(kr).astype(int)
    m = idx <= kmax
    tot = np.bincount(idx[m], P[m], minlength=kmax + 1)
    cnt = np.bincount(idx[m], minlength=kmax + 1)
    return tot / np.maximum(cnt, 1)


def ring_spectrum(frames, b=4):
    """Layered-front statistics of frames (T x N x N, or one N x N frame),
    coarse-grained to b x b blocks. Returns a dict:
      k_peak       ring radius in FFT index units (1 = one wave across N)
      wavelength   N / k_peak in cells
      prominence   S(k_peak) / median S over 1 .. n/4 (white noise ~ 1)
      contrast     S(k_peak) / min S over 1 .. k_peak - 1: how far the ring
                   rises above the dip on its low-k side. A reddened spectrum
                   (power falling from k = 1) has no ring: k_peak = 1 and
                   contrast nan. Read prominence only where contrast > ~2.
      layers       k_peak / FWHM of the ring, a crude count of ordered fronts
      isotropy     1 - |sum P e^{2i phi}| / sum P over the ring (ring ~ 1, plane wave 0)
      moving       ring power of the fluctuation about the time-mean field /
                   the ring power of the frames (nan for a single frame)
      S            the radial spectrum of the frames (index 0 .. n/2)
    """
    X = np.asarray(frames, float)
    if X.ndim == 2:
        X = X[None]
    N = X.shape[1]
    C = np.array([block_means(x, b) for x in X])
    n = C.shape[1]
    C = C - C.mean(axis=(1, 2), keepdims=True)
    kr, phi = _radial_index(n)
    kmax = n // 2
    P = (np.abs(np.fft.fft2(C)) ** 2).mean(axis=0)
    S = _radial_mean(P, kr, kmax)
    band = slice(1, max(n // 4, 3) + 1)
    kp = 1 + int(np.argmax(S[band]))
    half = S[kp] / 2
    lo = kp
    while lo > 1 and S[lo - 1] > half:
        lo -= 1
    hi = kp
    while hi < kmax and S[hi + 1] > half:
        hi += 1
    # linear interpolation of the half-max crossings, at least one bin wide
    left = lo - 1 + (half - S[lo - 1]) / (S[lo] - S[lo - 1]) if lo > 0 and S[lo] != S[lo - 1] else lo - 0.5
    right = hi + (S[hi] - half) / (S[hi] - S[hi + 1]) if hi < kmax and S[hi] != S[hi + 1] else hi + 0.5
    fwhm = max(right - left, 1.0)
    ring = np.abs(kr - kp) <= 0.5
    w = P[ring]
    iso = 1.0 - float(np.abs((w * np.exp(2j * phi[ring])).sum()) / w.sum()) if w.sum() > 0 else float('nan')
    if len(C) > 1:
        D = C - C.mean(axis=0)
        Pf = (np.abs(np.fft.fft2(D)) ** 2).mean(axis=0)
        moving = float(Pf[ring].sum() / P[ring].sum()) if P[ring].sum() > 0 else float('nan')
    else:
        moving = float('nan')
    med = float(np.median(S[band]))
    contrast = float(S[kp] / S[1:kp].min()) if kp > 1 and S[1:kp].min() > 0 else (
        float('inf') if kp > 1 else float('nan'))
    return {'k_peak': kp, 'wavelength': N / kp,
            'prominence': float(S[kp] / med) if med > 0 else float('inf'),
            'contrast': contrast,
            'layers': float(kp / fwhm), 'isotropy': iso, 'moving': moving, 'S': S}


def kw_spectrum(kymo, dt=1, b=1):
    """Moving-front statistics of a space-time plot (T x N; row t = step t*dt),
    coarse-grained to b cells along x. Returns a dict:
      moving       power at w != 0 / all power, over k != 0
      k_peak, w_peak   the strongest mode with k != 0 and w != 0 (index units)
      wavelength   N / k_peak in cells
      speed        phase speed of that mode in cells per step
      peak_share   its power (both signs) / all power at k != 0
    """
    X = np.asarray(kymo, float)
    T, N = X.shape
    if b > 1:
        X = X[:, :N // b * b].reshape(T, N // b, b).mean(axis=2)
    X = X - X.mean()
    P = np.abs(np.fft.fft2(X)) ** 2
    P[:, 0] = 0.0                                  # k = 0: bulk, not structure
    tot = P.sum()
    if tot == 0:
        return {'moving': float('nan'), 'k_peak': 0, 'w_peak': 0, 'wavelength': float('nan'),
                'speed': float('nan'), 'peak_share': float('nan')}
    moving = float(P[1:].sum() / tot)
    Q = P.copy()
    Q[0] = 0.0
    if Q.max() == 0:                               # nothing moves
        return {'moving': moving, 'k_peak': 0, 'w_peak': 0, 'wavelength': float('nan'),
                'speed': float('nan'), 'peak_share': 0.0}
    wi, ki = np.unravel_index(np.argmax(Q), Q.shape)
    w_idx = wi if wi <= T // 2 else wi - T
    k_idx = ki if ki <= X.shape[1] // 2 else ki - X.shape[1]
    share = float((Q[wi, ki] + Q[-wi % T, -ki % X.shape[1]]) / tot)
    speed = abs(w_idx / (T * dt)) / (abs(k_idx) / N)
    return {'moving': moving, 'k_peak': abs(int(k_idx)), 'w_peak': abs(int(w_idx)),
            'wavelength': N / abs(k_idx), 'speed': float(speed), 'peak_share': share}

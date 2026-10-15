"""Pure-Python reference implementation of the drift-shadow replicate
process (Docs/shadows.md §5.2), for bit-parity with the C core.

It ports, draw for draw, what C/evoca.c does per replicate per step
(drift_rep_step with genome_mutate / marker_mutate in replay mode and the
incremental tallies), including the xorshift32 per-replicate streams and
their seed derivation. Given the same founder genomes, seed, and event
stream (deaths D, births B, per-birth realized counts (nf, na, ne, nm)),
the tallies c[i] (allele-1 counts), e[key] (active-egene ternary-key
counts) and m[key] (marker key counts) must be bit-identical to the C
replicates at every step (tests/test_drift_shadow.py).

It is deliberately slow and literal: it exists to pin the C code's
semantics, not to run calibrations (those use the C synthetic hook).
Restricted mutation is not modeled (the synthetic hook requires it off).
"""
import numpy as np

LUT_BITS = 250
NEG = 8
M32 = 0xFFFFFFFF
M64 = 0xFFFFFFFFFFFFFFFF
DEFAULT_SEED = 0x12345678


def xorshift32(state):
    x = state
    x ^= (x << 13) & M32
    x ^= x >> 17
    x ^= (x << 5) & M32
    return x & M32


def seed_mix(z):
    z = (z + 0x9E3779B97F4A7C15) & M64
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & M64
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & M64
    z ^= z >> 31
    r = (z >> 32) & M32
    return r if r else 0x9E3779B9


def drift_seed_for(base, r):
    return seed_mix((base ^ 0x4452494654534844 ^ ((r & M32) * 0x9E3779B9)) & M64)


def ternary_key(value, mask):
    key, p = 0, 1
    for i in range(6):
        if (mask >> i) & 1:
            key += (1 if ((value >> i) & 1) == 0 else 2) * p
        p *= 3
    return key


def marker_key(m):
    return ternary_key(m & 0x3F, (m >> 8) & 0x3F)


class Replicate:
    """One shadow replicate: genome list + tallies + its own stream."""

    def __init__(self, lut, eg, egm, act, mk, seed, p_dup=1.0):
        self.lut = [bytearray(x) for x in lut]          # 32 bytes each
        self.eg = [bytearray(x) for x in eg]
        self.egm = [bytearray(x) for x in egm]
        self.act = [int(a) for a in act]
        self.mk = [int(m) for m in mk]
        self.s = seed
        self.thr = int(np.float32(p_dup) * np.float32(65536.0))
        self.recount()

    def next(self):
        self.s = xorshift32(self.s)
        return self.s

    # tallies
    def _tally(self, j, sign):
        L = self.lut[j]
        for bit in range(LUT_BITS):
            if (L[bit >> 3] >> (bit & 7)) & 1:
                self.c[bit] += sign
        a = self.act[j]
        for s in range(NEG):
            if (a >> s) & 1:
                self.e[ternary_key(self.eg[j][s] & 0x3F, self.egm[j][s] & 0x3F)] += sign
                self.T += sign
        self.m[marker_key(self.mk[j])] += sign

    def recount(self):
        self.c = np.zeros(LUT_BITS, np.int64)
        self.e = np.zeros(729, np.int64)
        self.m = np.zeros(729, np.int64)
        self.T = 0
        for j in range(len(self.lut)):
            self._tally(j, +1)

    def _copy(self, dst, src):
        vals = (bytearray(self.lut[src]), bytearray(self.eg[src]),
                bytearray(self.egm[src]), self.act[src], self.mk[src])
        if dst == len(self.lut):
            self.lut.append(vals[0]); self.eg.append(vals[1])
            self.egm.append(vals[2]); self.act.append(vals[3])
            self.mk.append(vals[4])
        else:
            (self.lut[dst], self.eg[dst], self.egm[dst], self.act[dst],
             self.mk[dst]) = vals

    def _mutate(self, j, nf, na, ne, nm):
        L, E, EM = self.lut[j], self.eg[j], self.egm[j]
        for _ in range(nf):
            bit = self.next() % LUT_BITS
            L[bit >> 3] ^= (1 << (bit & 7))
        before = self.act[j]
        for _ in range(na):
            b = self.next() & 7
            new = self.act[j] ^ (1 << b)
            if new != 0:
                self.act[j] = new
        newly = self.act[j] & ~before & 0xFF
        while newly:
            s = (newly & -newly).bit_length() - 1
            newly &= newly - 1
            r16 = self.next() & 0xFFFF
            if r16 < self.thr:
                n_act = bin(self.act[j]).count('1')
                pick = self.next() % n_act
                a2 = self.act[j]
                while pick > 0:
                    a2 &= a2 - 1
                    pick -= 1
                src = (a2 & -a2).bit_length() - 1
                E[s] = E[src]
                EM[s] = EM[src]
        for _ in range(ne):
            slot = self.next() & 7
            b12 = self.next() % 12
            if b12 < 6:
                E[slot] ^= (1 << b12)
            else:
                EM[slot] ^= (1 << (b12 - 6))
        for _ in range(nm):
            b = self.next() % 12
            self.mk[j] ^= (1 << b) if b < 6 else (1 << (b - 6 + 8))

    def step(self, D, recs):
        n = len(self.lut)
        D = min(D, n)
        for _ in range(D):
            v = self.next() % n
            self._tally(v, -1)
            if v != n - 1:
                self._copy(v, n - 1)
            for arr in (self.lut, self.eg, self.egm, self.act, self.mk):
                arr.pop()
            n -= 1
        n_surv = n
        if n_surv == 0:
            return
        for (nf, na, ne, nm) in recs:
            p = self.next() % n_surv
            self._copy(n, p)
            self._mutate(n, nf, na, ne, nm)
            self._tally(n, +1)
            n += 1

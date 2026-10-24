#ifndef CORE_H
#define CORE_H

/*
 * Shared core: pure helpers used by both models in the library
 * (EvoCA in evoca.c, ChemCA in chemca.c). Docs/2026-10-23_chem-model-v2.md §9b.
 *
 * Nothing here holds state. Each model owns its own PRNG state word and
 * passes a step function in where a stream is needed, so the two models
 * never share a random stream.
 */

#include <stdint.h>

/* xorshift32 PRNG step on an explicit state word. A valid state is never 0. */
static inline uint32_t core_xorshift32(uint32_t *s)
{
    uint32_t x = *s;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    return *s = x;
}

/* Periodic index wrap: ((i % N) + N) % N. */
static inline int core_wrap(int i, int N) { return ((i % N) + N) % N; }

/* Poisson sample (Knuth's algorithm) drawing 24-bit uniforms from the
 * stream whose step function is `next`. Returns 0 for lambda < 1e-6. */
int core_poisson(uint32_t (*next)(void), float lambda);

/* Diffusion of an N x N periodic field by strength d >= 0:
 * floor(d) full 3x3 box-blur passes, then one partial pass of strength
 * s = d - floor(d), F' = (1 - s) F + s * mean_3x3(F) (s = 1 is a full blur,
 * s = 0 is a no-op), so the knob is continuous below 1 as well as above.
 * Each pass writes into *tmp and swaps the two pointers, so after the call
 * *field points at the result and *tmp at scratch. d <= 0 does nothing. */
void core_diffuse(float **field, float **tmp, int N, float d);

#endif /* CORE_H */

#include "core.h"
#include <math.h>

/* Poisson sample (Knuth's algorithm). See core.h. */
int core_poisson(uint32_t (*next)(void), float lambda)
{
    if (lambda < 1e-6f) return 0;
    float L = expf(-lambda);
    float p = 1.0f;
    int k = 0;
    do {
        k++;
        p *= (float)(next() & 0xFFFFFF) / 16777216.0f;
    } while (p > L);
    return k - 1;
}

/* 3x3 periodic box-blur diffusion (CHEM-P1: the one copy, replacing the
 * food full pass, the food partial pass and the compound-field pass that
 * evoca.c used to carry). A full pass is the partial rule at s = 1, where
 * t = 0 and w9 = 1/9, so t*F + w9*sum == sum/9 bit for bit (0*F adds a
 * signed zero; IEEE x + (+-0) == x for x != 0, and +0 when x == 0).
 *
 * precise = 0 is EvoCA's float rule, kept bit-exact for its goldens. Its
 * float weights do not sum to 1 (t + 9*w9 != 1: 1/9 rounds up, 1 - s
 * rounds), a systematic mass bias of order 1e-8 per pass. precise = 1
 * (ChemCA, which books mass exactly) accumulates in double and writes
 * F + s*(sum/9 - F), so the only error is the unbiased final rounding. */
static void diffuse(float **field, float **tmp, int N, float d, int precise)
{
    if (!(d > 0.0f)) return;
    int n_full = (int)d;
    float frac = d - (float)n_full;
    for (int pass = 0; pass <= n_full; pass++) {
        float s = (pass < n_full) ? 1.0f : frac;
        if (s <= 0.0f) break;
        float t = 1.0f - s, w9 = s * (1.0f / 9.0f);
        float *src = *field, *dst = *tmp;
        for (int row = 0; row < N; row++) {
            for (int col = 0; col < N; col++) {
                if (precise) {
                    double sum = 0.0;
                    for (int di = -1; di <= 1; di++) {
                        int r = core_wrap(row + di, N);
                        for (int dj = -1; dj <= 1; dj++)
                            sum += src[r * N + core_wrap(col + dj, N)];
                    }
                    double f = src[row * N + col];
                    dst[row * N + col] = (float)(f + (double)s * (sum / 9.0 - f));
                    continue;
                }
                float sum = 0.0f;
                for (int di = -1; di <= 1; di++) {
                    int r = core_wrap(row + di, N);
                    for (int dj = -1; dj <= 1; dj++) {
                        int c = core_wrap(col + dj, N);
                        sum += src[r * N + c];
                    }
                }
                dst[row * N + col] = t * src[row * N + col] + w9 * sum;
            }
        }
        *tmp = src; *field = dst;
    }
}

void core_diffuse(float **field, float **tmp, int N, float d)
{
    diffuse(field, tmp, N, d, 0);
}

void core_diffuse_precise(float **field, float **tmp, int N, float d)
{
    diffuse(field, tmp, N, d, 1);
}

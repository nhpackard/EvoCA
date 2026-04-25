#include "evoca.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/* xorshift32 PRNG — stochastic tie-breaking in reproduction */
static uint32_t g_rng = 0x12345678u;
static inline uint32_t rng_next(void)
{
    g_rng ^= g_rng << 13;
    g_rng ^= g_rng >> 17;
    g_rng ^= g_rng << 5;
    return g_rng;
}

/*
 * Distance-ring classification for each (di,dj) offset, indexed [di+2][dj+2].
 *
 * ring_map: which of the three LUT rings (0-2) the cell belongs to
 *   0 : dist 1   (n1, 4 cells)
 *   1 : dist √2  (n2, 4 cells)
 *   2 : dist 2   (n3, 4 cells)
 *  -1 : centre or outside LUT neighbourhood (skipped)
 *
 * orbit_map: D4 orbit (0-5) for fiducial-pattern matching (full 5×5)
 *   0 : centre                    (1 cell)
 *   1 : dist-1 axis cells         (4 cells)
 *   2 : dist-2 axis cells         (4 cells)
 *   3 : dist-√2 diagonal cells    (4 cells)
 *   4 : dist-2√2 corner cells     (4 cells)
 *   5 : dist-√5 off-axis cells    (8 cells)
 */
static const int ring_map[5][5] = {
    {-1, -1,  2, -1, -1},
    {-1,  1,  0,  1, -1},
    { 2,  0, -1,  0,  2},
    {-1,  1,  0,  1, -1},
    {-1, -1,  2, -1, -1},
};

static const int orbit_map[5][5] = {
    {4, 5, 2, 5, 4},
    {5, 3, 1, 3, 5},
    {2, 1, 0, 1, 2},
    {5, 3, 1, 3, 5},
    {4, 5, 2, 5, 4},
};

/* ── Bit-packed LUT helpers ─────────────────────────────────────── */

static inline uint8_t lut_get(const uint8_t *b, int bit)
{
    return (b[bit >> 3] >> (bit & 7)) & 1;
}

static inline void lut_flip(uint8_t *b, int bit)
{
    b[bit >> 3] ^= (uint8_t)(1u << (bit & 7));
}

/* FNV-1a hash of a LUT → 32-bit value for genome coloring.
 * Hashes as uint32_t words (LUT_BYTES=32 = 8 words, no remainder). */
static uint32_t lut_hash_fn(const uint8_t *b)
{
    const uint32_t *w = (const uint32_t *)b;
    uint32_t h = 0x811c9dc5u;
    for (int i = 0; i < LUT_BYTES / 4; i++) {
        h ^= w[i];
        h *= 0x01000193u;
    }
    return h;
}

/* Map LUT hash → ARGB color.  Wild-type hash maps to white. */
static int32_t hash_to_color(uint32_t h, uint32_t wt)
{
    if (h == wt) return (int32_t)0xFFFFFFFFu;
    return (int32_t)(0xFF000000u | (h & 0x00FFFFFFu));
}

/* Poisson sample using Knuth's algorithm with the existing PRNG. */
static int poisson_sample(float lambda)
{
    if (lambda < 1e-6f) return 0;
    float L = expf(-lambda);
    float p = 1.0f;
    int k = 0;
    do {
        k++;
        p *= (float)(rng_next() & 0xFFFFFF) / 16777216.0f;
    } while (p > L);
    return k - 1;
}

/* ── Grid state ─────────────────────────────────────────────────── */

static int    gN          = 0;
static float  gfood_inc   = 0.0f;
static float  gm_scale    = 1.0f;
/* food_repro removed — reproduction threshold hardcoded to 1.0 */
static int    ggdiff      = 0;    /* diffusion passes per step */
static float  gmu_lut     = 0.0f; /* per-bit LUT mutation rate */
static float  gmu_egenome  = 0.0f; /* per-bit egenome mutation rate */
static float  gtax        = 0.0f; /* priv food decrement per step */
static int    grestricted_mu = 0; /* 0=random mutation, 1=restricted to active bits */
static int    g_diag        = 0; /* diagnostic prints */

/* Restricted mutation: tracks which LUT bit indices are queried each step */
static uint8_t lut_active[LUT_BYTES];   /* 250-bit mask */
static int     active_bits[LUT_BITS];   /* indices of set bits */
static int     n_active = 0;

static uint8_t *v_curr = NULL;   /* [N*N]            */
static uint8_t *v_next = NULL;   /* [N*N]            */
static uint8_t *lut    = NULL;   /* [N*N * LUT_BYTES] */
static uint8_t *egenome = NULL;   /* [N*N]            */
static float   *f_priv = NULL;   /* [N*N]            */
static float   *F_food = NULL;   /* [N*N]            */
static float   *F_temp = NULL;   /* [N*N] scratch for diffusion */
static uint8_t  *births    = NULL;   /* [N*N] birth events this step */
static uint32_t *lut_color = NULL;   /* [N*N] cached ARGB per cell */
static uint32_t *lut_hash_cache = NULL; /* [N*N] cached FNV-1a hash */
static uint32_t  wt_hash   = 0;     /* hash of wild-type LUT */
static uint32_t  g_step    = 0;     /* global step counter */
static uint8_t  *env_mask  = NULL;  /* [N*N] food regen mask; 1=regen, 0=no */
static uint8_t  *alive     = NULL;  /* [N*N] 1=alive organism, 0=dead slot */

/* Adaptive age colormap scale: tracks max observed age, decays slowly. */
static double    g_age_scale = 50.0;

/* ── Reproduction age histogram ────────────────────────────────── */

#define REPRO_AGE_MAX 1024          /* bins 0..1023; overflow in last bin */

static uint32_t *last_event_step = NULL;  /* [N*N] step of birth or repro */
static uint32_t  repro_age_hist[REPRO_AGE_MAX]; /* cumulative histogram */
static uint32_t  repro_age_t0 = 0;        /* start accumulating after this step */

/* ── Egenome activity (fixed-size, 64 entries) ──────────────────── */

#define EGENOME_COUNT 64

static uint64_t eg_act[EGENOME_COUNT];        /* cumulative activity */
static uint32_t eg_pop[EGENOME_COUNT];        /* current population */
static int32_t  eg_color[EGENOME_COUNT];      /* ARGB color per egenome */
static uint8_t  wt_egenome_val = 0;           /* wild-type egenome */
static int      eg_act_ymax = 2000;          /* Y-scale for egenome activity */

/* Precompute colors for all 64 egenomes.  Wild-type = white;
 * others = FNV-1a hash of the byte value → ARGB. */
static void eg_init_colors(uint8_t wt)
{
    wt_egenome_val = wt;
    for (int i = 0; i < EGENOME_COUNT; i++) {
        if ((uint8_t)i == wt) {
            eg_color[i] = (int32_t)0xFFFFFFFFu;
        } else {
            uint32_t h = 0x811c9dc5u ^ (uint32_t)i;
            h *= 0x01000193u;
            /* Ensure minimum brightness: set each channel to at least 0x40 */
            uint8_t r = (uint8_t)((h >> 16) & 0xFF); if (r < 0x40) r |= 0x80;
            uint8_t g = (uint8_t)((h >>  8) & 0xFF); if (g < 0x40) g |= 0x80;
            uint8_t b = (uint8_t)( h        & 0xFF); if (b < 0x40) b |= 0x80;
            eg_color[i] = (int32_t)(0xFF000000u | ((uint32_t)r << 16)
                                    | ((uint32_t)g << 8) | b);
        }
    }
}

/* ── Local-pattern activity & entropy ─────────────────────────── */

#define PAT_MAX_BITS 13
#define PAT_MAX_SIZE (1 << PAT_MAX_BITS)   /* 8192 */

static int      g_n_ent      = 2;          /* 1=VN, 2=Moore, 3=n3 */
static int      g_pat_bits   = 9;          /* bits in pattern key */
static int      g_pat_size   = 512;        /* 2^pat_bits */
static float    g_entropy    = 0.0f;       /* Shannon entropy (bits) */
static uint64_t pat_act[PAT_MAX_SIZE];     /* cumulative presence */
static uint32_t pat_pop[PAT_MAX_SIZE];     /* current step count */
static int32_t  pat_color[PAT_MAX_SIZE];   /* ARGB per pattern key */
static int      pat_act_ymax = 2000;

void evoca_set_n_ent(int n)
{
    if (n < 1) n = 1;
    if (n > 3) n = 3;
    g_n_ent    = n;
    g_pat_bits = (n == 1) ? 5 : (n == 2) ? 9 : 13;
    g_pat_size = 1 << g_pat_bits;
    memset(pat_act, 0, sizeof(pat_act));
    memset(pat_pop, 0, sizeof(pat_pop));
    /* Assign FNV-1a hash colors to all patterns */
    for (int i = 0; i < g_pat_size; i++) {
        uint32_t h = 0x811c9dc5u;
        h ^= (uint32_t)(i & 0xFF);    h *= 0x01000193u;
        h ^= (uint32_t)(i >> 8);      h *= 0x01000193u;
        uint8_t r = (uint8_t)((h >> 16) & 0xFF); if (r < 0x40) r |= 0x80;
        uint8_t g = (uint8_t)((h >>  8) & 0xFF); if (g < 0x40) g |= 0x80;
        uint8_t b = (uint8_t)( h        & 0xFF); if (b < 0x40) b |= 0x80;
        pat_color[i] = (int32_t)(0xFF000000u | ((uint32_t)r << 16)
                                 | ((uint32_t)g << 8) | b);
    }
}
int   evoca_get_n_ent(void)        { return g_n_ent; }
void  evoca_set_pat_act_ymax(int y){ pat_act_ymax = y > 1 ? y : 1; }
int   evoca_get_pat_act_ymax(void) { return pat_act_ymax; }

/* Extract local pattern key for cell (row,col).
 * Bit 0 = center; bits 1–4 = N,S,E,W; bits 5–8 = NE,NW,SE,SW;
 * bits 9–12 = (row-2,col),(row+2,col),(row,col+2),(row,col-2).
 * All accesses use periodic boundary. */
static uint16_t extract_pattern(int row, int col)
{
    int N = gN;
    uint16_t key = 0;
    int bit = 0;

    /* Center */
    key |= (uint16_t)(v_curr[row * N + col]) << bit++;

    /* Ring 1: N, S, E, W */
    int rN = (row - 1 + N) % N, rS = (row + 1) % N;
    int cE = (col + 1) % N,     cW = (col - 1 + N) % N;
    key |= (uint16_t)(v_curr[rN * N + col]) << bit++;
    key |= (uint16_t)(v_curr[rS * N + col]) << bit++;
    key |= (uint16_t)(v_curr[row * N + cE]) << bit++;
    key |= (uint16_t)(v_curr[row * N + cW]) << bit++;

    if (g_n_ent == 1) return key;

    /* Ring 2 diagonals: NE, NW, SE, SW */
    key |= (uint16_t)(v_curr[rN * N + cE]) << bit++;
    key |= (uint16_t)(v_curr[rN * N + cW]) << bit++;
    key |= (uint16_t)(v_curr[rS * N + cE]) << bit++;
    key |= (uint16_t)(v_curr[rS * N + cW]) << bit++;

    if (g_n_ent == 2) return key;

    /* Ring 3: (row±2,col) and (row,col±2) */
    int rN2 = (row - 2 + N) % N, rS2 = (row + 2) % N;
    int cE2 = (col + 2) % N,     cW2 = (col - 2 + N) % N;
    key |= (uint16_t)(v_curr[rN2 * N + col]) << bit++;
    key |= (uint16_t)(v_curr[rS2 * N + col]) << bit++;
    key |= (uint16_t)(v_curr[row * N + cE2]) << bit++;
    key |= (uint16_t)(v_curr[row * N + cW2]) << bit++;

    return key;
}

void evoca_pat_update(void)
{
    if (!v_curr) return;
    int N = gN;
    int sz = g_pat_size;

    /* Clear current-step counts */
    memset(pat_pop, 0, (size_t)sz * sizeof(uint32_t));

    /* Scan every site */
    for (int row = 0; row < N; row++)
        for (int col = 0; col < N; col++) {
            uint16_t key = extract_pattern(row, col);
            pat_pop[key]++;
            pat_act[key]++;
        }

    /* Shannon entropy H = -Σ p·log₂(p) */
    double H = 0.0;
    double total = (double)(N * N);
    for (int i = 0; i < sz; i++) {
        if (pat_pop[i] == 0) continue;
        double p = pat_pop[i] / total;
        H -= p * log2(p);
    }
    g_entropy = (float)H;
}

float evoca_get_entropy(void) { return g_entropy; }

void evoca_pat_activity_render_col(int32_t *col, int height)
{
    for (int y = 0; y < height; y++)
        col[y] = (int32_t)0xFF111111u;

    int sz = g_pat_size;
    uint64_t ymax = (uint64_t)pat_act_ymax;

    uint32_t ypop[height];
    memset(ypop, 0, (size_t)height * sizeof(uint32_t));

    /* Pass 1: extinct patterns — dimmed */
    for (int i = 0; i < sz; i++) {
        if (pat_act[i] == 0 || pat_pop[i] > 0) continue;
        uint64_t act = pat_act[i];
        int y = (height - 1) - (int)((uint64_t)(height - 1) * act / (act + ymax));
        if (y < 0) y = 0; if (y >= height) y = height - 1;
        uint32_t c = (uint32_t)pat_color[i];
        uint8_t r = (uint8_t)(((c >> 16) & 0xFF) * 15 / 100);
        uint8_t g = (uint8_t)(((c >>  8) & 0xFF) * 15 / 100);
        uint8_t b = (uint8_t)(( c        & 0xFF) * 15 / 100);
        col[y] = (int32_t)(0xFF000000u | ((uint32_t)r << 16)
                           | ((uint32_t)g << 8) | b);
    }

    /* Pass 2: alive patterns — full color, higher pop wins */
    for (int i = 0; i < sz; i++) {
        if (pat_pop[i] == 0) continue;
        uint64_t act = pat_act[i];
        int y = (height - 1) - (int)((uint64_t)(height - 1) * act / (act + ymax));
        if (y < 0) y = 0; if (y >= height) y = height - 1;
        if (pat_pop[i] >= ypop[y]) {
            col[y] = pat_color[i];
            ypop[y] = pat_pop[i];
        }
    }
}

/* ── Activity hash table ───────────────────────────────────────── */

#define ACT_INIT_CAP 4096
#define ACT_EMPTY    0        /* key=0 is the empty sentinel */

/* Per-bucket pop history ring buffer for the flux probe. */
#define FLUX_W_MAX   64

typedef struct {
    uint64_t activity;    /* cumulative count */
    uint32_t pop_count;   /* current population */
    int32_t  color;       /* ARGB (same as lut_color) */
} act_entry_t;

static uint32_t    *act_keys     = NULL;
static act_entry_t *act_vals     = NULL;
static uint32_t    *act_pop_hist = NULL;   /* [act_cap * FLUX_W_MAX] */
static int          act_cap      = 0;
static int          act_cnt      = 0;
static int          act_ymax     = 2000;   /* Y-scale for activity saturation */
static uint32_t     flux_t       = 0;      /* monotonic tick counter for pop_hist */

static void act_resize(void)
{
    int new_cap = act_cap * 2;
    if (g_diag)
        fprintf(stderr, "DIAG: act_resize %d -> %d (cnt=%d, step=%u)\n",
                act_cap, new_cap, act_cnt, g_step);
    uint32_t    *nk = calloc((size_t)new_cap, sizeof(uint32_t));
    act_entry_t *nv = calloc((size_t)new_cap, sizeof(act_entry_t));
    uint32_t    *nh = calloc((size_t)new_cap * FLUX_W_MAX, sizeof(uint32_t));
    for (int i = 0; i < act_cap; i++) {
        if (act_keys[i] == ACT_EMPTY) continue;
        uint32_t slot = act_keys[i] % (uint32_t)new_cap;
        while (nk[slot] != ACT_EMPTY)
            slot = (slot + 1) % (uint32_t)new_cap;
        nk[slot] = act_keys[i];
        nv[slot] = act_vals[i];
        memcpy(nh + (size_t)slot * FLUX_W_MAX,
               act_pop_hist + (size_t)i * FLUX_W_MAX,
               FLUX_W_MAX * sizeof(uint32_t));
    }
    free(act_keys); free(act_vals); free(act_pop_hist);
    act_keys = nk; act_vals = nv; act_pop_hist = nh; act_cap = new_cap;
}

/* Find or insert; returns pointer to value entry. */
static act_entry_t *act_find_or_insert(uint32_t key, int32_t color)
{
    if (act_cnt * 10 >= act_cap * 7) act_resize();
    uint32_t slot = key % (uint32_t)act_cap;
    while (act_keys[slot] != ACT_EMPTY) {
        if (act_keys[slot] == key) return &act_vals[slot];
        slot = (slot + 1) % (uint32_t)act_cap;
    }
    act_keys[slot] = key;
    act_vals[slot].activity  = 0;
    act_vals[slot].pop_count = 0;
    act_vals[slot].color     = color;
    act_cnt++;
    return &act_vals[slot];
}

/* Prune extinct genomes with low activity to keep the table small.
 * Keeps all alive genomes (pop_count > 0) and extinct genomes whose
 * activity >= threshold.  Already-rendered strip-chart columns preserve
 * the history; we only lose dim dots in *future* columns for pruned entries. */
static void act_compact(uint64_t threshold)
{
    int keep = 0;
    for (int i = 0; i < act_cap; i++) {
        if (act_keys[i] == ACT_EMPTY) continue;
        if (act_vals[i].pop_count > 0 || act_vals[i].activity >= threshold)
            keep++;
    }

    /* Pick capacity to keep load < 70% */
    int new_cap = ACT_INIT_CAP;
    while (new_cap * 7 < keep * 10 + 10) new_cap *= 2;

    uint32_t    *nk = calloc((size_t)new_cap, sizeof(uint32_t));
    act_entry_t *nv = calloc((size_t)new_cap, sizeof(act_entry_t));
    uint32_t    *nh = calloc((size_t)new_cap * FLUX_W_MAX, sizeof(uint32_t));

    for (int i = 0; i < act_cap; i++) {
        if (act_keys[i] == ACT_EMPTY) continue;
        if (act_vals[i].pop_count == 0 && act_vals[i].activity < threshold)
            continue;   /* prune */
        uint32_t slot = act_keys[i] % (uint32_t)new_cap;
        while (nk[slot] != ACT_EMPTY)
            slot = (slot + 1) % (uint32_t)new_cap;
        nk[slot] = act_keys[i];
        nv[slot] = act_vals[i];
        memcpy(nh + (size_t)slot * FLUX_W_MAX,
               act_pop_hist + (size_t)i * FLUX_W_MAX,
               FLUX_W_MAX * sizeof(uint32_t));
    }

    int old_cnt = act_cnt, old_cap = act_cap;
    free(act_keys); free(act_vals); free(act_pop_hist);
    act_keys = nk; act_vals = nv; act_pop_hist = nh;
    act_cap  = new_cap;
    act_cnt  = keep;

    if (g_diag)
        fprintf(stderr, "DIAG: act_compact %d/%d -> %d/%d thresh=%llu (step=%u)\n",
                old_cnt, old_cap, act_cnt, act_cap,
                (unsigned long long)threshold, g_step);
}

static void evoca_activity_init(void)
{
    act_cap = ACT_INIT_CAP;
    act_cnt = 0;
    act_keys     = calloc((size_t)act_cap, sizeof(uint32_t));
    act_vals     = calloc((size_t)act_cap, sizeof(act_entry_t));
    act_pop_hist = calloc((size_t)act_cap * FLUX_W_MAX, sizeof(uint32_t));
    flux_t       = 0;
}

static void evoca_activity_free(void)
{
    free(act_keys);     act_keys     = NULL;
    free(act_vals);     act_vals     = NULL;
    free(act_pop_hist); act_pop_hist = NULL;
    act_cap = 0;
    act_cnt = 0;
    flux_t  = 0;
}

/* ── Lifecycle ──────────────────────────────────────────────────── */

void evoca_init(int N, float food_inc, float m_scale)
{
    evoca_free();
    gN          = N;
    gfood_inc   = food_inc;
    gm_scale    = m_scale;

    size_t cells = (size_t)N * N;
    v_curr = calloc(cells,               sizeof(uint8_t));
    v_next = calloc(cells,               sizeof(uint8_t));
    lut    = calloc(cells * LUT_BYTES,   sizeof(uint8_t));
    egenome = calloc(cells,               sizeof(uint8_t));
    f_priv = calloc(cells,               sizeof(float));
    F_food = calloc(cells,               sizeof(float));
    F_temp = calloc(cells,               sizeof(float));
    births    = calloc(cells,               sizeof(uint8_t));
    lut_color = calloc(cells,               sizeof(uint32_t));
    lut_hash_cache = calloc(cells,          sizeof(uint32_t));
    last_event_step = calloc(cells,         sizeof(uint32_t));
    env_mask = malloc(cells * sizeof(uint8_t));
    memset(env_mask, 1, cells);   /* default: all sites regenerate */
    alive = malloc(cells * sizeof(uint8_t));
    memset(alive, 1, cells);      /* default: all cells alive */
    memset(repro_age_hist, 0, sizeof(repro_age_hist));
    g_step = 0;
    evoca_activity_init();
    memset(eg_act, 0, sizeof(eg_act));
    memset(eg_pop, 0, sizeof(eg_pop));
    evoca_set_n_ent(g_n_ent);   /* reset pattern arrays, keep current n_ent */
}

void evoca_free(void)
{
    free(v_curr); v_curr = NULL;
    free(v_next); v_next = NULL;
    free(lut);    lut    = NULL;
    free(egenome); egenome = NULL;
    free(f_priv); f_priv = NULL;
    free(F_food); F_food = NULL;
    free(F_temp); F_temp = NULL;
    free(births);    births    = NULL;
    free(lut_color); lut_color = NULL;
    free(lut_hash_cache); lut_hash_cache = NULL;
    free(last_event_step); last_event_step = NULL;
    free(env_mask);  env_mask  = NULL;
    free(alive);     alive     = NULL;
    evoca_activity_free();
    memset(eg_act,  0, sizeof(eg_act));
    memset(eg_pop,  0, sizeof(eg_pop));
    memset(pat_act, 0, sizeof(pat_act));
    memset(pat_pop, 0, sizeof(pat_pop));
    gN = 0;
}

/* ── Metaparam setters ──────────────────────────────────────────── */

void evoca_set_food_inc(float f)   { gfood_inc   = f; }
void evoca_set_m_scale(float m)    { gm_scale    = m; }
void evoca_set_gdiff(int d)        { ggdiff      = d; }
int  evoca_get_gdiff(void)         { return ggdiff;   }
void  evoca_set_mu_lut(float m)    { gmu_lut     = m; }
void  evoca_set_mu_egenome(float m) { gmu_egenome  = m; }
float evoca_get_mu_lut(void)       { return gmu_lut;    }
float evoca_get_mu_egenome(void)    { return gmu_egenome;  }
void  evoca_set_restricted_mu(int r) { grestricted_mu = r; }
int   evoca_get_restricted_mu(void)  { return grestricted_mu; }
void  evoca_set_tax(float t)      { gtax       = t; }
float evoca_get_tax(void)         { return gtax;      }
void  evoca_set_diag(int d)      { g_diag     = d; }
int   evoca_get_diag(void)       { return g_diag;    }

/* ── Bulk setters ───────────────────────────────────────────────── */

void evoca_set_v_all(const uint8_t *v, int len)
{
    memcpy(v_curr, v, (size_t)len);
}

void evoca_set_lut_all(const uint8_t *lb)
{
    size_t cells = (size_t)gN * gN;
    for (size_t i = 0; i < cells; i++)
        memcpy(lut + i * LUT_BYTES, lb, LUT_BYTES);
    wt_hash = lut_hash_fn(lb);
    int32_t wt_color = hash_to_color(wt_hash, wt_hash);  /* white */
    for (size_t i = 0; i < cells; i++) {
        lut_color[i] = (uint32_t)wt_color;
        lut_hash_cache[i] = wt_hash;
    }
}

void evoca_set_lut(int idx, const uint8_t *lb)
{
    memcpy(lut + (size_t)idx * LUT_BYTES, lb, LUT_BYTES);
    uint32_t h = lut_hash_fn(lb);
    lut_color[idx] = (uint32_t)hash_to_color(h, wt_hash);
    lut_hash_cache[idx] = h;
}

void evoca_set_egenome_all(uint8_t eg) {
    memset(egenome, eg, (size_t)gN * gN);
    eg_init_colors(eg);
}

void evoca_set_f_all(float f)
{
    for (size_t i = 0; i < (size_t)gN * gN; i++) f_priv[i] = f;
}

void evoca_set_F_all(float F)
{
    for (size_t i = 0; i < (size_t)gN * gN; i++) F_food[i] = F;
}

void evoca_set_env_mask(const uint8_t *mask)
{
    memcpy(env_mask, mask, (size_t)gN * gN);
}

uint8_t *evoca_get_env_mask(void) { return env_mask; }

/* ── Internal helpers ───────────────────────────────────────────── */

/*
 * Compute the LUT bit index for cell at (row,col).
 * Counts active cells in rings n1 (dist 1), n2 (dist √2), n3 (dist 2),
 * plus reads v_x (current cell state) from v_curr.
 */
static int compute_lut_bit(int row, int col)
{
    int N  = gN;
    int v_x = v_curr[row * N + col];
    int n[3] = {0, 0, 0};

    for (int di = -2; di <= 2; di++) {
        int r = ((row + di) % N + N) % N;
        for (int dj = -2; dj <= 2; dj++) {
            int ring = ring_map[di+2][dj+2];
            if (ring < 0) continue;
            int c = ((col + dj) % N + N) % N;
            if (v_curr[r * N + c]) n[ring]++;
        }
    }
    return LUT_IDX(v_x, n[0], n[1], n[2]);
}

/* Count matches between actual 5×5 config and fiducial pattern. */
static int fiducial_matches(int row, int col, uint8_t eg)
{
    int N = gN, matches = 0;
    for (int di = -2; di <= 2; di++) {
        int r = ((row + di) % N + N) % N;
        for (int dj = -2; dj <= 2; dj++) {
            int c      = ((col + dj) % N + N) % N;
            int orbit  = orbit_map[di+2][dj+2];
            int fid    = (eg >> orbit) & 1;
            if (v_curr[r * N + c] == fid) matches++;
        }
    }
    return matches;
}

/* ── Food diffusion (3×3 box blur, periodic) ───────────────────── */

static void diffuse_food_once(void)
{
    int N = gN;
    for (int row = 0; row < N; row++) {
        for (int col = 0; col < N; col++) {
            float sum = 0.0f;
            for (int di = -1; di <= 1; di++) {
                int r = ((row + di) % N + N) % N;
                for (int dj = -1; dj <= 1; dj++) {
                    int c = ((col + dj) % N + N) % N;
                    sum += F_food[r * N + c];
                }
            }
            F_temp[row * N + col] = sum * (1.0f / 9.0f);
        }
    }
    /* swap buffers */
    float *tmp = F_food; F_food = F_temp; F_temp = tmp;
}

/* ── Time step ──────────────────────────────────────────────────── */

void evoca_step(void)
{
    int    N     = gN;
    size_t cells = (size_t)N * N;
    g_step++;

    /* Phase 1: CA state update (double-buffered) */
    memset(lut_active, 0, LUT_BYTES);
    for (int row = 0; row < N; row++) {
        for (int col = 0; col < N; col++) {
            int idx = row * N + col;
            int bit = compute_lut_bit(row, col);
            if (alive[idx])
                lut_active[bit >> 3] |= (uint8_t)(1u << (bit & 7));
            v_next[idx] = lut_get(lut + (size_t)idx * LUT_BYTES, bit);
        }
    }
    uint8_t *tmp = v_curr; v_curr = v_next; v_next = tmp;

    /* Build active-bit index array for restricted mutation */
    if (grestricted_mu) {
        n_active = 0;
        for (int i = 0; i < LUT_BITS; i++)
            if ((lut_active[i >> 3] >> (i & 7)) & 1)
                active_bits[n_active++] = i;
    }

    /* Phase 2: Environmental food regeneration (respects env_mask) */
    for (size_t i = 0; i < cells; i++) {
        if (!env_mask[i]) continue;
        F_food[i] += gfood_inc;
        if (F_food[i] > 1.0f) F_food[i] = 1.0f;
    }

    /* Phase 2b: Food diffusion (gdiff passes of 3×3 box blur) */
    for (int d = 0; d < ggdiff; d++)
        diffuse_food_once();

    /* Phase 2c: Tax — decrement private food; death if depleted */
    if (gtax > 0.0f) {
        for (size_t i = 0; i < cells; i++) {
            if (!alive[i]) continue;
            f_priv[i] -= gtax;
            if (f_priv[i] <= 0.0f) {
                f_priv[i] = 0.0f;
                alive[i] = 0;
                /* Death: zero out LUT genome */
                memset(lut + i * LUT_BYTES, 0, LUT_BYTES);
                uint32_t dh = lut_hash_fn(lut + i * LUT_BYTES);
                lut_color[i] = (uint32_t)hash_to_color(dh, wt_hash);
                lut_hash_cache[i] = dh;
            }
        }
    }

    /* Phase 3: Eating (alive cells only) */
    for (int row = 0; row < N; row++) {
        for (int col = 0; col < N; col++) {
            int   idx      = row * N + col;
            if (!alive[idx]) continue;
            int   matches  = fiducial_matches(row, col, egenome[idx]);
            float mouthful = (gm_scale / 25.0f) * matches * F_food[idx];
            float headroom = 1.0f - f_priv[idx];
            if (mouthful > headroom) mouthful = headroom;
            F_food[idx] -= mouthful;
            f_priv[idx] += mouthful;
        }
    }

    /* Phase 4: Reproduction (alive parents only) */
    memset(births, 0, cells);
    for (int row = 0; row < N; row++) {
        for (int col = 0; col < N; col++) {
            int idx = row * N + col;
            if (!alive[idx]) continue;
            if (f_priv[idx] < 1.0f) continue;

            /* Reservoir-sample the neighbour with minimum f_priv.
             * Tie-breaking is uniformly random to avoid directional drift. */
            int   best_r = -1, best_c = -1;
            float best_f = 1e30f;
            int   n_best = 0;
            for (int di = -1; di <= 1; di++) {
                int r = ((row + di) % N + N) % N;
                for (int dj = -1; dj <= 1; dj++) {
                    if (di == 0 && dj == 0) continue;
                    int c = ((col + dj) % N + N) % N;
                    float f = f_priv[r * N + c];
                    if (f < best_f - 1e-7f) {
                        best_f = f; best_r = r; best_c = c; n_best = 1;
                    } else if (f < best_f + 1e-7f) {
                        if ((rng_next() % (uint32_t)++n_best) == 0) {
                            best_r = r; best_c = c;
                        }
                    }
                }
            }
            int child = best_r * N + best_c;
            alive[child] = 1;

            /* Record reproduction age and reset timestamps */
            last_event_step[child] = g_step;  /* child just born */
            if (g_step >= repro_age_t0 && last_event_step[idx] >= repro_age_t0) {
                uint32_t age = g_step - last_event_step[idx];
                if (age >= REPRO_AGE_MAX) age = REPRO_AGE_MAX - 1;
                repro_age_hist[age]++;
            }
            last_event_step[idx] = g_step;    /* parent just reproduced */

            memcpy(lut + (size_t)child * LUT_BYTES,
                   lut + (size_t)idx   * LUT_BYTES, LUT_BYTES);
            egenome[child] = egenome[idx];
            /* Mutate child's LUT */
            uint8_t *child_lut = lut + (size_t)child * LUT_BYTES;
            int nf;
            if (grestricted_mu && n_active > 0) {
                nf = poisson_sample(gmu_lut * n_active);
                for (int f = 0; f < nf; f++)
                    lut_flip(child_lut, active_bits[rng_next() % n_active]);
            } else {
                nf = poisson_sample(gmu_lut * LUT_BITS);
                for (int f = 0; f < nf; f++)
                    lut_flip(child_lut, rng_next() % LUT_BITS);
            }

            /* Mutate child's egenome */
            int nc = poisson_sample(gmu_egenome * 6);
            for (int f = 0; f < nc; f++)
                egenome[child] ^= (uint8_t)(1u << (rng_next() % 6));

            /* births: 1 = normal, 2 = mutant */
            births[child] = (nf > 0 || nc > 0) ? 2 : 1;

            /* Update cached genome color + hash (skip if unchanged) */
            if (nf > 0) {
                uint32_t ch = lut_hash_fn(child_lut);
                lut_color[child] = (uint32_t)hash_to_color(ch, wt_hash);
                lut_hash_cache[child] = ch;
            } else {
                lut_color[child] = lut_color[idx];
                lut_hash_cache[child] = lut_hash_cache[idx];
            }

            /* v_curr is dynamical state, not genome — do not copy */
            float half    = f_priv[idx] * 0.5f;
            f_priv[idx]   = half;
            f_priv[child] = half;
        }
    }
}

/* ── Activity tracking ─────────────────────────────────────────── */

void evoca_activity_update(void)
{
    if (!act_keys) return;
    size_t cells = (size_t)gN * gN;

    /* Clear all pop_counts */
    for (int i = 0; i < act_cap; i++)
        if (act_keys[i] != ACT_EMPTY)
            act_vals[i].pop_count = 0;

    /* Tally alive cells */
    for (size_t i = 0; i < cells; i++) {
        if (!alive[i]) continue;
        act_entry_t *e = act_find_or_insert(lut_hash_cache[i],
                                            (int32_t)lut_color[i]);
        e->pop_count++;
        e->activity++;
    }

    /* Compact when table gets large — prune extinct low-activity entries.
     * Threshold: keep extinct genomes with activity >= act_ymax/10, so they
     * remain visible on the chart.  Trigger at 50k entries. */
    if (act_cnt > 50000)
        act_compact((uint64_t)act_ymax / 10);

    /* Record per-bucket pop history for the flux probe. */
    int slot = (int)(flux_t % (uint32_t)FLUX_W_MAX);
    for (int i = 0; i < act_cap; i++) {
        act_pop_hist[(size_t)i * FLUX_W_MAX + slot] =
            (act_keys[i] == ACT_EMPTY) ? 0u : act_vals[i].pop_count;
    }
    flux_t++;
}

/* Activity flux probe: for each bucket, walk the past `window` ticks and
 * reconstruct the activity trajectory from stored pop history. A tick
 * contributes one slope sample = pop iff the activity increment interval
 * (A_prev, A_now] overlaps [a_lo, a_hi]. Writes up to max_n float slopes
 * and returns the count written. */
int evoca_activity_crossings(uint64_t a_lo, uint64_t a_hi, int window,
                             float *slopes_out, int max_n)
{
    if (!act_keys || window <= 0 || max_n <= 0) return 0;
    if (window > FLUX_W_MAX) window = FLUX_W_MAX;
    if ((uint32_t)window > flux_t) window = (int)flux_t;
    if (window <= 0) return 0;

    int n = 0;
    for (int i = 0; i < act_cap; i++) {
        if (act_keys[i] == ACT_EMPTY) continue;
        uint64_t A_after = act_vals[i].activity;
        size_t base = (size_t)i * FLUX_W_MAX;
        for (int k = 0; k < window; k++) {
            uint32_t slot = (flux_t - 1u - (uint32_t)k)
                            % (uint32_t)FLUX_W_MAX;
            uint32_t h = act_pop_hist[base + slot];
            uint64_t A_before = (A_after >= h) ? A_after - h : 0;
            if (h > 0 && A_before < a_hi && A_after >= a_lo) {
                if (n >= max_n) return n;
                slopes_out[n++] = (float)h;
            }
            if (A_before < a_lo) break;
            A_after = A_before;
        }
    }
    return n;
}

void evoca_activity_render_col(int32_t *col, int height)
{
    /* Fill with background */
    for (int y = 0; y < height; y++)
        col[y] = (int32_t)0xFF111111u;

    if (!act_keys || act_cnt == 0) return;

    /* Saturation formula: y = (H-1) - (H-1)*act/(act+ymax).
     * Waves rise toward the top as activity >> ymax. */
    uint64_t ymax = (uint64_t)act_ymax;

    /* Per-pixel pop tracking for priority (higher pop overwrites) */
    uint32_t ypop[height];
    memset(ypop, 0, (size_t)height * sizeof(uint32_t));

    /* Pass 1: extinct genomes (pop_count == 0) — dimmed color */
    for (int i = 0; i < act_cap; i++) {
        if (act_keys[i] == ACT_EMPTY) continue;
        if (act_vals[i].pop_count > 0) continue;
        uint64_t act = act_vals[i].activity;
        int y = (height - 1) - (int)((uint64_t)(height - 1) * act / (act + ymax));
        if (y < 0) y = 0;
        if (y >= height) y = height - 1;
        /* Dim: RGB × 0.15 */
        uint32_t c = (uint32_t)act_vals[i].color;
        uint8_t r = (uint8_t)(((c >> 16) & 0xFF) * 15 / 100);
        uint8_t g = (uint8_t)(((c >>  8) & 0xFF) * 15 / 100);
        uint8_t b = (uint8_t)(( c        & 0xFF) * 15 / 100);
        col[y] = (int32_t)(0xFF000000u | ((uint32_t)r << 16)
                           | ((uint32_t)g << 8) | b);
    }

    /* Pass 2: alive genomes — full color, higher pop wins */
    for (int i = 0; i < act_cap; i++) {
        if (act_keys[i] == ACT_EMPTY) continue;
        uint32_t pop = act_vals[i].pop_count;
        if (pop == 0) continue;
        uint64_t act = act_vals[i].activity;
        int y = (height - 1) - (int)((uint64_t)(height - 1) * act / (act + ymax));
        if (y < 0) y = 0;
        if (y >= height) y = height - 1;
        if (pop >= ypop[y]) {
            col[y] = act_vals[i].color;
            ypop[y] = pop;
        }
    }
}

int evoca_count_distinct_genomes(void)
{
    if (!act_keys) return 0;
    int n = 0;
    for (int i = 0; i < act_cap; i++) {
        if (act_keys[i] == ACT_EMPTY) continue;
        if (act_vals[i].pop_count > 0) n++;
    }
    return n;
}

int evoca_get_population(void)
{
    if (!alive) return 0;
    size_t cells = (size_t)gN * gN;
    int pop = 0;
    for (size_t i = 0; i < cells; i++) pop += alive[i];
    return pop;
}

void evoca_get_ages(int32_t *out)
{
    size_t cells = (size_t)gN * gN;
    if (!alive || !last_event_step) {
        for (size_t i = 0; i < cells; i++) out[i] = -1;
        return;
    }
    for (size_t i = 0; i < cells; i++) {
        if (!alive[i]) { out[i] = -1; continue; }
        uint32_t t0 = last_event_step[i];
        out[i] = (g_step >= t0) ? (int32_t)(g_step - t0) : 0;
    }
}

int evoca_activity_get(uint32_t *keys, uint64_t *activities,
                       uint32_t *pop_counts, int32_t *colors, int max_n)
{
    int n = 0;
    if (!act_keys) return 0;
    for (int i = 0; i < act_cap && n < max_n; i++) {
        if (act_keys[i] == ACT_EMPTY) continue;
        keys[n]       = act_keys[i];
        activities[n] = act_vals[i].activity;
        pop_counts[n] = act_vals[i].pop_count;
        colors[n]     = act_vals[i].color;
        n++;
    }
    return n;
}

/* ── Activity quantile (q_activity) probe ──────────────────────── */

static int _float_cmp(const void *a, const void *b)
{
    float fa = *(const float *)a, fb = *(const float *)b;
    return (fa > fb) - (fa < fb);
}

void evoca_q_activity_deciles(float *deciles_out)
{
    /* Collect normalized activities of alive genomes (pop_count > 0) */
    if (!act_keys || act_cnt == 0) {
        for (int i = 0; i < 9; i++) deciles_out[i] = 0.0f;
        return;
    }

    /* Count alive distinct genomes (diversity D) */
    int D = 0;
    for (int i = 0; i < act_cap; i++)
        if (act_keys[i] != ACT_EMPTY && act_vals[i].pop_count > 0)
            D++;

    if (D == 0) {
        for (int i = 0; i < 9; i++) deciles_out[i] = 0.0f;
        return;
    }

    /* Allocate temp buffer, collect and normalize */
    float *buf = (float *)malloc((size_t)D * sizeof(float));
    if (!buf) {
        for (int i = 0; i < 9; i++) deciles_out[i] = 0.0f;
        return;
    }
    int n = 0;
    float inv_D = 1.0f / (float)D;
    for (int i = 0; i < act_cap; i++) {
        if (act_keys[i] == ACT_EMPTY || act_vals[i].pop_count == 0) continue;
        buf[n++] = (float)act_vals[i].activity * inv_D;
    }

    /* Sort ascending */
    qsort(buf, (size_t)n, sizeof(float), _float_cmp);

    /* Extract deciles: p10, p20, ..., p90 */
    for (int d = 0; d < 9; d++) {
        int idx = (int)((float)(d + 1) * 0.1f * (float)n);
        if (idx >= n) idx = n - 1;
        deciles_out[d] = buf[idx];
    }
    free(buf);
}

/* ── Egenome activity tracking ──────────────────────────────────── */

void evoca_eg_activity_update(void)
{
    size_t cells = (size_t)gN * gN;
    memset(eg_pop, 0, sizeof(eg_pop));
    for (size_t i = 0; i < cells; i++) {
        if (!alive[i]) continue;
        uint8_t eg = egenome[i] & 0x3F;
        eg_pop[eg]++;
        eg_act[eg]++;
    }
}

void evoca_eg_activity_render_col(int32_t *col, int height)
{
    for (int y = 0; y < height; y++)
        col[y] = (int32_t)0xFF111111u;

    uint64_t ymax = (uint64_t)eg_act_ymax;

    uint32_t ypop[height];
    memset(ypop, 0, (size_t)height * sizeof(uint32_t));

    /* Pass 1: extinct egenomes — dimmed */
    for (int i = 0; i < EGENOME_COUNT; i++) {
        if (eg_act[i] == 0 || eg_pop[i] > 0) continue;
        uint64_t act = eg_act[i];
        int y = (height - 1) - (int)((uint64_t)(height - 1) * act / (act + ymax));
        if (y < 0) y = 0;
        if (y >= height) y = height - 1;
        uint32_t c = (uint32_t)eg_color[i];
        uint8_t r = (uint8_t)(((c >> 16) & 0xFF) * 15 / 100);
        uint8_t g = (uint8_t)(((c >>  8) & 0xFF) * 15 / 100);
        uint8_t b = (uint8_t)(( c        & 0xFF) * 15 / 100);
        col[y] = (int32_t)(0xFF000000u | ((uint32_t)r << 16)
                           | ((uint32_t)g << 8) | b);
    }

    /* Pass 2: alive egenomes — full color, higher pop wins */
    for (int i = 0; i < EGENOME_COUNT; i++) {
        if (eg_pop[i] == 0) continue;
        uint64_t act = eg_act[i];
        int y = (height - 1) - (int)((uint64_t)(height - 1) * act / (act + ymax));
        if (y < 0) y = 0;
        if (y >= height) y = height - 1;
        if (eg_pop[i] >= ypop[y]) {
            col[y] = eg_color[i];
            ypop[y] = eg_pop[i];
        }
    }
}

int evoca_eg_activity_get(uint64_t *activities, uint32_t *pop_counts,
                          int32_t *colors)
{
    for (int i = 0; i < EGENOME_COUNT; i++) {
        activities[i] = eg_act[i];
        pop_counts[i] = eg_pop[i];
        colors[i]     = eg_color[i];
    }
    return EGENOME_COUNT;
}

void evoca_set_eg_act_ymax(int y) { eg_act_ymax = y > 1 ? y : 1; }
int  evoca_get_eg_act_ymax(void)  { return eg_act_ymax; }

/* ── LUT complexity classification ─────────────────────────────── */

/* Return 1 (n1 only), 2 (n1+n2), or 3 (n1+n2+n3) for a given LUT. */
static int lut_complexity(const uint8_t *b)
{
    /* Check n3-dependence: for each (v_x,n1,n2), all 5 n3 entries equal? */
    for (int vx = 0; vx < 2; vx++)
        for (int n1 = 0; n1 < 5; n1++)
            for (int n2 = 0; n2 < 5; n2++) {
                int base = LUT_IDX(vx, n1, n2, 0);
                uint8_t first = lut_get(b, base);
                for (int n3 = 1; n3 < 5; n3++)
                    if (lut_get(b, base + n3) != first)
                        return 3;
            }

    /* n3-independent. Check n2-dependence. */
    for (int vx = 0; vx < 2; vx++)
        for (int n1 = 0; n1 < 5; n1++) {
            uint8_t first = lut_get(b, LUT_IDX(vx, n1, 0, 0));
            for (int n2 = 1; n2 < 5; n2++)
                if (lut_get(b, LUT_IDX(vx, n1, n2, 0)) != first)
                    return 2;
        }

    return 1;
}

void evoca_lut_complexity_counts(uint32_t *counts)
{
    counts[0] = counts[1] = counts[2] = 0;
    size_t cells = (size_t)gN * gN;
    for (size_t i = 0; i < cells; i++) {
        if (!alive[i]) continue;
        int lvl = lut_complexity(lut + i * LUT_BYTES);
        counts[lvl - 1]++;
    }
}

void evoca_lut_complexity_render_col(int32_t *col, int height)
{
    uint32_t counts[3];
    evoca_lut_complexity_counts(counts);
    uint32_t total = counts[0] + counts[1] + counts[2];

    /* Fill background */
    for (int y = 0; y < height; y++)
        col[y] = (int32_t)0xFF111111u;

    if (total == 0) return;

    /* Stacked from bottom: green (n1), yellow (n1+n2), red (full) */
    int h1 = (int)((uint64_t)counts[0] * height / total);
    int h2 = (int)((uint64_t)counts[1] * height / total);
    int h3 = height - h1 - h2;

    int y = height - 1;
    for (int i = 0; i < h1 && y >= 0; i++, y--)
        col[y] = (int32_t)0xFF00CC00u;   /* green: n1 only */
    for (int i = 0; i < h2 && y >= 0; i++, y--)
        col[y] = (int32_t)0xFFCCCC00u;   /* yellow: n1+n2 */
    for (int i = 0; i < h3 && y >= 0; i++, y--)
        col[y] = (int32_t)0xFFCC3300u;   /* red: full */
}

/* ── Egenome population banded chart ──────────────────────────── */

void evoca_eg_pop_render_col(int32_t *col, int height)
{
    /* Fill background */
    for (int y = 0; y < height; y++)
        col[y] = (int32_t)0xFF111111u;

    /* Sum total alive population across all egenomes */
    uint32_t total = 0;
    for (int i = 0; i < EGENOME_COUNT; i++)
        total += eg_pop[i];
    if (total == 0) return;

    /* Fixed centered order: odd indices descending, 0, even ascending.
       Wild-type (index 0) sits in the vertical centre; order never changes. */
    static const int order[EGENOME_COUNT] = {
        63,61,59,57,55,53,51,49,47,45,43,41,39,37,35,33,
        31,29,27,25,23,21,19,17,15,13,11, 9, 7, 5, 3, 1,
         0, 2, 4, 6, 8,10,12,14,16,18,20,22,24,26,28,30,
        32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62
    };

    /* Find last populated entry for remainder assignment */
    int last_k = -1;
    for (int k = EGENOME_COUNT - 1; k >= 0; k--)
        if (eg_pop[order[k]] > 0) { last_k = k; break; }

    /* Stack bands top to bottom in fixed order */
    int y = 0;
    int assigned = 0;
    for (int k = 0; k < EGENOME_COUNT && y < height; k++) {
        int eg = order[k];
        if (eg_pop[eg] == 0) continue;
        int h;
        if (k == last_k) {
            h = height - assigned;
        } else {
            h = (int)((uint64_t)eg_pop[eg] * height / total);
            if (h == 0) h = 1;  /* at least 1px if populated */
        }
        if (h <= 0) continue;
        for (int i = 0; i < h && y < height; i++, y++)
            col[y] = eg_color[eg];
        assigned += h;
    }
}

/* ── Visualisation ──────────────────────────────────────────────── */

void evoca_colorize(int32_t *pixels, int colormode)
{
    size_t cells = (size_t)gN * gN;
    switch (colormode) {
        case 0:
            for (size_t i = 0; i < cells; i++) {
                if (!alive[i])
                    pixels[i] = (int32_t)0xFF000000;
                else if (v_curr[i])
                    pixels[i] = (int32_t)lut_color[i];
                else
                    pixels[i] = (int32_t)0xFF333333;  /* alive, v=0: dark grey */
            }
            break;
        case 1:
            for (size_t i = 0; i < cells; i++) {
                float fv = F_food[i]; if (fv > 1.0f) fv = 1.0f;
                uint8_t g = (uint8_t)(fv * 255.0f);
                uint8_t r = alive[i] ? 180 : 0;
                pixels[i] = (int32_t)(0xFF000000u
                             | ((uint32_t)r << 16) | ((uint32_t)g << 8));
            }
            break;
        case 2:
            for (size_t i = 0; i < cells; i++) {
                if (!alive[i]) {
                    pixels[i] = (int32_t)0xFF000000;
                } else {
                    float fv = f_priv[i]; if (fv > 1.0f) fv = 1.0f;
                    uint8_t b = (uint8_t)(fv * 255.0f);
                    uint8_t r = v_curr[i] ? 180 : 0;
                    pixels[i] = (int32_t)(0xFF000000u
                                 | ((uint32_t)r << 16) | (uint32_t)b);
                }
            }
            break;
        case 3:
            for (size_t i = 0; i < cells; i++) {
                if (!alive[i]) {
                    pixels[i] = (int32_t)0xFF000000;
                } else if (births[i] == 2) {   /* mutant birth: bright magenta */
                    pixels[i] = v_curr[i] ? (int32_t)0xFFFF00FFu
                                          : (int32_t)0xFF800080u;
                } else if (births[i] == 1) {   /* normal birth: yellow */
                    pixels[i] = v_curr[i] ? (int32_t)0xFFFFFF00u
                                          : (int32_t)0xFF808000u;
                } else {
                    pixels[i] = v_curr[i] ? (int32_t)0xFF444444u
                                          : (int32_t)0xFF222222u;
                }
            }
            break;
        case 4: {
            /* Log-age gradient across the alive population. Tracks peak age
             * adaptively so hue range follows the distribution. */
            int max_age = 0;
            if (last_event_step) {
                for (size_t i = 0; i < cells; i++) {
                    if (!alive[i]) continue;
                    uint32_t t0 = last_event_step[i];
                    int a = (g_step >= t0) ? (int)(g_step - t0) : 0;
                    if (a > max_age) max_age = a;
                }
            }
            if ((double)max_age > g_age_scale) g_age_scale = (double)max_age;
            else                               g_age_scale *= 0.995;
            if (g_age_scale < 10.0) g_age_scale = 10.0;

            double denom = log1p(g_age_scale);
            for (size_t i = 0; i < cells; i++) {
                if (!alive[i] || !last_event_step) {
                    pixels[i] = (int32_t)0xFF000000;
                    continue;
                }
                uint32_t t0 = last_event_step[i];
                int a = (g_step >= t0) ? (int)(g_step - t0) : 0;
                float v = (float)(log1p((double)a) / denom);
                if (v < 0.0f) v = 0.0f;
                if (v > 1.0f) v = 1.0f;
                uint8_t R = (uint8_t)(80.0f  + 175.0f * v);
                uint8_t G = (uint8_t)(60.0f  + 560.0f * v * (1.0f - v));
                uint8_t B = (uint8_t)(40.0f  + 180.0f * (1.0f - v));
                pixels[i] = (int32_t)(0xFF000000u
                                       | ((uint32_t)R << 16)
                                       | ((uint32_t)G <<  8)
                                       |  (uint32_t)B);
            }
            break;
        }
        default:
            for (size_t i = 0; i < cells; i++)
                pixels[i] = (int32_t)0xFF000000;
    }
}

/* ── Accessors ──────────────────────────────────────────────────── */

uint8_t *evoca_get_v(void)      { return v_curr; }
float   *evoca_get_F(void)      { return F_food; }
float   *evoca_get_f(void)      { return f_priv; }
uint8_t *evoca_get_egenome(void) { return egenome; }
uint8_t *evoca_get_lut(void)    { return lut;    }
uint8_t *evoca_get_births(void) { return births; }
uint8_t *evoca_get_alive(void)  { return alive;  }
int      evoca_get_N(void)      { return gN;     }
int      evoca_get_cell_px(void) { return CELL_PX; }
void     evoca_set_act_ymax(int y) { act_ymax = y > 1 ? y : 1; }
int      evoca_get_act_ymax(void)  { return act_ymax; }
uint32_t *evoca_get_repro_age_hist(void) { return repro_age_hist; }
int      evoca_get_repro_age_max(void)   { return REPRO_AGE_MAX; }
uint32_t evoca_get_step(void)            { return g_step; }
void     evoca_set_repro_age_t0(uint32_t t) { repro_age_t0 = t; }
uint32_t evoca_get_repro_age_t0(void)       { return repro_age_t0; }
void     evoca_reset_repro_age_hist(void)   { memset(repro_age_hist, 0, sizeof(repro_age_hist)); }

/* ── Alive data plane ──────────────────────────────────────────── */

void evoca_set_alive(const uint8_t *arr)
{
    size_t cells = (size_t)gN * gN;
    memcpy(alive, arr, cells);
    /* Zero dead cells' data */
    for (size_t i = 0; i < cells; i++) {
        if (!alive[i]) {
            v_curr[i] = 0;
            f_priv[i] = 0.0f;
            memset(lut + i * LUT_BYTES, 0, LUT_BYTES);
            uint32_t dh = lut_hash_fn(lut + i * LUT_BYTES);
            lut_color[i] = (uint32_t)hash_to_color(dh, wt_hash);
            lut_hash_cache[i] = dh;
        }
    }
}

void evoca_set_alive_all(void)
{
    memset(alive, 1, (size_t)gN * gN);
}

void evoca_set_alive_fraction(float frac)
{
    size_t cells = (size_t)gN * gN;
    for (size_t i = 0; i < cells; i++)
        alive[i] = (rng_next() < (uint32_t)(frac * 4294967295.0f)) ? 1 : 0;
    /* Zero dead cells' data */
    for (size_t i = 0; i < cells; i++) {
        if (!alive[i]) {
            v_curr[i] = 0;
            f_priv[i] = 0.0f;
            memset(lut + i * LUT_BYTES, 0, LUT_BYTES);
            uint32_t dh = lut_hash_fn(lut + i * LUT_BYTES);
            lut_color[i] = (uint32_t)hash_to_color(dh, wt_hash);
            lut_hash_cache[i] = dh;
        }
    }
}

void evoca_set_alive_patch(int radius)
{
    int N = gN;
    int cx = N / 2, cy = N / 2;
    size_t cells = (size_t)N * N;
    memset(alive, 0, cells);
    for (int r = cy - radius; r < cy + radius; r++) {
        int rr = ((r % N) + N) % N;
        for (int c = cx - radius; c < cx + radius; c++) {
            int cc = ((c % N) + N) % N;
            alive[rr * N + cc] = 1;
        }
    }
    /* Zero dead cells' data */
    for (size_t i = 0; i < cells; i++) {
        if (!alive[i]) {
            v_curr[i] = 0;
            f_priv[i] = 0.0f;
            memset(lut + i * LUT_BYTES, 0, LUT_BYTES);
            uint32_t dh = lut_hash_fn(lut + i * LUT_BYTES);
            lut_color[i] = (uint32_t)hash_to_color(dh, wt_hash);
            lut_hash_cache[i] = dh;
        }
    }
}

void evoca_set_alive_halfplane(int axis)
{
    int N = gN;
    size_t cells = (size_t)N * N;
    memset(alive, 0, cells);
    for (int r = 0; r < N; r++) {
        for (int c = 0; c < N; c++) {
            int keep = (axis == 0) ? (c < N / 2) : (r < N / 2);
            if (keep) alive[r * N + c] = 1;
        }
    }
    /* Zero dead cells' data */
    for (size_t i = 0; i < cells; i++) {
        if (!alive[i]) {
            v_curr[i] = 0;
            f_priv[i] = 0.0f;
            memset(lut + i * LUT_BYTES, 0, LUT_BYTES);
            uint32_t dh = lut_hash_fn(lut + i * LUT_BYTES);
            lut_color[i] = (uint32_t)hash_to_color(dh, wt_hash);
            lut_hash_cache[i] = dh;
        }
    }
}

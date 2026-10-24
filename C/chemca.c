#include "chemca.h"
#include "core.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
 * ChemCA engine. See chemca.h for the model and Docs/plan-chemca.md for the
 * design choices. Step order (chemca_step):
 *   1. environment reactions (exponential integrator per site), then the
 *      open-mode source and sink;
 *   2. diffusion of each E_c (core_diffuse_precise);
 *   3. two-way transport for alive cells, d = u_c (E_c - a_c);
 *   4. agent reactions -> energy e;
 *   5. costs e -= m_e + c_g * popcount(active), death if they cannot be paid;
 *   6. toxicity: a private amount >= tox_c (> 0) kills;
 *   7. reproduction at e >= 1 into the Moore neighbor with the lowest e
 *      (dead = -inf), evicting an alive occupant; 50/50 split; mutation.
 *      With repro_gate = mass: body mass (u_c == 0 compounds) >= mass_repro,
 *      evicting the lowest body mass.
 * On every death the private stores return to the cell's site (mass is
 * conserved) and e is lost (booked as CHEMCA_EN_DEATH).
 */

/* ── PRNG (this model's own state; core.c holds only the step) ─────── */

#define CHEMCA_DEFAULT_SEED 0x2468ACE1u
static uint32_t cc_rng = CHEMCA_DEFAULT_SEED;
static uint32_t cc_rng_next(void) { return core_xorshift32(&cc_rng); }
/* Uniform in [0, 1) from 24 bits. */
static inline float cc_unif(void) { return (float)(cc_rng_next() & 0xFFFFFF) / 16777216.0f; }

void chemca_set_seed(uint32_t seed) { cc_rng = seed ? seed : CHEMCA_DEFAULT_SEED; }

/* ── State ─────────────────────────────────────────────────────────── */

static int      gN = 0, gK = 0;
static size_t   gcells = 0;
static long     g_step = 0;
static uint8_t *alive  = NULL;   /* [cells] */
/* Dense per-cell stores [cells * K], cell-major (store[i*K + c]). This
 * departs from the sparse slot layout of pinned decision A
 * (Docs/plan-chem.md §2), which targets universes of ~1000 compounds; for
 * K <= 8 dense is simpler and faster, and the layout is internal to this
 * file (Docs/2026-10-23_chem-model-v2.md §9b). */
static float   *store  = NULL;
static float   *e_en   = NULL;   /* [cells] energy */
static uint8_t *rg     = NULL;   /* [cells * NSLOT * GENE_BYTES] */
static uint8_t *ractive = NULL;  /* [cells] active-slot mask */
static uint8_t *births = NULL;   /* [cells] last step: 1 normal, 2 mutant */
static float   *E[CHEMCA_KMAX];  /* [cells] env field per compound */
static float   *E_tmp = NULL;    /* [cells] diffusion scratch */

/* Per-compound configuration. */
static float cG[CHEMCA_KMAX], cdiff[CHEMCA_KMAX], cu[CHEMCA_KMAX],
             ctox[CHEMCA_KMAX], cinc[CHEMCA_KMAX], csink[CHEMCA_KMAX];

/* Environment reactions. */
typedef struct { int A, B; float rate; float *field; } envrx_t;
static envrx_t envrx[CHEMCA_MAX_ENVRX];
static int     n_envrx = 0;

/* Metaparameters. */
static float g_k_min = 1e-3f, g_m_e = 0.0f, g_c_g = 0.0f;
static float g_mu_rgene = 0.0f, g_mu_ractive = 0.0f, g_p_dup = 0.5f;
static int   g_open = 0;
/* Reproduction gate: 0 energy (e >= 1), 1 mass (body mass >= g_mass_repro,
 * Kaneko style). Body compounds are those with transport rate u_c == 0. */
static int   g_repro_gate = 0;
static float g_mass_repro = 1.0f;
static float g_color_scale = 1.0f;
static double ktab[CHEMCA_NK];   /* k_min * 2^(k_idx/4) */

/* Bookkeeping. */
static double en_last[CHEMCA_EN_COUNT], en_cum[CHEMCA_EN_COUNT];
static double flux_in[CHEMCA_KMAX], flux_out[CHEMCA_KMAX];   /* cumulative */
static long n_births = 0, n_deaths = 0, n_tox = 0, n_edeath = 0, n_evict = 0;

#define GENE(i, s) (rg + ((size_t)(i) * CHEMCA_NSLOT + (s)) * CHEMCA_GENE_BYTES)

static void rebuild_ktab(void)
{
    for (int k = 0; k < CHEMCA_NK; k++)
        ktab[k] = (double)g_k_min * pow(2.0, k / 4.0);
}

/* ── Lifecycle ─────────────────────────────────────────────────────── */

void chemca_free(void)
{
    free(alive); free(store); free(e_en); free(rg); free(ractive); free(births);
    alive = NULL; store = NULL; e_en = NULL; rg = NULL; ractive = NULL; births = NULL;
    for (int c = 0; c < CHEMCA_KMAX; c++) { free(E[c]); E[c] = NULL; }
    free(E_tmp); E_tmp = NULL;
    chemca_env_rx_clear();
    gN = gK = 0; gcells = 0;
}

int chemca_init(int N, int K)
{
    if (N < 3 || K < 2 || K > CHEMCA_KMAX) return -1;
    chemca_free();
    gN = N; gK = K; gcells = (size_t)N * N;
    alive   = calloc(gcells, 1);
    store   = calloc(gcells * K, sizeof(float));
    e_en    = calloc(gcells, sizeof(float));
    rg      = calloc(gcells * CHEMCA_NSLOT * CHEMCA_GENE_BYTES, 1);
    ractive = calloc(gcells, 1);
    births  = calloc(gcells, 1);
    for (int c = 0; c < K; c++) E[c] = calloc(gcells, sizeof(float));
    E_tmp = calloc(gcells, sizeof(float));
    for (int c = 0; c < CHEMCA_KMAX; c++) {
        cG[c] = cdiff[c] = cu[c] = ctox[c] = cinc[c] = csink[c] = 0.0f;
        flux_in[c] = flux_out[c] = 0.0;
    }
    memset(en_last, 0, sizeof en_last);
    memset(en_cum, 0, sizeof en_cum);
    n_births = n_deaths = n_tox = n_edeath = n_evict = 0;
    g_step = 0;
    rebuild_ktab();
    return 0;
}

int  chemca_get_N(void)       { return gN; }
int  chemca_get_K(void)       { return gK; }
long chemca_get_step(void)    { return g_step; }
int  chemca_kmax(void)        { return CHEMCA_KMAX; }
int  chemca_nslot(void)       { return CHEMCA_NSLOT; }
int  chemca_gene_bytes(void)  { return CHEMCA_GENE_BYTES; }
int  chemca_nk(void)          { return CHEMCA_NK; }
int  chemca_energy_count(void){ return CHEMCA_EN_COUNT; }

/* ── Metaparameters ────────────────────────────────────────────────── */

static float nonneg(float v) { return v > 0.0f ? v : 0.0f; }
static float unit(float v)   { return v < 0.0f ? 0.0f : (v > 1.0f ? 1.0f : v); }

void  chemca_set_k_min(float v)      { g_k_min = nonneg(v); rebuild_ktab(); }
float chemca_get_k_min(void)         { return g_k_min; }
void  chemca_set_m_e(float v)        { g_m_e = nonneg(v); }
float chemca_get_m_e(void)           { return g_m_e; }
void  chemca_set_c_g(float v)        { g_c_g = nonneg(v); }
float chemca_get_c_g(void)           { return g_c_g; }
void  chemca_set_mu_rgene(float v)   { g_mu_rgene = nonneg(v); }
float chemca_get_mu_rgene(void)      { return g_mu_rgene; }
void  chemca_set_mu_ractive(float v) { g_mu_ractive = nonneg(v); }
float chemca_get_mu_ractive(void)    { return g_mu_ractive; }
void  chemca_set_p_dup(float v)      { g_p_dup = unit(v); }
float chemca_get_p_dup(void)         { return g_p_dup; }
void  chemca_set_open_mode(int on)   { g_open = on ? 1 : 0; }
int   chemca_get_open_mode(void)     { return g_open; }
void  chemca_set_repro_gate(int g)  { g_repro_gate = g ? 1 : 0; }
int   chemca_get_repro_gate(void)    { return g_repro_gate; }
void  chemca_set_mass_repro(float v) { g_mass_repro = v > 0.0f ? v : 1.0f; }
float chemca_get_mass_repro(void)    { return g_mass_repro; }
void  chemca_set_color_scale(float s){ g_color_scale = s > 0.0f ? s : 1.0f; }
float chemca_get_color_scale(void)   { return g_color_scale; }
float chemca_k_of(int k_idx)
{
    if (k_idx < 0 || k_idx >= CHEMCA_NK) return 0.0f;
    return (float)ktab[k_idx];
}

int chemca_set_compound(int c, float G, float diff, float u, float tox,
                        float inc, float sink)
{
    if (c < 0 || c >= gK) return -1;
    cG[c] = G; cdiff[c] = nonneg(diff); cu[c] = unit(u);
    ctox[c] = nonneg(tox); cinc[c] = nonneg(inc); csink[c] = nonneg(sink);
    return 0;
}

int chemca_get_compound(int c, float *out)
{
    if (c < 0 || c >= gK) return -1;
    out[0] = cG[c]; out[1] = cdiff[c]; out[2] = cu[c];
    out[3] = ctox[c]; out[4] = cinc[c]; out[5] = csink[c];
    return 0;
}

void chemca_env_rx_clear(void)
{
    for (int r = 0; r < n_envrx; r++) { free(envrx[r].field); envrx[r].field = NULL; }
    n_envrx = 0;
}

int chemca_env_rx_add(int A, int B, float rate)
{
    if (n_envrx >= CHEMCA_MAX_ENVRX || A < 0 || A >= gK || B < 0 || B >= gK || A == B)
        return -1;
    envrx[n_envrx] = (envrx_t){ A, B, nonneg(rate), NULL };
    return n_envrx++;
}

int chemca_env_rx_set_rate(int idx, float rate)
{
    if (idx < 0 || idx >= n_envrx) return -1;
    envrx[idx].rate = nonneg(rate);
    return 0;
}

int chemca_env_rx_set_field(int idx, const float *r)
{
    if (idx < 0 || idx >= n_envrx) return -1;
    if (!r) { free(envrx[idx].field); envrx[idx].field = NULL; return 0; }
    if (!envrx[idx].field) envrx[idx].field = malloc(gcells * sizeof(float));
    for (size_t i = 0; i < gcells; i++) envrx[idx].field[i] = nonneg(r[i]);
    return 0;
}

int chemca_env_rx_count(void) { return n_envrx; }

int chemca_env_rx_get(int idx, int *A, int *B, float *rate, int *has_field)
{
    if (idx < 0 || idx >= n_envrx) return -1;
    *A = envrx[idx].A; *B = envrx[idx].B; *rate = envrx[idx].rate;
    *has_field = envrx[idx].field != NULL;
    return 0;
}

/* ── Genes ─────────────────────────────────────────────────────────── */

static int gene_valid(const uint8_t *g)
{
    int A = g[CHEMCA_G_A], B = g[CHEMCA_G_B], C = g[CHEMCA_G_C], m = g[CHEMCA_G_M];
    if (A >= gK || B >= gK || A == B) return 0;
    if (g[CHEMCA_G_K] >= CHEMCA_NK || m > 2) return 0;
    if (C == CHEMCA_C_NONE) return m == 0;
    return C < gK && m >= 1;
}

/* (C, m) pairs: index 0 = (none, 0), then (c, 1), (c, 2) for each c. */
static int  cm_count(void) { return 1 + 2 * gK; }
static void cm_set(uint8_t *g, int j)
{
    if (j == 0) { g[CHEMCA_G_C] = CHEMCA_C_NONE; g[CHEMCA_G_M] = 0; return; }
    g[CHEMCA_G_C] = (uint8_t)((j - 1) / 2);
    g[CHEMCA_G_M] = (uint8_t)(1 + (j - 1) % 2);
}
static int cm_index(const uint8_t *g)
{
    if (g[CHEMCA_G_C] == CHEMCA_C_NONE) return 0;
    return 1 + 2 * g[CHEMCA_G_C] + (g[CHEMCA_G_M] - 1);
}

static void gene_random(uint8_t *g)
{
    g[CHEMCA_G_A] = (uint8_t)(cc_rng_next() % (uint32_t)gK);
    do g[CHEMCA_G_B] = (uint8_t)(cc_rng_next() % (uint32_t)gK);
    while (g[CHEMCA_G_B] == g[CHEMCA_G_A]);
    cm_set(g, (int)(cc_rng_next() % (uint32_t)cm_count()));
    g[CHEMCA_G_K] = (uint8_t)(cc_rng_next() % CHEMCA_NK);
}

static void clear_cell(size_t i)
{
    for (int c = 0; c < gK; c++) store[i * gK + c] = 0.0f;
    e_en[i] = 0.0f;
    memset(GENE(i, 0), 0, CHEMCA_NSLOT * CHEMCA_GENE_BYTES);
    ractive[i] = 0;
}

/* ── State setters ─────────────────────────────────────────────────── */

void chemca_set_alive(const uint8_t *a)
{
    for (size_t i = 0; i < gcells; i++) {
        alive[i] = a[i] ? 1 : 0;
        if (!alive[i]) clear_cell(i);
    }
}

void chemca_set_E(int c, const float *v)
{
    if (c < 0 || c >= gK) return;
    for (size_t i = 0; i < gcells; i++) E[c][i] = nonneg(v[i]);
}

void chemca_set_store(int c, const float *s)
{
    if (c < 0 || c >= gK) return;
    for (size_t i = 0; i < gcells; i++) store[i * gK + c] = alive[i] ? nonneg(s[i]) : 0.0f;
}

void chemca_set_e(const float *v)
{
    for (size_t i = 0; i < gcells; i++) e_en[i] = alive[i] ? nonneg(v[i]) : 0.0f;
}

int chemca_set_genes_all(const uint8_t *genes, uint8_t active)
{
    for (int s = 0; s < CHEMCA_NSLOT; s++)
        if (!gene_valid(genes + s * CHEMCA_GENE_BYTES)) return -1;
    if (!active) return -1;
    for (size_t i = 0; i < gcells; i++) {
        if (!alive[i]) continue;
        memcpy(GENE(i, 0), genes, CHEMCA_NSLOT * CHEMCA_GENE_BYTES);
        ractive[i] = active;
    }
    return 0;
}

int chemca_set_genes(const uint8_t *genes, const uint8_t *active)
{
    for (size_t i = 0; i < gcells; i++) {
        if (!alive[i]) continue;
        if (!active[i]) return -1;
        for (int s = 0; s < CHEMCA_NSLOT; s++)
            if (!gene_valid(genes + (i * CHEMCA_NSLOT + s) * CHEMCA_GENE_BYTES)) return -1;
    }
    for (size_t i = 0; i < gcells; i++) {
        if (!alive[i]) continue;
        memcpy(GENE(i, 0), genes + i * CHEMCA_NSLOT * CHEMCA_GENE_BYTES,
               CHEMCA_NSLOT * CHEMCA_GENE_BYTES);
        ractive[i] = active[i];
    }
    return 0;
}

void chemca_set_genes_random(void)
{
    for (size_t i = 0; i < gcells; i++) {
        if (!alive[i]) continue;
        for (int s = 0; s < CHEMCA_NSLOT; s++) gene_random(GENE(i, s));
        ractive[i] = (uint8_t)(1u << (cc_rng_next() & 7));
    }
}

/* ── Death ─────────────────────────────────────────────────────────── */

/* Every death: stores return to the site (mass conserved), e is lost. */
static void kill_cell(size_t i)
{
    for (int c = 0; c < gK; c++) E[c][i] += store[i * gK + c];
    en_last[CHEMCA_EN_DEATH] += (double)e_en[i];
    clear_cell(i);
    alive[i] = 0;
    n_deaths++;
}

/* ── Step phases ───────────────────────────────────────────────────── */

/* 1. Environment reactions, then open-mode source/sink. */
static void phase_env(void)
{
    for (int r = 0; r < n_envrx; r++) {
        envrx_t *x = &envrx[r];
        float *EA = E[x->A], *EB = E[x->B];
        double dG = (double)cG[x->B] - (double)cG[x->A], moved = 0.0;
        double f_scalar = -expm1(-(double)x->rate);
        for (size_t i = 0; i < gcells; i++) {
            double f = x->field ? -expm1(-(double)x->field[i]) : f_scalar;
            float conv = (float)((double)EA[i] * f);
            if (conv > EA[i]) conv = EA[i];
            EA[i] -= conv; EB[i] += conv;
            moved += conv;
        }
        en_last[CHEMCA_EN_LIGHT] += moved * dG;
    }
    if (!g_open) return;
    for (int c = 0; c < gK; c++) {
        if (cinc[c] <= 0.0f && csink[c] <= 0.0f) continue;
        double fs = -expm1(-(double)csink[c]), in = 0.0, out = 0.0;
        for (size_t i = 0; i < gcells; i++) {
            E[c][i] += cinc[c]; in += cinc[c];
            if (csink[c] > 0.0f) {
                float rm = (float)((double)E[c][i] * fs);
                if (rm > E[c][i]) rm = E[c][i];
                E[c][i] -= rm; out += rm;
            }
        }
        flux_in[c] += in; flux_out[c] += out;
        en_last[CHEMCA_EN_OPEN_IN]  += in  * cG[c];
        en_last[CHEMCA_EN_OPEN_OUT] += out * cG[c];
    }
}

/* 3. Two-way transport between each alive cell and its own site. */
static void phase_transport(void)
{
    for (int c = 0; c < gK; c++) {
        float u = cu[c];
        if (u <= 0.0f) continue;
        float *Ec = E[c];
        for (size_t i = 0; i < gcells; i++) {
            if (!alive[i]) continue;
            float *a = &store[i * gK + c];
            float d = u * (Ec[i] - *a);          /* > 0 uptake, < 0 export */
            if (d > Ec[i]) d = Ec[i];
            if (-d > *a)   d = -*a;
            *a += d; Ec[i] -= d;
            if (*a < 0.0f) { Ec[i] += *a; *a = 0.0f; }   /* rounding guard */
            if (Ec[i] < 0.0f) { *a += Ec[i]; Ec[i] = 0.0f; }
        }
    }
}

/* 4. Agent reactions. Rates come from the stores at the start of the
 * phase (a snapshot), so the result does not depend on slot order:
 * for each substrate A, lambda = sum_j k_j c_j^m_j over the active genes
 * using A, converted = a_A (1 - exp(-lambda)), split in proportion
 * lambda_j / lambda. Downhill genes (G_A >= G_B) apply first; uphill genes
 * then run in slot order while e covers their cost (clipped so e >= 0). */
static void phase_react(void)
{
    double yield = 0.0;
    int K = gK;
    for (size_t i = 0; i < gcells; i++) {
        if (!alive[i] || !ractive[i]) continue;
        double a0[CHEMCA_KMAX], a[CHEMCA_KMAX], lam_sum[CHEMCA_KMAX];
        double lam[CHEMCA_NSLOT];
        for (int c = 0; c < K; c++) { a0[c] = a[c] = store[i * K + c]; lam_sum[c] = 0.0; }
        uint8_t act = ractive[i];
        for (int s = 0; s < CHEMCA_NSLOT; s++) {
            lam[s] = 0.0;
            if (!(act >> s & 1)) continue;
            const uint8_t *g = GENE(i, s);
            double l = ktab[g[CHEMCA_G_K]];
            if (g[CHEMCA_G_C] != CHEMCA_C_NONE) {
                double cc = a0[g[CHEMCA_G_C]];
                l *= (g[CHEMCA_G_M] == 2) ? cc * cc : cc;
            }
            lam[s] = l;
            lam_sum[g[CHEMCA_G_A]] += l;
        }
        double conv_A[CHEMCA_KMAX];
        for (int c = 0; c < K; c++)
            conv_A[c] = lam_sum[c] > 0.0 ? a0[c] * -expm1(-lam_sum[c]) : 0.0;
        double e = e_en[i];
        for (int pass = 0; pass < 2; pass++) {       /* 0 downhill, 1 uphill */
            for (int s = 0; s < CHEMCA_NSLOT; s++) {
                if (lam[s] <= 0.0) continue;
                const uint8_t *g = GENE(i, s);
                int A = g[CHEMCA_G_A], B = g[CHEMCA_G_B];
                double dG = (double)cG[A] - (double)cG[B];
                if ((pass == 0) != (dG >= 0.0)) continue;
                double conv = conv_A[A] * (lam[s] / lam_sum[A]);
                if (dG < 0.0 && conv * -dG > e) conv = e / -dG;
                a[A] -= conv; a[B] += conv;
                e += conv * dG;
            }
        }
        for (int c = 0; c < K; c++) store[i * K + c] = a[c] > 0.0 ? (float)a[c] : 0.0f;
        /* Book the realized change of the float e, so the e identity holds
         * to double summation; the float rounding of e then shows only in
         * the compound free-energy identity, at float precision. */
        float e_new = e > 0.0 ? (float)e : 0.0f;
        yield += (double)e_new - (double)e_en[i];
        e_en[i] = e_new;
    }
    en_last[CHEMCA_EN_YIELD] += yield;
}

/* 5. Costs (only if m_e or c_g > 0): a cell that cannot pay dies; what it
 * had left counts as spent, so it dies with e = 0. 6. Toxicity. */
static void phase_costs_tox(void)
{
    int costs = (g_m_e > 0.0f || g_c_g > 0.0f);
    for (size_t i = 0; i < gcells; i++) {
        if (!alive[i]) continue;
        if (costs) {
            float cost = g_m_e + g_c_g * (float)__builtin_popcount(ractive[i]);
            if (e_en[i] >= cost) {
                float e_new = e_en[i] - cost;     /* realized, as for yield */
                en_last[CHEMCA_EN_COST] += (double)e_en[i] - (double)e_new;
                e_en[i] = e_new;
            } else {
                en_last[CHEMCA_EN_COST] += e_en[i];
                e_en[i] = 0.0f;
                kill_cell(i);
                n_edeath++;
                continue;
            }
        }
        for (int c = 0; c < gK; c++) {
            if (ctox[c] > 0.0f && store[i * gK + c] >= ctox[c]) {
                kill_cell(i);
                n_tox++;
                break;
            }
        }
    }
}

/* Mutation of a freshly copied child genome. Returns 1 if the genome
 * changed. Active-bit pass (rejecting an empty mask), duplication onto
 * newly active slots, then point mutations; as for egenes in evoca.c. */
static int mutate_child(size_t ch)
{
    uint8_t before[CHEMCA_NSLOT * CHEMCA_GENE_BYTES];
    memcpy(before, GENE(ch, 0), sizeof before);
    uint8_t act0 = ractive[ch];

    int na = core_poisson(cc_rng_next, g_mu_ractive * CHEMCA_NSLOT);
    for (int f = 0; f < na; f++) {
        uint8_t nw = ractive[ch] ^ (uint8_t)(1u << (cc_rng_next() & 7));
        if (nw) ractive[ch] = nw;
    }
    uint8_t newly = ractive[ch] & ~act0;
    while (newly) {
        int s = __builtin_ctz(newly);
        newly &= newly - 1;
        if (cc_unif() >= g_p_dup) continue;
        uint8_t others = ractive[ch] & ~(uint8_t)(1u << s);
        if (!others) continue;
        int pick = (int)(cc_rng_next() % (uint32_t)__builtin_popcount(others));
        while (pick-- > 0) others &= others - 1;
        memcpy(GENE(ch, s), GENE(ch, __builtin_ctz(others)), CHEMCA_GENE_BYTES);
    }

    int np = core_poisson(cc_rng_next, g_mu_rgene * CHEMCA_NSLOT);
    for (int f = 0; f < np; f++) {
        uint8_t *g = GENE(ch, cc_rng_next() & 7);
        float r = cc_unif();
        if (r < 0.7f) {                              /* k_idx +-1, clamped */
            int k = g[CHEMCA_G_K] + ((cc_rng_next() & 1) ? 1 : -1);
            g[CHEMCA_G_K] = (uint8_t)(k < 0 ? 0 : (k >= CHEMCA_NK ? CHEMCA_NK - 1 : k));
        } else if (r < 0.85f) {                      /* new (C, m) pair */
            int j = (int)(cc_rng_next() % (uint32_t)(cm_count() - 1));
            if (j >= cm_index(g)) j++;               /* skip the current pair */
            cm_set(g, j);
        } else {                                     /* A or B redrawn, A != B */
            int which = (cc_rng_next() & 1) ? CHEMCA_G_B : CHEMCA_G_A;
            int other = g[which == CHEMCA_G_A ? CHEMCA_G_B : CHEMCA_G_A];
            uint8_t v;
            do v = (uint8_t)(cc_rng_next() % (uint32_t)gK); while (v == other);
            g[which] = v;
        }
    }
    return ractive[ch] != act0 || memcmp(before, GENE(ch, 0), sizeof before) != 0;
}

/* Body mass: the stores of the impenetrable compounds (u_c == 0). */
static float body_mass(size_t i)
{
    float m = 0.0f;
    for (int c = 0; c < gK; c++) if (cu[c] == 0.0f) m += store[i * gK + c];
    return m;
}

/* The quantity that gates reproduction and ranks eviction targets. */
static inline float gate_value(size_t i)
{
    return g_repro_gate ? body_mass(i) : e_en[i];
}

/* 7. Reproduction, raster order; a child born this step does not itself
 * reproduce in the same step. Gate 'energy': e >= 1, evict the lowest e.
 * Gate 'mass': body mass >= mass_repro, evict the lowest body mass (with a
 * 50/50 split a newborn starts near mass_repro / 2, so this is division at
 * a doubling of birth mass). Dead neighbors count as -inf either way. */
static void phase_repro(void)
{
    int N = gN, K = gK;
    float thr = g_repro_gate ? g_mass_repro : 1.0f;
    for (int row = 0; row < N; row++) {
        for (int col = 0; col < N; col++) {
            size_t p = (size_t)row * N + col;
            if (!alive[p] || births[p] || gate_value(p) < thr) continue;
            /* Lowest-ranked Moore neighbor; dead counts as -inf; exact ties
             * broken uniformly by reservoir sampling. */
            float best = INFINITY; size_t child = 0; int n_best = 0;
            for (int di = -1; di <= 1; di++) {
                int r = core_wrap(row + di, N);
                for (int dj = -1; dj <= 1; dj++) {
                    if (!di && !dj) continue;
                    size_t q = (size_t)r * N + core_wrap(col + dj, N);
                    float v = alive[q] ? gate_value(q) : -INFINITY;
                    if (v < best) { best = v; child = q; n_best = 1; }
                    else if (v == best && cc_rng_next() % (uint32_t)++n_best == 0) child = q;
                }
            }
            if (alive[child]) { kill_cell(child); n_evict++; }
            alive[child] = 1;
            for (int c = 0; c < K; c++) {
                float h = store[p * K + c] * 0.5f;    /* exact: x/2 + (x - x/2) == x */
                store[child * K + c] = h;
                store[p * K + c] -= h;
            }
            float he = e_en[p] * 0.5f;
            e_en[child] = he; e_en[p] -= he;
            memcpy(GENE(child, 0), GENE(p, 0), CHEMCA_NSLOT * CHEMCA_GENE_BYTES);
            ractive[child] = ractive[p];
            births[child] = mutate_child(child) ? 2 : 1;
            n_births++;
        }
    }
}

void chemca_step(void)
{
    if (!gcells) return;
    memset(en_last, 0, sizeof en_last);
    memset(births, 0, gcells);
    n_births = n_deaths = n_tox = n_edeath = n_evict = 0;

    phase_env();
    for (int c = 0; c < gK; c++) core_diffuse_precise(&E[c], &E_tmp, gN, cdiff[c]);
    phase_transport();
    phase_react();
    phase_costs_tox();
    phase_repro();

    for (int k = 0; k < CHEMCA_EN_COUNT; k++) en_cum[k] += en_last[k];
    g_step++;
}

void chemca_step_n(int n) { for (int t = 0; t < n; t++) chemca_step(); }

/* ── Getters ───────────────────────────────────────────────────────── */

void chemca_get_alive(uint8_t *out)  { memcpy(out, alive, gcells); }
void chemca_get_births(uint8_t *out) { memcpy(out, births, gcells); }
void chemca_get_e(float *out)        { memcpy(out, e_en, gcells * sizeof(float)); }
void chemca_get_active(uint8_t *out) { memcpy(out, ractive, gcells); }
void chemca_get_genes(uint8_t *out)
{
    memcpy(out, rg, gcells * CHEMCA_NSLOT * CHEMCA_GENE_BYTES);
}
void chemca_get_E(int c, float *out)
{
    if (c >= 0 && c < gK) memcpy(out, E[c], gcells * sizeof(float));
}
void chemca_get_store(int c, float *out)
{
    if (c < 0 || c >= gK) return;
    for (size_t i = 0; i < gcells; i++) out[i] = store[i * gK + c];
}

void chemca_get_mass(double *me, double *ms)
{
    for (int c = 0; c < gK; c++) {
        double se = 0.0, ss = 0.0;
        for (size_t i = 0; i < gcells; i++) { se += E[c][i]; ss += store[i * gK + c]; }
        me[c] = se; ms[c] = ss;
    }
}

int chemca_get_mass_flux(int c, double *out)
{
    if (c < 0 || c >= gK) return -1;
    out[0] = flux_in[c]; out[1] = flux_out[c];
    return 0;
}

void chemca_get_energy_budget(double *last, double *cum)
{
    if (last) memcpy(last, en_last, sizeof en_last);
    if (cum)  memcpy(cum, en_cum, sizeof en_cum);
}

double chemca_get_e_total(void)
{
    double s = 0.0;
    for (size_t i = 0; i < gcells; i++) s += e_en[i];
    return s;
}

long chemca_get_population(void)
{
    long n = 0;
    for (size_t i = 0; i < gcells; i++) n += alive[i];
    return n;
}

long chemca_get_births_last(void)        { return n_births; }
long chemca_get_deaths_last(void)        { return n_deaths; }
long chemca_get_tox_deaths_last(void)    { return n_tox; }
long chemca_get_energy_deaths_last(void) { return n_edeath; }
long chemca_get_evictions_last(void)     { return n_evict; }

/* ── Display ───────────────────────────────────────────────────────── */

#define CHEMCA_NMODES 4
int chemca_colormode_count(void) { return CHEMCA_NMODES; }

static inline uint32_t chan(float v)
{
    float t = v / g_color_scale;
    return (uint32_t)(255.0f * (t < 0.0f ? 0.0f : (t > 1.0f ? 1.0f : t)));
}

static inline int32_t argb(uint32_t r, uint32_t g, uint32_t b)
{
    return (int32_t)(0xFF000000u | (r << 16) | (g << 8) | b);
}

/* Heat map t in [0, 1]: blue -> cyan -> yellow -> red. */
static int32_t heat(float t)
{
    t = t < 0.0f ? 0.0f : (t > 1.0f ? 1.0f : t);
    float r = t < 0.5f ? 0.0f : (t - 0.5f) * 2.0f;
    float g = t < 0.5f ? t * 2.0f : 1.0f - (t - 0.75f > 0.0f ? (t - 0.75f) * 4.0f : 0.0f);
    float b = t < 0.5f ? 1.0f - t * 2.0f : 0.0f;
    return argb((uint32_t)(255.0f * r), (uint32_t)(255.0f * g), (uint32_t)(255.0f * b));
}

void chemca_colorize(int32_t *buf, int mode)
{
    for (size_t i = 0; i < gcells; i++) {
        int32_t px = argb(0, 0, 0);
        switch (mode) {
        case 0:
            px = argb(chan(E[0][i]), gK > 1 ? chan(E[1][i]) : 0, gK > 2 ? chan(E[2][i]) : 0);
            break;
        case 1:
            if (alive[i]) px = argb(255, chan(store[i * gK]), chan(store[i * gK + 1]));
            break;
        case 2:
            if (alive[i]) px = heat(e_en[i]);
            break;
        case 3:
            if (alive[i] && ractive[i]) {
                int s = __builtin_ctz(ractive[i]);
                px = heat((float)GENE(i, s)[CHEMCA_G_K] / (CHEMCA_NK - 1));
            }
            break;
        default:
            break;
        }
        buf[i] = px;
    }
}

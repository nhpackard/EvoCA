#ifndef CHEMCA_H
#define CHEMCA_H

/*
 * ChemCA: a metabolism-based agent model (Docs/2026-10-23_chem-model-v2.md,
 * §4-10; Docs/plan-chemca.md). Lives in the same library as EvoCA but
 * shares only the pure helpers in core.c; all its state is private to
 * chemca.c, including its own PRNG state.
 *
 * World: an N x N periodic lattice; K <= CHEMCA_KMAX compounds, each with an
 * environmental field E_c(x), a free energy G_c, a diffusion strength, a
 * transport rate u_c, a toxicity threshold, and (open mode) a uniform source
 * and a proportional sink. Environment reactions A -> B (the "light") run at
 * every site. Alive cells carry private stores a_c, an energy e and
 * CHEMCA_NSLOT reaction genes with an active mask.
 *
 * Reaction gene: CHEMCA_GENE_BYTES bytes {A, B, C, m, k_idx}:
 *   A, B   substrate and product (A != B, both < K)
 *   C      catalyst compound, or CHEMCA_C_NONE
 *   m      catalyst order 0..2 (m == 0 iff C == CHEMCA_C_NONE)
 *   k_idx  0..CHEMCA_NK-1, rate k = k_min * 2^(k_idx/4)
 * Rate law: mass action lambda = k * c^m (c = private catalyst amount),
 * integrated exponentially; see chemca_step in chemca.c for the step order.
 */

#include <stdint.h>

#define CHEMCA_KMAX        8     /* max compounds */
#define CHEMCA_NSLOT       8     /* reaction-gene slots per cell */
#define CHEMCA_GENE_BYTES  5     /* A, B, C, m, k_idx */
#define CHEMCA_NK          64    /* k_idx range (6 bits): k_min..k_min*2^15.75,
                                  * 1e-3..55 at the default k_min */
#define CHEMCA_C_NONE      0xFF  /* "no catalyst" */
#define CHEMCA_MAX_ENVRX   8     /* max environment reactions */

/* Gene byte offsets within a slot. */
enum { CHEMCA_G_A = 0, CHEMCA_G_B = 1, CHEMCA_G_C = 2, CHEMCA_G_M = 3, CHEMCA_G_K = 4 };

/* Energy-budget entries (chemca_get_energy_budget). */
enum {
    CHEMCA_EN_LIGHT = 0,   /* env reactions: conversion * (G_B - G_A) */
    CHEMCA_EN_YIELD,       /* agent reactions: sum conversion * (G_A - G_B) -> e */
    CHEMCA_EN_COST,        /* spent on m_e + c_g * active genes */
    CHEMCA_EN_DEATH,       /* e lost at death (toxic, energy, eviction) */
    CHEMCA_EN_SPLIT,       /* reproduction overhead (0: the split is exact) */
    CHEMCA_EN_OPEN_IN,     /* open mode: source * G_c */
    CHEMCA_EN_OPEN_OUT,    /* open mode: sink * G_c */
    CHEMCA_EN_COUNT
};

/* ── Lifecycle ─────────────────────────────────────────────────────── */
int  chemca_init(int N, int K);         /* 0 on success; resets all state */
void chemca_free(void);
void chemca_set_seed(uint32_t seed);    /* this model's own PRNG */
int  chemca_get_N(void);
int  chemca_get_K(void);
long chemca_get_step(void);

/* Layout constants for the Python side. */
int chemca_kmax(void);
int chemca_nslot(void);
int chemca_gene_bytes(void);
int chemca_nk(void);

/* ── Metaparameters ────────────────────────────────────────────────── */
void  chemca_set_k_min(float v);      float chemca_get_k_min(void);
void  chemca_set_m_e(float v);        float chemca_get_m_e(void);
void  chemca_set_c_g(float v);        float chemca_get_c_g(void);
void  chemca_set_mu_rgene(float v);   float chemca_get_mu_rgene(void);
void  chemca_set_mu_ractive(float v); float chemca_get_mu_ractive(void);
void  chemca_set_p_dup(float v);      float chemca_get_p_dup(void);
void  chemca_set_open_mode(int on);   int   chemca_get_open_mode(void);
float chemca_k_of(int k_idx);         /* k_min * 2^(k_idx/4) */
/* Reproduction gate: 0 energy (e >= 1, evict lowest e), 1 mass (body mass
 * >= mass_repro, evict lowest body mass; body = compounds with u_c == 0). */
void  chemca_set_repro_gate(int g);   int   chemca_get_repro_gate(void);
void  chemca_set_mass_repro(float v); float chemca_get_mass_repro(void);

/* Per-compound configuration: free energy, diffusion (core_diffuse
 * strength), transport rate, toxicity threshold (0 = none), open-mode
 * uniform source per site and proportional sink rate. out[6] in that order. */
int  chemca_set_compound(int c, float G, float diff, float u, float tox,
                         float inc, float sink);
int  chemca_get_compound(int c, float *out6);

/* Environment reactions A -> B at a rate per step (scalar, or a per-site
 * field r(x) once set). Applied in the order added. */
void chemca_env_rx_clear(void);
int  chemca_env_rx_add(int A, int B, float rate);   /* index, or -1 */
int  chemca_env_rx_set_rate(int idx, float rate);
int  chemca_env_rx_set_field(int idx, const float *r);  /* NULL = scalar */
int  chemca_env_rx_count(void);
int  chemca_env_rx_get(int idx, int *A, int *B, float *rate, int *has_field);

/* ── State setters ─────────────────────────────────────────────────── */
/* Alive mask. Cells set dead have their stores, e and genes zeroed (an
 * init-time configuration, not a death: nothing is returned or booked). */
void chemca_set_alive(const uint8_t *a);
void chemca_set_E(int c, const float *E);
void chemca_set_store(int c, const float *s);       /* dead cells forced 0 */
void chemca_set_e(const float *e);                  /* dead cells forced 0 */
/* Every alive cell gets genes[CHEMCA_NSLOT * CHEMCA_GENE_BYTES] and active. */
int  chemca_set_genes_all(const uint8_t *genes, uint8_t active);
/* Full arrays [N*N*NSLOT*GENE_BYTES] and [N*N]; returns -1 on an invalid gene. */
int  chemca_set_genes(const uint8_t *genes, const uint8_t *active);
/* Every alive cell gets independent random valid genes, one active slot. */
void chemca_set_genes_random(void);

/* ── Dynamics ──────────────────────────────────────────────────────── */
void chemca_step(void);
void chemca_step_n(int n);

/* ── Getters (copy out into caller buffers) ───────────────────────── */
void chemca_get_alive(uint8_t *out);
void chemca_get_E(int c, float *out);
void chemca_get_store(int c, float *out);
void chemca_get_e(float *out);
void chemca_get_genes(uint8_t *out);
void chemca_get_active(uint8_t *out);
void chemca_get_births(uint8_t *out);     /* last step: 1 normal, 2 mutant, 0 none */
/* mass_env[K], mass_store[K]: double sums of E_c and stores. */
void chemca_get_mass(double *mass_env, double *mass_store);
/* Open-mode mass flux per compound, cumulative: out[2] = {in, out}. */
int  chemca_get_mass_flux(int c, double *out2);
/* Energy budget: last[CHEMCA_EN_COUNT] for the last step, cum[...] since
 * init. Either pointer may be NULL. */
void   chemca_get_energy_budget(double *last, double *cum);
int    chemca_energy_count(void);
double chemca_get_e_total(void);
long   chemca_get_population(void);
long   chemca_get_births_last(void);
long   chemca_get_deaths_last(void);       /* all deaths incl. eviction */
long   chemca_get_tox_deaths_last(void);
long   chemca_get_energy_deaths_last(void);
long   chemca_get_evictions_last(void);

/* ── Display ───────────────────────────────────────────────────────── */
/* Modes: 0 env fields of compounds 0..2 as R, G, B; 1 alive red with stores
 * of compounds 0, 1 as green, blue (dead black); 2 energy e / 1 (dead
 * black); 3 k_idx of the first active gene as a heat map (dead black).
 * Amounts are scaled by the color scale (default 1) and clipped. */
void chemca_colorize(int32_t *buf, int mode);
int  chemca_colormode_count(void);
void  chemca_set_color_scale(float s);
float chemca_get_color_scale(void);

#endif /* CHEMCA_H */

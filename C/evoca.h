#ifndef EVOCA_H
#define EVOCA_H

#include <stdint.h>

/*
 * EvoCA — Evolutionary Cellular Automata core
 *
 * Rule space: outer-totalistic on 3-ring neighbourhood (n1, n2, n3).
 *
 * The LUT is indexed by (v_x, n1, n2, n3) where:
 *   v_x ∈ {0,1}   — current cell state
 *   n1  ∈ {0..4}  — active cells at distance 1  (orthogonal Moore nbrs)
 *   n2  ∈ {0..4}  — active cells at distance √2 (diagonal Moore nbrs)
 *   n3  ∈ {0..4}  — active cells at distance 2
 *
 * GoL is exactly encodable: Moore count = n1+n2; GoL ignores n3.
 *
 * LUT size: 2 × 5 × 5 × 5 = 250 bits  → 32 bytes (bit-packed)
 * Flat bit index:  v_x*125 + n1*25 + n2*5 + n3
 *
 * Fiducial pattern c(x): D4-symmetric binary pattern on the 5×5 grid.
 * 6 independent bits (one per D4 orbit) in the lower 6 bits of egenome.
 * (The fiducial still uses the full 5×5 neighbourhood for eating.)
 */

#define LUT_BITS   250   /* 2*5*5*5 */
#define LUT_BYTES   32   /* ceil(250/8) */

/* Display scale: screen pixels per simulation cell.  Change and recompile. */
#define CELL_PX  2

/* Flat bit index from per-ring counts. */
#define LUT_IDX(vx,n1,n2,n3) \
    ((vx)*125 + (n1)*25 + (n2)*5 + (n3))

/* ── Lifecycle ─────────────────────────────────────────────────────── */

void evoca_init(int N, float food_inc, float m_scale);
void evoca_free(void);

/* ── Metaparam setters ──────────────────────────────────────────────── */

void evoca_set_food_inc(float f);
void evoca_set_m_scale(float m);
void evoca_set_gdiff(float d);
void evoca_set_mu_lut(float m);
void evoca_set_mu_egenome(float m);
void evoca_set_tax(float t);
void evoca_set_restricted_mu(int r);
int  evoca_get_restricted_mu(void);
void evoca_set_diag(int d);
int  evoca_get_diag(void);

/* ── Bulk setters ───────────────────────────────────────────────────── */

void evoca_set_v_all(const uint8_t *v, int len);

/* Set all cells' LUT from a LUT_BYTES-length bit-packed byte array. */
void evoca_set_lut_all(const uint8_t *lut_bytes);

/* Set one cell's LUT. */
void evoca_set_lut(int idx, const uint8_t *lut_bytes);

void evoca_set_egenome_all(uint8_t eg);
void evoca_set_f_all(float f);
void evoca_set_F_all(float F);
void     evoca_set_env_mask(const uint8_t *mask);
uint8_t *evoca_get_env_mask(void);

/* ── Update ─────────────────────────────────────────────────────────── */

void evoca_step(void);

/* ── Activity tracking ─────────────────────────────────────────────── */

void evoca_activity_update(void);
void evoca_activity_render_col(int32_t *col, int height);
int  evoca_activity_get(uint32_t *keys, uint64_t *activities,
                        uint32_t *pop_counts, int32_t *colors, int max_n);

/* ── Activity quantile probe ──────────────────────────────────────── */

void evoca_q_activity_deciles(float *deciles_out);  /* 9 floats: p10..p90 */

/* ── Per-step demography counters ──────────────────────────────────── */

int evoca_get_births_last(void);
int evoca_get_deaths_last(void);

/* ── Neutral shadow population (Channon-style activity calibration) ────
 *
 * Shadow population mirroring the real run's demography under random
 * selection. Same alive count 1:1 with real, same reproduction
 * procedure (LUT bit flips at the same gmu_lut), but births pick a
 * uniform-random surviving shadow parent and deaths pick a uniform-
 * random shadow member. Initialised by copying the current alive cells'
 * LUTs at enable time; mirrored each evoca_step thereafter.
 *
 * N-activity is bucketised by the same FNV-1a LUT-content hash as
 * G-activity, so the two distributions live in the same magnitude
 * space and can be directly compared. */

void evoca_neutral_enable(void);
void evoca_neutral_disable(void);
int  evoca_neutral_is_enabled(void);
int  evoca_neutral_get_population(void);

void evoca_n_activity_update(void);
void evoca_n_activity_render_col(int32_t *col, int height);
int  evoca_n_activity_get(uint32_t *keys, uint64_t *activities,
                          uint32_t *pop_counts, int32_t *colors, int max_n);
void evoca_set_n_act_ymax(int y);
int  evoca_get_n_act_ymax(void);

/* N-overlay alpha (percent, 0..100). 100 = opaque white. */
void   evoca_set_n_overlay_alpha(int pct);
int    evoca_get_n_overlay_alpha(void);

/* Latest N-activity top-decile value, computed each
 * evoca_n_activity_render_col() call. Useful for displaying the
 * threshold-line value in the UI. 0 if no live shadow buckets yet. */
double evoca_get_n_p90(void);

/* Nq-activity: 9 deciles of N-activity (p10..p90). */
void evoca_nq_activity_deciles(float *deciles_out);

/* Activity flux probe: for each live-genome bucket, walk the past `window`
 * ticks (capped at a small compile-time bound) and reconstruct the activity
 * trajectory from stored pop history. A tick contributes one slope sample
 * = pop iff the activity increment interval (A_prev, A_now] overlaps
 * [a_lo, a_hi]. Writes up to max_n float slopes and returns the count
 * written. Use evoca_q_activity_deciles + live diversity count to pick
 * [a_lo, a_hi] bounds from the current distribution. */
int evoca_activity_crossings(uint64_t a_lo, uint64_t a_hi, int window,
                             float *slopes_out, int max_n);

/* Count distinct LUT content hashes across live cells. Useful for driving
 * band bounds in the flux probe and as a diversity trace. */
int  evoca_count_distinct_genomes(void);

/* Total alive cells. */
int  evoca_get_population(void);

/* Copy per-cell age (g_step - last_event_step[x]) into out[N*N].
 * Dead cells get -1. */
void evoca_get_ages(int32_t *out);

/* ── Egenome activity tracking ─────────────────────────────────────── */

void evoca_eg_activity_update(void);
void evoca_eg_activity_render_col(int32_t *col, int height);
int  evoca_eg_activity_get(uint64_t *activities, uint32_t *pop_counts,
                           int32_t *colors);
void evoca_set_eg_act_ymax(int y);
int  evoca_get_eg_act_ymax(void);

/* ── LUT complexity ───────────────────────────────────────────────── */

void evoca_lut_complexity_counts(uint32_t *counts);  /* counts[3]: n1, n1+n2, full */
void evoca_lut_complexity_render_col(int32_t *col, int height);
void evoca_eg_pop_render_col(int32_t *col, int height);  /* stacked egenome pop */

/* ── Local-pattern activity & entropy ─────────────────────────────── */

void  evoca_set_n_ent(int n);         /* 1=VN, 2=Moore, 3=n3; resets history */
int   evoca_get_n_ent(void);
void  evoca_pat_update(void);         /* scan lattice; updates pat_pop/act + entropy */
float evoca_get_entropy(void);        /* Shannon entropy from last evoca_pat_update() */
void  evoca_pat_activity_render_col(int32_t *col, int height);
void  evoca_set_pat_act_ymax(int y);
int   evoca_get_pat_act_ymax(void);

/* ── Visualisation ──────────────────────────────────────────────────── */

/* Fill pixels[N*N] with int32 ARGB values.
   colormode 0: cell state with LUT-hashed genome color (wild-type=white)
   colormode 1: env food F as green; alive cells tinted red
   colormode 2: private food f as blue; alive cells tinted red
   colormode 3: birth events (yellow=birth, magenta=mutant birth, dim=alive, black=dead)
   colormode 4: cell age (cool→hot log gradient; black=dead) */
void evoca_colorize(int32_t *pixels, int colormode);

/* ── Accessors ──────────────────────────────────────────────────────── */

uint8_t *evoca_get_v(void);
float   *evoca_get_F(void);
float   *evoca_get_f(void);
uint8_t *evoca_get_egenome(void);
uint8_t *evoca_get_lut(void);    /* [N*N * LUT_BYTES] */
uint8_t *evoca_get_births(void); /* [N*N] 0=none, 1=birth, 2=mutant birth */
uint8_t *evoca_get_alive(void);  /* [N*N] 1=alive organism, 0=dead slot */
int      evoca_get_N(void);
int      evoca_get_cell_px(void);
float    evoca_get_gdiff(void);
float    evoca_get_mu_lut(void);
float    evoca_get_mu_egenome(void);
float    evoca_get_tax(void);
void     evoca_set_act_ymax(int y);
int      evoca_get_act_ymax(void);
uint32_t *evoca_get_repro_age_hist(void);
int      evoca_get_repro_age_max(void);
uint32_t evoca_get_step(void);
void     evoca_set_repro_age_t0(uint32_t t);
uint32_t evoca_get_repro_age_t0(void);
void     evoca_reset_repro_age_hist(void);

/* ── Alive data plane ──────────────────────────────────────────────── */

void     evoca_set_alive(const uint8_t *arr);  /* set alive array; zeroes dead cells' data */
void     evoca_set_alive_all(void);            /* all cells alive */
void     evoca_set_alive_fraction(float frac); /* random fraction alive */
void     evoca_set_alive_patch(int radius);    /* square patch at center */
void     evoca_set_alive_halfplane(int axis);  /* 0=left, 1=top */

#endif /* EVOCA_H */

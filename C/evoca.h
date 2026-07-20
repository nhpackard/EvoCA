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

/* Multi-egene cell genome: each cell carries up to NEGENOME_MAX egenes
 * (each a 6-bit fiducial pattern), with a per-cell 8-bit `active` mask
 * marking which slots are live. Negene = popcount(active). Inactive
 * slots' egene bits drift under mutation (pseudogenes). */
#define NEGENOME_MAX 8

/* Display scale: screen pixels per simulation cell.  Change and recompile. */
#define CELL_PX  2

/* Flat bit index from per-ring counts. */
#define LUT_IDX(vx,n1,n2,n3) \
    ((vx)*125 + (n1)*25 + (n2)*5 + (n3))

/* ── Lifecycle ─────────────────────────────────────────────────────── */

void evoca_init(int N, float food_inc, float m_scale);
void evoca_free(void);

/* Deterministically reset the PRNG state (seed==0 remapped internally). */
void evoca_set_seed(uint32_t seed);

/* ── Metaparam setters ──────────────────────────────────────────────── */

void evoca_set_food_inc(float f);
void evoca_set_m_scale(float m);
void evoca_set_gdiff(float d);
void evoca_set_mu_lut(float m);
void evoca_set_mu_egene(float m);
void evoca_set_mu_egenome(float m);
void evoca_set_p_dup_egene(float p);
void evoca_set_tax(float t);
void evoca_set_tax_per_egene(float t);
void evoca_set_tax_lut(float t);
void evoca_set_tax_ring(float t);
void evoca_set_restricted_mu(int r);
int  evoca_get_restricted_mu(void);
void evoca_set_food_per(int p);    /* env food-regen period (>=1); default 1 */
int  evoca_get_food_per(void);
void evoca_set_diag(int d);
int  evoca_get_diag(void);

/* ── Bulk setters ───────────────────────────────────────────────────── */

void evoca_set_v_all(const uint8_t *v, int len);

/* Set all cells' LUT from a LUT_BYTES-length bit-packed byte array. */
void evoca_set_lut_all(const uint8_t *lut_bytes);

/* Set one cell's LUT. */
void evoca_set_lut(int idx, const uint8_t *lut_bytes);

/* Set every cell to: all NEGENOME_MAX slots = `eg`, only slot 0 active.
 * Negene = 1 everywhere. Used by Python set_egenome_all path. */
void evoca_set_egenome_all(uint8_t eg);

/* Set every cell to: all NEGENOME_MAX slots = independent random 6-bit
 * values, exactly one random slot active per cell. Negene = 1
 * everywhere; the inactive slots act as a pseudogene reservoir. */
void evoca_set_egenome_random_all(uint8_t wt);

/* Set every cell to a single active egene with the given (value, mask)
 * ternary pair. Use this to seed runs at any cognitive complexity. */
void evoca_set_egenome_pair_all(uint8_t value, uint8_t mask);
void evoca_set_f_all(float f);
void evoca_set_F_all(float F);
/* Per-cell agent mutation-rate gene (mutation-of-mutation Wave-1 storage
 * scaffold). Inherited on reproduction but read nowhere that affects
 * dynamics, and not itself mutated, so runs stay bit-for-bit identical to
 * pre-gene behaviour. Seed a scalar with _all, or an N*N row-major float
 * array for the distribution-seeding assay. */
void  evoca_set_mu_lut_gene_all(float g);
void  evoca_set_mu_lut_gene(const float *g);
void     evoca_set_env_mask(const uint8_t *mask);
uint8_t *evoca_get_env_mask(void);

/* ── Update ─────────────────────────────────────────────────────────── */

void evoca_step(void);

/* ── Activity tracking ─────────────────────────────────────────────── */

void evoca_activity_update(void);
void evoca_activity_render_col(int32_t *col, int height);
int  evoca_activity_get(uint32_t *keys, uint64_t *activities,
                        uint32_t *pop_counts, int32_t *colors, int max_n);

/* ── Effective-rule activity (eff_activity) ─────────────────────────
 * Genome-species activity whose species key is the cell's LUT MASKED to
 * the entries the POPULATION has exercised (eff_union_mask = OR of the
 * per-cell used-masks accumulated during Phase 1), folded with the same
 * egene stream as the full-genome key. Two genomes differing only in
 * entries never exercised by any cell collapse to one effective-rule
 * species, so silent LUT churn does not mint a new species (unlike
 * full-genome `activity`). Additive no-op to dynamics. The alive
 * eff-species count is <= the alive full-genome `activity` count each
 * step (the effective-rule map is a coarsening of genome space).
 *
 * NOTE: the key uses the population-union mask, not each cell's own
 * per-cell mask. A per-cell key would SPLIT one genome across cells with
 * divergent neighbourhood histories, making eff-species larger than
 * full-species and defeating the collapse goal. The per-cell masks are
 * still accumulated and exposed (get_eff_used_mask) as the diagnostic
 * substrate the union is built from. */
void evoca_eff_activity_update(void);
void evoca_eff_activity_render_col(int32_t *col, int height);
int  evoca_eff_activity_get(uint32_t *keys, uint64_t *activities,
                            uint32_t *pop_counts, int32_t *colors, int max_n);
int  evoca_eff_count_distinct(void);              /* alive eff-species count */
void evoca_get_eff_hashes(uint32_t *out);         /* out[N*N] per-cell eff hash */
void evoca_get_eff_used_mask(uint8_t *out);       /* out[N*N*LUT_BYTES] per-cell */
void evoca_get_eff_union_mask(uint8_t *out);      /* out[LUT_BYTES] population union */
void evoca_set_eff_act_ymax(int y);
int  evoca_get_eff_act_ymax(void);

/* ── Activity quantile probe ──────────────────────────────────────── */

void evoca_q_activity_deciles(float *deciles_out);  /* 9 floats: p10..p90 */

/* ── Per-step demography counters ──────────────────────────────────── */

int evoca_get_births_last(void);
int evoca_get_deaths_last(void);

/* ── Realised per-component mutation flux (fixed-space shadows) ──────
 * Bit-flip counts actually applied to children in the most recent
 * evoca_step (measured, not the nominal mu_* metaparameter). Used to
 * drive the eg_activity / dyn_activity neutral baselines. */
long evoca_get_eg_bitflips_last(void);
long evoca_get_lut_bitflips_last(void);
long evoca_get_repro_last(void);

/* ── Fixed-space neutral baselines (eg_activity / dyn_activity) ──────
 *
 * Closed-form Wright–Fisher drift shadows over the two FIXED finite
 * bucket spaces (729 ternary egene keys; 500 LUT (input,output)
 * buckets), driven by the realised per-component mutation flux. They
 * close the gap left by the LUT-hash Channon shadow, which models LUT
 * bytes only while get_activity hashes the full genome — so excess was
 * undefined for egene-driven runs.
 *
 * evoca_eg_excess_pc / evoca_dyn_excess_pc return the per-component-
 * normalised cumulative excess (observed ΣG/D_G − shadow ΣN/D_N),
 * consistent with python/evoca_explore.py's excess_pc_*. ≈ 0 under
 * neutral drift, clearly > 0 under selection. Advanced automatically
 * inside evoca_step(); the caller must keep eg_pop/dyn_pop current via
 * evoca_eg_activity_update() (dyn is updated by Phase 1 every step). */
double evoca_eg_excess_pc(void);
double evoca_dyn_excess_pc(void);
double evoca_eg_shadow_total(void);
double evoca_dyn_shadow_total(void);

/* Per-bucket z-scores — the per-component ("per-gene") neutral shadow,
 * well-posed in these FIXED spaces (ill-posed in the open LUT-hash
 * space, do not attempt there). For each bucket k:
 *     z[k] = (observed cumulative activity[k] − shadow drift mean[k])
 *            / sqrt(closed-form WF drift variance[k])
 * z ≈ 0 for every bucket under a null/frozen (zero-flux) arm; z ≫ 0 on
 * buckets the real run concentrates on faster than drift (selection),
 * z < 0 where drift stays warmer than selection. Buckets with no
 * accumulated drift variance yield z=0. Fill caller buffers of length
 * evoca_eg_nbuckets() (729) / evoca_dyn_nbuckets() (500). Advanced
 * automatically inside evoca_step() alongside the excess accessors.
 * NOTE: dyn z-scores are interpretable only where mu_lut>0. */
void evoca_eg_zscores(double *out);
void evoca_dyn_zscores(double *out);
int  evoca_eg_nbuckets(void);
int  evoca_dyn_nbuckets(void);

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

/* ── Opt-in per-cell lineage record (research-directions §1) ──────────
 * Default OFF. Call evoca_enable_lineage(1) BEFORE evoca_init to
 * allocate the arrays and have Phase-4 reproduction record, per child:
 * the parent's full-genome FNV-1a hash (same hash as the activity
 * tracker), a monotonic birth id, and the parent's birth id (for
 * parent->child chain walking). When OFF there is zero added work in
 * the hot path and no memory used. Additive only — no effect on the
 * CA or any existing metric. */
void evoca_enable_lineage(int on);
int  evoca_lineage_enabled(void);
void evoca_get_lineage_parent_hash(uint32_t *out);  /* out[N*N] */
void evoca_get_lineage_birth_id(uint64_t *out);     /* out[N*N] */
void evoca_get_lineage_parent_id(uint64_t *out);    /* out[N*N] */

/* ── Per-lineage net metabolic return reducer (checkpoint §D2) ────────
 * Also gated by evoca_enable_lineage(); zero-cost when lineage is OFF.
 * Net return of a lineage (keyed by full-genome FNV-1a hash) = cumulative
 * intake by its members − cumulative tax paid by its members, accumulated
 * over the whole run and retained after extinction. */
void   evoca_lineage_return_update(void);   /* refresh pop_count/last_step */
int    evoca_lineage_return_count(void);    /* # distinct lineages         */
int    evoca_get_lineage_return(uint32_t *hash, double *intake, double *tax,
                                uint64_t *births, uint32_t *pop,
                                uint32_t *first_step, uint32_t *last_step);
double evoca_get_lineage_intake_total(void);
double evoca_get_lineage_tax_total(void);
void   evoca_get_lineage_intake(double *out); /* out[N*N] per-cell, this life */
void   evoca_get_lineage_tax(double *out);    /* out[N*N] per-cell, this life */
void evoca_get_lut_hashes(uint32_t *out);           /* out[N*N] LUT-only FNV-1a */

/* ── Egenome activity tracking ─────────────────────────────────────── */

void evoca_eg_activity_update(void);
void evoca_eg_activity_render_col(int32_t *col, int height);
/* Egene food intake probe: cumulative food obtained per egene byte
 * value (mouthful split equally across max-match-tied winners). */
void evoca_eg_food_render_col(int32_t *col, int height);
int  evoca_eg_food_get(uint64_t *food_out, uint32_t *pop_counts,
                        int32_t *colors);
void evoca_set_eg_food_ymax(int y);
int  evoca_get_eg_food_ymax(void);
/* Egenome scalar stats: fills out[0..4] = [mean Negene, std Negene,
 * distinct egene values (out of 64), mean max-match, frac at max]. */
void evoca_egenome_stats(float *out);

/* Egene cognitive stats: fills out[0..2] = [mean cognitive specificity
 * (non-wildcard cell-positions per active egene, range 0..25),
 * mean cognitive load (per-cell sum of cell-positions across active
 * egenes, range 0..200), mean food intake per alive eater this step
 * (range 0..1)]. */
void evoca_egene_stats(float *out);

/* Dissipation / energy-flux accounting (checkpoint §L, action L1).
 * Fills out[0..15] with the source→sink free-energy through-flux:
 *   [0] source_step [1] diss_step [2] tax_step [3] death_step
 *   [4] overflow_step [5] source_cum [6] diss_cum [7] tax_cum
 *   [8] death_cum [9] overflow_cum [10] alive [11] mouthful_step
 *   [12] phi (=diss_step/alive) [13] storage_env [14] storage_priv
 *   [15] step. All values in food units. See evoca.c for full doc. */
#define EVOCA_DISS_N 16
void evoca_get_dissipation(double *out);

/* Dyn-activity probe: 500 buckets indexed by (input * 2 + output),
 * where input is the 250-value LUT bit index from
 * LUT_IDX(v_x, n1, n2, n3) and output is the new cell state.
 * Cumulative activity tally + per-step population (count of alive
 * cells whose Phase-1 transition fell in this bucket). */
void evoca_dyn_activity_render_col(int32_t *col, int height);
int  evoca_dyn_activity_get(uint64_t *activities, uint32_t *pop_counts,
                             int32_t *colors);
void evoca_set_dyn_act_ymax(int y);
int  evoca_get_dyn_act_ymax(void);
/* Count of (input, output) buckets ever observed (range [0, 500]). */
int  evoca_get_dyn_distinct(void);
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
float   *evoca_get_mu_lut_gene(void);  /* [N*N] per-cell LUT mutation rate */
/* Returns pointer to a [N*N] scratch array holding, for each cell, the
 * lowest-index active slot's egene byte (0 if no slot is active). The
 * scratch is rebuilt on each call. Read-only; for Python display use. */
uint8_t *evoca_get_egenome(void);

/* Returns pointer to the live [N*N*NEGENOME_MAX] egene-value-bytes array. */
uint8_t *evoca_get_egenes(void);

/* Returns pointer to the live [N*N*NEGENOME_MAX] egene-mask-bytes array.
 * Bit i in mask byte = 1 means orbit i is non-wildcard (consults the
 * value bit); bit i = 0 means orbit i is wildcard (no contribution to
 * eating, no per-position tax). */
uint8_t *evoca_get_egenes_mask(void);

/* Returns pointer to the live [N*N] active-mask byte array (one byte
 * per cell, bit i = slot i is active). */
uint8_t *evoca_get_active(void);

uint8_t *evoca_get_lut(void);    /* [N*N * LUT_BYTES] */
uint8_t *evoca_get_births(void); /* [N*N] 0=none, 1=birth, 2=mutant birth */
uint8_t *evoca_get_alive(void);  /* [N*N] 1=alive organism, 0=dead slot */
int      evoca_get_N(void);
int      evoca_get_cell_px(void);
float    evoca_get_gdiff(void);
float    evoca_get_mu_lut(void);
float    evoca_get_mu_egene(void);
float    evoca_get_mu_egenome(void);
float    evoca_get_p_dup_egene(void);
float    evoca_get_tax(void);
float    evoca_get_tax_per_egene(void);
float    evoca_get_tax_lut(void);
float    evoca_get_tax_ring(void);
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

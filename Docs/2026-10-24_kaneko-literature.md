# Kaneko-group cell models with catalytic networks, and how chemca relates (2026-10-24)

For N. Literature check on Kaneko's simulations of cells containing autocatalytic networks, plus
the closest spatial relatives (Hogeweg's group), against the chemca design
(`Docs/2026-10-23_chem-model-v2.md`). PDFs are in `Docs/references/kaneko/` and
`Docs/references/hogeweg/`; bibliographic entries are in `Docs/references/references.md`. I read
the model sections of the key papers and the abstracts of the rest. Where I state a detail, it is
from the paper. Interpretations are marked as mine.

## Summary

- **Yes, chemca is heading toward a 2D-lattice relative of the Kaneko cell models.** It is not a copy, though. The Kaneko line differs from chemca in four ways that matter:
  1. **No genome.** The reaction network is fixed per run, and a cell inherits only its chemical composition (split almost equally at division). Differentiation into cell types is a dynamical instability ("isologous diversification"), not genetic change.
  2. **Evolution, where present, is outside the simulation.** Networks are mutated and selected in a generational loop: n mothers, m mutants each, keep the fastest-dividing (Furusawa & Kaneko 2005; 2000 compares random networks rather than evolving them in place).
  3. **Growth is mass, not energy.** A cell grows by taking up penetrable nutrients and converting them into impenetrable chemicals, and divides when its volume (total content) doubles. There is no free energy and no waste toxicity. Supply comes from a chemical bath at the boundary.
  4. **Space is a continuum medium.** Cells sit in a diffusing medium (a 2D grid in the 1998 *Artificial Life* paper, a 1D chain in the 2000 PRL) and push neighbors aside on division. There is no eviction.
- **chemca would combine both kinds of inheritance in situ**: heritable reaction rates (rgenes), plus compositional inheritance (stores split 50/50), with selection happening on the lattice by eviction rather than in an external loop. As far as I found, that combination (Kaneko-style compositional dynamics, genetic rates and in-situ spatial selection) is not what the Kaneko papers do. That is my reading of the papers I checked, not an exhaustive search.
- **Two results transfer directly as hypotheses for chemca:**
  - Isologous diversification: genetically identical cells can split into distinct compositional types through intracellular instability plus coupling through the medium. chemca could show this with `mu_rgene = 0`, and it needs a readout (clustering cells by store composition and dynamics, independent of genes).
  - Minority-molecule control (Kaneko & Yomo 2002; Kamimura & Kaneko 2010): a slowly replicated, scarce species becomes the carrier of heredity and controls the cell's state. In gs3, the activator A, inherited at division and needed for metabolism, is a candidate minority-control molecule.
- **The rate law matches.** Kaneko's 1998/2000 cells use catalysis quadratic in the catalyst, `c_substrate · c_catalyst^α` with α = 2 ("in most simulations"), which is the gs3 order-2 gene (F + 2A → 3A). Takagi & Kaneko (2001) extend Gray-Scott itself to many chemicals and find spots that replicate and differentiate. That is the RD analog of the evolutionary Gray-Scott direction.
- **Hogeweg's group is the closest spatial and evolutionary relative, and it has spirals.** Boerlijst & Hogeweg (1991) showed that hypercycles on a lattice (a CA) self-organize into **spiral waves**, which make them resistant to parasites: selection acts at the level of spirals. Takeuchi & Hogeweg (2009) compare spatial self-organization (traveling waves) with explicit compartments as two routes to multilevel selection. Given N's interest in spirals and our CHEM-N6/N7 wave results, this line probably matters more for chemca's evolutionary questions than the Kaneko cell papers do.

## The papers

| Paper (file) | Model | Space | Heredity / evolution | Key result |
|---|---|---|---|---|
| Kaneko & Yomo 1997, *Bull. Math. Biol.* (adap-org/9606002) | cells with catalytic networks; divide after chemical consumption | well-mixed medium | composition only | **isologous diversification**: synchronized oscillations lose synchrony as cells multiply, forming distinct, heritable cell types |
| Furusawa & Kaneko 1998, *Artificial Life* 4:79 (`..._ALife.pdf`) | k chemicals, random fixed network Con(i,j,ℓ), rate `e1·c_m·c_j²`; penetrable vs impenetrable chemicals; growth = uptake; divide at volume doubling; state-dependent adhesion | **2D grid medium**, diffusion De, bath boundary | composition only | ring pattern of differentiated cells; regeneration after damage; with adhesion, release of seed cells, so a multicellular life cycle |
| Furusawa & Kaneko 1998, adap-org/9802002 | same family | | composition | hierarchical differentiation from a "stem" type; stable ensembles |
| Kaneko & Furusawa 1999, chao-dyn/9912005; Furusawa & Kaneko 2000, nlin/0009007 | same | | composition | chaotic intracellular dynamics give stem cells (stochastic proliferate-or-differentiate) |
| Furusawa & Kaneko 2000, *PRL* 84:6130 (nlin/0009008) | k = 20 chemicals, 6 paths each, 3 penetrable, α = 2, 800 random networks | 1D chain in a medium (2D said consistent) | compares random networks, no in-run evolution | about 5% of networks give **exponential** growth; those are the ones with differentiation and chaotic dynamics ("cooperative" resource use); the rest grow linearly (only edge cells divide) |
| Kaneko & Yomo 2002 (nlin/0203039); Kaneko 2003 (cond-mat/0302313) | stochastic molecules, division at a molecule-count threshold, replication errors | single cell lineages | composition, kinetics | three phases (no recursive production / itinerancy / fixed); **minority molecules carry heredity**; evolvability in the itinerant phase |
| Furusawa & Kaneko 2003, *PRL* 90:088102 (physics/0209103) | random catalytic network `i + ℓ → j + ℓ`, nutrient uptake by membrane diffusion D, division at N_max | single cell | composition | at the critical uptake D_c (maximal growth), abundances follow **Zipf's law**, with a hierarchy of catalysis |
| Furusawa & Kaneko 2005 (q-bio/0504032) | same, uptake via transporters | single cell | **network evolution by external selection** (n × m mutants, keep the fastest) | selection drives cells to the critical point; the power law becomes embedded in the network topology (scale-free) |
| Kaneko 2005 review, *Adv. Chem. Phys.* (q-bio/0504016) | survey: recursive production, minority control, constructive biology | | | the program's overview; §3.2 is the modeling strategy |
| Kamimura & Kaneko 2010, *PRL* 105:268103 (1005.1142); 2017 (1711.08988) | two-species hypercycle, molecule-level, crowding | protocell | composition | a protocell reproduces when one species is slow (a minority); resists parasites; 2017 has a 2D RD appendix |
| Takagi & Kaneko 2001 (nlin/0104058) | **Gray-Scott extended to many chemicals** | 2D RD, no cells | none | self-replicating spots; homogeneous, oscillatory and "multipotent" spot types; spot differentiation |
| Boerlijst & Hogeweg 1991, *Physica D* 48:17 (`hogeweg/`) | hypercycle of self-replicating molecules on a CA lattice | 2D CA | molecules replicate with mutation (parasites) | **spiral waves** emerge and make the hypercycle resistant to parasites; selection at the level of spirals |
| Boerlijst & Hogeweg 1991, *Artificial Life II* (`hogeweg/`) | same, self-structuring | 2D CA | | companion paper |
| Takeuchi & Hogeweg 2009, *PLoS Comput. Biol.* 5:e1000542 (`hogeweg/`) | RNA-like replicators | 2D CA vs compartments | in situ | spatial self-organization (traveling waves) vs compartments as multilevel selection |

Not downloaded (scripted download refused): Takeuchi, Kaneko & Hogeweg 2016, *Proc. R. Soc. B*
283:20153109, "Evolutionarily stable disequilibrium"
(https://dspace.library.uu.nl/bitstream/1874/347072/1/Evo.pdf); von der Dunk, Colizzi & Hogeweg
2017, *Life* 7:43 (https://www.mdpi.com/2075-1729/7/4/43). Both open; they need a browser download.

## What this means for chemca (my reading)

1. **Reproduction gate: energy (chemca's default) or mass (Kaneko)?**
   - Kaneko's growth-by-uptake-and-division-at-doubling is §4's option (iii), biomass, which I had deferred.
   - The Kaneko results show it alone suffices for rich behavior (differentiation, Zipf, life cycles), with no thermodynamics.
   - chemca's energy currency adds what Kaneko lacks: waste, toxicity, and a free-energy budget, which is where our waves come from.
   - Worth offering both as a gate option: `repro_gate = 'energy' | 'mass'`. With mass, impenetrable compounds (transport u = 0) are the body, and division happens at doubling.
   - This is a design question for N, not something to build silently. The engine agent is building the energy gate.
2. **Compositional heredity is already in chemca**, because stores split 50/50. So chemca can test isologous diversification in situ by switching genes off (`mu_rgene = 0`, identical genes) and asking whether cells split into compositional types. That needs a readout: clustering of per-cell store vectors and their time series. It is cheap to add to `run_chemca`.
3. **Network size.** Kaneko's interesting cases need k of about 20 species and several paths each, with only about 5% of random networks giving exponential growth. chemca's stage 3 (`n_chem > 2`, (A, B, C, k) genes) is where this lives. In chemca the network is per-cell and heritable, so it can evolve in situ rather than in an external loop. That is a genuine difference and a candidate research question: does in-situ spatial selection find Kaneko's "cooperative, differentiating" networks?
4. **Spirals and multilevel selection** (Boerlijst & Hogeweg). Our food/waste waves (CHEM-N6) and the CHEM-N7 result (genotype lifetimes collapse in the wave regime, yet the egene count rises) look like the beginning of selection mediated by spatial pattern. The hypercycle literature says spirals can turn selection into selection among spirals. That is a concrete thing to look for in vanilla food/waste (Track A) as well as in chemca.
5. **Rate law and toxicity.** The α = 2 catalysis is consistent with gs3. Kaneko's models have no toxicity, so "need some, die of too much" (§6) is our addition. In Kaneko's cells the analogous tension is the critical nutrient uptake D_c: too much uptake leaves no room for catalysts, too little starves.

## Questions for N

1. Add a mass-based reproduction gate (Kaneko style) as an option next to the energy gate?
2. Add an isologous-diversification readout (compositional clustering with genes frozen) to the chemca test plan?
3. Pursue the Hogeweg spiral and multilevel-selection angle in vanilla food/waste now (Track A), using the N6/N7 machinery?

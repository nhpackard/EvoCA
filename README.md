# EvoCA 

Evolutionary Cellular Automata: a spatially inhomogeneous binary 2D CA where every cell is regarded as an organism, carrying a genome governing its local rule.  There is a resource field (food) modeled as a spatial field over a lattice the same size as the CA lattice.  Each organism at each lattice site can interact with the food field by eating a mouthful of food, transferring a fraction of the food from the food field to the organism's private stash of food.  The eating algorithm is controlled by a local pattern match using an additional piece of genetic data carried by each organism.  Every time step, the organism is taxed; its private food stash is diminished by a fixed amount.  The survival of an organism is a competition between the tax and the ability of the organism to evolve a rule genome and an eating genome to get more food than is lost by the tax.  If an organism's food stash reaches 1.0, it reproduces with genetic mutation; the offspring replaces the neighbor with the least food.  If an organism's food stash reaches 0, it dies: its  LUT is set to all zeroes. It no longer eats, its CA state is zero.  Dead cells can only come alive through a reproductive event.

# Implementation notes

Find detailed implementation notes in [model.md](./Docs/model.md)

The initial implementation of ths model was a collaboration between N. Packard and Claude Code

## brief summary

CA dynamics are implemented in C for speed.  Control is via a python API running inside a jupyter notebook.  The API enables the launch of SDL graphics windows with widgets controling the simulation.  Launch of an example run:
```
N = 512
params = {'N':N,
          'food_inc':0.12,
          'm_scale':0.4,
          'mu_lut':0.001,
          'tax':0.05,
          'restricted_mu':True}
sim = EvoCA()
sim.init(**params)
#sim.set_lut_all(gol_lut)
sim.set_lut_random(n_init=1)
sim.set_egenome_all(0b000011)
sim.set_v(rng2.integers(0, 2, (N, N), dtype=np.uint8))
sim.set_f_all(0.1)
sim.set_F_all(0.5)

run_with_controls(sim, probes={'lut_complexity': True, 'activity': True},diag=True)
```

Metaparameters may be changed mid-stream with widgets.  Initial and current metaparams can be exported with one of the widgets, and then re-imported: a run repeated with 
```
sim,kw = import_run(recipefile)
run_with_controls(sim, **kw)
```

# License

The license is an MIT license, found in the file [LICENSE](./LICENSE).
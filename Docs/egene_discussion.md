
This doc is to sharpen the notion of egene evolution, and the aims of its study.

# overall goals

The main goal is to study evolution of cognition, and to have a concrete model that enables the observation of evolution of increased cognitive ability through an evolutionary process driven by selection of the fittest in a population of agents that have evolving cognitive aparatus. 

The EvoCA platform:  provides a cognitive mechanism in the form of egenes.  Tax on egenes makes cognition costly; survival is a balance between cognition enabling greater food harvesting to overcome the cost of the cognition.

The prior:  we should see a flow of resources (food) to support evolution of cognitive ability.  This depends on the cognitive aparatus being capable of increasing food harvest of agents.

# sharpening of egene / feeding

To sharpen the egenes functionality I propose some changes:

- rather than have egenes code for a pattern of 1's and 0's that get matched to a local pattern, we change to gene values of 0, 1, or *, where * means 'don't care'

- rather than simply summing up matches to establish the food harvesting mouthful, over all non-* genes sum +1 for match and -1 for non-match.  Penalty for non-match increases value of the egene template actually matching.

- etax should be proportional to number of non-* elements of egene.

- we need to establish a proper normalization of eating capability given matching of genome to local pattern.  Open to suggestion.

# review measures for progress

We have had a discussion of measuring evolutionary progress of agent's eating capabilities.  I sould like to review these, and implement probes that reflect the sharpest measurements.

# evolutionary test

This section is not particular to egenes and evolution of cognition, but a feature for all EvoCA evolutionary experiments:

I propose to be able to extract a patch of the population (e.g. middle square of lattice, side length 1/3 of lattice), then insert that into another population, either initialized or itself evolved, to then observe whether the patch grows or shrinks.  growth => "evolutionary success", shrink => "evolutionary failure"

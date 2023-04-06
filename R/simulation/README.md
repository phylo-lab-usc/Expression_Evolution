This folder contains code used for simulations.

function.R contains functions that are repeatedly used in the simulations.

sim_basic.R is for simulation along replicate lineages with stabilizing selection on the protein level or under neutrality.

sim_basic_auto.R is for simulation along replicate lineages under different combination of effective population size and strength of selection.

sim_tree_big.R is for simulation along a phylogenetic tree with stabilizing selection on the protein level or under neutrality.

sim_tree_big_rescale.R is for simulations along lambda transformations of a given phylogenetic tree.

sim_interaction_2g_auto.R is for simulation along replicate lineages for two interacting genes, including a target gene that is subject to stabilizing selection and a regulator that is not directly under stabilizing selection.

sim_interaction_2g_auto2.R is for simulation along replicate lineages for two interacting genes that are both subject to direct selection.

sim_interaction_3g.R is for simulation along replicate lineages for three interacting genes.

sim_basic_dir.R is for simulation along replicate lineages with directional selection on the protein level.

sim_basic_dir2.R is for simulation along replicate lineages with directional selection on multiple genes with different protein levels.

tr_50.txt is the phylogenetic tree with 50 species used in the study. It is an input file used by sim_tree_big_rescale.R and sim_tree_big_rescale.R.

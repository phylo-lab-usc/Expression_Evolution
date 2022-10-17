This folder contains simulations for 2 interacting genes.<br />

Output for each parameter combination would be compressed into a package.<br />

File names are in "width_C12_C21" or "width1_width2_C12_C21" format.<br />
If only gene 2's protein level is under selection, name of the package would be in "width_C12_C21" format, where "width" is the SD of the fitness function for gene 2 (set to be 0 under neutrality), C12 is gene 1's effect on gene 2, and C21 is gene 2's effect on gene 1.<br />
If both gene's protein levels are under selection, "width1" and "width2" represent SDs of the fitness functions for gene 1 and 2, respectively.<br />

Package content: <br />
End-point genotypic values (gt_end.txt) and phenotypes (pt_end.txt) for all lineages<br />
Correlation matrix for end-point genotypic values and phenotypes (cor_mat.txt)<br />
Variances of genotypic values and phenotypes through time (var_all.txt)<br />
Folder for each replicate lineage's genotypic values (gt_all.txt) and phenotypes (pt_all.txt) through time<br />

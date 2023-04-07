This folder contains results of simulations along a phylogenetic tree.

1_out_tr.txt: End-point phenotypes at all tips when there is stabilizing selection. Each row represents an independent simulation. Columns of the table correspond to: mRNA level of species 1, translation rate of species 1, protein level of species 1, mRNA level of species 2, translation rate of species 2, protein level of species 2, ..., mRNA level of species 50, translation rate of species 50, protein level of species 50.

1_rmat_tr.txt: Average evolutionary variance-covariance (VCV) matrix across simulations under stabilizing selection.

0_out_tr.txt: End-point phenotypes at all tips under neutrality. Each row represents an independent simulation. Columns of the table correspond to: mRNA level of species 1, translation rate of species 1, protein level of species 1, mRNA level of species 2, translation rate of species 2, protein level of species 2, ..., mRNA level of species 50, translation rate of species 50, protein level of species 50.

0_rmat_tr.txt: Average evolutionary variance-covariance (VCV) matrix across simulations under neutrality.

1_cor_error_out.txt: Evolutionary correlations between traits in the presence of stabilizing selection. Columns of the table correspond to: SD of measurement error, evolutionary correlation between the mRNA level and the translation rate, evolutionary correlation between the mRNA level and the protein level.

0_cor_error_out.txt: Evolutionary correlations between traits under neutrality. Columns of the table correspond to: SD of measurement error, evolutionary correlation between the mRNA level and the translation rate, evolutionary correlation between the mRNA level and the protein level.

out_tr_rescale_final.txt: Results of simulations along transformed trees in the presence of stabilizing selection. Columns of the table correspond to: Pagel's Lambda, evolutionary correlation between the mRNA level and the translation rate, evolutionary correlation between the mRNA level and the protein level.

Directory ./error contains data subject to measurement error.

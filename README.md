# ModelTeller
A machine learning algorithm for phylogenetic model selection.

ModelTeller extracts features from the input nucleotide multiple sequence alignment and predicts the best substitution model for phylogenetic reconstruction, with emphasis on branch-lengths accuracy!
Possible nucleotide models: JC, F81, K2P/K80, HKY, SYM, and GTR -- all with an additional parameter for the proportion of invariable sites (+I), rate-heterogeneity across sites sampled from the gamma distribution (+G), both (+I+G), or none.

Output:
1. The best predicted model and the ranking of alternative models, in case that you cannot use the best model.
2. The maximum-likelihood tree, reconstructed using the best model.

Prerequisites:
PhyML3 from 
Python-3
modules: h2o.ai, ete3, numpy, scipy, biopython

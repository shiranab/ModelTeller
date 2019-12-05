# ModelTeller
A machine learning algorithm for phylogenetic model selection.

ModelTeller extracts features from the input nucleotide multiple sequence alignment and predicts the best substitution model for phylogenetic reconstruction, with emphasis on branch-lengths accuracy!
Possible nucleotide models: JC, F81, K2P/K80, HKY, SYM, and GTR -- all with an additional parameter for the proportion of invariable sites (+I), rate-heterogeneity across sites sampled from the gamma distribution (+G), both (+I+G), or none.

# Output:
The output files will be located in the directory in which you execute ModelTeller.
1. The best predicted model and the ranking of alternative models, in case that you cannot use the best model.
2. The maximum-likelihood tree, reconstructed using the best model.

# Prerequisites:
1. Python-3.5 or newer
2. python modules: h2o.ai, ete3, argparse, numpy, scipy, biopython
3. Phyml (can be found in our code, in directory phyml_exe).

# How to install?
1. download the ModelTeller code with all its scripts and directories.
2. download the cached ModelTeller model from: https://modelteller.tau.ac.il/modelteller_drf_models.tar.gz
and unzip it in the location of the script.
3. The definitions.py file contains the paths to the ModelTeller trained model (step 2 above) and the phyml executables.
Verify that: 
- The Phyml exe that is suitable for your OS is defined
- The ModelTeller trained model is located in the directory "drf_models"
4. Run modelTeller!

# How to run?
## The -m <msa_file> parameter:
A mandatory input for ModelTeller is an MSA file in a nucleotide format, in one of the following formats: fasta, phylip (interleaved or sequential), clustel, emboss, nexus, Ig, nexus, mauve, and Stockholm.
## Without any additional parameters
this option will rapidly select a single model for maximum-likelihood phylogeny reconstruction. The maximum-likelihood phylogeny, i.e., the model parameters, the branch-lengths, and the topology will be optimized for you after the model is predicted.
## The -g parameter:
this option will first optimized the GTR+I+G phylogeny and then predict the optimal model for branch-length estimation. The maximum-likelihood phylogeny will be optimized given the GTR+I+G fixed topology. Prediction of the best model might take a little while due to the preliminary GTR+I+G phylogeny reconstruction, but the final tree will be computed rapidly. This procedure was demonstrated to perform best both for topologies and branch-length estimation accuracies on complex simulation models.
## The -u <user_tree> parameter:
if you have a good, validated topology for your data, please provide it and ModelTeller will predict the best model for branch-length estimation. The maximum-likelihood phylogeny will be computed for you given your fixed topology.

# Examples:
python modelteller.py -m example/test_msa.phy

python modelteller.py -m example/test_msa.phy -g

python modelteller.py -m example/test_msa.phy -u example/test_tree.txt

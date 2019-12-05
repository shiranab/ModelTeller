import argparse, os, sys, re, logging, itertools, shutil, math, copy
import pandas as pd
import numpy as np
from Bio import AlignIO
from Bio.Align import AlignInfo
from ete3 import Tree
import ete3.coretype.tree
import h2o
from collections import Counter


script_dir = os.path.dirname(__file__)
MODELTELLER_DRF_MODEL = os.path.join(script_dir, 'drf_models','ModelTeller_model')
MODELTELLERg_DRF_MODEL = os.path.join(script_dir, 'drf_models','ModelTellerG_model')
PHYML_SCRIPT = os.path.join(script_dir, "phyml_exe", "PhyML_3.0_linux64")


BASE_MODELS = ["JC", "F81", "K80", "HKY", "SYM", "GTR"]
MODELS_TAGS = ["", "+I", "+G", "+I+G"]
ALL_PHYML_MODELS = [base_model + tag for base_model in BASE_MODELS for tag in MODELS_TAGS]

FEATURES_TO_INCLUDE = ["ac_subs", "ag_subs", "aln_entropy", "at_subs", "base_freqs_entropy", "bollback_multinomial",
                       "cg_subs", "ct_subs", "freq_A", "freq_C", "freq_G", "freq_T", "gt_subs", "GTR+I+G_entropy_bl",
                       "GTR+I+G_entropy_diam", "GTR+I+G_entropy_diam_cnt", "GTR+I+G_frac_cherries", "GTR+I+G_gamma",
                       "GTR+I+G_logL", "GTR+I+G_max_bl", "GTR+I+G_max_diam", "GTR+I+G_max_diam_cnt", "GTR+I+G_mean_bl",
                       "GTR+I+G_mean_diam", "GTR+I+G_mean_diam_cnt", "GTR+I+G_min_bl", "GTR+I+G_min_diam",
                       "GTR+I+G_mu_rate", "GTR+I+G_parsimony", "GTR+I+G_pInv", "GTR+I+G_std_bl", "GTR+I+G_std_diam",
                       "GTR+I+G_std_diam_cnt", "GTR+I+G_tree_size", "n_unique_sites", "nchars", "ntaxa",
                       "pinv_sites_100p", "rmsa_bollback_multinomial", "rmsa_aln_entropy", "rmsa_pinv_sites_100p",
                       "rmsa_n_unique_sites", "sop_score", "transition_avg", "transversion_avg",
                       "GTR+I+G_stemminess85_idx", "GTR+I+G_stemminess90_idx",
                       "model_I", "model_G", "model_F", "model_matrix"]

FEATURE_NAMES_MAPPING = {'ac_subs': 'A-C substitution rate', 'ag_subs': 'A-G substitution rate',
                         'aln_entropy': 'MSA entropy', 'at_subs': 'A-T substitution rate',
                         'base_freqs_entropy': 'base frequencies entropy',
                         'bollback_multinomial': 'Multinomial test statistic', 'cg_subs': 'C-G substitution rate',
                         'ct_subs': 'C-T substitution rate', 'freq_A': 'A frequency', 'freq_C': 'C frequency',
                         'freq_G': 'G frequency', 'freq_T': 'T frequency', 'gt_subs': 'G-T substitution rate',
                         'GTR+I+G_entropy_bl': 'Branch lengths (entropy)',
                         'GTR+I+G_entropy_diam': 'p-distances (entropy)',
                         'GTR+I+G_entropy_diam_cnt': '# Branches in p-distance (entropy)',
                         'GTR+I+G_frac_cherries': 'Fraction of cherries', 'GTR+I+G_gamma': 'Gamma parameter',
                         'GTR+I+G_logL': 'logL', 'GTR+I+G_max_bl': 'Branch lengths (max)',
                         'GTR+I+G_max_diam': 'p-distances (max)',
                         'GTR+I+G_max_diam_cnt': '# Branches in p-distance (max)',
                         'GTR+I+G_mean_bl': 'Branch lengths (avg)', 'GTR+I+G_mean_diam': 'p-distances (avg)',
                         'GTR+I+G_mean_diam_cnt': '# Branches in p-distance (avg)',
                         'GTR+I+G_min_bl': 'Branch lengths (min)', 'GTR+I+G_min_diam': 'p-distances (min)',
                         'GTR+I+G_mu_rate': '# Substitutions per unit time', 'GTR+I+G_parsimony': 'Parsimony',
                         'GTR+I+G_pInv': 'pInv parameter', 'GTR+I+G_std_bl': 'Branch lengths (std)',
                         'GTR+I+G_std_diam': 'p-distances (std)',
                         'GTR+I+G_std_diam_cnt': '# Branches in p-distance (std)',
                         'GTR+I+G_tree_size': 'Total branch lengths', 'n_unique_sites': '# Different site-patterns',
                         'nchars': '# MSA sites', 'ntaxa': '# MSA sequences',
                         'pinv_sites_100p': '% Fully conserved sites',
                         'rmsa_bollback_multinomial': 'Subgroup MSA - multinomial test statistic',
                         'rmsa_aln_entropy': 'Subgroup MSA - MSA entropy',
                         'rmsa_pinv_sites_100p': 'Subgroup MSA - %fully conserved sites',
                         'rmsa_n_unique_sites': 'Subgroup MSA - # Different site-patterns', 'sop_score': 'SOP score',
                         'transition_avg': 'Transitions (avg)', 'transversion_avg': 'Transversions (avg)',
                         'GTR+I+G_stemminess85_idx': 'Cumulative stemminess index',
                         'GTR+I+G_stemminess90_idx': 'Noncumulative stemminess index'}
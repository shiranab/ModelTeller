import traceback

from definitions import *
import compute_features
import tree_functions
from utils import *
import phyml

def validate_input(msa_file, user_tree_file):
	"""
	:param msa_file: the path to an MSA file, one of biopython's formats
	:param user_tree_file: (optional) the path to a user tree file, if fixed tree was desired
	:return: a biopython object of the msa and an ete3 object of the tree if exists
	"""

	# identify format and retrieve all MSAs
	for aln_format in ["clustal", "emboss", "fasta", "fasta-m10", "ig", "maf", "mauve", "nexus", "phylip-relaxed", "phylip-sequential", "stockholm"]:
		try:
			msa_obj = AlignIO.read(msa_file, format=aln_format)
			logger.info("The MSA file is format: " + aln_format)
			break
		except Exception:
			msa_obj = None
	if msa_obj is None:
		logger.error("Error occured: the input file is not a valid alignmnet in a supported format.\n"
		             "Please verify that all sequences are at the same length and that the input format is correct.")

	# validate MSA characters
	msa_info = AlignInfo.SummaryInfo(msa_obj)
	aln_letters = msa_info._get_all_letters()
	for let in aln_letters:
		if not (let.lower() in "acgt-"):
			logger.warning("There are characters that are not nucleotides or gaps in your input MSA.")
			break

	# validate tree file in Newick format and suits the msa
	tree_obj = None
	if user_tree_file:
		try:
			with open(user_tree_file) as fpr:
				tree_obj = tree_functions.get_newick_tree(fpr.read().strip())
		except:
			logger.error("Tree file is invalid. Please verify that it's in Newick format.")

		# assert that the tree matches the corresponding MSA
		leaves = sorted([node.name for node in tree_obj.get_leaves()])
		seq_names = sorted([rec.id for rec in msa_obj])
		if len(leaves) != len(seq_names) or (not all(x == y for x,y  in zip(seq_names,leaves))):
			logger.error("The tips of the tree and the MSA sequences names do not match")

	return msa_obj


def predict_h2o(ext_df, drf_model_path):

	try:
		h2o.init(max_mem_size="8G")
		covtype_df = h2o.H2OFrame(ext_df[FEATURES_TO_INCLUDE])
		drf_model = h2o.load_model(drf_model_path)

		ext_df["pred_Bs"] = drf_model.predict(covtype_df).as_data_frame()
	except:
		exc_type, exc_value, exc_traceback = sys.exc_info()
		traceback.print_exception(exc_type, exc_value, exc_traceback, file=sys.stdout)
	finally:
		h2o.cluster().shutdown(prompt=False)

	return


def main(msa_obj, msa_filepath, GTRIG_topology, user_tree_file):
	"""
	:param msa_obj: a biopython.AlignIO obj of the input MSA
	:param GTRIG_topology: True - compute GTR+I+G ml tree and fix the topology for ModelTeller computation, else --
	:param user_tree_file: if GTR+I+G is False, use the given topology for ModelTeller computation, if None --
	If both GTRIG_topology and user_tree_file topology are empty, compute a ml tree for a single model
	:return:
	"""
	ext_df, features_tree_file = compute_features.prepare_features_df(msa_obj, msa_filepath, GTRIG_topology, user_tree_file)
	if not GTRIG_topology:
		drf_model_path = MODELTELLER_DRF_MODEL
	else:
		drf_model_path = MODELTELLERg_DRF_MODEL

	predict_h2o(ext_df, drf_model_path)

	probs_df = ext_df.pivot(index='index', columns='model', values='pred_Bs')[ALL_PHYML_MODELS]
	ranked_df = pd.DataFrame.rank(probs_df, axis=1, method="min")
	ext_df["model_rank"] = ranked_df.stack().values

	# save features nicely
	ext_df.drop(["model_matrix", "model_F", "model_I", "model_G"], inplace=True, axis=1)
	ext_df.rename(mapper=FEATURE_NAMES_MAPPING, axis="columns", inplace=True)
	ext_df.to_csv("features_with_models_rankings.csv")

	selected_model = ext_df.loc[ext_df["model_rank"]==1, "model"].to_list()[0] # in case there multiple minimals, take the first
	logger.info("Success: ModelTeller selected model is: " + selected_model)

	logger.info("Now computing the final phylogeny... Please wait until PhyML is done.")

	fixed_tree = None
	if GTRIG_topology:
		fixed_tree = features_tree_file
	elif user_tree_file:
		fixed_tree = user_tree_file

	#reconstruct maximum-likelihood tree (with fixed topology if selected)
	_, opt_phyml_tree_filepath = phyml.run_phyml(msa_filepath, selected_model,
	                                             topology="fixed" if fixed_tree else "ml",
	                                             tree_file=fixed_tree)

	logger.info("Done. ML tree is in: " + opt_phyml_tree_filepath)


if __name__ == '__main__':
	logger = logging.getLogger('ModelTeller main script')
	init_commandline_logger(logger)

	parser = argparse.ArgumentParser(description='ModelTeller running')
	parser.add_argument('--msa_filepath', '-m', default=None,
						help='A file with an alignment or several ones of the same format.')
	parser.add_argument('--GTRIG_topology', '-g', action='store_true',
						help="Reconstruct a maximum-likelihood tree using GTR+I+G model and use this as a fixed topology.")
	parser.add_argument('--user_tree_file', '-u', default=None,
						help="Specify your tree file in Newick format and use this tree as a fixed topology.")  # if p=2
	args = parser.parse_args()

	GTRIG_topology = args.GTRIG_topology
	user_tree_file = args.user_tree_file
	msa_filepath = args.msa_filepath

	assert bool(GTRIG_topology) != bool(user_tree_file) or not bool(user_tree_file), \
		"Please select either a GTR+I+G tree or a user-defined topology. ModelTeller cannot accept both"

	msa_obj = validate_input(msa_filepath, user_tree_file)
	main(msa_obj, msa_filepath, GTRIG_topology, user_tree_file)


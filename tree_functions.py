from definitions import *
from utils import compute_entropy, lists_diff


def get_newick_tree(tree):
	"""
	:param tree: newick tree string or txt file containing one tree
	:return:	tree: a string of the tree in ete3.Tree format
	"""
	if type(tree) == str:
		if os.path.exists(tree):
			with open(tree, 'r') as tree_fpr:
				tree = tree_fpr.read().strip()
		tree = Tree(tree, format=1)
	return tree


def reroot_tree(tree, outgroup_name):
	tree = get_newick_tree(tree)
	tree.set_outgroup(tree & outgroup_name)
	return tree


def get_frac_of_cherries(tree):
	"""
	McKenzie, Andy, and Mike Steel. "Distributions of cherries for two models of trees."
	 Mathematical biosciences 164.1 (2000): 81-92.
	:param tree:
	:return:
	"""
	tree = get_newick_tree(tree)
	tree_root = tree.get_tree_root()
	leaves = list(tree_root.iter_leaves())
	cherries_cnt = 0
	for leaf1, leaf2 in itertools.combinations(leaves, 2):
		if leaf1.up is leaf2.up:
			cherries_cnt += 1
	return 2*cherries_cnt/len(leaves)


def get_leaves_branches(tree):
	"""
	:param tree:
	:return: a list of pendant edges lengths, i.e., the brnches that lead to leaves
	"""
	tree = get_newick_tree(tree)
	tree_root = tree.get_tree_root()
	leaves_bl = []
	for node in tree_root.iter_leaves():
		leaves_bl.append(node.dist)
	return leaves_bl


def get_stemminess_indexes(tree):
	"""
		:param tree: ete3 tree; not rerooting!
		:return: cumulative stemminess index Fiala and Sokal 1985
				 noncumulative stemminess index Rohlf 1990
		formula cumulative stemminess: https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1558-5646.1985.tb00398.x
		formula noncumulative stemminess: https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1558-5646.1990.tb03855.x
		"""
	subtree_blsum_dict = {}
	nodes_height_dict = {}
	stem85_index_lst = []
	stem90_index_lst = []
	for node in tree.traverse(strategy="postorder"):
		if node.is_leaf():
			subtree_blsum_dict[node] = 0
			nodes_height_dict[node] = 0
		elif node.is_root():
			continue
		else:
			subtree_blsum_dict[node] = subtree_blsum_dict[node.children[0]] + subtree_blsum_dict[node.children[1]] + \
			                           node.children[0].dist + node.children[1].dist
			nodes_height_dict[node] = max(nodes_height_dict[node.children[0]] + node.children[0].dist,
			                              nodes_height_dict[node.children[1]] + node.children[1].dist)
			stem85_index_lst.append(node.dist/(subtree_blsum_dict[node] + node.dist))
			stem90_index_lst.append(node.dist/(nodes_height_dict[node]) + node.dist)

	return np.mean(stem85_index_lst), np.mean(stem90_index_lst)


def get_branch_lengths(tree):
	"""
	:param tree: Tree node or tree file or newick tree string;
	:return: total branch lengths
	"""
	# TBL
	tree = get_newick_tree(tree)
	tree_root = tree.get_tree_root()
	branches = []
	for node in tree_root.iter_descendants(): # the root dist is 1.0, we don't want it
		branches.append(node.dist)
	return branches


def get_branch_lengths_estimates(tree):
	"""
	:param tree: Tree node or tree file or newick tree string;
	:return:
	"""
	# TBL
	branches = get_branch_lengths(tree)
	entropy = compute_entropy(branches)

	return max(branches), min(branches), np.mean(branches), np.std(branches), entropy


def get_diameters_estimates(tree_filepath, actual_bl=True):
	"""
	if not actual_bl - function changes the tree! send only filepath
	:param tree_filepath: tree file or newick tree string;
	:param actual_bl: True to sum actual dists, False for num of branches
	:return: min, max, mean, and std of tree diameters
	"""
	# tree = copy.deepcopy(get_newick_tree(tree)) # do not deepcopy! when trees are large it exceeds recursion depth
	if not actual_bl:
		assert isinstance(tree_filepath, str)
	tree = get_newick_tree(tree_filepath)
	tree_root = tree.get_tree_root()
	if not actual_bl:
		for node in tree_root.iter_descendants():
			node.dist = 1.0
	tree_diams = []
	leaves = list(tree_root.iter_leaves())
	for leaf1, leaf2 in itertools.combinations(leaves, 2):
		tree_diams.append(leaf1.get_distance(leaf2))
	entropy = compute_entropy(tree_diams)

	return max(tree_diams), min(tree_diams), np.mean(tree_diams), np.std(tree_diams), entropy


def get_internal_and_external_leaves_relative_to_subroot(tree_root, subroot):
	all_leaves = tree_root.get_leaves()
	subtree_leaves = subroot.get_leaves()
	other_leaves = lists_diff(all_leaves, subtree_leaves)
	return subtree_leaves, other_leaves


def get_largest_branch(tree):
	tree = get_newick_tree(tree)
	tree_nodes = list(tree.traverse("levelorder"))
	max_bl_node = max(tree_nodes, key=lambda node: node.dist)
	return max_bl_node

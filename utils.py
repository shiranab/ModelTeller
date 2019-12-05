from definitions import *


def init_commandline_logger(logger):
	logger.setLevel(logging.DEBUG)
	# create console handler and set level to debug
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.DEBUG)
	# create formatter
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	# add formatter to ch
	ch.setFormatter(formatter)
	# add ch to logger
	logger.addHandler(ch)


def get_avg(l):
	avg = sum(l)/len(l)
	return avg


def get_var(l):
	avg = get_avg(l)
	dist_from_avg = list(map(lambda x: (x - avg)**2, l))
	var = sum(dist_from_avg)/len(dist_from_avg)
	return var


def get_std(l):
	var = get_var(l)
	std = math.sqrt(var)
	return std


def is_file_empty(filepath):
	"""
	:param filepath:
	:return: True if filepath doesn't exist or is empty
	"""
	if os.path.exists(filepath):
		with open(filepath) as fpr:
			if re.search("\S", fpr.read()):
				return False
	return True


def compute_entropy(lst, epsilon=0.000001):
	if np.sum(lst) != 0:
		lst_norm = np.array(lst)/np.sum(lst)
	else:
		lst_norm = np.array(lst) + epsilon
	entropy = -1*sum(np.log2(lst_norm)*lst_norm)
	if np.isnan(entropy):
		lst_norm += epsilon
		entropy = -1*sum(np.log2(lst_norm)*lst_norm)
	return entropy


def lists_diff(li1, li2):
	"""
	:param li1, li2: two lists
	:return: diff between lists
	"""
	return (list(set(li1) - set(li2)))
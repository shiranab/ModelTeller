from definitions import *
from utils import *

def remove_gaps_from_sequence(seq):
	return re.sub("[^agctAGCT]+", "", seq, re.I)


def count_fully_conserved_fraction(msa):
	"""
	:param msa:
	:param thresholds: a list of percentages
	:return:  a list - for every percentage, how many sites above this conservation thresholds
	"""
	msa_length = msa.get_alignment_length()
	invariant_sites = 0
	for col_i in range(0, msa_length):
		col_gapless = remove_gaps_from_sequence(msa[:,col_i])
		# if something that is not gaps but not AGCT appears in the column
		if len(col_gapless) == 0:
			continue
		if col_gapless.count(col_gapless[0]) == len(col_gapless):
			invariant_sites += 1

	return invariant_sites/msa_length


def calculate_column_entropy(col_string):
	# column_entropy = - sum(for every nucleotide x) {count(x)*log2(Prob(nuc x in col i))}
	col_gapless = remove_gaps_from_sequence(col_string).upper()
	col_entropy = 0
	for x in ['A', 'G', 'C', 'T']:
		count_x = str.count(col_gapless, x)
		if count_x == 0:
			entropy_x = 0
		else:
			prob_x = count_x/len(col_gapless)
			entropy_x = count_x*math.log2(prob_x)
		col_entropy += entropy_x

	return -col_entropy


def get_msa_avg_entropy(msa):
	msa_length = msa.get_alignment_length()
	sum_entropy = 0
	for col_i in range(0, msa_length):
		sum_entropy += calculate_column_entropy(msa[:,col_i])

	return sum_entropy/msa_length


def calculate_bollback_multinomial(msa):
	msa_length = msa.get_alignment_length()
	cols = []
	for col_i in range(0, msa_length):
		cols.append(msa[:,col_i])
	counts = Counter(cols)

	multinomial = 0
	for k in counts:
		c = counts[k]
		multinomial += c*math.log(c)
	multinomial -= msa_length*math.log(msa_length)
	return multinomial, len(counts), len(counts)/msa_length


def infer_pairwise_substitution_matrix(seq1, seq2):
	substitution_count_dictionary = {"AC": 0, "AG": 0, "AT": 0, "CG": 0, "CT": 0, "GT": 0, "AA": 0, "GG": 0, "CC": 0,
	                                 "TT": 0, "1s": 0, "2s": 0} #1s for one space vs nucleotide, 2s for 2 spaces
	pa_length = 0
	for i in range(0, len(seq1)):
		ch1 = min(seq1[i].upper(), seq2[i].upper())
		ch2 = max(seq1[i].upper(), seq2[i].upper())
		if ch1 not in ["A", "G", "C", "T"] and ch2 not in ["A", "G", "C", "T"]:
			substitution_count_dictionary["2s"] +=1
		elif ch1 in ["A", "G", "C", "T"] and ch2 in ["A", "G", "C", "T"]: #both nucleotides
			substitution_count_dictionary[ch1+ch2] += 1
			pa_length +=1
		else: #1 space
			substitution_count_dictionary["1s"] += 1

	return substitution_count_dictionary, pa_length


def compute_pairwise_substitution_rates(seq1, seq2):
	MATCH_SCORE = 1
	MISMATCH_SCORE = -1
	GAP_SCORE = -1

	transition_rate, transversion_rate, match, mismatch, gap = 0, 0, 0, 0, 0

	substitution_count_dictionary, pa_length = infer_pairwise_substitution_matrix(seq1, seq2)
	if pa_length != 0:
		transition_rate = float(substitution_count_dictionary["AG"] + substitution_count_dictionary["CT"]) / pa_length
		transversion_rate = float(substitution_count_dictionary["AC"] + substitution_count_dictionary["AT"] +
								  substitution_count_dictionary["CG"] + substitution_count_dictionary["GT"]) / pa_length
		match = (substitution_count_dictionary["AA"]+substitution_count_dictionary["CC"]+substitution_count_dictionary["GG"]+substitution_count_dictionary["TT"])*MATCH_SCORE
		mismatch = (substitution_count_dictionary["AC"]+substitution_count_dictionary["AG"]+substitution_count_dictionary["AT"] +
				   substitution_count_dictionary["CG"]+substitution_count_dictionary["CT"]+substitution_count_dictionary["GT"]+substitution_count_dictionary["1s"])*MISMATCH_SCORE
	gap = substitution_count_dictionary["1s"]*GAP_SCORE

	return transition_rate, transversion_rate, match+mismatch+gap,\
			substitution_count_dictionary["AC"], substitution_count_dictionary["AG"], substitution_count_dictionary["AT"], \
			substitution_count_dictionary["CG"], substitution_count_dictionary["CT"], \
			substitution_count_dictionary["GT"]


def compute_base_frequencies(msa):
	freqs = []
	seqs_msa = list(msa)
	allchars = "".join([rec._seq._data for rec in seqs_msa])

	for nuc in "ACGT":
		freqs.append(len(re.findall(nuc, allchars, re.I)))
	freqs = {"freq_" + nuc : freq / sum(freqs) for nuc, freq in zip("ACGT", freqs)}
	return freqs


def calculate_substitution_rates(msa):
	seqs_msa = list(msa)
	transition_rates = []
	transversion_rates = []
	sop_score = 0
	ac_cnt = ag_cnt = at_cnt = cg_cnt = ct_cnt = gt_cnt = 0
	for i in range(0, len(seqs_msa)-1):
		for j in range(i+1,len(seqs_msa)):
			seq1 = seqs_msa[i]._seq._data
			seq2 = seqs_msa[j]._seq._data
			transition_rate, transversion_rate, pair_sop, ac, ag, at, cg, ct, gt = compute_pairwise_substitution_rates(seq1, seq2)
			transition_rates.append(transition_rate)
			transversion_rates.append(transversion_rate)
			sop_score += pair_sop
			ac_cnt, ag_cnt, at_cnt, cg_cnt, ct_cnt, gt_cnt = \
				ac_cnt + ac, ag_cnt + ag, at_cnt + at, cg_cnt + cg, ct_cnt + ct, gt_cnt + gt

	subs_sum = sum([ac_cnt, ag_cnt, at_cnt, cg_cnt, ct_cnt, gt_cnt])
	return {"transition_avg": get_avg(transition_rates),
			"transversion_avg": get_avg(transversion_rates),
			"sop_score": sop_score},\
		   {c: x/subs_sum if subs_sum != 0 else 0 for c,x in
			zip(["ac_subs", "ag_subs", "at_subs", 'cg_subs', 'ct_subs', 'gt_subs'],
				[ac_cnt, ag_cnt, at_cnt, cg_cnt, ct_cnt, gt_cnt])}


def get_msa_properties(msa):
	"""
	:param msa: bio.AlignIO format or path to msa file
	:return:
	"""
	ntaxa = len(msa)
	nchars = msa.get_alignment_length()

	return ntaxa, nchars


def remove_nonACTGU_sites(msa):
	"""
	removes sites that contain non ACGT characters (i.e., N's, gaps etc.)
	:param msa: bio.AlignIO format or path to msa file
	:return:
	"""
	nchars = msa.get_alignment_length()
	new_msa = copy.deepcopy(msa)
	for col_i in range(nchars - 1, -1, -1):
		col = msa[:, col_i]
		if not re.search("[ACGTUacgtu]", col, re.IGNORECASE):
			new_msa = new_msa[:, :col_i] + new_msa[:, col_i + 1:]
	return new_msa


def reduce_msa_to_seqs_by_name(msa, keep_names_lst):
	new_msa = []
	all_names = [rec.id for rec in list(msa)]
	for name in keep_names_lst:
		new_msa.append(msa[all_names.index(name), :])
	#remove positions that are just gaps after removal of sequences
	new_msa = remove_nonACTGU_sites(AlignIO.MultipleSeqAlignment(new_msa))
	return new_msa


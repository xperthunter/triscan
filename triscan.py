#!/usr/bin/python3

from itertools import product
import json
import math
import sys

#import biopython
import cppyy
import numpy as np
import scipy.stats

import msalib

q = [
	'A','C','D','E','F',
	'G','H','I','K','L',
	'M','N','P','Q','R',
	'S','T','V','W','Y',
	'.'
]

kd_scale = {
	"I":4.5,
	"V":4.2,
	"L":3.8,
	"F":2.8,
	"C":2.5,
	"M":1.9,
	"A":1.8,
	"G":-0.4,
	"T":-0.7,
	"S":-0.8,
	"W":-0.9,
	"Y":-1.3,
	"P":-1.6,
	"H":-3.2,
	"E":-3.5,
	"Q":-3.5,
	"D":-3.5,
	"N":-3.5,
	"K":-3.9,
	"R":-4.5
}

pI = {
	"I":6.02,
	"A":6.00,
	"L":5.98,
	"R":10.76,
	"K":9.74,
	"N":5.41,
	"M":5.74,
	"D":2.77,
	"F":5.48,
	"C":5.05,
	"P":6.30,
	"Q":5.65,
	"S":5.68,
	"E":3.22,
	"T":5.66,
	"G":5.97,
	"W":5.89,
	"H":7.59,
	"Y":5.66,
	"V":5.96
}

norm_kd = {k:v/4.5 for k,v in kd_scale.items()}
norm_pI = {k:v/10.76 for k,v in pI.items()}
q_squared = list(product(q, q))


def index_mapper(lid, seq):

	lid_info = lid.split('/')
	ends = lid_info[-1].split('-')

	beg, end = int(ends[0]), int(ends[1])

	j = 0
	mapper = {}
	for i, sym in enumerate(seq):
		if sym == '.':
			continue
		j += 1
		mapper[int(i)] = beg + j - 1

	return mapper


def one_body_potentials(col1, col2):
	cijs = []
	for ei, ej in zip(col1, col2):
		cij = 0
		if ei in norm_kd and ej in norm_kd:
			#cij = norm_kd[ei]*norm_kd[ej] + 0.1*norm_pI[ei]*norm_pI[ej]
			#cij = norm_kd[ei] + norm_kd[ej]
			cij = norm_kd[ei]*norm_kd[ej]
			cijs.append(cij)

	if len(cijs) == 0:
		return None

	cijs = np.array(cijs)

	#snr = abs(np.mean(cijs) / np.std(cijs))
	snr = np.mean(cijs)
	return snr


def m_eff(similarity_scores):
	m_eff = 0
	for v in similarity_scores.values():
		#print(v)
		m_eff += 1/v
	return m_eff


def make_ma(seqs, cutoff):
	sim_scores = dict()
	for i, seq in enumerate(seqs):
		sim_scores[i] = m_a(seq, seqs, cutoff)
	
	return sim_scores


def m_a(seq_i, seqs, cutoff):
	
	m_a = 0
	L = len(seq_i)
	max_mismatch = L * (1 - cutoff)
	for i, seq_j in enumerate(seqs):
		assert(len(seq_j) == len(seq_i))
		mismatches = 0
		similar = True
		for j, (ei, ej) in enumerate(zip(seq_i, seq_j)):
			if ei != ej:
				mismatches += 1

			if mismatches > max_mismatch:
				similar = False
				break

			if mismatches + (L-j ) < max_mismatch:
				break

		if similar:
			m_a += 1
	
	m_a = m_a
	return m_a


def f_i(msa, ma_scores, ll):
	
	meff = m_eff(ma_scores)
	fi = dict()
	
	for i in range(msa.length):
		if msa.cons[i] == '.': continue
		fi[i] = {aa: ((ll*meff) / len(q)) for aa in q}
		
		col = msa.column(i)
		
		for n, elm in enumerate(col):
			if elm not in q: continue
			fi[i][elm] += 1 / ma_scores[n]
		
		for k in fi[i].keys(): fi[i][k] = fi[i][k] * (1 / ((1+ll)*meff))
	
	return fi


def f_ij(msa, ma_scores, ll):
	
	meff = m_eff(ma_scores)
	ll_rescaled = (ll**2 / (1 + 2* ll)) * meff
	fij = dict()
	
	for i in range(msa.length):
		if msa.cons[i] == '.': continue
		for j in range(i+1, msa.length):
			if msa.cons[j] == '.': continue
			fij[(i,j)] = {pair: (ll*meff) / len(q_squared) for pair in q_squared}

			#fij[(i,j)] = dict()

			col_i = msa.column(i)
			col_j = msa.column(j)
		
			for n, (ei, ej) in enumerate(zip(col_i, col_j)):
				if ei not in q or ej not in q: continue
				#if (ei,ej) not in fij[(i,j)]: fij[(i,j)][(ei,ej)] = ll_rescaled / (len(q_squared))
				fij[(i,j)][(ei,ej)] += 1 / ma_scores[n]
			
			for k in fij[(i,j)].keys(): fij[(i,j)][k] = fij[(i,j)][k] / (meff + ll_rescaled)
	
	return fij
	

def mutual_information_ij(f1, f2):
	
	mutual = dict()
	
	for (i,j) in f2.keys():
		info = 0
		
		#for (ai, aj) in f2[(i,j)].keys():
		for (ai, aj) in q_squared:
			info += f2[(i,j)][(ai,aj)] * math.log2(f2[(i,j)][(ai,aj)] / (f1[i][ai] * f1[j][aj]))
		
		mutual[(i,j)] = info
	
	return mutual


counter = 0
for msa in msalib.read_stockholm(sys.argv[1]):
	print(msa.identifier)
	print(msa.accession)
	print(f"num of seqs: {len(msa.seqs)}")
	print(f"msa width: {len(msa.seqs[0])}")

	# Compute Mutual Information
	max_similarity = 0.8
	#gap_count = msa.cons.count('.')
	#L = len(msa.seqs[0]) - gap_count
	#mismatch_max = math.ceil(L * (1 - max_similarity))
	results = np.zeros(msa.depth, dtype=np.int32)
	lens = np.array(msa.lens, dtype=np.intc)
	cppyy.gbl.get_ma(msa.seqs, lens, msa.depth, max_similarity, results)
	#cppyy.gbl.get_ma(msa.seqs, msa.lens, msa.depth, max_similarity, results)
	ma = {k:v for k, v in enumerate(results)}
	meff = m_eff(ma)

	print(f"m_eff: {meff:.2f} # of seqs: {len(msa.seqs)}")
	print(f"msa width: {len(msa.seqs[0])}")

	# Compute Single Point Correlations -- amino acid preferences per column
	single_corr = f_i(msa, ma, 1.25) # msa obj, seq. sim scaling, pseudo-counts
	
	# Compute Two Point Correlations -- how often does aa_A and aa_B appear in columns i,j
	two_corr = f_ij(msa, ma, 1.25)
	
	# Mutual information -- info. gain when assuming two columns occur together, versus separate
	mutual = mutual_information_ij(single_corr, two_corr)
	
	# Display mutual
	for k,v in sorted(mutual.items(), key=lambda x: x[1], reverse=False):
		print(k,v, msa.cons[k[0]], msa.cons[k[1]])

	sys.exit()
	
	
	# gather column frequencies
	col_aafreqs = []
	aafreqs     = []
	col_freqs   = []
	for i in range(msa.length):
		if msa.cons[i] == '.': continue
		c = msa.column(i)
		aafreqs = {}
		for aa in c:
			if aa not in aafreqs: aafreqs[aa] = 0
			
			aafreqs[aa] += 1 / msa.depth
		
		col_aafreqs.append(aafreqs)
		aafreqs = {int(ii):v for ii, (k,v) in enumerate(sorted(aafreqs.items(), key = lambda x: x[1], reverse=True))}
		col_freqs.append(aafreqs)
	
	# compute entropy per column
	entropys = []
	s = dict()
	for i, c in enumerate(col_freqs):
		entropy = 0
		for f in c.values():
			entropy += f*math.log2(f)
		
		entropys.append(-1.0*entropy)
		s[i] = -1.0*entropy
	
	entropys = np.array(entropys)
	print(f"avg S: {np.mean(entropys):.2f} 0.1-q: {np.quantile(entropys,0.1):.2f} 0.9-q: {np.quantile(entropys,0.9):.2f}")

	# compute frequency distr. dis. on pairs of columns
	distances = {}
	cij_scores = {}
	mij = {}
	lower_b = np.quantile(entropys,0.25)
	upper_b = np.quantile(entropys,0.75)
	for i, ci in enumerate(col_freqs):
		
		if s[i] < lower_b or s[i] > upper_b:
			continue
		
		if msa.cons[i] == '.':
			continue

		for j in range(i+5, msa.length):
			if msa.cons[j] == '.': continue
			if j not in s: continue
			if s[j] < lower_b or s[j] > upper_b:
				continue

			coupling = one_body_potentials(msa.column(i), msa.column(j))
			#print(coupling)
			if coupling is None:
				continue

			cij_scores[(i,j)] = coupling
			
			dis = 0
			cj = col_freqs[j]
			for ii in range(21):
				if ii in ci and ii in cj:
					dis += abs(ci[ii] - cj[ii])
				elif ii in ci and ii not in cj:
					dis += ci[ii]
				elif ii not in ci and ii in cj:
					dis += cj[ii]
				elif ii not in ci and ii not in cj:
					break
				else:
					print("some error")
					sys.exit()
			
			distances[(i,j)] = dis
			mij[(i,j)] = mutual[(i,j)]

			"""
			if coupling > 7:
				print()
				print(coupling)
				print(i,j)
				print(dis)
				print(mutual[(i,j)])
				print(' '.join(msa.column(i)[:50]))
				print(' '.join(msa.column(j)[:50]))
				print()
			"""

	
	avg_mutual = np.mean(np.array(list(mij.values())))
	std_mutual = np.std(np.array(list(mij.values())))
	mutual_norm = {k:((v - avg_mutual)/std_mutual) for k,v in mij.items()}

	norm_distances = {k:(2.0 - v) for k,v in distances.items()}
	mean_dis = np.mean(np.array(list(norm_distances.values())))
	std_dis  = np.std(np.array(list(norm_distances.values())))
	norm_distances = {k:((v - mean_dis) / std_dis) for k,v in norm_distances.items()}

	mean_cij  = np.mean(np.array(list(cij_scores.values())))
	std_cij   = np.std(np.array(list(cij_scores.values())))
	norm_cijs = {k:((v - mean_cij) / std_cij) for k,v in cij_scores.items()}
	#print(mean_cij, std_cij)
	#sys.exit()
	sum_zscores = {k:(mutual_norm[k] + norm_distances[k] + norm_cijs[k]) for k in norm_cijs.keys()}

	#imapper = index_mapper(msa.lids[0], msa.seqs[0])
	for k,v in sorted(sum_zscores.items(), key = lambda x: x[1]):
		#try:
		#l = (msa.resindices[0][int(k[0])], msa.resindices[0][int(k[1])])
		print(f"pair: {k} sum of z-scores: {v:.4f} m_ij: {mutual_norm[k]:.4f} d_ij: {norm_distances[k]:.4f} c_ij: {norm_cijs[k]:.4f}")
		#except:
		#	continue

	mutual_norm    = {k:v for k,v in sorted(mutual_norm.items(), key=lambda x: x[1])}
	norm_distances = {k:norm_distances[k] for k in mutual_norm.keys()}
	norm_cijs      = {k:norm_cijs[k] for k in mutual_norm.keys()}

	data = np.array(
		[
			np.array(list(mutual_norm.values())),
			np.array(list(norm_distances.values())),
			np.array(list(norm_cijs.values()))
		]
	)

	corrs = np.corrcoef(data)
	print(corrs)

	print(msa.identifier)
	print(msa.accession)
	print(len(msa.seqs[0]))


	sys.exit()

	"""
		mins = []
	mutuals = []
	for ii, (k,v) in enumerate(sorted(distances.items(), key=lambda x: x[1], reverse=True)):
		print(f'pair: {k} dis: {v:.4f} m_ij: {mutual[k]:.4f} max: {max_mutual:.4f} min: {min_mutual:.4f}')

		mins.append(v)
		mutuals.append(mutual[k])

		#if ii == 15: break



	cijs = []
	for e22, e133 in zip(msa.column(22), msa.column(133)):
		cij = None
		if e22 in norm_kd and e133 in norm_kd:
			cij = -1.0*norm_kd[e22]*norm_kd[e133] + 0.1*norm_pI[e22]*norm_pI[e133]
			cijs.append(cij)
		print(e22, e133, cij)
	
	cijs = np.array(cijs)
	cijbar = np.mean(cijs)
	cijstd = np.std(cijs)
	snr = cijbar/cijstd

	print(f"cij info: mean: {cijbar:.4f} std: {cijstd:.4f} snr: {snr:.4f}")
	#sys.exit()
	#print(json.dumps({k:v for k,v in sorted(col_aafreqs[22].items(), key = lambda x: x[1])},indent=2))
	#print()
	#print(json.dumps({k:v for k,v in sorted(col_aafreqs[133].items(), key = lambda x: x[1])},indent=2))
	#print(dict(k:v for k,v in sorted(col_aafreqs[111].items(), lambda x: x[1]))
	
	
	print(len(distances))
	print(msa.length ** 2)

	for i in range(msa.length):
		for j in range(i+1, msa.length):



	counter += 1
	sys.exit()
	if counter == 10: sys.exit()
	"""

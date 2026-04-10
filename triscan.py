#!/usr/bin/python3

from itertools import product
import json
import math
import sys

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


def one_body_potentials(col1, col2):

	cijs = []

	for ei, ej in zip(col1, col2):
		cij = 0
	if e22 in norm_kd and e133 in norm_kd:
			cij = -1.0*norm_kd[e22]*norm_kd[e133] + 0.1*norm_pI[e22]*norm_pI[e133]


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
	for i, seq_j in enumerate(seqs):
		assert(len(seq_j) == len(seq_i))
		similarity = 0
		for ei, ej in zip(seq_i, seq_j):
			if ei == ej:
				similarity += 1 / len(seq_i)
		
		#print(similarity)
		if similarity > cutoff: m_a += 1
	
	m_a = m_a
	return m_a

def f_i(msa, ma_scores, ll):
	
	meff = m_eff(ma_scores)
	fi = dict()
	
	for i in range(msa.length):
		fi[i] = {aa: ((ll*meff) / len(q)) for aa in q}
		
		col = msa.column(i)
		
		for n, elm in enumerate(col):
			fi[i][elm] += 1 / ma_scores[n]
		
		for k in fi[i].keys(): fi[i][k] = fi[i][k] * (1 / ((1+ll)*meff))
	
	return fi

def f_ij(msa, ma_scores, ll):
	
	meff = m_eff(ma_scores)
	fij = dict()
	
	for i in range(msa.length):
		for j in range(i+1, msa.length):
			fij[(i,j)] = {pair: (ll*meff) / len(q_squared) for pair in q_squared}
		
			col_i = msa.column(i)
			col_j = msa.column(j)
		
			for n, (ei, ej) in enumerate(zip(col_i, col_j)):
				fij[(i,j)][(ei,ej)] += 1 / ma_scores[n]
			
			for k in fij[(i,j)].keys(): fij[(i,j)][k] = fij[(i,j)][k] * (1 / ((1+ll)*meff))
	
	return fij	
	

def mutual_information_ij(f1, f2):
	
	mutual = dict()
	
	for (i,j) in f2.keys():
		info = 0
		
		for (ai, aj) in q_squared:
			info += f2[(i,j)][(ai,aj)] * math.log2(f2[(i,j)][(ai,aj)] / (f1[i][ai] * f1[j][aj]))
		
		mutual[(i,j)] = info
	
	return mutual

counter = 0
for msa in msalib.read_stockholm(sys.argv[1]):
	print(msa.identifier)
	print(msa.accession)
	
	col_freqs = []
	col_aafreqs = []
	
	max_similarity = 0.7
	ma = make_ma(msa.seqs, max_similarity)
	meff = m_eff(ma)
	
	print(meff)
	print(len(msa.seqs))
	
	single_corr = f_i(msa, ma, 1.25)
	#print(json.dumps(single_point_correlations,indent=2))
	
	two_corr = f_ij(msa, ma, 1.25)
	
	mutual = mutual_information_ij(single_corr, two_corr)
	#for k,v in sorted(mutual.items(), key=lambda x: x[1]):
	#	print(k,v)
	#sys.exit()
	
	for i in range(msa.length):
		c = msa.column(i)
		aafreqs = {}
		for aa in c:
			if aa not in aafreqs: aafreqs[aa] = 0
			
			aafreqs[aa] += 1 / msa.depth
		
		col_aafreqs.append(aafreqs)
		aafreqs = {int(ii):v for ii, (k,v) in enumerate(sorted(aafreqs.items(), key = lambda x: x[1], reverse=True))}
		col_freqs.append(aafreqs)
	
	entropys = []
	s = dict()
	for i, c in enumerate(col_freqs):
		entropy = 0
		for f in c.values():
			entropy += f*math.log2(f)
		
		entropys.append(-1.0*entropy)
		s[i] = -1.0*entropy
	
	entropys = np.array(entropys)
	print(f"\tavg S: {np.mean(entropys):.2f} 0.1-q: {np.quantile(entropys,0.1):.2f} 0.9-q: {np.quantile(entropys,0.9):.2f}")
	
	distances = {}
	cij_scores = {}
	lower_b = np.quantile(entropys,0.25)
	upper_b = np.quantile(entropys,0.75)
	for i, ci in enumerate(col_freqs):
		
		if s[i] < lower_b or s[i] > upper_b:
			continue
		
		for j in range(i+1, msa.length):
			
			if s[j] < lower_b or s[j] > upper_b:
				continue
			
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

			for e1, e2

	
	mutuals_sorted = list(sorted(mutual.values()))
	max_mutual = mutuals_sorted[-1]
	min_mutual = mutuals_sorted[0]
	
	mins = []
	mutuals = []
	for ii, (k,v) in enumerate(sorted(distances.items(), key=lambda x: x[1], reverse=True)):
		print(f'pair: {k} dis: {v:.4f} m_ij: {mutual[k]:.4f} max: {max_mutual:.4f} min: {min_mutual:.4f}')

		mins.append(v)
		mutuals.append(mutual[k])

		#if ii == 15: break
	
	r, p = scipy.stats.pearsonr(mins, mutuals)
	print(f"corr: {r:.4f} p: {p:1.2E}")

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

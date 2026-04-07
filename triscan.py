#!/usr/bin/python3

import json
import math
import sys

import numpy as np

import msalib

for msa in msalib.read_stockholm(sys.argv[1]):
	print(msa.identifier)
	print(msa.accession)
	
	col_freqs = []
	col_aafreqs = []
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
	print(f"\tavg S: {np.mean(entropys):.2f} 0.25-q: {np.quantile(entropys,0.25):.2f} 0.75-q: {np.quantile(entropys,0.75):.2f}")
	
	distances = {}
	lower_b = np.quantile(entropys,0.30)
	upper_b = np.quantile(entropys,0.70)
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
	
	for k,v in sorted(distances.items(), key=lambda x: x[1], reverse=True):
		print(f'pair: {k} dis: {v:.4f}')
	
	print(json.dumps({k:v for k,v in sorted(col_aafreqs[74].items(), key = lambda x: x[1])},indent=2))
	print()
	print(json.dumps({k:v for k,v in sorted(col_aafreqs[75].items(), key = lambda x: x[1])},indent=2))
	#print(dict(k:v for k,v in sorted(col_aafreqs[111].items(), lambda x: x[1]))
	
	
	print(len(distances))
	print(msa.depth ** 2)
	
	sys.exit()
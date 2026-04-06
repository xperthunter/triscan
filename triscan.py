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
	for i in range(msa.length):
		c = msa.column(i)
		aafreqs = {}
		for aa in c:
			if aa not in aafreqs: aafreqs[aa] = 0
			
			aafreqs[aa] += 1 / msa.depth
		
		aafreqs = {int(ii):v for ii, (k,v) in enumerate(sorted(aafreqs.items(), key = lambda x: x[1], reverse=True))}
		#print(json.dumps(aafreqs,indent=2))
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
	print(f"avg S: {np.mean(entropys):.2f} 0.25-q: {np.quantile(entropys,0.25):.2f} 0.75-q: {np.quantile(entropys,0.75):.2f}")
	
	distances = {}
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
	
	for k,v in sorted(distances.items(), key=lambda x: x[1], reverse=True):
		if k[0] != 22: continue
		print(k,v)
	
	print(len(distances))
	print(msa.depth ** 2)
	
	for e23, e33 in zip(msa.column(22), msa.column(133)):
		print(e23, e33)
	
	print(msa.length)
	sys.exit()
			
			
			
			
			
		
		
	
	
	"""
	for seq in msa.seqs:
		print(seq)
		
	c = msa.column(11)
	print(c)
	
	for i in range(msa.length):
		for j in range(i, msa.length):
			ci = msa.column(i)
			cj = msa.column(j)
			
			for a,b in zip(ci, cj):
	"""
	sys.exit()
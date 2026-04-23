import gzip
import math
import statistics
import sys


AA = { # Amino Acid Frequencies (from Pfam 38.1)
	'A': 0.0846,
	'C': 0.0140,
	'D': 0.0563,
	'E': 0.0647,
	'F': 0.0421,
	'G': 0.0680,
	'H': 0.0221,
	'I': 0.0585,
	'K': 0.0527,
	'L': 0.1022,
	'M': 0.0209,
	'N': 0.0393,
	'P': 0.0442,
	'Q': 0.0377,
	'R': 0.0569,
	'S': 0.0647,
	'T': 0.0541,
	'V': 0.0702,
	'W': 0.0142,
	'Y': 0.0325,
}


B62 = { # BLOSUM62 scoring matrix
	'A':{'A':4,'R':-1,'N':-2,'D':-2,'C':0,'Q':-1,'E':-1,'G':0,'H':-2,'I':-1,'L':-1,'K':-1,'M':-1,'F':-2,'P':-1,'S':1,'T':0,'W':-3,'Y':-2,'V':0},
	'R':{'A':-1,'R':5,'N':0,'D':-2,'C':-3,'Q':1,'E':0,'G':-2,'H':0,'I':-3,'L':-2,'K':2,'M':-1,'F':-3,'P':-2,'S':-1,'T':-1,'W':-3,'Y':-2,'V':-3},
	'N':{'A':-2,'R':0,'N':6,'D':1,'C':-3,'Q':0,'E':0,'G':0,'H':1,'I':-3,'L':-3,'K':0,'M':-2,'F':-3,'P':-2,'S':1,'T':0,'W':-4,'Y':-2,'V':-3},
	'D':{'A':-2,'R':-2,'N':1,'D':6,'C':-3,'Q':0,'E':2,'G':-1,'H':-1,'I':-3,'L':-4,'K':-1,'M':-3,'F':-3,'P':-1,'S':0,'T':-1,'W':-4,'Y':-3,'V':-3},
	'C':{'A':0,'R':-3,'N':-3,'D':-3,'C':9,'Q':-3,'E':-4,'G':-3,'H':-3,'I':-1,'L':-1,'K':-3,'M':-1,'F':-2,'P':-3,'S':-1,'T':-1,'W':-2,'Y':-2,'V':-1},
	'Q':{'A':-1,'R':1,'N':0,'D':0,'C':-3,'Q':5,'E':2,'G':-2,'H':0,'I':-3,'L':-2,'K':1,'M':0,'F':-3,'P':-1,'S':0,'T':-1,'W':-2,'Y':-1,'V':-2},
	'E':{'A':-1,'R':0,'N':0,'D':2,'C':-4,'Q':2,'E':5,'G':-2,'H':0,'I':-3,'L':-3,'K':1,'M':-2,'F':-3,'P':-1,'S':0,'T':-1,'W':-3,'Y':-2,'V':-2},
	'G':{'A':0,'R':-2,'N':0,'D':-1,'C':-3,'Q':-2,'E':-2,'G':6,'H':-2,'I':-4,'L':-4,'K':-2,'M':-3,'F':-3,'P':-2,'S':0,'T':-2,'W':-2,'Y':-3,'V':-3},
	'H':{'A':-2,'R':0,'N':1,'D':-1,'C':-3,'Q':0,'E':0,'G':-2,'H':8,'I':-3,'L':-3,'K':-1,'M':-2,'F':-1,'P':-2,'S':-1,'T':-2,'W':-2,'Y':2,'V':-3},
	'I':{'A':-1,'R':-3,'N':-3,'D':-3,'C':-1,'Q':-3,'E':-3,'G':-4,'H':-3,'I':4,'L':2,'K':-3,'M':1,'F':0,'P':-3,'S':-2,'T':-1,'W':-3,'Y':-1,'V':3},
	'L':{'A':-1,'R':-2,'N':-3,'D':-4,'C':-1,'Q':-2,'E':-3,'G':-4,'H':-3,'I':2,'L':4,'K':-2,'M':2,'F':0,'P':-3,'S':-2,'T':-1,'W':-2,'Y':-1,'V':1},
	'K':{'A':-1,'R':2,'N':0,'D':-1,'C':-3,'Q':1,'E':1,'G':-2,'H':-1,'I':-3,'L':-2,'K':5,'M':-1,'F':-3,'P':-1,'S':0,'T':-1,'W':-3,'Y':-2,'V':-2},
	'M':{'A':-1,'R':-1,'N':-2,'D':-3,'C':-1,'Q':0,'E':-2,'G':-3,'H':-2,'I':1,'L':2,'K':-1,'M':5,'F':0,'P':-2,'S':-1,'T':-1,'W':-1,'Y':-1,'V':1},
	'F':{'A':-2,'R':-3,'N':-3,'D':-3,'C':-2,'Q':-3,'E':-3,'G':-3,'H':-1,'I':0,'L':0,'K':-3,'M':0,'F':6,'P':-4,'S':-2,'T':-2,'W':1,'Y':3,'V':-1},
	'P':{'A':-1,'R':-2,'N':-2,'D':-1,'C':-3,'Q':-1,'E':-1,'G':-2,'H':-2,'I':-3,'L':-3,'K':-1,'M':-2,'F':-4,'P':7,'S':-1,'T':-1,'W':-4,'Y':-3,'V':-2},
	'S':{'A':1,'R':-1,'N':1,'D':0,'C':-1,'Q':0,'E':0,'G':0,'H':-1,'I':-2,'L':-2,'K':0,'M':-1,'F':-2,'P':-1,'S':4,'T':1,'W':-3,'Y':-2,'V':-2},
	'T':{'A':0,'R':-1,'N':0,'D':-1,'C':-1,'Q':-1,'E':-1,'G':-2,'H':-2,'I':-1,'L':-1,'K':-1,'M':-1,'F':-2,'P':-1,'S':1,'T':5,'W':-2,'Y':-2,'V':0},
	'W':{'A':-3,'R':-3,'N':-4,'D':-4,'C':-2,'Q':-2,'E':-3,'G':-2,'H':-2,'I':-3,'L':-2,'K':-3,'M':-1,'F':1,'P':-4,'S':-3,'T':-2,'W':11,'Y':2,'V':-3},
	'Y':{'A':-2,'R':-2,'N':-2,'D':-3,'C':-2,'Q':-1,'E':-2,'G':-3,'H':2,'I':-1,'L':-1,'K':-2,'M':-1,'F':3,'P':-3,'S':-2,'T':-2,'W':2,'Y':7,'V':-1},
	'V':{'A':0,'R':-3,'N':-3,'D':-3,'C':-1,'Q':-2,'E':-2,'G':-3,'H':-3,'I':3,'L':1,'K':-2,'M':1,'F':-1,'P':-2,'S':-2,'T':0,'W':-3,'Y':-1,'V':4},
}

def get_fp(filename):
	"""Returns a file pointer for reading based on file name"""
	if   filename.endswith('.gz'): return gzip.open(filename, 'rt')
	elif filename == '-':          return sys.stdin
	else:                          return open(filename)


def read_fasta(filename):
	"""Simple fasta file iterator: yields defline, seq"""
	name = None
	seqs = []
	fp = get_fp(filename)
	while True:
		line = fp.readline()
		if line == '': break
		line = line.rstrip()
		if line.startswith('>'):
			if len(seqs) > 0:
				seq = ''.join(seqs)
				yield name, seq
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:
			seqs.append(line)
	yield name, ''.join(seqs)
	fp.close()


class MSA:
	"""Simple class for MSAs (coming from Pfam seeds)"""
	def __init__(self, lines):
		self.identifier = None # release identifier
		self.accession = None  # stable identifier
		self.description = []  # to be joined later
		self.type = None       # Domain, Family, Repeat, Colied-Coil, Disorderd, Motif
		self.length = None
		self.depth = None
		self.uids = []         # short identifier
		self.lids = []         # long identifier (shows subsequence)
		self.sub_seqs = []
		self.seqs = []
		self.lens = []
		self.uid_index = {}
		self.lid_index = {}
		self.cons = None       # Probably not useful
		self.resindices = dict()

		if not lines[0].startswith('# STOCKHOLM 1.0'):
			sys.exit('MSA constructor error')

		for line in lines[1:]:
			if line.startswith('#=GF'):
				tag = line[5:7]
				val = line[10:]
				match tag:
					case 'ID': self.identifier = val
					case 'AC': self.accession = val
					case 'DE': self.description.append(val)
					case 'TP': self.type = val
			elif line.startswith('#=GS'):
				foo, lid, ac, uid = line.split()
				self.uids.append(uid)
				self.lids.append(lid)
				ll = lid.split('/')[-1].split('-')
				ll = [int(l) for l in ll]
				self.sub_seqs.append(tuple(ll))
				seqlen = ll[1] - ll[0]
				assert(seqlen > 0)
				self.lens.append(seqlen)
			elif line.startswith('#=GC'):
				self.cons = line.split()[2]
			elif line.startswith('#'):
				if line.startswith('#=GR'): pass
				else: sys.exit('unrecognized line')
			else:
				lid, seq = line.split()
				self.seqs.append(seq.upper())

		self.description = ' '.join(self.description)
		self.length = len(self.seqs[0])
		self.depth = len(self.seqs)
		
		for k, (ends, lid, uid, seq) in enumerate(zip(self.sub_seqs, self.lids, self.uids, self.seqs)):
			self.uid_index[uid] = seq
			self.lid_index[lid] = seq

			beg, end = ends
			
			j = 0
			self.resindices[k] = dict()
			for i, sym in enumerate(seq):
				if sym.upper() not in AA: continue

				j += 1
				self.resindices[k][i] = beg + j - 1


	def write(self, fp):
		print('# STOCKHOLM 1.0', file=fp)
		print('#=GF ID', self.identifier, file=fp)
		print('#=GF AC', self.accession, file=fp)
		print('#=GF DE', self.description, file=fp)
		print('#=GF TP', self.type, file=fp)
		for lid, uid in zip(self.lids, self.uids):
			print('#=GS', lid, 'AC', uid, file=fp)
		for lid, seq in zip(self.lids, self.seqs):
			print(lid, seq, sep='\t', file=fp)
		print('//', file=fp)


	def column(self, n):
		letters = []
		for seq in self.seqs:
			letters.append(seq[n])
		return ''.join(letters)


def read_stockholm(filename):
	"""Stockholm file iterator: yields msa objects"""
	fp = get_fp(filename)
	lines = []
	for line in fp:
		if line.startswith('//'):
			yield MSA(lines)
			lines = []
		else:
			lines.append(line.rstrip())
	fp.close()


def column_discretizer(col):
	gap_count = col.count('-')
	if gap_count / len(col) > 0.5: return 9 # gap code

	# get average score among all pairwise comparisons
	scores = []
	for i, a in enumerate(col):
		if a not in B62: continue
		for b in col[i+1:]:
			if b not in B62: continue
			scores.append(B62[a][b])

	x = round(statistics.mean(scores))
	if x < -2: x = -2
	if x > 6: x = 6
	return x + 2


#############
## C area ##
###########

import cppyy
import numpy as np

cppyy.cppdef("""
extern "C" {
#include <string.h>
#include <stdio.h>

void get_ma(char **seqs, int *lens, int size, float max_similarity, int *results) {
	int max_mismatch = 0;
	for (int i = 0; i < size; i++) {
		char *s1 = seqs[i];
		int slen = strlen(s1);
		results[i] = 1;
		max_mismatch = (int) (lens[i] * (1 - max_similarity) + 1.0);
		for (int j = 0; j < size; j++) {
			if (i == j) continue;
			char *s2 = seqs[j];
			int mismatch = 0;
			for (int k = 0; k < slen; k++){
				if (s1[k] == '.') continue;
				if (s1[k] == '-') continue;
				if (s2[k] == '.') continue;
				if (s2[k] == '-') continue;
				if (s1[k] != s2[k]) mismatch++;
				if (mismatch >= max_mismatch) break;
			}
			if (mismatch < max_mismatch) results[i] = results[i] + 1;
		}
	}
}}""")



cppyy.cppdef("""
extern "C" {
#include <string.h>
#include <stdio.h>

// this could be 2x faster by mirroring the half matrix
// could also use a thread-pool
// could also be outside python FFS
void get_similarities(char **seqs, int size, float *results) {
	for (int i = 0; i < size-1; i++) {
		double sum = 0;
		for (int j = i+1; j < size; j++) {
			char *s1 = seqs[i];
			char *s2 = seqs[j];
			int slen = strlen(s1);
			int match = 0;
			int total = 0;
			for (int k = 0; k < slen; k++) {
				if (s1[k] == '.') continue;
				if (s1[k] == '-') continue;
				if (s2[k] == '.') continue;
				if (s2[k] == '-') continue;
				if (s1[k] == s2[k]) match++;
				total++;
			}
			sum += (double)match/(double)total;
		}
		results[i] = sum/(double)(size -1);
	}
}}""")

if __name__ == '__main__':
	import argparse
	import math
	
	parser = argparse.ArgumentParser()
	parser.add_argument('stockholm')
	parser.add_argument('--verbose', action='store_true');
	arg = parser.parse_args();
	for msa in read_stockholm(arg.stockholm):
		print(msa.accession, msa.depth, msa.description)
		# results = np.zeros(msa.depth, dtype=np.float32)
		# cppyy.gbl.get_similarities(msa.seqs, msa.depth, results)
		# for uid, dis in zip(msa.uids, results):
		# 	if arg.verbose: print(uid, dis)
		
		# sequence identity clustering
		id_threshold = float(0.7)
		results = np.zeros(msa.depth, dtype=np.intc)
		#L = len(msa.seqs[0]) - msa.cons.count('.')
		#lens = [L for i in range(msa.depth)]
		#lens = np.array(lens, dtype=np.intc)
		lens = np.array(msa.lens, dtype=np.intc)
		cppyy.gbl.get_ma(msa.seqs, lens, msa.depth, id_threshold, results)
		#sys.exit()
		M_eff = 0
		for seq, sims in zip(msa.seqs, results):
			#print("\n"+seq)
			#print(sims)
			#print("\n")
			M_eff += 1/sims
		
		print(f"M_eff: {M_eff:.2f}")
		sys.exit()

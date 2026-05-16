#!/usr/bin/python3


import argparse
from itertools import product
import json
import math
import sys

from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1 as three_to_one
import cppyy
import numpy as np
import scipy.stats

from msalib import MSA
from msalib import read_stockholm


parser = argparse.ArgumentParser(description="measure mutual information with new psuedo counts")

parser.add_argument('--msa', '-m', required=True, type=str, metavar='<str>', help="path to a msa in stockholm format")
parser.add_argument('--pdb', '-p', required=True, type=str, metavar='<str>', help="path to a pdb file for a member of the msa")
parser.add_argument('--pid', '-d', required=True, type=str, metavar='<str>', help="id for protein member to score against")
parser.add_argument('--sim', '-s', default=0.8, type=float, metavar='<float>', help="max sequence identity for sequence clustering")
parser.add_argument('--psuedo', '-l', default=1.0, type=float, metavar='<float>', help="pseudo count scale. psuedo counts added as l * n_eff")
parser.add_argument('--rescale', '-r', action='store_true', help="rescale psuedo counts (optional)")
parser.add_argument('--cutoff', '-c', default=8.0, type=float, metavar='<float>', help="distance cutoff to call contacts")

args = parser.parse_args()

for msa in read_stockholm(args.msa):
	#print(msa.identifier)
	#print(msa.accession)
	#print(f"num of seqs: {len(msa.seqs)}")
	#print(f"msa width: {len(msa.seqs[0])}")
	print()
	pdbser = PDBParser(PERMISSIVE=1)
	structure = pdbser.get_structure(args.pid, args.pdb)
	
	msa.measure_mij(similarity_cutoff=args.sim, psuedo=args.psuedo, rescale=args.rescale)
	
	msa.score_mij(test_id=args.pid, pdb=structure, cutoff=args.cutoff)
	
	for k,v in msa.cumulative_scores.items():
		if (k+1) % 10 != 0: continue
		
		print(f"rank {k+1:>4} correct: {v:>4}   acc: {v / (k+1):6.4f}")
		
		if (k+1) == 100: break
	
	print()
	
	
	
	
	
	sys.exit()
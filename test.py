
from os import system, walk
from main import load_matrix, fastaread
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
#format_alignment


# HelicoverpaArmigera-cMyc.fasta  
# LarimichthysCrocea-cMyc.fasta
# HomoSapiens-SHH.fasta
# RatusNorvegicus-SHH.fasta

FASTAS_PATH = "./ex1/fastas/" 
def gen_fastas_list():
	seqs = [ ]
	for root, dirs, files in walk(FASTAS_PATH):
		for fatfile in filter(lambda _filename: _filename[0] != ".", files):
			seqs.append( f"{FASTAS_PATH}/{fatfile}")
	return seqs
	# seqs.append(  )


import subprocess
def sanity():
	# result = 
	# result.stdout
	seqs = gen_fastas_list()
	
	for testlen in [-1] :#[ 1 , 2, 3, 4, 5, 6, 10, 20, 40, -1]:
		system(f"python3 ./main.py --align_type=golbal --testlen={testlen} --score=./ex1/score_matrix.tsv {seqs[0]} {seqs[1]}")
		matrix = load_matrix( "./ex1/score_matrix.tsv")
		# matrix = load_matrix( "./my_score.tsv")  
		print(seqs[0], seqs[1])
		seq_a, seq_b  = fastaread( seqs[0] ).__next__()[1], fastaread(seqs[1]).__next__()[1]
		for a in pairwise2.align.globalds(seq_a[:testlen], seq_b[:testlen], matrix, -8, -8)[:1]:
			print(pairwise2.format_alignment(*a))


def sanity_local():
	# result = 
	# result.stdout
	seqs = gen_fastas_list()
	
	for testlen in [-1]: #[ 1 , 2, 3, 4, 5, 6, 10, 20, 40]:
		system(f"python3 ./main.py --align_type=local --testlen={testlen} --score=./ex1/score_matrix.tsv {seqs[0]} {seqs[1]}")
		# matrix = load_matrix( "./ex1/score_matrix.tsv") 
		matrix = load_matrix( "./my_score.tsv") 
		seq_a, seq_b  = fastaread( seqs[0] ).__next__()[1], fastaread(seqs[1]).__next__()[1]
		for a in pairwise2.align.localds(seq_a[:testlen], seq_b[:testlen], matrix, -8, -8)[:1]:
			print(pairwise2.format_alignment(*a))


def sanity_overlap():
	# result = 
	# result.stdout
	seqs = gen_fastas_list()
	
	for testlen in [-1] :#, 6, 10]:
		system(f"python3 ./main.py --align_type=overlap --testlen={testlen} --score=my_score.tsv {seqs[0]} {seqs[1]}")
		matrix = load_matrix( "./my_score.tsv") 
		# matrix = load_matrix( "./ex1/score_matrix.tsv") 
		seq_a, seq_b  = fastaread( seqs[0] ).__next__()[1], fastaread(seqs[1]).__next__()[1]
		for a in pairwise2.align.globalds(seq_a[:testlen], seq_b[:testlen], matrix, -8, -8)[:1]:
			print(pairwise2.format_alignment(*a))


import sys

CMP_DIR = 'cmp'
FILES = [
    'ex1/fastas/HelicoverpaArmigera-cMyc.fasta',
    'ex1/fastas/HomoSapiens-SHH.fasta',
    'ex1/fastas/LarimichthysCrocea-cMyc.fasta',
    'ex1/fastas/RatusNorvegicus-SHH.fasta',
]

if __name__ == "__main__":
	# clean_score()
	# sanity_overlap()
	# sanity()
	sanity_local()
	# gen_fastas_list()
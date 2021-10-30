
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


def sanity():
	seqs = gen_fastas_list()
	testlen = 100
	system(f"python3 ./main.py --align_type=golbal --testlen={testlen} --score=./ex1/score_matrix.tsv {seqs[0]} {seqs[1]}")	
	matrix = load_matrix( "./ex1/score_matrix.tsv") 
	seq_a, seq_b  = fastaread( seqs[0] ).__next__()[1], fastaread(seqs[1]).__next__()[1]
	for a in pairwise2.align.globaldx(seq_a[:testlen], seq_b[:testlen], matrix):
		print(pairwise2.format_alignment(*a))
		break
	

if __name__ == "__main__":
	sanity()
	# gen_fastas_list()
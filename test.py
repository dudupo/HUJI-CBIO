
from os import system, walk




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
	system(f"python3 ./main.py --align_type=golbal --score=./ex1/score_matrix.tsv {seqs[0]} {seqs[1]} ")	


if __name__ == "__main__":
	sanity()
	# gen_fastas_list()
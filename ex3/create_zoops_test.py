from motif_find import generate_tau, sample
from ZOOPS_EM import write_motif_tsv

import argparse

def parse_args():
    """
    Parse the command line arguments.
    :return: The parsed args.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='File path with list of sequences (e.g. seqs_ATTA.fasta)')
    parser.add_argument('seed', help='Guess for the motif (e.g. ATTA)')
    parser.add_argument('p', type=float, help='Initial guess for the p transition probability (e.g. 0.01)')
    parser.add_argument('q', type=float, help='Initial guess for the q transition probability (e.g. 0.9)')
    parser.add_argument('alpha', type=float, help='Softening parameter for the initial profile (e.g. 0.1)')
    parser.add_argument('convergenceThr', type=float, help='ll improvement threshold for the stopping condition'
                                                           ' (e.g. 0.1)')
    return parser.parse_args()

def main():
    args = parse_args()
    write_motif_tsv("my_motif.tsv", args.seed, args.alpha )    
    p,q = args.p, args.q
    tau, emission, p, q = generate_tau("my_motif.tsv", p, q)
    
    _filename = f"cases/seqs_{args.seed}_{args.alpha}a_{args.p}p_{args.q}q.fasta"
    _file = open(_filename, "w+")
    _sample = sample(emission, tau ,q)
    _file.write(f">seq1\n{_sample}")
    for _ in range(2,200):
        _sample = sample(emission, tau ,q)
        _file.write(f"\n>seq{_}\n{_sample}")
    
    from random import random
    def addnoise(prob):
        return 0.5* ( prob + random())
    
    pp, qq = addnoise(p), addnoise(q)
    
    command = f"python3 -W ignore ZOOPS_EM.py {_filename} {args.seed} {pp} {qq} {args.alpha} {args.convergenceThr}\n"
    open("generated_tests.sh", "a+").write(command)

if __name__ == '__main__':
    main()
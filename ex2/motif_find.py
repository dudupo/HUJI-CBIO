import argparse
import numpy as np
from typing import List
from itertools import groupby



def load_matrix(path: str):
    """
    this function reads the matrix file and returns a python dict object such that every 2-tuple of letters from the
    alphabet is mapped to a float that is the value for the pair according to the file
    :param path: the path to the file
    :return: dict such that every 2-tuple of letters from the
    alphabet is mapped to a float that is the value for the pair according to the file
    """

    mat = [ dict( { 'A' :0.25, 'T' :0.25, 'C' :0.25, 'G' :0.25} ) ]
    file = open(path, 'r')
    alphabet = file.readline().split()
    for row in file.readlines():
        mat.append(dict())
        for char, val in zip(alphabet, row.split()):
            mat[-1][char] =  float(val)
    mat.append( dict( { 'A' :0.25, 'T' :0.25, 'C' :0.25, 'G' :0.25} ))
    return mat

def fastaread(fasta_name):
    """
    Read a fasta file. For each sequence in the file, yield the header and the actual sequence.
    In Ex1 you may assume the fasta files contain only one sequence.
    You may keep this function, edit it, or delete it and implement your own reader.
    """
    f = open(fasta_name)
    faiter = (x[1] for x in groupby(f, lambda line: line.startswith(">")))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq



def viterbi(emission, tau):
    pass

def forward(X, emission, tau):
    n, m = len(X), len(tau)
    F = np.zeros( shape=(n,m))
    # for l in range(m):
    F[0][0] =  1
    print(tau)
    for i in range(1,n):
        for l in range(m):
            F[i][l] += emission[l][X[i]] * (F[i-1].T @ tau[l])
        # print(F)
        # print()
        print()
    return F
def backward(X, emission, tau):
    n, m = len(X), len(emission)
    B = np.zeros( shape=(n,m))
    for i in reversed(range(n)):
        for l in range(m):
            B[i][l] =  tau[l].T @ (emission[l][X[i+1]]* (B[i+1]))
    return B

def posterior(k, X, emission, tau, i):
    F = forward(X, emission, tau)
    B = backward(X, emission, tau)
    Px = sum(F[-1])
    return F[i][k] * B[i][k]/Px

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--alg', help='Algorithm (e.g. viterbi)', required=True)
    parser.add_argument('seq', help='A sequence over the alphabet [A,C,G,T] (e.g. ACTGGACTACGTCATGCA)')
    parser.add_argument('initial_emission', help='Path to emission table (e.g. initial_emission.tsv)')
    parser.add_argument('p', help='transition probability p (e.g. 0.01)', type=float)
    parser.add_argument('q', help='transition probability q (e.g. 0.5)', type=float)
    args = parser.parse_args()
    
    X, p, q = str(args.seq), float(args.p), float(args.q)

    emission = load_matrix( args.initial_emission )

    tau = np.zeros( (len(emission), len(emission)))
    tau[0][0], tau[0][1] = 1-p, p
    tau[-1][-2], tau[-1][-1] = 1-q, q

    # tau = [ [ 1-p , p ] + [ 0 ] * (len(emission)-2)  ]
    for i in range(1,len(tau)-1): 
        tau[i][i+1] = 1
    # tau.append( [ 0 ] * (len(emission)-2)  + [ 1-q, q ] )
    tau = np.array(tau)
    



    if args.alg == 'viterbi':
        pass
        raise NotImplementedError

    elif args.alg == 'forward':
        ret =  forward(X, emission, tau)
        print(ret)

    elif args.alg == 'backward':
        pass
        raise NotImplementedError

    elif args.alg == 'posterior':
        pass
        raise NotImplementedError


if __name__ == '__main__':
    main()

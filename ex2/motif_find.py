import argparse
import numpy as np
from typing import List
from itertools import groupby

from scipy.special import logsumexp



def load_matrix(path: str):
    """
    this function reads the matrix file and returns a python dict object such that every 2-tuple of letters from the
    alphabet is mapped to a float that is the value for the pair according to the file
    :param path: the path to the file
    :return: dict such that every 2-tuple of letters from the
    alphabet is mapped to a float that is the value for the pair according to the file
    """

    logquad = np.log(0.25)

    mat = [ dict( { 'A' :logquad, 'T' :logquad, 'C' :logquad, 'G' :logquad} ) ]
    file = open(path, 'r')
    alphabet = file.readline().split()
    for row in file.readlines():
        mat.append(dict())
        for char, val in zip(alphabet, row.split()):
            mat[-1][char] =  np.log(float(val))
    mat.append( dict( { 'A' :logquad, 'T' :logquad, 'C' :logquad, 'G' :logquad} ))
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
    F[0][0] =  0
    print(tau)
    for i in range(1,n):
        for l in range(m):
            F[i][l] = emission[l][X[i]] + ( logsumexp( F[i-1] + tau[l]))
    return np.exp(F)

def backward(X, emission, tau):
    n, m = len(X), len(tau)
    B = np.zeros( shape=(n,m))
    
    B[n-1][m-1] = 0
    for i in reversed(range(n-1)):
        for l in range(m):
            B[i][l] =  logsumexp(tau[l] + (emission[l][X[i+1]]+ (B[i+1])))
    return np.exp(B)

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

    tau = np.ones( (len(emission), len(emission))) * -np.inf
    tau[0][0], tau[0][1] = np.log(1-p), np.log(p)
    tau[-1][-2], tau[-1][-1] = np.log(1-q), np.log(q)

    # tau = [ [ 1-p , p ] + [ 0 ] * (len(emission)-2)  ]
    for i in range(1,len(tau)-1): 
        tau[i][i+1] = 0
    # tau.append( [ 0 ] * (len(emission)-2)  + [ 1-q, q ] )
    tau = np.array(tau)
    



    if args.alg == 'viterbi':
        pass
        raise NotImplementedError

    elif args.alg == 'forward':
        ret =  forward(X, emission, tau)
        print(ret)

    elif args.alg == 'backward':
        ret = backward(X, emission, tau)
        print(ret)

    elif args.alg == 'posterior':
        pass
        raise NotImplementedError


if __name__ == '__main__':
    main()

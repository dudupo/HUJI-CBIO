import argparse
import numpy as np
from typing import List
from itertools import groupby

from scipy.special import logsumexp

from copy import deepcopy

def load_matrix(path: str):
    """
    this function reads the matrix file and returns a python dict object such that every 2-tuple of letters from the
    alphabet is mapped to a float that is the value for the pair according to the file
    :param path: the path to the file
    :return: dict such that every 2-tuple of letters from the
    alphabet is mapped to a float that is the value for the pair according to the file
    """

    logquad = np.log(0.25)
    mat = [ dict( { '^' : -np.inf, '$' : -np.inf,\
         'A' :logquad, 'T' :logquad, 'C' :logquad, 'G' :logquad} ) ]
    file = open(path, 'r')
    alphabet = file.readline().split()
    for row in file.readlines():
        mat.append(dict())
        for char, val in zip(alphabet, row.split()):
            mat[-1][char] =  np.log(float(val))
        for char in ['^', '$']:
            mat[-1][char] = -np.inf
    mat.append(  deepcopy(mat[0]) )
    begin, end = ( { key : -np.inf for key in mat[0].keys() } for _ in range(2))
    begin['$'], end['^'] = 0, 0
    return [begin] + mat + [end]

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
        seq = "$" + "".join(s.strip() for s in next(faiter)) + "^"
        yield header, seq

def forward(X, emission, tau):
    n, m = len(X), len(tau)    
    F = np.ones( shape=(n,m) ) * -np.inf
    F[0][0] = emission[0][X[0]]    
    for i in range(1,n):
        for l in range(m):
            F[i][l] = emission[l][X[i]] + ( logsumexp( F[i-1] + (tau.T)[l]))
                
    return np.exp(F)

def backward(X, emission, tau):
    n, m = len(X), len(tau)
    B = np.ones(shape=(n,m)) * -np.inf   
    B[n-1][m-1] = 0 #emission[m-1][X[n-1]] 
    for i in reversed(range(1, n)):
        for l in range(m):
            emission_k = np.array( [ emission[k][X[i]] for k in range(m) ]) 
            B[i-1][l] = logsumexp(emission_k + tau[l] + B[i])

    return np.exp(B)

def printHiddens(X, emission, tau, states):
    ret = ""
    # print(states)
    for state in states:
        if (state > 1) and (state < (len(tau)-2)): 
            ret += "M"
        else:
            ret += "B"
    
    def split_to_lines(_str):
        chunks, chunk_size = len(_str), 50
        return [_str[i:i + chunk_size] for i in range(0, chunks, chunk_size)]

    for upper_line, bottom_line in \
        zip(split_to_lines(ret[1:-1]), split_to_lines(X[1:-1])):
        print(upper_line)
        print(bottom_line)
        print()


def posterior(X, emission, tau):
    F = np.log(forward(X, emission, tau))
    B = np.log(backward(X, emission, tau))
    states = np.argmax((F + B), axis=1)
    printHiddens(X, emission, tau, states)

def viterbi(X, emission, tau):
    def vitforward():
        n, m = len(X), len(tau)    
        V = np.ones( shape=(n,m) ) * -np.inf
        P = np.zeros(shape=(n,m))
        V[0][0] = emission[0][X[0]]    
        for i in range(1,n):
            for l in range(m):
                V[i][l] = emission[l][X[i]] + ( np.max( V[i-1] + (tau.T)[l]))                    
                P[i][l] = int( np.argmax( V[i-1] + (tau.T)[l]) )
        return V, P

    V, P = vitforward()
    states = [len(tau)-1]
    for i in range(1,len(X)):
        states=  [ int(P[-i][states[0]]) ] + states  
    printHiddens(X, emission, tau, states)
    
def main():

    import warnings
    warnings.filterwarnings("ignore")

    parser = argparse.ArgumentParser()
    parser.add_argument('--alg', help='Algorithm (e.g. viterbi)', required=True)
    parser.add_argument('seq', help='A sequence over the alphabet [A,C,G,T] (e.g. ACTGGACTACGTCATGCA)')
    parser.add_argument('initial_emission', help='Path to emission table (e.g. initial_emission.tsv)')
    parser.add_argument('p', help='transition probability p (e.g. 0.01)', type=float)
    parser.add_argument('q', help='transition probability q (e.g. 0.5)', type=float)
    args = parser.parse_args()
    
    X, p, q = str(args.seq), float(args.p), float(args.q)
    X = "$" + X + "^"
    emission = load_matrix( args.initial_emission )

    # #############################
    # TODO: add transition from the first to the end with q probability.
    # #############################
    tau = np.ones( (len(emission), len(emission))) * -np.inf
    tau[1][1], tau[1][2] = np.log(1-p), np.log(p)
    tau[-2][-2] = np.log(1-p) #np.log(p)
    tau[-2][-1] = np.log(p)
    tau[0][1]   = np.log(q)
    tau[0][-2]  = np.log(1-q)
    
    for i in range(2,len(tau)-2): 
        tau[i][i+1] = 0

    tau = np.array(tau)


    if args.alg == 'viterbi':
        viterbi(X, emission, tau)
        
    elif args.alg == 'forward':
        ret =  forward(X, emission, tau)
        print(np.log(ret[-1][-1]))

    elif args.alg == 'backward':
        ret = backward(X, emission, tau)
        print(np.log(ret[0][0]))

    elif args.alg == 'posterior':
        posterior(X, emission, tau)
        


if __name__ == '__main__':
    main() 
    
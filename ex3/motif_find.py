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

def forward(X, emission, tau, q):
    n, m = len(X), len(tau)    
    F = np.ones( shape=(n+2,m) ) * -np.inf
    
    # for l in range(m):
    F[0][0] = np.log(q)     #+  emission[0][X[0]]
    F[0][-1] = np.log(1- q)   #+ emission[-1][X[0]]
    
    for i in range(1,n+1):
        for l in range(m):
            F[i][l] = emission[l][X[i-1]] + ( logsumexp( F[i-1] + (tau.T)[l])) 
    
    F[-1][-1] = F[-2][-1] +  tau[0][1] - tau[0][0]
    # F[0][0] =  logsumexp(np.array([np.log(q) + F[1][0], np.log(1- q) + F[1][-1]]))
    return np.exp(F)

def backward(X, emission, tau, q):
    n, m = len(X), len(tau)
    B = np.ones(shape=(n+2,m)) * -np.inf 
    
    B[n+1][m-1] = tau[0][1] - tau[0][0] 
    
    for i in reversed(range(1, n+2)):
        for l in range(m):
            emission_k = np.array( [ emission[k][X[i-2]] for k in range(m) ]) 
            B[i-1][l] =  logsumexp(tau[l] + emission_k + B[i])

    B[0][0] =  logsumexp(np.array([np.log(q) + B[1][0], np.log(1- q) + B[1][-1]])) 
    return np.exp(B)

def printHiddens(X, emission, tau, states):
    ret = ""
    # print(states)
    for state in states:
        if (state > 0) and (state < (len(tau)-1)): 
            ret += "M"
        else:
            ret += "B"
    
    def split_to_lines(_str):
        chunks, chunk_size = len(_str), 50
        return [_str[i:i + chunk_size] for i in range(0, chunks, chunk_size)]

    for upper_line, bottom_line in \
        zip(split_to_lines(ret), split_to_lines(X)):
        print(upper_line)
        print(bottom_line)
        print()


def AGCT(char):
    return int({'A':0, 'G':1, 'C':2, 'T':3}[char])

from itertools import product
def transition_event(X, emission, tau, q):
    F = np.log(forward(X, emission, tau, q))
    B = np.log(backward(X, emission, tau, q))
    stats = np.zeros(tau.shape)
    
    emissions = np.array([[emission[l][x] for x in X]\
         for l in range(tau.shape[0])]) 
    
    for k in range(tau.shape[0]):
        for l in range(tau.shape[0]): 
            stats[k][l] = logsumexp((F[0:-2,k] + B[2:,l] - F[-1][-1])\
                 + tau[k,l] + emissions[l])
    
    
    stats_emis = np.zeros((tau.shape[0]-2,4))
    for k in range(tau.shape[0]-2):
        for neckloied in ['A', 'G', 'C', 'T']:
            for i,char in enumerate(X):
                if char == neckloied:
                    stats_emis[k][AGCT(char)] = logsumexp([stats_emis[k][AGCT(char)],\
                         F[i+1,k+1] + B[i+1,k+1] - F[-1][-1]])
    
    Nnq  = np.log(1- q) + B[1][-1] - F[-1][-1]
    Nq = np.log(q) + B[1][0] - F[-1][-1]

    Np = logsumexp([ np.log(2) + stats[0][1],  stats[0][-1]  ]) #, tau[0][1] - F[-1][-1] 
    Nnp = logsumexp([ stats[0][0], stats[-1][-1]]) #  
    S = np.array([
        Nq, 
        Nnq,
        Np,
        Nnp
    ])   
    ret = np.hstack((S, stats_emis.flatten()))
    # print(ret)
    return S, stats_emis 

def posterior(X, emission, tau, q):
    F = np.log(forward(X, emission, tau, q))
    B = np.log(backward(X, emission, tau, q))
    states = np.argmax((F[0:-1,:] + B[1:,:])[1:], axis=1)
    printHiddens(X, emission, tau, states)

def viterbi(X, emission, tau, q):
    def vitforward():
        n, m = len(X), len(tau)    
        F = np.ones(shape=(n,m) ) * -np.inf
        P = np.zeros(shape=(n,m) )
        
        F[0][0] = np.log(q) + emission[0][X[0]]
        F[0][-1] = np.log(1- q) + emission[-1][X[0]]

        for i in range(1,n):
            for l in range(m):
                F[i][l] =  np.max( emission[l][X[i]] + F[i-1] + (tau.T)[l])
                P[i][l] = int(np.argmax( emission[l][X[i]] + F[i-1] + (tau.T)[l]))
        F[-1][-1] += tau[0][1]
        return np.exp(F), P

    F, P = vitforward()
    states = [len(tau)-1]
    for i in range(len(X)-1):
        states=  [ int(P[-i-1][states[0]]) ] + states    
    printHiddens(X, emission, tau, states)

def _viterbi(X, emission, tau, q):
    def vitforward():
        n, m = len(X), len(tau)    
        F = np.ones(shape=(n,m) ) * -np.inf
        P = np.zeros(shape=(n,m) )
        
        F[0][0] = np.log(q) + emission[0][X[0]]
        F[0][-1] = np.log(1- q) + emission[-1][X[0]]

        for i in range(1,n):
            for l in range(m):
                F[i][l] =  np.max( emission[l][X[i]] + F[i-1] + (tau.T)[l])
                P[i][l] = int(np.argmax( emission[l][X[i]] + F[i-1] + (tau.T)[l]))
        F[-1][-1] += tau[0][1]
        return np.exp(F), P

    F, P = vitforward()
    states = [len(tau)-1]
    for i in range(len(X)-1):
        states=  [ int(P[-i-1][states[0]]) ] + states 

    for i,state in enumerate(states):
        if (state > 0) and (state < (len(tau)-1)):
            return i               
    return -1
def generate_tau(initial_emission, p, q):
    emission = load_matrix( initial_emission )

    tau = np.ones( (len(emission), len(emission))) * -np.inf
    tau[0][0], tau[0][1] = np.log(1-p), np.log(p)
    tau[-1][-1] = np.log(1-p) #np.log(p)

    for i in range(1,len(tau)-1): 
        tau[i][i+1] = 0
    tau = np.array(tau)
    return tau, emission, p ,q


from random import random, choice
def sample(emission, tau, q):
    _tau = np.exp(tau)

    def sample_unif():
        return choice([ ['A'] , ['G'], ['C'], ['T'] ])
    def sample_motif():
        motif = []
        for _ in range(1, tau.shape[0] - 1):
            uni_sample, current_prob = random(), 0
            for char, prob in emission[_].items():
                current_prob += np.exp(prob)
                if uni_sample < current_prob: 
                    motif.append(char)
                    break
        return motif
   
    p = _tau[0][1] 
    ret = []
    if (random() < q):
        ret += sample_unif()
        while random() < (1-p):
            ret += sample_unif()
        ret += sample_motif()
    
    ret += sample_unif()
    while random() < (1-p):
        ret += sample_unif()
    
    return "".join(ret)

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

    emission = load_matrix( args.initial_emission )

    tau = np.ones( (len(emission), len(emission))) * -np.inf
    tau[0][0], tau[0][1] = np.log(1-p), np.log(p)
    tau[-1][-1] = np.log(1-p) #np.log(p)

    for i in range(1,len(tau)-1): 
        tau[i][i+1] = 0

    tau = np.array(tau)


    if args.alg == 'viterbi':
        viterbi(X, emission, tau, q)
        
    elif args.alg == 'forward':
        ret =  forward(X, emission, tau, q)
        print(np.log(ret[-1][-1]))

    elif args.alg == 'backward':
        ret = backward(X, emission, tau, q)
        print(np.log(ret[0][0]))

    elif args.alg == 'posterior':
        posterior(X, emission, tau, q)
        


if __name__ == '__main__':
    main() 
    
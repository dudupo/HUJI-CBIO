import argparse
import numpy as np

def viterbi(emission, tau):
    pass

def forward(X, emission, tau):
    n, m = len(X), len(emission)
    F = np.zeros( shape=(n,m))
    for i in range(n):
        for l in range(m):
            F[i][l] =  emission[l][X[i]] * (F[i-1].T @ tau[l])
    return F
def backward(X, emission, tau):
    n, m = len(X), len(emission)
    B = np.zeros( shape=(n,m))
    for i in reversed(range(n)):
        for l in range(m):
            B[i][l] =  tau[l].T @ (emission[X[i+1]]* (B[i+1]))
    return B

def posterior(X, emission, tau, i):
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

    if args.alg == 'viterbi':
        raise NotImplementedError

    elif args.alg == 'forward':
        raise NotImplementedError

    elif args.alg == 'backward':
        raise NotImplementedError

    elif args.alg == 'posterior':
        raise NotImplementedError


if __name__ == '__main__':
    main()

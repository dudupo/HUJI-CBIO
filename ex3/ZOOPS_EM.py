# import emission as emission

from motif_find import transition_event, generate_tau, sample , forward, _viterbi, AGCT,backward
import argparse
import numpy as np
from scipy.special import logsumexp


def expectaion(seqs, emission, tau ,q, parity=True):
    '''given distributions returns the expectaion of the stats'''
    p = np.exp(tau[0][1])

    S, emit =  transition_event(seqs[0], emission, tau, q)
    S, emit = np.exp(S), np.exp(emit)
    emit -= 1
    for seq in seqs[1:]:
        _S, _emit =  transition_event(seq, emission, tau, q)
        S += np.exp(_S)
        emit += np.exp(_emit)
        emit -= 1

    # if True : #parity:
    q = S[0]/(S[0] + S[1])
    p = S[2]/(S[2] + S[3])
    retemissions = emission
    # print("q:{0}, p:{1}".format(q, p))

    emissiontilde = emit 
    temp =  np.sum(emissiontilde, axis=1)
    for j in range(len(emissiontilde)):
        emissiontilde[j] = np.divide(emissiontilde[j], temp[j],\
            out=np.zeros_like(emissiontilde[j]), where=temp[j]  != 0)

    retemissions = [ ]
    retemissions.append({ a : np.log(0.25) for a in ['A','G', 'C', 'T'] })
    for k in range(len(tau)-2):
        retemissions.append({ a : np.log(emissiontilde[k][AGCT(a)]) for a in ['A','G', 'C', 'T'] })
    retemissions.append({ a : np.log(0.25) for a in ['A','G', 'C', 'T'] })
    return p,q, retemissions

def compute_log_likelihood_for_sequences(sequences,tau,p,q,emissiones):
    """
    given a list of sequences transition probability matrix and probability to go in and out of a motif and emissiones
    matrix returns the log likelihod for the sequences
    :param sequences: the sequences
    :param tau: transition probability matrix
    :param p: probability to go in and out of a motif
    :param q: probability to go out of a motif
    :return: log likelihood of the sequences given tau p q
    """
    return np.sum([ np.log(forward(seq,emissiones,tau,q))[-1][-1] for seq in sequences ])




def BaumWelch(seqs, emission, tau, q, convergenceThr):

    result_history = []
    result_history.append(
        compute_log_likelihood_for_sequences(
            seqs, tau, np.exp(tau[0][1]), q, emission))

    p, q, emission = expectaion(seqs, emission, tau, q)    
    tau[0][0], tau[0][1] = np.log(1 - p), np.log(p)
    tau[-1][-1] = np.log(1 - p)
    result_history.append(
        compute_log_likelihood_for_sequences(
            seqs, tau, p, q, emission))
    # parity = not parity

    while (result_history[-1]-result_history[-2]) > convergenceThr:
        # retemission, retp, retq = emission, np.exp(tau[0][1]), q
        p,q, emission = expectaion(seqs, emission, tau ,q)
        tau[0][0], tau[0][1] = np.log(1-p), np.log(p)
        tau[-1][-1] = np.log(1-p)
        result_history.append(
            compute_log_likelihood_for_sequences(
                seqs,tau,p,q,emission))
        # parity = not parity
        
    return tau, p, q, result_history, emission

def dump_results(result_history, emissiones, p, q, sequences,tau):
    """write results"""
    # shaked
    with open("ll_history.txt","w+") as history_file:
        history_file.write(f"{str(result_history[0],)}")
        for i in result_history[1:]:
            history_file.write(f"\n{str(i,)}")
    with open("motif_profile.txt","w+") as motif_profile_file:
        for base in ["A","C","G","T"]:
            for s in emissiones[1:-2]:
                motif_profile_file.write(f"{str(round(np.exp(s[base]),2))}\t")
            motif_profile_file.write(f"{str(round(np.exp(emissiones[-2][base]),2))}")
            motif_profile_file.write("\n")
        motif_profile_file.write(f"{str(round(q,4))}\n")
        motif_profile_file.write(f"{str(round(p, 4))}")
    with open("motif_positions.txt","w+") as motif_positions_file:
        for seq in sequences:
            motif_positions_file.write(f"{_viterbi(seq,emissiones,tau,q)}\n")
        











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

def write_motif_tsv(_file_path, seed, alpha):
    beta = 1 - 3 * alpha
    with open(_file_path, 'w+') as _file:
        probabilities = { "A" : alpha, "C" : alpha, "G" : alpha, "T": alpha }
        _file.write("A C G T")
        for charter in seed:
            probabilities[charter] = beta
            _file.write("\n{0} {1} {2} {3}".format(*probabilities.values()))
            probabilities[charter] = alpha

def readseqs( _file_path ):
    return [s[:-1] for s in open(_file_path, 'r').readlines()[1::2] if len(s) > 1 ]

def main():
    args = parse_args()
    write_motif_tsv("my_motif.tsv", args.seed, args.alpha )    
    p,q = args.p, args.q
    tau, emission, p, q = generate_tau("my_motif.tsv", p, q)
    
    
    seqs = readseqs( args.fasta )    
    tau, p, q ,log_likelihood_history, emission= BaumWelch(seqs, emission, tau ,q, args.convergenceThr)
    dump_results(log_likelihood_history, emission, p, q,seqs,tau)

if __name__ == "__main__":
    main()


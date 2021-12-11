from motif_find import transition_event, generate_tau
import argparse
import numpy as np
# from motif_find import...

def transitions ():
    # shaked
    pass
def emissions ():
    # shaked
    pass

def maximize(stats):
    ''' estimate the emissions and the transitions tables '''
    # shaked
    pass

def expectaion(seqs, emission, tau ,q):
    '''given distributions returns the expectaion of the stats'''
    stats = np.array(list(map(\
        lambda x: transition_event(x, emission, tau, q), seqs)))
    stats = np.sum(stats, axis=0) #/ len(seqs) # Np,Nq
    # print(stats)
    r = stats[0] + stats[1]
    h = stats[2] + stats[3]
    p,q = stats[0] / r , stats[2] / h
    print(p,q)
    return p,q 

def BaumWelch(seqs, emission, tau, q, convergenceThr):
    for j in range(100): #convergenceThr
        p,q = expectaion(seqs, emission, tau ,q)
        tau[0][0], tau[0][1] = np.log(1-p), np.log(p)
        tau[-1][-1] = np.log(1-p) 
    return tau, p, q

def dump_results():
    # shaked
    pass

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
    _file = open(_file_path, 'w+')
    probabilities = { "A" : alpha, "C" : alpha, "G" : alpha, "T": alpha }
    _file.write("A C G T")
    for charter in seed:
        probabilities[charter] = beta
        _file.write("\n{0} {1} {2} {3}".format(*probabilities.values()))
        probabilities[charter] = alpha

def readseqs( _file_path ):
    print(_file_path)
    return [s[:-1] for s in open(_file_path, 'r').readlines()[1::2]]

def main():
    args = parse_args()
    write_motif_tsv("my_motif.tsv", args.seed, args.alpha )    
    p,q = args.p, args.q
    tau, emission, p, q = generate_tau("my_motif.tsv", p, q)
    seqs = readseqs( args.fasta )    
    print(seqs)
    tau, p, q = BaumWelch(seqs, emission, tau ,q, args.convergenceThr)
    print(p,q)


    # build transitions 
    # build emissions 

    # load fasta

    # run Baum-Welch

    # dump results


if __name__ == "__main__":
    main()


import argparse
from itertools import groupby

import numpy as np


def load_matrix(path: str):
    """
    this function reads the matrix file and returns a python dict object such that every 2-tuple of letters from the
    alphabet is mapped to a float that is the value for the pair according to the file
    :param path: the path to the file
    :return: dict such that every 2-tuple of letters from the
    alphabet is mapped to a float that is the value for the pair according to the file
    """

    mat = dict()
    file = open(path, 'r')
    alphabet = file.readline().split()
    for row in alphabet:
        line = file.readline().split()
        char = line[0]
        for val, col in zip(line[1:], alphabet):
            mat[(char, col)] = float(val)
    return mat


def local_base_case(seq1: str, seq2: str, mat: dict)
    """
    given the sequences for the program return a matrix field according to the base case and the
    score matrix for the local alignment
    :param seq1:the first sequence
    :param seq1:the second sequence
    :param mat:the score matrix that was returned from read_matrix function
    :return:ndArray with sides the size of seq1 and seq2 and the values in the leftmost column and upper most row will 
    be 0 
    """
    shape = (len(seq1), len(seq2))
    # todo is this faster then using np.empty and iterating with a python loop over the uppermost row an leftmost column?
    table = np.zeros(shape, dtype=float)
    return table


def global_base_case(seq1: str, seq2: str, mat: dict):
    """
    given the sequences for the program return a matrix field according to the base case and the
    score matrix for the global alignment
    :param seq1:the first sequence
    :param seq1:the second sequence
    :param mat:the score matrix that was returned from read_matrix function
    :return:ndArray with sides the size of seq1 and seq2 and the values in the leftmost column would be acording to the
    alignment of the beginning of seq1 with '-'
    """
    shape = (len(seq1), len(seq2))
    table = np.empty(shape, dtype=float)
    table[0, 0] = mat[(seq1[0], '-')]
    for row, char in enumerate(seq1[1:]):
        table[row, 0] = mat[(char, '-')] + table[row - 1, 0]
    return table


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


def main():
    a = fastaread("ex1/fastas/HomoSapiens-SHH.fasta").__next__()[1]
    parser = argparse.ArgumentParser()
    parser.add_argument('seq_a', help='Path to first FASTA file (e.g. fastas/HomoSapiens-SHH.fasta)')
    parser.add_argument('seq_b', help='Path to second FASTA file')
    parser.add_argument('--align_type', help='Alignment type (e.g. local)', required=True)
    parser.add_argument('--score', help='Score matrix in.tsv format (default is score_matrix.tsv) ',
                        default='score_matrix.tsv')
    command_args = parser.parse_args()
    if command_args.align_type == 'global':
        raise NotImplementedError
    elif command_args.align_type == 'local':
        raise NotImplementedError
    elif command_args.align_type == 'overlap':
        raise NotImplementedError
    # print the best alignment and score


if __name__ == '__main__':
    main()

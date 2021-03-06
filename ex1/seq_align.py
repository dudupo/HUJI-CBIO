import argparse
from itertools import groupby
from typing import List
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


def local_base_case(seq1: str, seq2: str, mat: dict):
    """
    given the sequences for the program return a matrix field according to the base case and the
    score matrix for the local alignment
    :param seq1:the first sequence
    :param seq1:the second sequence
    :param mat:the score matrix that was returned from read_matrix function
    :return:ndArray with sides the size of seq1 and seq2 and the values in the leftmost column and upper most row will 
    be 0 
    """
    shape = (len(seq1) + 1, len(seq2) + 1)
    table = np.zeros(shape) 
    table[:, 0] = 0
    table[0, :] = 0
    trace = np.ones(shape) * -1 # np.ones(shape, dtype=np.int8) * 3
    trace, table = fill_tables_for_local(seq1, seq2, mat, table, trace)
    return table, trace


def fill_cell_for_local(seq1, seq2, mat, table, trace, i, j):
    options = [table[i - 1, j] + mat[(seq1[i - 1], '-')], table[i - 1, j - 1] + mat[(seq1[i - 1], seq2[j - 1])],
               table[i, j - 1] + mat[('-', seq2[j - 1])], 0]
    table_val, trace_val = max_argmax(options)
    table[i, j] = table_val
    trace[i, j] = trace_val
    


def fill_tables_for_local(seq1: str, seq2: str, mat: dict, table, trace):
    """
    fill the whole table for the local aliment using fill_cell_for_local function and returns the trace and table for
    seq1 seq2
    :param seq1:the first sequence
    :param seq2:the second sequence
    :param mat: the score matrix
    :param table: the dynamic programming table for the best scores of partial alignments
    :param trace:the matrix we use to indicate what is the partial aliment we used in order to get to this alignment
    :return:trace and table
    """
    for col in range(1, table.shape[1]):
        for row in range(1, table.shape[0]):
            fill_cell_for_local(seq1, seq2, mat, table, trace, row, col)
    return trace, table


def extract_solution_local(seq1: str, seq2: str, mat: dict, table, trace):
    """
    extract the optimal local alignment for sqe1 seq2 based on table, trace
    :param seq1:the first sequence
    :param seq2:the second sequence
    :param mat: the score matrix
    :param table: the dynamic programming table for the best scores of partial alignments
    :param trace:the matrix we use to indicate what is the partial aliment we used in order to get to this alignment
    :return:trace and table
    :return: two strings that are sqe1 and seq2 with '-' in them according to the optimal global alignment
    """
    i, j = np.unravel_index(table.argmax(), table.shape)
    score = table[i, j]
    str1, str2 = "", ""
    while i > 0 and j > 0:
        if trace[i, j] == 2:
            str1 += '-'
            str2 += seq2[j - 1]
            j -= 1
            continue
        elif trace[i, j] == 1:
            str1 += seq1[i - 1]
            str2 += seq2[j - 1]
            j -= 1
            i -= 1
            continue
        elif trace[i, j] == 0:
            str1 += seq1[i - 1]
            str2 += '-'
            i -= 1
            continue
        else: 
            break
    return str1[::-1], str2[::-1], score


def overlap_base_case(seq1: str, seq2: str, mat: dict):
    """
    given the sequences for the program return a matrix field according to the base case and the
    score matrix for the overlap alignment
    :param seq1:the first sequence
    :param seq1:the second sequence
    :param mat:the score matrix that was returned from read_matrix function
    :return:ndArray with sides the size of seq1 and seq2 and the values in the leftmost column will be 0 
    """

    shape = (len(seq1)+1, len(seq2)+1)
    table = np.zeros(shape)
    trace = np.zeros(shape)
    table[:, 0] = 0
    trace, table = fill_tables_for_overlap(seq1, seq2, mat, table, trace)
    return table, trace



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
    shape = (len(seq1) + 1, len(seq2) + 1)
    trace = np.zeros(shape)  
    trace[:, 0] = 0
    cost = np.zeros(shape)  
    cost[0, 0] = 0  
    for row, char in enumerate(seq1):
        cost[row + 1, 0] = mat[(char, '-')] + cost[row, 0]
    trace, cost = fill_tables_for_global(seq1, seq2, mat, cost, trace)
    return cost, trace


def max_argmax(options: List[float]) -> tuple:
    """
    retuen the min value and its index from a array
    :param options: a list with numbers
    :return: min,argmin of the list
    """
    argmax = np.argmax(options)
    return options[argmax], argmax


def fill_cell_for_global(seq1: str, seq2: str, mat: dict, table, trace, i: int, j: int):
    """
    fill one cell of the matrix according to the global alignment algorithm
    :param seq1:the first sequence
    :param seq2: the second sequence
    :param mat: the score matrix
    :param table: the dynamic programming table for the best scores of partial alignments
    :param trace:the matrix we use to indicate what is the partial aliment we used in order to get to this alignment
    :param i:the row of the cell we wish to calculate
    :param j:the column of the cell we wish to calculate
    :return:none
    """
    options = [table[i - 1, j] + mat[(seq1[i - 1], '-')], table[i - 1, j - 1] + mat[(seq1[i - 1], seq2[j - 1])],
               table[i, j - 1] + mat[('-', seq2[j - 1])]]
    table_val, trace_val = max_argmax(options)
    table[i, j] = table_val
    trace[i, j] = trace_val


def fill_tables_for_global(seq1: str, seq2: str, mat: dict, table, trace):
    """
    fill the whole table for the global aliment using fill_cell_for_global function and returns the trace and table for
    seq1 seq2
    :param seq1:the first sequence
    :param seq2:the second sequence
    :param mat: the score matrix
    :param table: the dynamic programming table for the best scores of partial alignments
    :param trace:the matrix we use to indicate what is the partial aliment we used in order to get to this alignment
    :return:trace and table
    """
    for col in range(1, table.shape[1]):
        trace[0, col] = 2
        table[0, col] = table[0, col - 1] + mat[(seq2[col - 1], '-')]
        for row in range(1, table.shape[0]):
            fill_cell_for_global(seq1, seq2, mat, table, trace, row, col)
    return trace, table


def fill_cell_for_ovarlap(seq1, seq2, mat, table, trace, i, j):
    if i!=table.shape[0]-1:
        options = [table[i - 1, j] + mat[(seq1[i - 1], '-')], table[i - 1, j - 1] + mat[(seq1[i - 1], seq2[j - 1])],
                table[i, j - 1] + mat[('-', seq2[j - 1])]]
    else:
        options = [table[i - 1, j] + mat[(seq1[i - 1], '-')], table[i - 1, j - 1] + mat[(seq1[i - 1], seq2[j - 1])],
                   table[i, j - 1]]
    table_val, trace_val = max_argmax(options)
    table[i, j] = table_val
    trace[i, j] = trace_val




def fill_tables_for_overlap(seq1: str, seq2: str, mat: dict, table, trace):
    """
    fill the whole table for the overlap aliment using fill_cell_for_global function and returns the trace and table for
    seq1 seq2
    :param seq1:the first sequence
    :param seq2:the second sequence
    :param mat: the score matrix
    :param table: the dynamic programming table for the best scores of partial alignments
    :param trace:the matrix we use to indicate what is the partial aliment we used in order to get to this alignment
    :return:trace and table
    """
    for col in range(1, table.shape[1]):
        trace[0, col] = 2
        table[0, col] = table[0, col - 1] + mat[(seq2[col - 1], '-')]
        for row in range(1, table.shape[0]):
            fill_cell_for_ovarlap(seq1, seq2, mat, table, trace, row, col)
    return trace, table


# todo finish
def extract_solution_global(seq1: str, seq2: str, mat: dict, table, trace):
    """
    extract the optimal global alignment for sqe1 seq2 based on table, trace
    :param seq1:the first sequence
    :param seq2:the second sequence
    :param mat: the score matrix
    :param table: the dynamic programming table for the best scores of partial alignments
    :param trace:the matrix we use to indicate what is the partial aliment we used in order to get to this alignment
    :return:trace and table
    :return: two strings that are sqe1 and seq2 with '-' in them according to the optimal global alignment
    """
    i, j = table.shape[0] - 1, table.shape[1] - 1
    str1, str2 = "", ""
    while i > 0 and j > 0:
        if trace[i, j] == 2:
            str1 += '-'
            str2 += seq2[j - 1]
            j -= 1
            continue
        # change
        elif trace[i, j] == 1:
            str1 += seq1[i - 1]
            str2 += seq2[j - 1]
            j -= 1
            i -= 1
            continue
        else:
            str1 += seq1[i - 1]
            str2 += '-'
            i -= 1
            continue
    return str1[::-1], str2[::-1], table[-1,-1]


def extract_solution_overlap(seq1: str, seq2: str, mat: dict, table, trace):
    

    i, j = table.shape[0] - 1, np.argmax(table[table.shape[0] - 1])
    score = table[i,j]
    str1, str2 = "", ""    
    k = table.shape[1] - 1

    while k > j:
        str2 += seq2[k - 1]
        str1 += '-'
        k -= 1


    while j > 0:
        if trace[i, j] == 2:
            str1 += '-'
            str2 += seq2[j - 1]
            j -= 1
            continue
        # change
        elif trace[i, j] == 1:
            str1 += seq1[i - 1]
            str2 += seq2[j - 1]
            j -= 1
            i -= 1
            continue
        else:
            str1 += seq1[i - 1]
            str2 += '-'
            i -= 1
            continue
    
    while i > 0:
        str1 += seq1[i-1]
        str2 += '-'
        i -= 1

    return str1[::-1], str2[::-1], score


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


def general_alignment(seq_a, seq_b, mat, base_case_func, extract_func):
    def split_to_lines(_str):
        chunks, chunk_size = len(_str), 50
        return [_str[i:i + chunk_size] for i in range(0, chunks, chunk_size)]

    table, trace = base_case_func(seq_a, seq_b, mat)
    ret_seq1, ret_seq2, score = extract_func(seq_a, seq_b, mat, table, trace)
    for upper_line, bottom_line in \
            zip(split_to_lines(ret_seq1), split_to_lines(ret_seq2)):
        print(upper_line)
        print(bottom_line)
        print()
    return table, int(score)


def global_alignment(seq_a, seq_b, mat):
    table, score = general_alignment(seq_a, seq_b, mat, \
                              global_base_case, extract_solution_global)
    print("global:{0}".format(score))
def local_alignment(seq_a, seq_b, mat):
    table, score = general_alignment(seq_a, seq_b, mat, \
                              local_base_case, extract_solution_local)
    print("local:{0}".format(score))
def overlap_alignment(seq_a, seq_b, mat):
    table, score = general_alignment(seq_a, seq_b, mat, \
                              overlap_base_case, extract_solution_overlap)
    
    print("overlap:{0}".format(score))

def test_overlap():
    matrix = load_matrix( "./ex1/score_matrix.tsv")

    seqs = [ ("AGG" , "GAC" ), ("AGGCTTTTTTT" , "GGCCAGCCTTTCCC" )]
    for (x,y) in seqs:
        print(x,y)
        overlap_alignment(x,y, matrix)



import sys

sys.setrecursionlimit(50000)


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('seq_a', help='Path to first FASTA file (e.g. fastas/HomoSapiens-SHH.fasta)')
    parser.add_argument('seq_b', help='Path to second FASTA file')
    parser.add_argument('--align_type', help='Alignment type (e.g. local)', required=True)
    parser.add_argument('--score', help='Score matrix in.tsv format (default is score_matrix.tsv) ',
                        default='score_matrix.tsv')
    command_args = parser.parse_args()
    
    seq_a, seq_b = fastaread(command_args.seq_a).__next__()[1], fastaread(command_args.seq_b).__next__()[1]
    seq_a, seq_b = seq_a[:int(command_args.testlen)], seq_b[:int(command_args.testlen)]
    mat = load_matrix(command_args.score)

    if command_args.align_type == 'global':
        alignment = global_alignment(seq_a, seq_b, mat)
    elif command_args.align_type == 'local':
        alignment = local_alignment(seq_a, seq_b, mat)
    elif command_args.align_type == 'overlap':
        alignment = overlap_alignment(seq_a, seq_b, mat)

if __name__ == '__main__':
    main()


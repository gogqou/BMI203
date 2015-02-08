'''
Created on Feb 7, 2015

@author: gogqou
'''

import sub_matrix_mk_dict as subMdict
import numpy as np
import sys
import read_fasta as rFasta
import readwrite_file as RWFile
import SW_run as SW
import os

def SW_one_round(seq1, seq2, sub_Matrix):
    
    [sim_Matrix, pointers] = SW.similarity_matrix(seq1,seq2, sub_Matrix)
    
    [aligned_sequence, fitted_seq1, fitted_seq2, sim_Matrix, score] = SW.trace_aligned_seq(seq1, seq2, sim_Matrix, pointers)
    print aligned_sequence
    print fitted_seq1
    print fitted_seq2
    print 'score = ', score
    
    return score
def list_of_sequences_from_txt(sequence_list_file):
    seq_list = []
    lines = RWFile.readtxt(sequence_list_file)
    for line in lines:
        seq_list.append(line.split(' '))
    return seq_list
def main():
    if len(sys.argv)>4:
        print 'provide sequences to align, directory, and substitution matrix '
        sys.exit()
    
    home = sys.argv[2]
    subMatrixFile = sys.argv[3]
    
    sub_Matrix = subMdict.mk_dict(subMatrixFile)
    
    seq_list = list_of_sequences_from_txt(sys.argv[1])
    print seq_list
    for i in range(0,len(seq_list)):
        seq_locations = seq_list[i]
        seqfile1=seq_locations[0]
        seqfile2=seq_locations[1]
        seq1=rFasta.read_fa(home+'/'+seqfile1)
        print seq1, len(seq1)
        seq2=rFasta.read_fa(home+'/'+seqfile2)
        print seq2, len(seq2)
        score = SW_one_round(seq1, seq2, sub_Matrix)
        

    return 'done'
if __name__ == '__main__':
    main()
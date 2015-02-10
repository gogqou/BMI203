'''
Created on Feb 9, 2015

@author: gogqou
'''
import sub_matrix_mk_dict as subMdict
import numpy as np
import sys
import read_fasta as rFasta
import readwrite_file as RWFile
import SW_run as SW
import os
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d


def SW_one_round(seq1, seq2, sub_Matrix, origSubMatrix, gap_init, gap_ext):
    
    [sim_Matrix, pointers] = SW.similarity_matrix(seq1,seq2, sub_Matrix, gap_init, gap_ext)
    [aligned_sequence, fitted_seq1, fitted_seq2, sim_Matrix, score, count_array] = SW.trace_aligned_seq(seq1, seq2, sim_Matrix, pointers, origSubMatrix, sub_Matrix)
    print sim_Matrix
    np.savetxt('sim_Matrix.txt', sim_Matrix, fmt='%.1d', delimiter = ' ')
    print count_array
    print aligned_sequence
    print fitted_seq1
    print fitted_seq2
    print 'score = ', score
    
    #print count_array
    
    return score, count_array

def main():
    np.set_printoptions(threshold=1000, linewidth=1000, precision = 5, suppress = False)
    
    if len(sys.argv)>5:
        print 'provide positive and negative pairs of sequences to align, directory, and substitution matrix choice '
        sys.exit()
    #seq1 = sys.argv[1]
    #seq2 = sys.argv[2]
    seqfile1 = sys.argv[1]
    seqfile2= sys.argv[2]
    
    home = sys.argv[3]
    seq1=rFasta.read_fa(home+seqfile1)
    seq2=rFasta.read_fa(home+seqfile2)
    subMatrixFile = sys.argv[4]
    
    gap_init = 13
    gap_ext = 3
    
    [sub_Matrixdict, origSubMatrix, AAlist] = subMdict.mk_dict(home+subMatrixFile)
    print origSubMatrix
    [score, count_array] = SW_one_round(seq1, seq2, sub_Matrixdict,origSubMatrix, gap_init, gap_ext)
    
    
    np.savetxt(home+'test_count_array.txt', count_array, fmt ='%.1d', delimiter = '  ', header = ' A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *')
if __name__ == '__main__':
    main()
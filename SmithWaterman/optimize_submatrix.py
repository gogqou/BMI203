'''
Created on Feb 8, 2015

@author: gogqou


start with a static alignment with a given matrix and optimize it against objective function :
sum TP for FP rates of 0, .1, .2, .3


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
import seq_list_SW as slSW
import pylab


def main():
    if len(sys.argv)>5:
        print 'provide positive and negative pairs of sequences to align, directory, and substitution matrix choice '
        sys.exit()
    pos_seq_list_name = sys.argv[1]
    neg_seq_list_name = sys.argv[2]
    home = sys.argv[3]
    pos_seq_list_file = home+pos_seq_list_name+'.txt'
    neg_seq_list_file = home+neg_seq_list_name+'.txt'
    subMatrixFile = sys.argv[4]
    gap_init = 11
    gap_ext = 2
    [sub_Matrix, origSubMatrix] = subMdict.mk_dict(home+subMatrixFile)

    [pos_scores, pos_align_array] = slSW.scores_from_seq_list(home, pos_seq_list_file, sub_Matrix, origSubMatrix, gap_init, gap_ext)
    [neg_scores, neg_align_array]= slSW.scores_from_seq_list(home, neg_seq_list_file, sub_Matrix, origSubMatrix, gap_init, gap_ext)

    print 'done'
    
    return 1

if __name__ == '__main__':
    main()
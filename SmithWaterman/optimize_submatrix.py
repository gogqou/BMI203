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
import ROC_curves as ROC
def obj_function(pos_scores, neg_scores):
    obj_fnc_val=0
    for FPrate in range(0,40,10):
        TPrate = ROC.true_pos_rate(pos_scores, neg_scores, FPrate)
        obj_fnc_val = obj_fnc_val + TPrate/100
    return obj_fnc_val
def calc_new_scores(pos_align_array,neg_align_array, subMatrix):
    print len(pos_align_array)
    new_pos_scores=np.zeros(len(pos_align_array))
    new_neg_scores=np.zeros(len(neg_align_array))
    for i in range(0,len(pos_align_array)):
        new_pos_scores [i]= np.sum(pos_align_array[i]*subMatrix)
    for i in range(0,len(neg_align_array)):
        new_neg_scores [i]= np.sum(neg_align_array[i]*subMatrix)
    return new_pos_scores, new_neg_scores


def generate_newsubMatrix(origsubMatrix):
    
    
    newsubMatrix = 0
    return newsubMatrix
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
    
    [sub_Matrixdict, origSubMatrix] = subMdict.mk_dict(home+subMatrixFile)

    #[pos_scores, pos_align_array] = slSW.scores_from_seq_list(home, pos_seq_list_file, sub_Matrixdict, origSubMatrix, gap_init, gap_ext)
    #[neg_scores, neg_align_array]= slSW.scores_from_seq_list(home, neg_seq_list_file, sub_Matrixdict, origSubMatrix, gap_init, gap_ext)
    #RWFile.writetxt(pos_align_array, home, 'pos_align_array.txt')
    #RWFile.writetxt(neg_align_array, home, 'neg_align_array.txt')
    
    neg_align_array=np.genfromtxt(home+'neg_align_array.txt', comments="#",  unpack=False)
    pos_align_array=np.genfromtxt(home+'pos_align_array.txt', comments="#",  unpack=False)
    print pos_align_array
    [new_pos_scores, new_neg_scores] = calc_new_scores(pos_align_array, neg_align_array, origSubMatrix)
    print new_pos_scores
    print new_neg_scores
    
    print 'done'
    
    return 1

if __name__ == '__main__':
    main()
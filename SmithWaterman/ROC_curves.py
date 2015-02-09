'''
Created on Feb 8, 2015

@author: gogqou

calculate FP and TP rates and graphs them 
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

def false_pos_rate(pos_score_list, neg_score_list, TPrate):
    percentile_pos_scores = np.percentile(pos_score_list, 100-TPrate)
    above_threshold= [k for k in neg_score_list if k>=percentile_pos_scores]
    false_pos_rate = len(above_threshold)/float(len(neg_score_list))
    return false_pos_rate

def true_pos_rate(pos_score_list, neg_score_list, FPrate):
    percentile_neg_scores = np.percentile(neg_score_list, 100-FPrate)
    above_threshold= [k for k in pos_score_list if k>=percentile_neg_scores]
    true_pos_rate = len(above_threshold)/float(len(pos_score_list))
    return true_pos_rate
def ROC_graph(pos_scores, neg_scores):
    ROC_array = np.zeros([50,2])
    i = 0
    for TPrate in range(0, 100, 2):
        FPrate = false_pos_rate(pos_scores,neg_scores, TPrate)
        ROC_array[i]=[TPrate, FPrate]
        i= i+1  
    return ROC_array
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
    subMatrixFile_list = ['BLOSUM50', 'BLOSUM62', 'MATIO', 'PAM100', 'PAM250']
    gap_init = 13
    gap_ext = 2
    [sub_Matrix, origSubMatrix, AAlist] = subMdict.mk_dict(home+subMatrixFile)
    
    
    # for testing different scoring matrices

    for i in range(0,len(subMatrixFile_list)):
        [sub_Matrix, origsubMatrix, AAlist] = subMdict.mk_dict(home+subMatrixFile_list[i]) 
        [pos_scores, pos_align_array] = slSW.scores_from_seq_list(home, pos_seq_list_file, sub_Matrix,origSubMatrix, gap_init, gap_ext)
        [neg_scores, neg_align_array]= slSW.scores_from_seq_list(home, neg_seq_list_file, sub_Matrix, origSubMatrix, gap_init, gap_ext)
        ROC_array = ROC_graph(pos_scores, neg_scores)
        pylab.plot(ROC_array[:,1], ROC_array[:,0]/100, label = subMatrixFile_list[i])

#     [pos_scores, pos_align_array] = slSW.scores_from_seq_list(home, pos_seq_list_file, sub_Matrix, origSubMatrix, gap_init, gap_ext)
#     [neg_scores, neg_align_array]= slSW.scores_from_seq_list(home, neg_seq_list_file, sub_Matrix, origSubMatrix, gap_init, gap_ext)
#     ROC_array = ROC_graph(pos_scores, neg_scores,)
#     pylab.plot(ROC_array[:,1], ROC_array[:,0]/100, label = subMatrixFile)
    pylab.axis([0,1,0,1])
    pylab.legend(loc = 'lower right')
    pylab.ylabel('True Positive Rate')
    pylab.xlabel('False Positive Rate')
    pylab.savefig('ROC.png')
    print 'done'
    
    return 1

if __name__ == '__main__':
    main()
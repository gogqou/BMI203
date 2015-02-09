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
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d

def SW_one_round(seq1, seq2, sub_Matrix, origSubMatrix, gap_init, gap_ext):
    
    [sim_Matrix, pointers] = SW.similarity_matrix(seq1,seq2, sub_Matrix, gap_init, gap_ext)
    
    [aligned_sequence, fitted_seq1, fitted_seq2, sim_Matrix, score, count_array] = SW.trace_aligned_seq(seq1, seq2, sim_Matrix, pointers, origSubMatrix, sub_Matrix)

    #print aligned_sequence
    #print fitted_seq1
    #print fitted_seq2
    #print 'score = ', score
    #print count_array
    
    return score, count_array
def list_of_sequences_from_txt(sequence_list_file):
    seq_list = []
    lines = RWFile.readtxt(sequence_list_file)
    for line in lines:
        seq_list.append(line.split(' '))
    return seq_list

def scores_from_seq_list(home,seq_list_file,sub_Matrix, origSubMatrix, gap_init, gap_ext):
    seq_list = list_of_sequences_from_txt(seq_list_file)
    print seq_list
    score_list = np.zeros([len(seq_list),1])
    #has the count array for each pairing in the sequence list
    alignment_array = np.zeros([len(seq_list),24,24])
    for i in range(0,len(seq_list)):
        seq_locations = seq_list[i]
        seqfile1=seq_locations[0]
        seqfile2=seq_locations[1]
        seq1=rFasta.read_fa(home+seqfile1)
        seq2=rFasta.read_fa(home+seqfile2)
        [score, count_array] = SW_one_round(seq1, seq2, sub_Matrix,origSubMatrix, gap_init, gap_ext)
        score_list[i] = score
        alignment_array[i] = count_array
        
    return score_list, alignment_array


def false_pos_rate(pos_score_list, neg_score_list):
    
    seventy_percentile_pos_scores = np.percentile(pos_score_list, 30)
    above_threshold= [k for k in neg_score_list if k>=seventy_percentile_pos_scores]
    false_pos_rate = len(above_threshold)/float(len(neg_score_list))
    return false_pos_rate


def graph_cost_array():
    
    arrayfile = sys.argv[3]+'FPrateBLOSUM50.txt'
    array = RWFile.readcsv(arrayfile)
    fig=plt.figure()
    ax= fig.add_subplot(111, projection = '3d')
    x= array[:,0]
    y= array[:,1]
    z= array[:,2]
    cset=ax.plot(x,y,z, ls= 'None', marker = 'o')
    ax.set_xlabel('Gap Initiation Cost')
    ax.set_ylabel('Gap Extension Cost')
    ax.set_zlabel('False Positive Rate')
    plt.show()
    return 1


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
    #FP_array = np.zeros([120,3])
    FP_array = []
    i=0
    #gap_init_cost, gap_ext_cost
    [sub_Matrix, origSubMatrix] = subMdict.mk_dict(home+subMatrixFile)
    
    '''
    #for testing gap costs
    
    for gap_init in range(1,21):
        for gap_ext in range(1,6):
            
            [pos_score_list, pos_alignment_array] = scores_from_seq_list(home,pos_seq_list_file, sub_Matrix,origSubMatrix, gap_init, gap_ext)
            [neg_score_list, neg_alignment_array] = scores_from_seq_list(home,neg_seq_list_file, sub_Matrix,origSubMatrix, gap_init, gap_ext)
            
            false_pos_rt= false_pos_rate(pos_score_list, neg_score_list)
            FP_array[i] = [gap_init, gap_ext, false_pos_rt]
            i=i+1
        
    print FP_array
    RWFile.writetxt(FP_array, home, 'FPrate SubMatrices.txt')
    graph_cost_array()
    '''
    
    # for testing input matrices
    for i in range(0,len(subMatrixFile_list)):
        [sub_Matrix, origsubMatrix] = subMdict.mk_dict(home+subMatrixFile_list[i])
    
        [pos_score_list, pos_alignment_array] = scores_from_seq_list(home,pos_seq_list_file, sub_Matrix, origSubMatrix,13, 2)
        [neg_score_list, neg_alignment_array] = scores_from_seq_list(home,neg_seq_list_file, sub_Matrix, origSubMatrix, 13, 2)
        false_pos_rt= false_pos_rate(pos_score_list, neg_score_list)
        FP_array.append([subMatrixFile_list[i], false_pos_rt])
        
    print FP_array
  
            
          
            
            
    print 'done'
if __name__ == '__main__':
    main()

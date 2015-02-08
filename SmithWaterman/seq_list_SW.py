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

def SW_one_round(seq1, seq2, sub_Matrix, gap_init, gap_ext):
    
    [sim_Matrix, pointers] = SW.similarity_matrix(seq1,seq2, sub_Matrix, gap_init, gap_ext)
    
    [aligned_sequence, fitted_seq1, fitted_seq2, sim_Matrix, score] = SW.trace_aligned_seq(seq1, seq2, sim_Matrix, pointers)
    #print aligned_sequence
    #print fitted_seq1
    #print fitted_seq2
    #print 'score = ', score
    
    return score
def list_of_sequences_from_txt(sequence_list_file):
    seq_list = []
    lines = RWFile.readtxt(sequence_list_file)
    for line in lines:
        seq_list.append(line.split(' '))
    return seq_list

def scores_from_seq_list(home,seq_list_file,sub_Matrix, gap_init, gap_ext):
    seq_list = list_of_sequences_from_txt(seq_list_file)
    print seq_list
    score_list = np.zeros([len(seq_list),1])
    for i in range(0,len(seq_list)):
        seq_locations = seq_list[i]
        seqfile1=seq_locations[0]
        seqfile2=seq_locations[1]
        seq1=rFasta.read_fa(home+seqfile1)
        seq2=rFasta.read_fa(home+seqfile2)
        score = SW_one_round(seq1, seq2, sub_Matrix, gap_init, gap_ext)
        score_list[i] = score
        
    return score_list


def false_pos_rate(pos_score_list, neg_score_list):
    
    seventy_percentile_pos_scores = np.percentile(pos_score_list, 30)
    above_threshold= [k for k in neg_score_list if k>=seventy_percentile_pos_scores]
    false_pos_rate = len(above_threshold)/float(len(neg_score_list))
    return false_pos_rate
    
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
    FP_array = np.zeros([100,3])
    i=0
    #gap_init_cost, gap_ext_cost
    sub_Matrix = subMdict.mk_dict(home+subMatrixFile)
    for gap_init in range(1,21):
        for gap_ext in range(1,6):
            
            pos_score_list = scores_from_seq_list(home,pos_seq_list_file, sub_Matrix, gap_init, gap_ext)
            neg_score_list = scores_from_seq_list(home,neg_seq_list_file, sub_Matrix, gap_init, gap_ext)
            
            false_pos_rt= false_pos_rate(pos_score_list, neg_score_list)
            FP_array[i] = [gap_init, gap_ext, false_pos_rt]
            i=i+1
    print FP_array
    #print np.percentile(score_list, 70)
    #RWFile.writecsv(score_list, home, seq_list_name+subMatrixFile+'score_list.txt')
    return 'done'
if __name__ == '__main__':
    main()
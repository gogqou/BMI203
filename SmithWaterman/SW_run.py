'''
Created on Feb 4, 2015

@author: gogqou

Implements Smith Waterman Algorithm for local sequence alignment



Implement the Smith-Waterman algorithm and instrument the code such
that it can use any scoring matrix provided (i.e. will read it in from
a separate file). 

Scoring matrix includes gap-scoring scheme and substitution matrix; needs to be able to read in from outside file.



Take input: seq1 and seq2
Substitution matrix file--should call a parsing function to read in the file
(includes gap costs)

Output: filled out similarity matrix H and best sequence alignment and indices from seq1 and seq2


'''

import sub_matrix_mk_dict as subMdict
import numpy as np
import sys
import read_fasta as rFasta
def similarity_matrix(seq1,seq2,substitution_Matrix_dictionary, gap_init, gap_ext):
    
    similarity_matrix = np.zeros([len(seq1), len(seq2)]) #initializes similarity matrix with zeros everywhere; will replace all but first row and column
    H = similarity_matrix.copy()
    C = substitution_Matrix_dictionary # for easier typing later on
    pointers = {} #dictionary to keep pointers of where each space got its value
    Gappointers = {}
    for i in range(0,len(seq1)):
        Gappointers[(i,0)] = False
    for j in range(0,len(seq2)):
        Gappointers[(0,j)] = False   
    
    for j in range(1,len(seq2)):
        for i in range(1,len(seq1)):
        
            #calculates the value for each of the possible moves so it's easier to compare them later
            if Gappointers[(i,j-1)] is False:
                left = H[i, j-1]-gap_init 
            else:  
                left = H[i, j-1] -gap_ext
            if Gappointers[(i-1,j)] is False:
                up = H[i-1, j]-gap_init
            else:  
                up = H[i-1, j]-gap_ext

            #gap in seq1; if gap=1, then there was already a gap started, no gap_init cost
            #if gap = 0, then choosing either left or up requires starting a gap, multiplying -1 by the gap_init cost means subtracting it
            matchscore_dict= C[seq1[i]+seq2[j]] #necessary because the value associated with this key includes the score as well as the index in the original scoring matrix
            diagonal = H[i-1,j-1]+matchscore_dict[0] #match
            H[i,j] = max(left, up, diagonal, 0) 
            if H[i,j]== left:
                pointers[(i,j)] = 'left'
                Gappointers[(i,j)] = True
            elif H[i,j]== up:
                pointers[(i,j)] = 'up'
                Gappointers[(i,j)] = True
            elif H[i,j]== diagonal:
                pointers[(i,j)] = 'diagonal'
                Gappointers[(i,j)] = False
            else:
                pointers[(i,j)]='0'
                Gappointers[(i,j)] = False
    return H, pointers

def trace_aligned_seq(seq1, seq2, similarity_matrix, pointers, origsubMatrix, subMatrixdict):
    H=similarity_matrix.copy()
    C = subMatrixdict # for easier typing later on
    seq = ''
    newseq1 = ''
    newseq2= ''
    indices = np.unravel_index(H.argmax(), H.shape)
    count_array = np.zeros(origsubMatrix.shape)
    i = indices[0]
    j = indices[1]
    start=i
    end= j
    #seq=seq+seq2[j]
    #newseq1=newseq1+seq2[j]
    #newseq2=newseq2+seq2[j]
    score = H[i,j]
    gapcount = 0
    gap_open = 0
    altscore = 0
    while H[i,j]>0 and pointers[(i,j)] is not '0':

        if pointers[(i,j)] is 'left':
            #print 'left'
            seq=seq+seq2[j]
            newseq2 = newseq2 + seq2[j]
            newseq1 = newseq1+'-'
            j = j-1
            if H[i,j]>0 and pointers[(i,j)] is 'diagonal':           
                gap_open = gap_open +1
            else:
                gapcount = gapcount+1
        elif pointers[(i,j)] is 'up':
            #print 'up'
            seq=seq+seq1[i]
            newseq1 = newseq1 + seq1[i]
            newseq2 = newseq2+'-'            
            i = i-1
            if H[i,j]>0 and pointers[(i,j)] is 'diagonal':            
                gap_open = gap_open +1
            else:
                gapcount = gapcount+1
                
        elif pointers[(i,j)] is 'diagonal': 
            #print 'diag'
            seq=seq+seq1[i]
            newseq1=newseq1+seq1[i]
            newseq2=newseq2+seq2[j]
            #finds the AAtoAA comparison in the dictionary
            #2nd and 3rd entries are the indices in the original substitution matrix
            #only necessary in a match because in other cases, the cost is just gap initiation or extension
            #at these same indices, increment count of the times that that score was used
   
            matchscore_dict= C[seq1[i]+seq2[j]]
            #print matchscore_dict[0]
            altscore = altscore + matchscore_dict[0]
            count_array[matchscore_dict[1],matchscore_dict[2]]= count_array[matchscore_dict[1],matchscore_dict[2]] + 1
            i = i-1
            j = j-1
    
    seq=seq[::-1]
    newseq1=newseq1[::-1]
    newseq2=newseq2[::-1]
    altscore = altscore - gap_open*13 - gapcount*3
    #print 'altscore', altscore
    count_array[0,23] = gap_open
    count_array[23,23]=gapcount
    #H[start,end]=0
    # to test ROC, average by minimum len of the compared pair
    #score = score/min(len(seq1), len(seq2))
    return seq, newseq1, newseq2, H, score, count_array


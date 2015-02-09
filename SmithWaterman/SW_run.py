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
    H = similarity_matrix
    C = substitution_Matrix_dictionary # for easier typing later on
    pointers = {} #dictionary to keep pointers of where each space got its value
    gap = 0
    for i in range(1,len(seq1)):
        for j in range(1,len(seq2)):
            #calculates the value for each of the possible moves so it's easier to compare them later
            
            #left = H[i-1,j]+C[seq1[i]+'*'] +(gap-1)*gap_init - gap*gap_ext #gap in seq2
            #up=H[i,j-1]+C['*'+seq2[j]] +(gap-1)*gap_init - gap*gap_ext 
            left = H[i-1,j]+(gap-1)*gap_init - gap*gap_ext #gap in seq2
            up=H[i,j-1] +(gap-1)*gap_init - gap*gap_ext 
            #gap in seq1; if gap=1, then there was already a gap started, no gap_init cost
                       
            #if gap = 0, then choosing either left or up requires starting a gap, multiplying -1 by the gap_init cost means subtracting it
            matchscore_dict= C[seq1[i]+seq2[j]] #necessary because the value associated with this key includes the score as well as the index in the original scoring matrix
            diagonal = H[i-1,j-1]+matchscore_dict[0] #match
            
            H[i,j] = max(left, up, diagonal, 0) 
            
            if H[i,j]== left:
                pointers[str(i)+str(j)] = 'left'
                gap = 1
            elif H[i,j]== up:
                pointers[str(i)+str(j)] = 'up'
                gap = 1
            else:
                pointers[str(i)+str(j)] = 'diagonal'
                gap = 0
            j=j+1 #increment across the matrix
        i=i+1        #then increment down the matrix, so you know all values you need will be available
            
    return H, pointers

def trace_aligned_seq(seq1, seq2, similarity_matrix, pointers, origsubMatrix, subMatrixdict):
    H=similarity_matrix
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
    seq=seq+seq2[j]
    newseq1=newseq1+seq2[j]
    newseq2=newseq2+seq2[j]
    score = H[i,j]
    
    #finds the AAtoAA comparison in the dictionary
    #2nd and 3rd entries are the indices in the original substitution matrix
    matchscore_dict= C[seq1[i]+seq2[j]]
    #at these same indices, increment count of the times that that score was used
    count_array[matchscore_dict[1],matchscore_dict[2]]= count_array[matchscore_dict[1],matchscore_dict[2]] + 1
    while H[i,j]>0:
        ind1=str(i)
        ind2=str(j)
        
        if pointers[ind1+ind2] is 'left':
            j = j-1
            seq=seq+seq2[j]
            newseq2 = newseq2 + seq2[j]
            newseq1 = newseq1+'-'
            score = score+ H[i,j]
        elif pointers[ind1+ind2] is 'up':
            i = i-1
            seq=seq+seq1[i]
            newseq1 = newseq1 + seq1[i]
            newseq2 = newseq2+'-'
            score = score+ H[i,j]
        else: 
            i = i-1
            j = j-1
            seq=seq+seq1[i]
            newseq1=newseq1+seq1[i]
            newseq2=newseq2+seq2[j]
            score = score+ H[i,j]
            #only necessary in a match because in other cases, the cost is just gap initiation or extension
            matchscore_dict= C[seq1[i]+seq2[j]]
            count_array[matchscore_dict[1],matchscore_dict[2]]= count_array[matchscore_dict[1],matchscore_dict[2]] + 1
    seq=seq[::-1]
    newseq1=newseq1[::-1]
    newseq2=newseq2[::-1]
    H[start,end]=0
    # to test ROC, average by minimum len of the compared pair
    #score = score/min(len(seq1), len(seq2))
    return seq, newseq1, newseq2, H, score, count_array


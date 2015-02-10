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
    gap = 0
    for i in range(1,len(seq1)):
        for j in range(1,len(seq2)):
            #calculates the value for each of the possible moves so it's easier to compare them later
            
            #left = H[i-1,j]+C[seq1[i]+'*'] +(gap-1)*gap_init - gap*gap_ext #gap in seq2
            #up=H[i,j-1]+C['*'+seq2[j]] +(gap-1)*gap_init - gap*gap_ext 
            #left = H[i-1,j]+(gap-1)*gap_init - gap*gap_ext #gap in seq2
            left = max(H[i-1,j]- gap_ext, H[i-1, j-1]-gap_init- gap_ext)
            #up=H[i,j-1] +(gap-1)*gap_init - gap*gap_ext 
            
            up = max(H[i, j-1]-gap_ext, H[i-1,j-1]-gap_init-gap_ext)
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
            elif H[i,j]== diagonal:
                pointers[str(i)+str(j)] = 'diagonal'
                gap = 0
            else:
                pointers[str(i)+str(j)]='0'
            j=j+1 #increment across the matrix
        i=i+1        #then increment down the matrix, so you know all values you need will be available
            
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
    while H[i,j]>0 and pointers[str(i)+str(j)] is not '0':
        ind1=str(i)
        ind2=str(j)
        if pointers[ind1+ind2] is 'left':
            print 'left'
            seq=seq+seq2[j]
            newseq2 = newseq2 + seq2[j]
            newseq1 = newseq1+'-'
            matchscore_dict= C[seq2[j]+seq2[j]]
            matchscore_dict= C[seq1[i]+seq1[i]]
            print seq1[i]
            print seq2[j]
            print H[i,j]
            #print H[i,j]-H[i, j-1]
            print matchscore_dict[0]
            altscore = altscore + matchscore_dict[0]
            count_array[matchscore_dict[1],matchscore_dict[2]]= count_array[matchscore_dict[1],matchscore_dict[2]] + 1
            j = j-1
            if H[i,j]>0 and pointers[ind1+str(j)] is 'diagonal':           
                gap_open = gap_open +1
            else:
                gapcount = gapcount+1
        elif pointers[ind1+ind2] is 'up':
            print 'up'
            seq=seq+seq1[i]
            newseq1 = newseq1 + seq1[i]
            newseq2 = newseq2+'-'
            print seq1[i]
            print seq2[j]
            
            matchscore_dict= C[seq1[i]+seq1[i]]
            print H[i,j]
            #print H[i,j]-H[i-1, j]
            print matchscore_dict[0]
            altscore = altscore + matchscore_dict[0]
            count_array[matchscore_dict[1],matchscore_dict[2]]= count_array[matchscore_dict[1],matchscore_dict[2]] + 1
            i = i-1
            if H[i,j]>0 and pointers[str(i)+ind2] is 'diagonal':            
                gap_open = gap_open +1
            else:
                gapcount = gapcount+1
                
        elif pointers[ind1+ind2] is 'diagonal': 
            print 'diag'
            seq=seq+seq1[i]
            newseq1=newseq1+seq1[i]
            newseq2=newseq2+seq2[j]
            print seq1[i]
            print seq2[j]
            #finds the AAtoAA comparison in the dictionary
            #2nd and 3rd entries are the indices in the original substitution matrix
            #only necessary in a match because in other cases, the cost is just gap initiation or extension
            #at these same indices, increment count of the times that that score was used
            print H[i,j]
            matchscore_dict= C[seq1[i]+seq2[j]]
            print matchscore_dict[0]
            altscore = altscore + matchscore_dict[0]
            count_array[matchscore_dict[1],matchscore_dict[2]]= count_array[matchscore_dict[1],matchscore_dict[2]] + 1
            i = i-1
            j = j-1
        else: 
            print 'new alignment'
            break
    
    matchscore_dict= C[seq1[i]+seq2[j]]
    altscore = altscore + matchscore_dict[0]
    count_array[matchscore_dict[1],matchscore_dict[2]]= count_array[matchscore_dict[1],matchscore_dict[2]] + 1
    seq=seq+seq1[i]
    newseq1=newseq1+seq1[i]
    newseq2=newseq2+seq2[j]
    seq=seq[::-1]
    newseq1=newseq1[::-1]
    newseq2=newseq2[::-1]
    
    print 'altscore', altscore
    count_array[0,23] = gap_open
    count_array[23,23]=gapcount
    #H[start,end]=0
    # to test ROC, average by minimum len of the compared pair
    #score = score/min(len(seq1), len(seq2))
    return seq, newseq1, newseq2, H, score, count_array


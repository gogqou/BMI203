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
    score = 0
    for i in range(1,len(seq1)):
        for j in range(1,len(seq2)):
            #calculates the value for each of the possible moves so it's easier to compare them later
            
            
            
            # NEED TO ADD IN A CHECK FOR IF THIS IS A GAP INIT OR EXTENSION
            
            
            
            
            
            left = H[i-1,j]+C[seq1[i]+'*'] #gap in seq2
            up=H[i,j-1]+C['*'+seq2[j]] #gap in seq1
            diagonal = H[i-1,j-1]+C[seq1[i]+seq2[j]] #match
            
            H[i,j] = max(left, up, diagonal, 0) 
            
            if H[i,j]== left:
                pointers[str(i)+str(j)] = 'left'
            elif H[i,j]== up:
                pointers[str(i)+str(j)] = 'up'
            else:
                pointers[str(i)+str(j)] = 'diagonal'
            j=j+1 #increment across the matrix
        i=i+1        #then increment down the matrix, so you know all values you need will be available
            
    return H, pointers

def trace_aligned_seq(seq1, seq2, similarity_matrix, pointers):
    H=similarity_matrix
    seq = ''
    newseq1 = ''
    newseq2= ''
    indices = np.unravel_index(H.argmax(), H.shape)
    i = indices[0]
    j = indices[1]
    start=i
    end= j
    seq=seq+seq2[j]
    newseq1=newseq1+seq2[j]
    newseq2=newseq2+seq2[j]
    score = H[i,j]
    while H[i,j]>0:
        ind1=str(i)
        ind2=str(j)
        
        score = score+ H[i,j]
        if pointers[ind1+ind2] is 'left':
            j = j-1
            seq=seq+seq2[j]
            newseq2 = newseq2 + seq2[j]
            newseq1 = newseq1+'-'
        elif pointers[ind1+ind2] is 'up':
            i = i-1
            seq=seq+seq1[i]
            newseq1 = newseq1 + seq1[i]
            newseq2 = newseq2+'-'
        else: 
            i = i-1
            j = j-1
            seq=seq+seq1[i]
            newseq1=newseq1+seq1[i]
            newseq2=newseq2+seq2[j]
    seq=seq[::-1]
    newseq1=newseq1[::-1]
    newseq2=newseq2[::-1]
    H[start,end]=0
    return seq, newseq1, newseq2, H, score


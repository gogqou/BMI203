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
def similarity_matrix(seq1,seq2,substitution_Matrix_dictionary):
    
    similarity_matrix = np.zeros([len(seq1), len(seq2)]) #initializes similarity matrix with zeros everywhere; will replace all but first row and column
    H = similarity_matrix
    C = substitution_Matrix_dictionary # for easier typing later on
    pointers = {} #dictionary to keep pointers of where each space got its value
    for i in range(1,len(seq1)):
        for j in range(1,len(seq2)):
            #calculates the value for each of the possible moves so it's easier to compare them later
            left = H[i-1,j]+C[seq1[i]+'*'] #gap in seq2
            up=H[i,j-1]+C['*'+seq2[j]] #gap in seq1
            diagonal = H[i-1,j-1]+C[seq1[i]+seq2[j]] #match
            
            H[i,j] = max(left, up, diagonal, 0) #prevents negative values by taking the max with 0
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
    indices = np.unravel_index(H.argmax(), H.shape)
    i = indices[0]
    j = indices[1]
    print i,j
    print seq2[j]
    end_seq1_index = i
    end_seq2_index = j
    seq=seq+seq2[j]
    while H[i,j]>0:
        ind1=str(i)
        ind2=str(j)
        if pointers[ind1+ind2] is 'left':
            j = j-1
            seq=seq+seq2[j]
        elif pointers[ind1+ind2] is 'up':
            i = i-1
            seq=seq+seq1[i]
        else: 
            i = i-1
            j = j-1
            seq=seq+seq1[i]
    seq=seq[::-1]
    start_seq1_index = i
    start_seq2_index = j
    return seq, start_seq1_index,start_seq2_index,end_seq1_index,end_seq2_index
def main():
    if len(sys.argv)>4:
        print 'provide sequences to align and substitution matrix '
        sys.exit()
    seq1=rFasta.read_fa(sys.argv[1])
    print seq1, len(seq1)
    seq2=rFasta.read_fa(sys.argv[2])
    print seq2, len(seq2)
    subMatrixFile = sys.argv[3]
    
    sub_Matrix = subMdict.mk_dict(subMatrixFile)
    
    [sim_Matrix, pointers] = similarity_matrix(seq1,seq2, sub_Matrix)
    print sim_Matrix
    [aligned_sequence, seq1a,seq2a, seq1b, seq2b] = trace_aligned_seq(seq1, seq2, sim_Matrix, pointers)
    print aligned_sequence
    start = min(seq1a, seq2a)
    end = max(seq1b, seq2b)
    print seq1[start:end]
    print seq2[start:end]
    return 'done'
if __name__ == '__main__':
    main()
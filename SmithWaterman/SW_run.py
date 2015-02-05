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
def similarity_matrix(seq1,seq2,substitution_Matrix_dictionary):
    
    similarity_matrix = np.zeros([len(seq1), len(seq2)])
    print similarity_matrix
    H = similarity_matrix
    C = substitution_Matrix_dictionary
    pointers = {}

    print len(seq1)
    for i in range(1,len(seq1)):
        for j in range(1,len(seq2)):
            left = H[i-1,j]+C[seq1[i]+'*']
            up=H[i,j-1]+C['*'+seq2[j]]
            diagonal = H[i-1,j-1]+C[seq1[i]+seq2[j]]
            
            H[i,j] = max(left, up, diagonal, 0)
            if H[i,j]== left:
                pointers[str(i)+str(j)] = 'left'
            elif H[i,j]== up:
                pointers[str(i)+str(j)] = 'up'
            else:
                pointers[str(i)+str(j)] = 'diagonal'
            j=j+1
        i=i+1        
        
    print H    
    print pointers 
    
    return similarity_matrix
def main():
    if len(sys.argv)>4:
        print 'provide sequences to align and substitution matrix '
        sys.exit()
    seq1=sys.argv[1]
    seq2=sys.argv[2]
    subMatrixFile = sys.argv[3]
    
    sub_Matrix = subMdict.mk_dict(subMatrixFile)
    
    sim_Matrix = similarity_matrix(seq1,seq2, sub_Matrix)
    
    return 'done'
if __name__ == '__main__':
    main()
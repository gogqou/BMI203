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

import sub_matrix_mk_dict as subMatrixdict

import sys

def main():
    if len(sys.argv)>3:
        print 'provide sequences to align and substitution matrix '
        sys.exit()
    subMatrixFile = sys.argv[1]
    subMatrixdict.submatrix_dict(subMatrixFile)
    
    
    return 'done'
if __name__ == '__main__':
    main()
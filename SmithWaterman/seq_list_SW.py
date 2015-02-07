'''
Created on Feb 7, 2015

@author: gogqou
'''

import sub_matrix_mk_dict as subMdict
import numpy as np
import sys
import read_fasta as rFasta

import SW_run as SW
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
    [sim_Matrix, pointers] = SW.similarity_matrix(seq1,seq2, sub_Matrix)
            
    [aligned_sequence, fitted_seq1, fitted_seq2, sim_Matrix, score] = SW.trace_aligned_seq(seq1, seq2, sim_Matrix, pointers)
    print aligned_sequence
    print fitted_seq1
    print fitted_seq2
    print 'score = ', score
    return 'done'
if __name__ == '__main__':
    main()
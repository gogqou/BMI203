'''
Created on Feb 4, 2015

@author: gogqou

reads in substitution matrix and gives back a dictionary of [a,b] keys that are associated with substitution values

'''

import readwrite_file as rwF
def submatrix_dict(inputFile):
    
    sub_dict = {}
    lines=rwF.readtxt(inputFile)
    for line in lines:
        print line
    return sub_dict
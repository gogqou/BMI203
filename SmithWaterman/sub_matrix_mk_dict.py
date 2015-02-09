'''
Created on Feb 4, 2015

@author: gogqou

reads in substitution matrix and gives back a dictionary of [a,b] keys that are associated with substitution values

'''

import readwrite_file as rwF

from numpy import genfromtxt
import re
import numpy as np
def mk_dict(inputFile):
    
    sub_dict = {}
    temp = []
    lines =rwF.readtxt(inputFile)
    for line in lines:
        if '#' not in line:
            if 'A' in line:
                names = line
    names= re.sub(r" ", "", names)
    AAlist=list(names)
    submatrixVals=genfromtxt(inputFile, comments="#",  unpack=False)
    submatrixVals = submatrixVals[1:]
    AArange=len(AAlist)
    for i in range(AArange):
        for j in range(AArange):
            sub_dict[AAlist[i]+AAlist[j]] = [submatrixVals[i,j], i, j]
    return sub_dict, submatrixVals, AAlist

def mk_dict_np(inputFile, AAlist):
    AArange = len(AAlist)
    optimized_Matrix = np.load(inputFile)
    sub_dict = {}
    for i in range(AArange):
        for j in range(AArange):
            sub_dict[AAlist[i]+AAlist[j]] = [optimized_Matrix[i,j], i, j]
    return sub_dict, optimized_Matrix, AAlist
    
    
'''
Created on Feb 4, 2015

@author: gogqou

reads in substitution matrix and gives back a dictionary of [a,b] keys that are associated with substitution values

'''

import readwrite_file as rwF

from numpy import genfromtxt
import re
def submatrix_dict(inputFile):
    
    sub_dict = {}
    temp = []
    i= 0
    lines =rwF.readtxt(inputFile)
    for line in lines:
        if '#' not in line:
            if 'A' in line:
                names = line
    names= re.sub(r"\W", "", names)
    print names
    a=list(names)
    print a
    a=genfromtxt(inputFile, comments="#", skip_header=0, unpack=True)
    #lines=rwF.readtxt(inputFile)
    print a[1]
            
    return sub_dict
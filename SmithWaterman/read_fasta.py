'''
Created on Feb 5, 2015

@author: gogqou

reads in .fa file and outputs a string for the sequence
'''

from numpy import genfromtxt
import readwrite_file as rwF
def read_fa(inputFile):
    sequence = ''
    lines =rwF.readtxt(inputFile)
    
    for line in lines:
        if '>' not in line:
            sequence=sequence+line
    return sequence
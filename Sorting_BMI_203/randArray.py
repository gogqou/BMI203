'''
Created on Jan 7, 2015

@author: Guanqing Ou


allows generation of random array of given size
'''

import numpy as np
import random
def randArray(array_range,size):
    
    rArray = random.sample(range(1,array_range),size)
    
    
    return rArray

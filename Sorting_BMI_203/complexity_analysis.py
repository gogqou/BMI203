'''
Created on Jan 7, 2015

@author: gogqou

Run sort algorithms multiple times for arrays of multiple sizes to allow analysis of algorithm complexity

Fits output of array size n vs. time required to run and polynomial or log curve fits

Then plots both raw data and fitted curve data 
'''
import sys
import bubbleSort
import randArray
import numpy as np
import time
from pylab import plot, show, ylim, yticks
from datetime import datetime
import quickSort
from numpy.polynomial import polynomial as P
import math
from scipy.optimize import curve_fit
from pylab import *

def curvefit_choice(data, fitfunc):
    xdata=data[:,0]
    ydata=data[:,1]
    x = np.array(xdata)
    y = np.array(ydata)
    if fitfunc is 'log':
        
        def func(x,a):
            return a*x*np.log(x)
        popt, pcov = curve_fit(func, x, y)
        print popt
        print pcov
        xx = np.linspace(100,1000,1000)
        yy = func(xx, *popt)
        plot(xdata,ydata, 'ko')
        plot(xx,yy)
        savefig('log_fit_quicksort.png')
    elif fitfunc is 'polynomial':
        
    
        coef,stats = P.polyfit(xdata,ydata,deg=4, full=True)
        print coef
        print stats
        coef,stats = P.polyfit(xdata,ydata,deg=3, full=True)
        print coef
        print stats
        coef,stats = P.polyfit(xdata,ydata,deg=2, full=True)
        print coef
        print stats
        coef,stats = P.polyfit(xdata,ydata,deg=1, full=True)
        print coef
        print stats
    else: 
        print 'unexpected function name for fit, try polynomial or log'
    return 1
def main():
    complex_array= np.zeros((10000,2)) #allocate space for the array of n, and time to sort array with n elements
    #print complex_array
    index = 0
    array_range = 10000
    for array_size in range(100,1000,100): #iterate over the range (100,1000) in increments of 100
        for arraynum in range(2): #randomly generate 100 arrays of the given array_size
            array=randArray.randArray(array_range,array_size)
            #[sortedArray, runTime, assign_count, cond_count]=bubbleSort.bubbleSort(array) #bubble sort and get the sorted array, the execution time, assignment and conditionals count
            [sortedArray, runTime, assign_count, cond_count]=quickSort.quickSort(array)            
           # print assign_count, cond_count
            complex_array[index,:]= [array_size, runTime]
            index = index + 1
    complex_array = complex_array[0:index]
    print complex_array
    xdata=complex_array[:,0]
    ydata=complex_array[:,1]
    curvefit_choice(complex_array, 'log')
    return complex_array
if __name__ == '__main__':
    main()
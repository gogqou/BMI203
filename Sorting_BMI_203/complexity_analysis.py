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
from datetime import datetime
import quickSort
from numpy.polynomial import polynomial as P
import math
from scipy.optimize import curve_fit
import pylab

def curvefit_choice(data, fitfunc):
    #curve fit based on whether expected is log or polynomial time
        
    xdata=data[:,0]
    ydata=data[:,1]
    x = np.array(xdata)
    y = np.array(ydata)
    
    if fitfunc is 'log':
        
        def func(x,a):
            return a*x*np.log(x)
        popt, pcov = curve_fit(func, x, y)
        #print popt
        #print pcov
        xx = np.linspace(100,1000,1000)
        yy = func(xx, *popt)
        pylab.plot(xdata,ydata, 'ko', label = 'quickSort')
        pylab.plot(xx,yy, label = 'axlogx() log fit')
        pylab.legend()
        pylab.ylabel('Runtime (s)')
        pylab.xlabel('Array Size')
        pylab.savefig('log_fit_quicksort.png')
        pylab.plt.clf()
        
    elif fitfunc is 'polynomial':
        def func(x,a,b,c):
            return a*x**2 + b*x + c

        popt, pcov = curve_fit(func, x, y)
        xx = np.linspace(100,1000,1000)
        yy = func(xx, *popt)
        
        pylab.plot(xdata,ydata, 'ko', label = 'bubbleSort')
        pylab.plot(xx,yy, label='ax^2 + bx + c polynomial fit')
        pylab.legend()
        pylab.ylabel('Runtime (s)')
        pylab.xlabel('Array Size')
        pylab.savefig('polyfit_bubbleSort.png')
        pylab.plt.clf()
    else: 
        print 'unexpected function name for fit, try polynomial or log'
    return 1

def graph(array, name):
    # graphs a given array and saves it based on input name
    
    array_x = array[:,0]
    array_y = array[:,1]
    pylab.plot(array_x,array_y, 'o', label = name)
    print array_y
    pylab.ylabel('Count')
    pylab.xlabel('Array Size')
    pylab.legend()
    pylab.savefig(name + '.png')
    
    
def main():
    if len(sys.argv)>2:
        print 'provide method input; bubble = 1, quickSort = 2'
        sys.exit()
    if len(sys.argv)<2:
        print 'provide method input; bubble = 1, quickSort = 2'
        sys.exit()
        
    method=int (sys.argv[1])
    
    complex_array= np.zeros((10000,2)) #allocate space for the array of n, and time to sort array with n elements
    #print complex_array
    index = 0
    array_range = 10000
    assign_count_array =np.zeros((10000,2)) #allocate space for the array of n, and num of asignments to sort array with n elements
    cond_count_array = np.zeros((10000,2)) #allocate space for the array of n, and num of conditional comparisons to sort array with n elements
    
    if method == 1:
        print 'bubbleSort'
        
        for array_size in range(100,1000,100): #iterate over the range (100,1000) in increments of 100
            for arraynum in range(100): #randomly generate 100 arrays of the given array_size
                array=randArray.randArray(array_range,array_size)
                [sortedArray, runTime, assign_count, cond_count]=bubbleSort.bubbleSort(array) #bubble sort and get the sorted array, the execution time, assignment and conditionals count          
                #print assign_count, cond_count
                complex_array[index,:]= [array_size, runTime]
                assign_count_array[index,:] = [array_size, assign_count]           
                cond_count_array[index,:] = [array_size, cond_count]
                #print cond_count_array[index,:]
                index = index + 1
    
    elif method == 2:
        print 'quickSort'
        for array_size in range(100,1000,100): #iterate over the range (100,1000) in increments of 100
            for arraynum in range(100): #randomly generate 100 arrays of the given array_size
                array=randArray.randArray(array_range,array_size)
                [sortedArray, runTime, assign_count, cond_count]=quickSort.quickSort(array)            
                #print assign_count, cond_count
                complex_array[index,:]= [array_size, runTime]
                assign_count_array[index,:] = [array_size, assign_count]
                cond_count_array[index,:] = [array_size, cond_count]
                print cond_count_array[index,:]
                index = index + 1
    complex_array = complex_array[0:index]
    assign_count_array = assign_count_array[0:index]
    cond_count_array = cond_count_array[0:index]
    #print complex_array
    #print cond_count_array
    xdata=complex_array[:,0]
    ydata=complex_array[:,1]
    if method == 1:
        curvefit_choice(complex_array, 'polynomial')
    if method == 2:
        curvefit_choice(complex_array, 'log')
    graph(assign_count_array, 'Assignments Count')
    graph(cond_count_array, 'Conditionals Count')
    return complex_array
if __name__ == '__main__':
    main()
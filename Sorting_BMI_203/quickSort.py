'''
Created on Jan 7, 2015

@author: Guanqing Ou

quick Sort: 

The steps are:

    Pick an element, called a pivot, from the array.
    Reorder the array so that all elements with values less than the pivot come before the pivot, while all elements with values greater than the pivot come after it (equal values can go either way). After this partitioning, the pivot is in its final position. This is called the partition operation.
    Recursively apply the above steps to the sub-array of elements with smaller values and separately to the sub-array of elements with greater values.

'''

import random
import randArray
from datetime import datetime

def quickSort(array):    
    assign = 0
    cond = 0
    start_time=datetime.now()
    [array, assign_count, cond_count]=quickSort_implement (array, 0, len(array)-1, assign, cond)    
    print assign_count, cond_count
    runTime = datetime.now()-start_time
    runTime = runTime.total_seconds()
    return array, runTime, assign_count, cond_count

def quickSort_implement(array,startIndex,endIndex, assign, cond):
    if startIndex < endIndex:
        pivot_randIndex = random.randint(startIndex, endIndex)
        print pivot_randIndex, array[pivot_randIndex]
        [pivot_sortedIndex, assign, cond] = partition(array, startIndex, endIndex, pivot_randIndex, assign, cond)  
        print pivot_sortedIndex, array[pivot_sortedIndex]
        [array1, assign, cond] = quickSort_implement(array, startIndex, pivot_sortedIndex-1, assign, cond)
        [array2, assign, cond] = quickSort_implement(array, pivot_sortedIndex+1, endIndex, assign, cond)
    
    return array, assign, cond
def partition(array, start, end, pivot_randIndex, assign, cond):
    left = start
    assign=assign+1
    right = end-1
    assign=assign+1
    pivot = array[pivot_randIndex]
    assign=assign+1
    while True:
        cond=cond+1
        while left <= right and array[left]< pivot:
            cond=cond+1
            left= left+1
            assign=assign+1
        while left <= right and array[right]>= pivot:
            cond=cond+1
            right = right-1
            assign=assign+1
        if left > right:
            cond=cond+1 
            break
        array[left], array[right] = array[right], array[left]
        assign=assign+1

    array[pivot_randIndex]=array[left]
    assign=assign+1
    array[end]=array[left]
    assign=assign+1
    array[left]=pivot       
    assign=assign+1 
    return left, assign, cond

def main():
    array=randArray.randArray(100,20)
    print array
    sorted_array = quickSort(array)
    print sorted_array
if __name__ == '__main__':
    main()
    
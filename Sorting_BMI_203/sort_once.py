'''
Created on Jan 19, 2015

@author: gogqou

runs code once to show functionality
'''
import sys
import randArray
import quickSort
import bubbleSort

def main():
    if len(sys.argv)>4:
        print 'too many inputs; provide method (bubble = 1 or quick = 2), size of random array, and range of random array values'
        sys.exit()
    elif len(sys.argv)<3:
        print 'too few inputs; provide method (bubble = 1 or quick = 2), size of random array, and range of random array values'
        sys.exit()
        
    method=int (sys.argv[1])
    print method
    size = int(sys.argv[2])
    range= int(sys.argv[3])
    
    array=randArray.randArray(range,size)
    print array
    if method ==1:
        [sortedArray, runTime, assign_count, cond_count]=bubbleSort.bubbleSort(array) #bubble sort and get the sorted array, the execution time, assignment and conditionals count
    
    elif method == 2:
        [sortedArray, runTime, assign_count, cond_count]=quickSort.quickSort(array)            
    else:
        print 'unrecognized method, either bubble = 1 or quick = 2'
        sys.exit()       
    print sortedArray
    return sortedArray
if __name__ == '__main__':
    main()
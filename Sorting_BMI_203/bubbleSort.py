'''
Created on Jan 7, 2015

@author: Guanqing Ou

Bubble Sort implementation: 
list of n numbers

check misordered adjacent pairs
swap any adjacent pair where second val < first

repeat procedure "n" times

'''

from datetime import datetime
def bubbleSort(array):
    assign_count = 0
    cond_count = 0
    start_time=datetime.now()
    for passnum in range(len(array)-1,0,-1):
        for i in range(passnum):
            if array[i]>array[i+1]:
                cond_count = cond_count +1
                temp2=array[i+1]
                assign_count = assign_count +1
                temp1=array[i]
                assign_count = assign_count +1
                array[i]=temp2
                assign_count = assign_count +1
                array[i+1]= temp1
                assign_count = assign_count +1
    runTime= datetime.now()- start_time
    runTime = runTime.total_seconds()
    return array, runTime, assign_count, cond_count   


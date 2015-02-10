'''
Created on Feb 8, 2015

@author: gogqou


start with a static alignment with a given matrix and optimize it against objective function :
sum TP for FP rates of 0, .1, .2, .3

'''
import sub_matrix_mk_dict as subMdict
import numpy as np
import sys
import SW_run as SW
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
import seq_list_SW as slSW
import pylab
import ROC_curves as ROC

def calc_new_scores(pos_align_array,neg_align_array, subMatrix):
    #multiple alignment arrays by substitution matrix to get the score for each aligned pair
    new_pos_scores=np.zeros(len(pos_align_array))
    new_neg_scores=np.zeros(len(neg_align_array))
    posgaps = np.zeros([len(pos_align_array), 2])
    neggaps = np.zeros([len(neg_align_array), 2])
    for i in range(0,len(pos_align_array)):
        new_pos_scores [i]= np.sum(pos_align_array[i]*subMatrix)
        array1 = pos_align_array[i].copy()
        posgaps[i] = [array1[0,23],array1[23,23]] 
    for j in range(0,len(neg_align_array)):
        new_neg_scores [j]= np.sum(neg_align_array[j]*subMatrix)
        array2 = neg_align_array[j].copy()
        neggaps[j] = [array2[0,23],array2[23,23]] 
        
    return new_pos_scores, new_neg_scores, posgaps, neggaps

def obj_function(pos_scores, neg_scores):
    obj_fnc_val=0
    for FPrate in range(0,40,10):
        TPrate = ROC.true_pos_rate(pos_scores, neg_scores, FPrate)
        obj_fnc_val = obj_fnc_val + TPrate
    return obj_fnc_val
def generate_newsubMatrix(origsubMatrix, gap_init, gap_ext):
    
    newsubMatrix = origsubMatrix.copy()
    newsubMatrix[:,23]=-gap_init
    newsubMatrix[23,:]=-gap_init
    newsubMatrix [23,23] = -gap_ext
    return newsubMatrix

def priority_pos_neg(pos_align_array, neg_align_array):
    sum_pos = np.zeros([24,24])
    sum_neg = np.zeros([24,24])
    #sum all the alignment arrays for each pairing together to get an average sense of which values in the substitution matrix are used 
    for i in range(0,len(pos_align_array)):
        sum_pos = sum_pos+pos_align_array [i]
    for j in range(0,len(neg_align_array)):
        sum_neg = sum_neg+neg_align_array [j]
    #find the values that are most differentially used between positive and negative pairings
    #takes absolute value so that both large negative and positive values are noted
    cumul =sum_pos-sum_neg
    abs_cumul = np.absolute(cumul)
    #print abs_cumul
    priority_list = []
    max=1
    #while there is a nonzero value left, keep finding the max, adding its index to the list, then setting it to zero to find a new max
    while max>0: 
        max_indices = np.unravel_index(abs_cumul.argmax(), abs_cumul.shape)
        priority_list.append(max_indices)
        max = abs_cumul[max_indices]
        abs_cumul[max_indices]=0
    delete_list = []
    for i in range(0,24):
        delete_list.append((23,i))
        delete_list.append((i,23))
    new_priority_list =[x for x in priority_list if x not in delete_list]
    return new_priority_list

def optimization(priority_list, origsubMatrix, pos_align, neg_align, delta, epsilon, max_itr):
    #delta = amount to vary values in subMatrix
    #epsillon = deviation from maximum value of objective function at which the optimization may stop
    #max_itr = maximum iterations to run if epsilon condition never reached NOT USED RIGHT NOW
    bestsubMatrix = origsubMatrix
    [new_pos_scores, new_neg_scores, posgaps, neggaps] = calc_new_scores(pos_align, neg_align, origsubMatrix)
    best_obj_func = obj_function(new_pos_scores, new_neg_scores)
    itr = 0
    for i in range(0,len(priority_list)):
        for newval in range(-50,50, 4):
            newMatrix= newsubMatrix(newval, priority_list[i], bestsubMatrix)
            
            [new_pos_scores, new_neg_scores,posgaps, neggaps] = calc_new_scores(pos_align, neg_align, newMatrix)
            old_best_obj_func = best_obj_func
            best_obj_func = max(best_obj_func, obj_function(new_pos_scores, new_neg_scores))
            eps = best_obj_func-old_best_obj_func
            if eps>0.0:
                'yes'
                bestsubMatrix = newMatrix.copy()
    
    return bestsubMatrix, best_obj_func



def newsubMatrix(newval, index, inputMatrix):
    newMatrix = inputMatrix.copy()
    newMatrix[index[0], index[1]] = newval
    newMatrix[index[1], index[0]] = newval
    return newMatrix
def graph_ROC(array, graphlabel):
    pylab.plot(array[:,1], array[:,0]/100, label = graphlabel)
    
    pylab.axis([0,1,0,1])
    pylab.legend(loc = 'lower right')
    pylab.ylabel('True Positive Rate')
    pylab.xlabel('False Positive Rate')
    pylab.savefig(graphlabel+'_ROC.png')
def main():
    np.set_printoptions(threshold=1000, linewidth=1000, precision = 5, suppress = False)
    if len(sys.argv)>5:
        print 'provide positive and negative pairs of sequences to align, directory, and substitution matrix choice '
        sys.exit()
    pos_seq_list_name = sys.argv[1]
    neg_seq_list_name = sys.argv[2]
    home = sys.argv[3]
    
    pos_seq_list_file = home+pos_seq_list_name+'.txt'
    neg_seq_list_file = home+neg_seq_list_name+'.txt'
    
    subMatrixFile = sys.argv[4]
    
    gap_init = 7
    gap_ext = 3
    
    [sub_Matrixdict, origSubMatrix, AAlist] = subMdict.mk_dict(home+subMatrixFile)
    '''
    [pos_scores, pos_align_array] = slSW.scores_from_seq_list(home, pos_seq_list_file, sub_Matrixdict, origSubMatrix, gap_init, gap_ext)
    
    np.save(home+'pos_align_array', pos_align_array)
    np.save(home+'pos_scores', pos_scores )
    
    [neg_scores, neg_align_array]= slSW.scores_from_seq_list(home, neg_seq_list_file, sub_Matrixdict, origSubMatrix, gap_init, gap_ext)
    
    np.save(home+'neg_align_array', neg_align_array )
    np.save(home+'neg_scores', neg_scores )
    '''
    origSubMatrix = generate_newsubMatrix(origSubMatrix, gap_init, gap_ext)
    
    
    #if already saved versions and made no changes to SW_run and seq_list_SW: just load from saved npy files
    neg_align_array=np.load(home+'neg_align_array.npy')
    neg_scores=np.load(home+'neg_scores.npy')
    pos_align_array=np.load(home+'pos_align_array.npy')
    pos_scores=np.load(home+'pos_scores.npy')
    
    obj_func = obj_function(pos_scores, neg_scores)
    print obj_func
    [new_pos_scores, new_neg_scores, posgaps, neggaps] = calc_new_scores(pos_align_array, neg_align_array, origSubMatrix)

    #generate priority of values in the scoring matrix to change
    priority_list = priority_pos_neg(pos_align_array, neg_align_array)
    
    #optimize using obj func calculation and output best performing matrix
    [bestMatrix, best_obj_func] = optimization(priority_list, origSubMatrix, pos_align_array, neg_align_array, 1, .3, 10)
    
    [opt_pos_scores, opt_neg_scores, posgaps, neggaps] = calc_new_scores(pos_align_array, neg_align_array, bestMatrix)
    print bestMatrix
    print best_obj_func
    '''
    np.save(home+'bestMatrix',bestMatrix)
    
    np.savetxt(home+'bestMatrix.txt', bestMatrix, fmt ='%+.2d', delimiter = '   ', header = ' A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V    B    Z    X    *')
    '''
    #try out new optimized matrix
    
    optimized_Matrix = np.load(home+'bestMatrix.npy')
    print optimized_Matrix
    [optimized_subMatrix_dict, newSubMatrix, AAlist] = subMdict.mk_dict_np(home+'bestMatrix.npy', AAlist)
    #print optimized_subMatrix_dict
    #print len(optimized_subMatrix_dict.keys())
    #print np.transpose(optimized_subMatrix_dict.keys())
    
    
    [new_opt_pos_scores, opt_pos_align_array] = slSW.scores_from_seq_list(home, pos_seq_list_file, optimized_subMatrix_dict, newSubMatrix, gap_init, gap_ext)
    [new_opt_neg_scores,opt_neg_align_array] = slSW.scores_from_seq_list(home, neg_seq_list_file, optimized_subMatrix_dict, newSubMatrix, gap_init, gap_ext)

    #[opt_pos_scores, opt_neg_scores, posgaps, neggaps] = calc_new_scores(pos_align_array, neg_align_array, optimized_Matrix)
    obj_func = obj_function(new_opt_pos_scores, new_opt_neg_scores)
    
    print obj_func
   
    orig_ROC_array = ROC.ROC_graph(new_pos_scores, new_neg_scores)
    pylab.plot(orig_ROC_array[:,1], orig_ROC_array[:,0]/100, label = 'Original Alignment')  
    optimized_ROC_array = ROC.ROC_graph(opt_pos_scores, opt_neg_scores)
    pylab.plot(optimized_ROC_array[:,1], optimized_ROC_array[:,0]/100, label = 'Optimized Matrix')
    realigned_optimized_ROC_array = ROC.ROC_graph(new_opt_pos_scores, new_opt_neg_scores)
    pylab.plot(realigned_optimized_ROC_array[:,1], realigned_optimized_ROC_array[:,0]/100, label = 'Realigned')
    pylab.axis([0,1,0,1])
    pylab.legend(loc = 'lower right')
    pylab.ylabel('True Positive Rate')
    pylab.xlabel('False Positive Rate')
    pylab.savefig('Compare_Optimized_ROC_Realigned.png')


    
    print 'done'
    
    return 1

if __name__ == '__main__':
    main()
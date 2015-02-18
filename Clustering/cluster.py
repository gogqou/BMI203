'''
Created on Feb 16, 2015

@author: gogqou

modified Scott Pegg's initial version


clustering algorithms:


hierarchical: 
bottom up

partitioning:
k medioids 
'''
#!C:/Python22/python.exe

###############################################################################
###############################################################################
##                                                                           ##
##  MIS-203                                                                  ##
##  Clustering Lecture   (April 12, 2005)                                    ##
##  Programming Assignment                                                   ##
##                                                                           ##
##  Author: Scott Pegg                                                       ##
##                                                                           ##
##                                                                           ##
##  The abbreviated instructions:                                            ##
##    You will be given a set of enzyme active sites in PDB format           ##
##    (1) Implement a similarity metric                                      ##
##    (2) Implement a clustering method based on a partitioning algorithm    ##
##    (3) Implement a clustering method based on a hierarchical algorithm    ##
##    (4) Answer the questions given in the homework assignment              ##
##                                                                           ##
##  Please read the full instructions from the course website _before_ you   ##
##  start this assignment!                                                   ##
##                                                                           ##
###############################################################################
###############################################################################



from string import *
from math import *
import sys, os
import glob 

import numpy as np


###############################################################################
#                                                                             #
# A simple class for an atom                                                  #

class Atom:

    def __init__(self, type):
        self.type = type
        self.coords = [0.0, 0.0, 0.0]

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.type

###############################################################################



###############################################################################
#                                                                             #
# A simple class for an amino acid residue                                    #

class Residue:

    def __init__(self, type, number):
        self.type = type
        self.number = number
        self.atoms = []
        self.sum_coords=[0.0, 0.0, 0.0]
        self.avg_coords=[0.0, 0.0, 0.0]
        self.max_atom=[]
        self.min_atom = []

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.type + " " + self.number

###############################################################################



###############################################################################
#                                                                             #
# A simple class for an active site                                           #

class ActiveSite:

    def __init__(self, name):
        self.name = name
        self.residues = []
        

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.name

###############################################################################



###############################################################################
#                                                                             #
# Read in all of the active sites from the given directory.                   #
#                                                                             #
# Input: directory                                                            #
# Output: list of ActiveSite instances                                        #

def read_active_sites(dir):
    files = glob.glob(dir + '/*.pdb')

    active_sites = []
    for file in files:
        if os.name == 'nt':
            name = splitfields(file, '\\')[-1][:-4]
        else:
            name = splitfields(file, '/')[-1][:-4]
        active_site = ActiveSite(name)
        
        l = open(file, "r").readlines()
        r_num = 0
        for i in range(len(l)):
            t = split(l[i])
            if t[0] != 'TER':
    
                # read in an atom
                atom_type = l[i][13:17]
                x_coord = float(l[i][30:38])
                y_coord = float(l[i][38:46])
                z_coord = float(l[i][46:54])
                atom = Atom(atom_type)
                atom.coords = [x_coord, y_coord, z_coord]
                  
                residue_type = l[i][17:20]
                residue_number = l[i][23:26]
            
                # make a new residue if needed
                if residue_number != r_num:
                    residue = Residue(residue_type, residue_number)
                    r_num = residue_number
    
                # add the atom to the residue
                residue.atoms.append(atom)
                residue.sum_coords = np.add(residue.sum_coords, atom.coords)
          
            else:  # I've reached a TER card
                active_site.residues.append(residue)       
            
            temp_atoms = residue.atoms 
            residue.avg_coords = residue.sum_coords/len(temp_atoms)
        active_sites.append(active_site)
    
    print "Read in %d active sites" % len(active_sites)
    
    return active_sites

###############################################################################
#
#calculate the major axes within an active site                               #
# basically figure out the largest distance between major carbons             #
#
#
#

def compute_major_axes(site):
    length = 0
    axes = 0
    
    return axes, length


###############################################################################
#                                                                             #
# Compute the similarity between two given ActiveSite instances.              #
#                                                                             #
# Input: two ActiveSite instances                                             #
# Output: the similarity between them (a floating point number)               #
#                                                                             #

def compute_similarity(site_A, site_B):

    similarity = 0.0
    
    #Fill in your code here!
    #sum of all the distances between primary carbons?

    return similarity

###############################################################################



###############################################################################
#                                                                             #
# Cluster a given set of ActiveSite instances using a partitioning method.    #
#                                                                             #
# Input: a list of ActiveSite instances                                       #
# Output: a clustering of ActiveSite instances                                #
#         (this is really a list of clusters, each of which is list of        #
#         ActiveSite instances)                                               #

def cluster_by_partitioning(active_sites):


  # Fill in your code here!

    active_sites = read_active_sites('/home/gogqou/Documents/Classes/bmi-203-hw3/active sites')
    print type(active_sites)
    residues= active_sites[2].residues
    print residues
    
    for i in range(len(residues)):
        print residues[i].avg_coords
        atoms = residues[i].atoms
        print atoms
    print active_sites
    return []

###############################################################################



###############################################################################
#                                                                             #
# Cluster the given set of ActiveSite instances using a hierarchical          #
# algorithm.                                                                  #
#                                                                             #
# Input: a list of ActiveSite instances                                       #
# Output: a list of clusterings                                               #
#         (each clustering is a list of lists of Sequence objects)            #

def cluster_hierarchically(active_sites):


  # Fill in your code here!


    return []

###############################################################################



###############################################################################
#                                                                             #
# Write the clustered ActiveSite instances out to a file.                     #
#                                                                             #
# Input: a filename and a clustering of ActiveSite instances                  #
# Output: none                                                                #

def write_clustering(filename, clusters):

    out = open(filename, 'w')

    for i in range(len(clusters)):
        out.write("\nCluster %d\n--------------\n" % i)
    for j in range(len(clusters[i])):
        out.write("%s\n" % clusters[i][j])

    out.close()

###############################################################################



###############################################################################
#                                                                             #
# Write a series of clusterings of ActiveSite instances out to a file.        #
#                                                                             #
# Input: a filename and a list of clusterings of ActiveSite instances         #
# Output: none                                                                #

def write_mult_clusterings(filename, clusterings):

    out = open(filename, 'w')

    for i in range(len(clusterings)):
        clusters = clusterings[i]

    for j in range(len(clusters)):
        out.write("\nCluster %d\n------------\n" % j)
    for k in range(len(clusters[j])):
        out.write("%s\n" % clusters[j][k])

    out.close()

###############################################################################






def main():
    # Some quick stuff to make sure the program is called correctly
    if len(sys.argv) < 4:
        print "Usage: cluster.py [-P| -H] <pdb directory> <output file>"
        sys.exit(0)
    
    ###############################################################################
#                                                                             #
# Top Level                                                                   #

    active_sites = read_active_sites(sys.argv[2])
    
    # Choose clustering algorithm
    if sys.argv[1][0:2] == '-P':
        print "Clustering using Partitioning method"
        clustering = cluster_by_partitioning(active_sites)
        #write_clustering(sys.argv[3], clustering)
    
    if sys.argv[1][0:2] == '-H':
        print "Clustering using hierarchical method"
        clusterings = cluster_hierarchically(active_sites)
        write_mult_clusterings(sys.argv[3], clusterings)


###############################################################################
    
    return 1
if __name__ == '__main__':
    main()
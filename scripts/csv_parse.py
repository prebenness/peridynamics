# -*- coding: utf-8 -*-
"""

For storing the output csv files

Created on Thu Oct 10 15:47:06 2019

@author: Ben Boys
"""

import os, csv, math
import numpy as np
from scipy import sparse
from scipy.constants import pi as PI
from geom_utils import inv_materials, materials, mat_cons
from phys_props import Geometry

class csvParse:
    def __init__(self, geom, verbose=False):
        self.v = verbose
        self.geom = geom
    def save_results(out_dir, time_tag, geom, load_tag, n_num_nodes):
        """ Input: root output directory
                    load_tag = TEST_1A_LOAD400000 for example
                    time_tag = "150"
                    
            Output: None
        """
        # Initiate output file path
        out_dir = out_dir + "/" + load_tag + "/" + str(n_num_nodes) + "-nodes/"
        
        
        # Attempt to create the directory
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        
        out_path = out_dir + "coors" + str(time_tag) + ".csv" # set out_path to file name
        
        # Attempt to write coordinates to output file
        with open(out_path, 'w', newline = '') as myfile:
            wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
            wr.writerow(['X','Y','Z'])
            
            for i in range(len(geom.coors)):
                X = geom.coors[i][0]
                Y = geom.coors[i][1]
                Z = geom.coors[i][2]
                mylist = [X, Y, Z]
                
                wr.writerow(mylist)
                
        
        # Store a "point" for each bond and associated strain and health value
        out_path = out_dir + "bonds" + str(time_tag) + ".csv" # set out_path to file name
        
        # only store bonds that are not broken?
        # for the geom.connect matrix
        
        
        with open(out_path, 'w', newline = '') as myfile:
            wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
            wr.writerow(['X','Y','Z', 'strain', 'health'])
            
            conn_lower_triangular = sparse.tril(geom.conn) # only iterate through the bonds once, i.e. lower triangular is all that is needed
            for k in range(len(np.nonzero(conn_lower_triangular)[0])): # only iterate over non zero values of the sparse matrix
                
                i = np.nonzero(conn_lower_triangular)[0][k]
                j = np.nonzero(conn_lower_triangular)[1][k]
                
                
# =============================================================================
#                 TODO: assertion
#                 assert geom.conn[i][j] == True,"geom.conn had no value for this index!"
# =============================================================================
                
                X = (geom.coors[i][0] + geom.coors[j][0])*0.5
                Y = (geom.coors[i][1] + geom.coors[j][1])*0.5
                Z = (geom.coors[i][2] + geom.coors[j][2])*0.5
                
                strain = geom.bond_stretches[i, j]
                health = 1.0 - (abs(strain) / geom.fail_stretches[i, j])
                     
                # Health cannot be less than 0
                if health <= 0.0:
                    health = 0.0
                    strain = geom.fail_stretches[i, j]
                         
                mylist = [X, Y, Z, strain, health]
                wr.writerow(mylist)
            
# =============================================================================
#             # NAIVE IMPLEMENTATION, ABOUT 2X SLOWER
#             for j in range(geom.conn_0.shape[1]):
#                 for i in range(j, geom.conn_0.shape[0]):
#                     if geom.conn[i, j] == True:
#                         # i is one node, j is other node, below is the average position of the bond
#                         X = (geom.coors[i][0] + geom.coors[j][0])*0.5
#                         Y = (geom.coors[i][1] + geom.coors[j][1])*0.5
#                         Z = (geom.coors[i][2] + geom.coors[j][2])*0.5
#                 
#                         strain = geom.bond_stretches[i, j]
#                         health = 1.0 - (abs(strain) / geom.fail_stretches[i, j])
#                     
#                         # Health cannot be less than 0
#                         if health <= 0.0:
#                             health = 0.0
#                             strain = geom.fail_stretches[i, j]
#                         
#                         mylist = [X, Y, Z, strain, health]
#                         wr.writerow(mylist)
# =============================================================================
                    
                
                
        
                
                
        
                
                

        
            
        
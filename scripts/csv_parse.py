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
        out_path = out_dir + "/" + load_tag + "/" + str(n_num_nodes) + "-nodes/"
        
        
        # Attempt to create the directory
        if not os.path.exists(out_path):
            os.makedirs(out_path)
        
        out_path = out_path + str(time_tag) + "coors.csv" # set out_path to file name
        
        # Attempt to write coordinates to output file
        with open(out_path, 'w', newline = '') as myfile:
            wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
            wr.writerow(['X','Y','Z'])
            
            print('LENGTH OF LIST IS {} columns, {} rows'.format(len(geom.coors[0]), len(geom.coors)))
            for i in range(len(geom.coors)):
                X = geom.coors[i][0]
                Y = geom.coors[i][1]
                Z = geom.coors[i][2]
                mylist = [X, Y, Z]
                
                wr.writerow(mylist)
        # Store a "point" for each bond and associated strain value
        
        #out_path = out_dir + "/" + test_tag + "/" + str(n_num_nodes) + "-nodes/" + "bonds.csv" + time_tag
                
                
                

        
            
        
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
    def __init__(self, geom, settings, verbose=False):
        self.v = verbose
        self.geom = geom
        self.settings = settings
# =============================================================================
#         self.func_order = [
#             'construct_H',
#             'calc_fails',
#             'calc_node_forces',
#             'calc_accs',
#             'calc_vels',
#             'calc_disps',
#             'calc_bond_stretches',
#             'total',
#             ]
#         self.t_rec = [timedelta(microseconds=0)]*len(self.func_order)
#         self.construct_H(init=True)
# ============================================================================= 
    def save_results(out_dir, time_tag, geom, load_tag, test_tag, n_num_nodes):
        """ Input: root output directory
                    test_tag = TEST_1A_LOAD400000 for example
                    time_tag = "150"
                    
            Output: None
        """
        # Initiate output file path
        out_path = out_dir + "/" + test_tag + "/" + str(n_num_nodes) + "-nodes/" + load_tag + "/" + "coors.csv" + time_tag
        
        # Attempt to write coordinates to output file
        with open('out_path', 'w') as myfile:
            wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
            wr.writerow(['X','Y','Z'])
            for i in range(len(geom.coors[0])):
                X = geom.coors[0][i]
                Y = geom.coors[1][i]
                Z = geom.coors[2][i]
                
                wr.writerow([X, Y, Z])
        # Store a "point" for each bond and associated strain value
        out_path = out_dir + "/" + test_tag + "/" + str(n_num_nodes) + "-nodes/" + load_tag + '/' + "bonds.csv" + time_tag
                
                
                

        
            
        
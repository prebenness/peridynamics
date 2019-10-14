#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 19:33:12 2019

@author: preben
"""
import os, csv, math
import numpy as np
from scipy import sparse
from scipy.constants import pi as PI
from geom_utils import inv_materials, materials, mat_cons
class Material:

    def __init__(self, flag):
        if flag == 'concrete':
            self.E_mod = 22e9
            self.v_ratio = 0.2
            self.G_mod = 8.8e9
            self.dens = 2400.0
            self.crit_ts = 0.000533	# ORIGINALLY 0.000533

            self.calc_eff_mod()
        elif flag == 'steel':
            self.E_mod = 210e9
            self.v_ratio = 0.3
            self.G_mod = 78e9
            self.dens = 8000.0
            self.crit_ts = 0.01     # ORIGINALLY 0.01

            self.calc_eff_mod()
        else:
            raise NotImplementedError('Does not recognise matflag {}'.format(flag))
        self.name = flag
        
    def calc_eff_mod(self):
        self.eff_mod = self.E_mod / ((1 - 2 * self.v_ratio)*(1 + self.v_ratio))
    
    def set_c_bond(self, delta):
        self.c_bond = ( 12.0 * self.E_mod ) / ( PI * math.pow(delta, 4.0) )       
        
class Geometry:
    def __init__(self, N, verbose=False):
        ## Scalars
        self.num_nodes = N
        self.v = verbose
        
        ## Sparse matrices for bonds
        self.conn = sparse.csr_matrix( (N, N), dtype=np.bool_ )
        self.conn_0 = sparse.csr_matrix( (N, N), dtype=np.bool_ )
        self.vol_coors = sparse.csr_matrix( (N, N), dtype=np.float64 )
        self.bond_types = sparse.csr_matrix( (N, N), dtype=np.int8 )
        self.stiff_facs = sparse.csr_matrix( (N, N), dtype=np.float64 )
        self.init_bond_lens = sparse.csr_matrix( (N, N), dtype=np.float64 )
        self.no_fails = sparse.csr_matrix( (N, N), dtype=np.bool_ )
        self.bond_stretches = sparse.csr_matrix( (N, N), dtype=np.float64 )
        self.fail_stretches = sparse.csr_matrix( (N, N), dtype=np.float64 )
         
        ## Dense matrices for nodes
        self.mat_flags = np.zeros( (N, 1), dtype=np.int8 )
        self.dens = np.zeros( (N, 1), dtype=np.float64 )
        self.coors = np.zeros( (N, 3), dtype=np.float64 )
        self.coors_0 = np.zeros(( N, 3), dtype=np.float64 )
        self.constraints = np.zeros( (N, 3), dtype=np.bool_ )
        self.body_forces = np.zeros( (N, 3), dtype=np.float64 )
        self.forces = np.zeros( (N, 3), dtype=np.float64 )
        self.accs = np.zeros( (N, 3), dtype=np.float64 )
        self.vels = np.zeros( (N, 3), dtype=np.float64 )
        self.disps = np.zeros( (N, 3), dtype=np.float64 )

    def load_geom(self, path):
        '''
        path should be path to coors.csv file
        '''
        N = self.num_nodes
        root_name = '_'.join( os.path.split(path)[-1].split('_')[1:] )
        in_dir = os.path.abspath( os.path.split(path)[0] )
        
        with open(os.path.join(in_dir, 'coors_' + root_name), 'r') as csv_file:
            reader = csv.reader(csv_file)
            i = 0
            for row in reader:
                if 'x' in row:
                    continue
                try:
                    self.coors_0[i, :] = np.array([ float(s) for s in row ])
                    i += 1
                except ValueError:
                    print( 'Skipping line in coors.csv:\n{}'.format(row) )
        print('Finished loading node coordinates')
        
        with open(os.path.join(in_dir, 'neighbours_' + root_name), 'r') as r:
            n = 0
            row = -1
            for l in r:
                try:
                    cols = np.array( [float(s) for s in l.split(',') ])
                    rows = np.full((len(cols)), row)
                    vals = np.full((len(cols)), 1)
                    self.conn_0 += sparse.csr_matrix( (vals,(rows, cols)), shape=(N, N) )
                    n += len(cols)
                except ValueError:
                    print('Skipping line in neighbours.csv:\n{}'.format(l))
                row += 1
                
                if row % 1000 == 0 and self.v:
                    print('Loaded {} out of {} lines'.format(row, N))
        print('Finished loading connectivity matrix')
        self.num_bonds = int(n/2)
        
        with open(os.path.join(in_dir, 'd_xyz_' + root_name), 'r') as r:
            for l in r:
                try:
                    self.d_xyz = [ float(s) for s in l.split(',')]
                except ValueError:
                    print('Skipping line in d_xyz:\n{}'.format(l))
        
        with open(os.path.join(in_dir, 'mat_flags_' + root_name), 'r') as r:
            i = 0
            for l in r:
                try:
                    self.mat_flags[i] = inv_materials[ l.strip() ]
                    mat = Material( flag=l.strip() )
                    self.dens[i] = mat.dens
                    i += 1
                except KeyError:
                    print('Skipping line in mat_flags:\n{}'.format(l))
        print('Finished loading material flags')
        
    def set_body_forces(self, settings, tag, scale=1.0):
        
        z_force = -1 * settings.MAX_REAC * scale / (self.num_nodes * settings.volume)
        base_vec = np.array([ 0, 0, z_force ])
        if tag == '1A':
            
            x_min_1_1A = 0.20625 - 0.040
            x_min_2_1A = x_min_1_1A + 0.4125
            x_min_3_1A = x_min_2_1A + 0.4125
            x_min_4_1A = x_min_3_1A + 0.4125
            x_min_5_1A = x_min_4_1A + 0.4125
            x_min_6_1A = x_min_5_1A + 0.4125
            x_min_7_1A = x_min_6_1A + 0.4125
            x_min_8_1A = x_min_7_1A + 0.4125
        
            x_max_1_1A = 0.20625 + 0.040
            x_max_2_1A = x_max_1_1A + 0.4125
            x_max_3_1A = x_max_2_1A + 0.4125
            x_max_4_1A = x_max_3_1A + 0.4125
            x_max_5_1A = x_max_4_1A + 0.4125
            x_max_6_1A = x_max_5_1A + 0.4125
            x_max_7_1A = x_max_6_1A + 0.4125
            x_max_8_1A = x_max_7_1A + 0.4125                   
    
            x_mins_1A = [
                    x_min_1_1A,
                    x_min_2_1A, 
                    x_min_3_1A, 
                    x_min_4_1A, 
                    x_min_5_1A, 
                    x_min_6_1A,
                    x_min_7_1A,
                    x_min_8_1A,
                    ]
            
            x_maxs_1A = [
                    x_max_1_1A,
                    x_max_2_1A,
                    x_max_3_1A,
                    x_max_4_1A,
                    x_max_5_1A,
                    x_max_6_1A,
                    x_max_7_1A,
                    x_max_8_1A,
                    ]
    
            y_min_1A = 0.125 - 0.080
            y_max_1A = 0.125 + 0.080
            
            z_min_1A = 0.600 - 2.0 * settings.delta
            z_max_1A = 0.600
            
            i = 0
            num_loaded = 0
            for node in self.coors_0:
                x_check = any([ ( (node[0]>mi) and (node[0]<ma) ) for mi, ma in zip(x_mins_1A, x_maxs_1A) ])
                y_check = ( (node[1]>y_min_1A) and (node[1]<y_max_1A) )
                z_check = ( (node[2]>z_min_1A) and (node[2]<z_max_1A) )
                
                if all([ x_check, y_check, z_check ]):
                    self.body_forces[i, :] = base_vec
                    num_loaded += 1
                i += 1
            
            self.num_loaded = num_loaded
            
        elif tag == '1B':
            
            x_min_1B = 1.550
            x_max_1B = 1.750
            y_min_1B = 0.000
            y_max_1B = 0.250
            z_min_1B = 0.600 - 2.0 * settings.delta
            z_max_1B = 0.600
            
            i = 0
            num_loaded = 0
            for node in self.coors_0:
                x_check = ( (node[0]>x_min_1B) and (node[0]<x_max_1B) )
                x_check = ( (node[1]>y_min_1B) and (node[1]<y_max_1B) )
                x_check = ( (node[2]>z_min_1B) and (node[2]<z_max_1B) )
                
                if all([ x_check, y_check, z_check ]):
                    self.body_forces[i, :] = base_vec
                    num_loaded += 1
                i += 1
                
            self.num_loaded = num_loaded
            
        elif tag == 'debug':
            self.body_forces = np.tile(base_vec, (self.num_nodes, 1))
        
        else:
            raise NotImplementedError('Did not recognise model tag {}'.format(tag))
            
        self.body_forces *= (self.num_nodes)/(self.num_loaded)
        assert( np.sum(self.body_forces != 0.0) )
        print('Finished loading {} out of {} nodes'.format(self.num_loaded, self.num_nodes))
        
    def set_constraints(self, settings, const_tag):
        if const_tag == 'cantilever':
            i = 0
            for node in self.coors_0:
                if node[0] <= settings.CLAMP_DIST:
                    self.constraints[i, :] = np.array([0.0, 0.0, 0.0])
        else:
            raise NotImplementedError('Did not recognise constraint type {}'.format(const_tag))
        print('Finished setting node constraints')
    
    def calc_vol_corrs(self, settings):
        n = 0
        cx = self.conn_0.tocoo()
        tmp_vol_coors = self.vol_coors.tolil()
        for i, j in zip(cx.row, cx.col):
            if i <= j:
                n += 1
                continue
            
            len_ij = np.linalg.norm(self.coors_0[i, :] - self.coors_0[j, :])
            
            if len_ij <= settings.delta - settings.radij:
                tmp_fac = 1.0
            elif (len_ij > settings.delta - settings.radij) and (len_ij <= settings.delta):
                tmp_fac = (settings.delta - settings.radij - len_ij) / (2.0 * settings.radij)
            else:
                tmp_fac = 0.0
            n += 1
            tmp_vol_coors[i, j] = tmp_fac
            
            if n % 1000 == 0 and self.v:
                print('Finished calculating {} of {} volume corrections ({:4.2f}%)'.format(n, self.num_bonds, 100*n/self.num_bonds))

        print('Finished calculating all volume corrections')
        self.vol_coors = tmp_vol_coors.tocsr()
    
    def set_bond_types_stiff_facs(self, settings, bulk_mat, rebar_mat):
        cx = self.conn_0.tocoo()
        n = 0
        tmp_stiff_facs = self.stiff_facs.tolil()
        tmp_bond_types = self.bond_types.tolil()
        tmp_fail_stretches = self.fail_stretches.tolil()
        
        for i, j in zip(cx.row, cx.col):
            if i <= j:
                continue
            
            if (materials[self.mat_flags[i,0]] == bulk_mat.name) and (materials[self.mat_flags[j,0]] == bulk_mat.name):
                c_tmp = bulk_mat.c_bond
                tmp_bond_types[i, j] = mat_cons['btb']
                tmp_fail_stretches[i, j] = bulk_mat.crit_ts
            elif (materials[self.mat_flags[i,0]] == bulk_mat.name) and (materials[self.mat_flags[j,0]] == rebar_mat.name):
                c_tmp = bulk_mat.c_bond
                tmp_bond_types[i, j] = mat_cons['btr']
                tmp_fail_stretches[i, j] = bulk_mat.crit_ts * 3.0
            elif (materials[self.mat_flags[i,0]] == rebar_mat.name) and (materials[self.mat_flags[j,0]] == bulk_mat.name):
                c_tmp = bulk_mat.c_bond
                tmp_bond_types[i, j] = mat_cons['rtb']
                tmp_fail_stretches[i, j] = bulk_mat.crit_ts * 3.0
            elif (materials[self.mat_flags[i,0]] == rebar_mat.name) and (materials[self.mat_flags[j,0]] == rebar_mat.name):
                c_tmp = rebar_mat.c_bond
                tmp_bond_types[i, j] = mat_cons['rtr']
                tmp_fail_stretches[i, j] = bulk_mat.crit_ts
            else:
                raise RuntimeError('Invalid material flags')
                
            neigh_vol_node = self.conn_0.getrow(i).nnz
            neigh_vol_neigh = self.conn_0.getrow(i).nnz
            stiff_fac = ( 2.0 * settings.neigh_vol ) * c_tmp / ( neigh_vol_node + neigh_vol_neigh )
            tmp_stiff_facs[i, j] = stiff_fac
            n += 1
        
            if n % 1000 == 0 and self.v:
                print('Finished calculating {} out of {} bond types stiff facs ({:4.2f}%)'.format(n, self.num_bonds, 100*n/self.num_bonds))

        print('Finished calculating all bond type stiff facs')
        self.stiff_facs = tmp_stiff_facs.tocsr()
        self.bond_types = tmp_bond_types.tocsr()
        self.fail_stretches = tmp_fail_stretches.tocsr()
        
    def init_sim(self):
        self.conn = self.conn_0.copy()
        self.coors = self.coors_0.copy()

class SimulationSettings:
    def __init__(self, geom, bulk_mat, rebar_mat, saf_fac=1.5, nt=2e4):
        self.MAX_REAC = 4e5
        self.CLAMP_DIST = 0.05
        
        self.build_up = 500
        self.damping = 2.5e6
        self.delta = np.mean(geom.d_xyz) * PI
        self.neigh_vol = 4.0/3.0 * PI * math.pow(self.delta, 3.0)
        self.volume = np.prod(geom.d_xyz)
        self.radij = (1.0 / 2.0) * np.sum(geom.d_xyz) / 3.0
        
        bulk_mat.set_c_bond(self.delta)
        rebar_mat.set_c_bond(self.delta)
        
        self.saf_fac = saf_fac
        self.nt = nt
        self.calc_dt(bulk_mat)
        
    def calc_dt(self, bulk_mat):
        dx = self.delta / PI
        self.dt = ( 0.8 * math.pow( 2.0 * bulk_mat.dens * dx / ( PI * pow(self.delta, 2.0) * dx * bulk_mat.c_bond), 0.5)) / self.saf_fac


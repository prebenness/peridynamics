#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 17:37:55 2019

@author: preben
"""

import numpy as np
from datetime import datetime, timedelta
from scipy import sparse
from csv_parse import csvParse

class Simulation:
    
    def __init__(self, geom, settings, verbose=False):
        self.v = verbose
        self.geom = geom
        self.settings = settings
        self.func_order = [
            'construct_H',
            'calc_fails',
            'calc_node_forces',
            'calc_accs',
            'calc_vels',
            'calc_disps',
            'calc_bond_stretches',
            'total',
            ]
        self.t_rec = [timedelta(microseconds=0)]*len(self.func_order)
        self.construct_H(init=True)
        
    def prune_mats(self):
        self.geom.conn.prune()
        self.geom.vol_coors.prune()
        self.geom.bond_types.prune()
        self.geom.stiff_facs.prune()
        self.geom.no_fails.prune()
        self.geom.bond_stretches.prune()
        self.geom.fail_stretches.prune()
        
    def run_time_integration(self, bulk_mat, rebar_mat, args, load_tag, n_num_nodes, verbose=False):
        old = False
        tol = 1e-5
        prev_perc_damage = 0.0
        
        save_every = 50
        max_its = self.settings.nt
        
        # Stretches -> fails -> forces -> accs -> vels -> disps -> stretches
        it = 0
        t = list(range(len(self.func_order)))
        while True:
            gst = datetime.now()
            if old:
                version = 'NAIVE FOR LOOP IMPLEMENTATION'
                t[0] = self.construct_H_old(init=False)
                t[1] = self.calc_fails_old()
                t[2] = self.calc_node_forces_old()
            else:
                version = 'EFFICIENT MATRIX IMPLEMENTATION'
                t[0] = self.construct_H(init=False)
                t[1] = self.calc_fails()
                t[2] = self.calc_node_forces()
                
            t[3] = self.calc_accs()
            t[4] = self.calc_vels()
            t[5] = self.calc_disps()
            self.geom.coors += self.geom.disps
            assert( not np.array_equal(self.geom.coors, self.geom.coors_0) )
            
            if old:
                t[6] = self.calc_bond_stretches_old()
            else:
                t[6] = self.calc_bond_stretches()
                
            it += 1
            rep = '\nCompleted iteration {}    {}\n{} out of {} bonds remaining ({:4.2f}%)'.format(it, version, int(self.geom.conn.nnz/2),\
                                       self.geom.num_bonds, 100*self.geom.conn.nnz/self.geom.num_bonds)
            t[7] = datetime.now() - gst
            print(rep)
#            self.t_rec = [ tt/(it+1) + tc for tt, tc in zip(t, self.t_rec) ]
            self.t_rec = t
            self.print_times()
            self.prune_mats()
            
            if it % 50 ==0:
                csvParse.save_results(args.out_dir, it, self.geom, load_tag, n_num_nodes)
                print('SAVING NOW')
            
            
            
    def calc_bond_stretches(self):
        st = datetime.now()
        shape = (self.geom.num_nodes, self.geom.num_nodes)
        H_x, H_y, H_z = self.H_x, self.H_y, self.H_z
        H_x0, H_y0, H_z0 = self.H_x0, self.H_y0, self.H_z0
        
        diff_x = H_x - H_x0
        diff_y = H_y - H_y0
        diff_z = H_z - H_z0
        
        tmp_x = sparse.csr_matrix( diff_x[H_x0.nonzero()] / H_x0[H_x0.nonzero()] )
        tmp_x.resize(shape)
        
        tmp_y = sparse.csr_matrix( diff_y[H_y0.nonzero()] / H_y0[H_y0.nonzero()] )
        tmp_y.resize(shape)
        
        tmp_z = sparse.csr_matrix( diff_z[H_z0.nonzero()] / H_z0[H_z0.nonzero()] )
        tmp_z.resize(shape)
        
        tmp = tmp_x.power(2) + tmp_y.power(2) + tmp_z.power(2)
        self.geom.bond_stretches = tmp.sqrt()
        return datetime.now() - st
    
    def construct_H(self, init):
        st = datetime.now()
        shape = (self.geom.num_nodes, self.geom.num_nodes)
        
        if init:
            cs = sparse.csr_matrix(self.geom.coors)
#            coors = self.geom.coors
        else:
            cs = sparse.csr_matrix(self.geom.coors_0)
#            coors = self.geom.coors_0
            
#        print('loaded conn {}'.format(datetime.now()-st))
        
#        cxo = self.geom.conn.tocoo()
#        N = self.geom.num_bonds
        cols, rows, data_x, data_y, data_z = np.array((1,1)), np.array((1,1)), np.array((1,1)), np.array((1,1)), np.array((1,1))
        for i in range(shape[0]):
            row = self.geom.conn.getrow(i)
            np.append(rows, row.indices)
            
#            np.append(data_x, coors[:, 0])
#            np.append(data_y, coors[:, 1])
#            np.append(data_z, coors[:, 2])
            
            np.append(data_x, cs.getrow(0).data)
            np.append(data_y, cs.getrow(1).data)
            np.append(data_z, cs.getrow(2).data)

            np.append(cols, np.full([1, row.nnz], i))
            
#            data_x.extend( list(cs.getrow(0).data) )
#            data_y.extend( list(cs.getrow(1).data) )
#            data_z.extend( list(cs.getrow(2).data) )
#            cols.append( [i] * row.nnz )
#            n += row.nnz
            
        lam_x = sparse.csr_matrix( (data_x, (rows, cols)), shape=shape)
        lam_y = sparse.csr_matrix( (data_y, (rows, cols)), shape=shape)
        lam_z = sparse.csr_matrix( (data_z, (rows, cols)), shape=shape)

#        lam_x = sparse.lil_matrix(shape)
#        lam_y = sparse.lil_matrix(shape)
#        lam_z = sparse.lil_matrix(shape)
#        x_col, y_col, z_col = cs.getcol(0), cs.getcol(1), cs.getcol(2)
#        for i in range(shape[0]):
#            con_col = self.geom.conn.getrow(i)
#            lam_x[:, i] = x_col.multiply(con_col.transpose())
#            lam_y[:, i] = y_col.multiply(con_col.transpose())
#            lam_z[:, i] = z_col.multiply(con_col.transpose())
#            i += 1

#        print('computed lam {}'.format(datetime.now()-st))

#        lam_x, lam_y, lam_z = lam_x.tocsr(), lam_y.tocsr(), lam_z.tocsr()
            
        H_x = lam_x - lam_x.transpose()
        H_y = lam_y - lam_y.transpose()
        H_z = lam_z - lam_z.transpose()
        
#        print('computed H {}'.format(datetime.now()-st))
        
        if init:
            self.H_x0 = H_x
            self.H_y0 = H_y
            self.H_z0 = H_z
        else:
            self.H_x = H_x
            self.H_y = H_y
            self.H_z = H_z
            norms = H_x.power(2) + H_y.power(2) + H_z.power(2)
            self.norms = norms.sqrt()
        
        print('Stored H {}'.format(datetime.now()-st))
        return st - datetime.now()
        
    def construct_H_old(self, init):
        st = datetime.now()
        
        if init:
            cs = self.geom.coors
        else:
            cs = self.geom.coors_0
            
        V_x = cs[:, 0]
        lam_x = sparse.csr_matrix(np.tile(V_x, (self.geom.num_nodes, 1)))
        del V_x
        
        V_y = cs[:, 1]
        lam_y = sparse.csr_matrix(np.tile(V_y, (self.geom.num_nodes, 1)))
        del V_y
        
        V_z = cs[:, 2]
        lam_z = sparse.csr_matrix(np.tile(V_z, (self.geom.num_nodes, 1)))
        del V_z
        
        H_x = sparse.csr_matrix(lam_x.multiply(self.geom.conn) - lam_x.transpose().multiply(self.geom.conn))
        H_y = sparse.csr_matrix(lam_y.multiply(self.geom.conn) - lam_y.transpose().multiply(self.geom.conn))
        H_z = sparse.csr_matrix(lam_z.multiply(self.geom.conn) - lam_z.transpose().multiply(self.geom.conn))
        
        if init:
            self.H_x0 = H_x
            self.H_y0 = H_y
            self.H_z0 = H_z
        else:
            self.H_x = H_x
            self.H_y = H_y
            self.H_z = H_z
            norms = H_x.power(2) + H_y.power(2) + H_z.power(2)
            self.norms = norms.sqrt()
        
        return st - datetime.now()
    
    def calc_bond_stretches_old(self):
        st = datetime.now()
        tmp_bond_stretches = self.geom.bond_stretches.tolil()
        cx = self.geom.conn_0.tocoo()
        
        for i, j in zip(cx.row, cx.col):
            if i <= j:
                continue

            init = np.linalg.norm(self.geom.coors_0[i, :] - self.geom.coors_0[j, :])
            cur = np.linalg.norm(self.geom.coors[i, :] - self.geom.coors[j, :])
            s = ( cur - init ) / init

            if init == 0:
                raise RuntimeError('DIV ZERO ERROR:\ni: {} j: {}\nCoor i: {} Coor j: {}'.format(i, j, self.geom.coors[i, :], self.geom.coors[j, :]))
            
            tmp_bond_stretches[i, j] = s
            
        if self.v:
            print('Finished calculating bond stretches')
        
        self.geom.bond_stretches = tmp_bond_stretches.tocsr()
        return datetime.now() - st
    
    def calc_fails(self):
        st = datetime.now()
        
        bond_healths = self.geom.fail_stretches - self.geom.bond_stretches.sign().multiply(self.geom.bond_stretches)
        bond_healths = bond_healths + self.geom.bond_stretches
        bond_healths = bond_healths.multiply(0.5)
        
        self.geom.conn = bond_healths.astype(bool)
        self.geom.conn.eliminate_zeros()
        self.geom.conn.prune()
        
        return datetime.now() - st
    
    def calc_fails_old(self):
        st = datetime.now()
        tmp_conn = self.geom.conn.tolil()
        
        bond_healths = self.geom.fail_stretches - self.geom.bond_stretches.sign().multiply(self.geom.bond_stretches)

        cx = bond_healths.tocoo()
        for i, j, val in zip(cx.row, cx.col, cx.data):
            assert(val)
            if val <= 0:
                tmp_conn[i, j] = False
                
        self.geom.conn = tmp_conn.tocsr()
        
        return datetime.now() - st
    
    def calc_node_forces(self):
        st = datetime.now()
        shape = (self.geom.num_nodes, self.geom.num_nodes)
        
        forces = self.geom.body_forces
        
        force_mags = self.geom.bond_types * self.geom.stiff_facs * self.geom.bond_stretches * self.settings.volume *\
                    self.geom.vol_coors

        force_mags = self.geom.bond_types.multiply(self.geom.stiff_facs)
        force_mags = force_mags.multiply(self.geom.bond_stretches)
        force_mags = force_mags.multiply(self.geom.vol_coors)
        force_mags = force_mags.multiply(self.settings.volume)

        H_x, H_y, H_z = self.H_x, self.H_y, self.H_z
        
        norm_force_mags = sparse.csr_matrix( force_mags[self.norms.nonzero()] / self.norms[self.norms.nonzero()] )
        norm_force_mags.resize(shape)
        
        cxx = sparse.coo_matrix( H_x.multiply(force_mags) )
        cxy = sparse.coo_matrix( H_y.multiply(force_mags) )
        cxz = sparse.coo_matrix( H_z.multiply(force_mags) )
        
        forces[cxx.row, 0] += cxx.data
        forces[cxy.row, 1] += cxy.data
        forces[cxz.row, 2] += cxz.data

        if self.v:
            print('Finished calculating node forces')
        
        self.geom.forces = forces + self.geom.vels * self.settings.damping
        return datetime.now() - st
    
    def calc_node_forces_old(self):
        st = datetime.now()
        forces = self.geom.body_forces
        
        cx = self.geom.conn.tocoo()
        for i, j, val in zip(cx.row, cx.col, cx.data):
            assert(val)
            if i <= j:
                continue
            
            bf_mult = self.geom.bond_types[i, j]
                
            tmp_force = bf_mult * self.geom.stiff_facs[i, j] * self.geom.bond_stretches[i, j] * self.settings.volume * self.geom.vol_coors[i, j]
            vec = self.geom.coors[j, :] - self.geom.coors[i, :]
            s = np.linalg.norm(vec)
            if s == 0:
                raise RuntimeError('DIV ZERO ERROR\n i: {} j: {}\nCoor i: {} Coor j: {}'.format(i, j, self.geom.coors[i, :], self.geom.coors[j, :]))
                s = 1.0
            vec = vec / s
            
            force_vec = tmp_force * vec
            
            forces[i, :] += force_vec
            forces[j, :] -= force_vec

        if self.v:        
            print('Finished calculating node forces')
        
        self.geom.forces = forces + self.geom.vels * self.settings.damping
        return datetime.now() - st
    
    def calc_accs(self):
        st = datetime.now()
        self.geom.accs = self.geom.forces / np.tile(self.geom.dens, (1, 3))
        return datetime.now() - st
    
    def calc_vels(self):
        st = datetime.now()
        self.geom.vels = self.geom.vels + self.geom.accs * self.settings.dt
        return datetime.now() - st
    
    def calc_disps(self):
        st = datetime.now()
        self.geom.disps = self.geom.disps + self.geom.vels * self.settings.dt
        return datetime.now() - st
    
    def print_times(self):
        tot = sum(tt.microseconds for tt in self.t_rec)
        tot = self.t_rec[-1].microseconds
        rep = 'Timing report\n'
        for f, tt in zip(self.func_order[:-1], self.t_rec[:-1]):
            rep += 'Func: {:20}    Time[ms]: {:8d}    Of total: {:4.2f}%\n'.format(f, tt.microseconds, 100*tt.microseconds/tot)
#        rep += 'Total: {:>40}'.format(tot)
        rep += '\n'
        
        print(rep)
        
        
        
        
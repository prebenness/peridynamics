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
# =============================================================================
#         TODO:
#         DONE Use stretches to update fails
#         DONE Use stretch to calc forces on all nodes
#         ??? if in build-up phase, scale forces accordingly:
#         DONE from forces get accs, vels, disps for all nodes
#         DONE? use disps and old positions to get new positions
#         save results and check for convergence every N steps
#             check if has converged or spent iteration budget
#             average absolute displacement over loaded nodes
#         incriment iteration counter
# =============================================================================
        
        # Define local variables
        
        it = 0 # No. of iterations
        save_every = 50 # Save and check for convergence every 50 timesteps
        
        old = False # Old or new implementation
        
        tol = 1e-5 # Tolerance
        perc_damage = 0.0 # percent damage
        prev_perc_damage = 0.0 # previous iteration percent damage
        max_its = self.settings.nt # max no. of iterations before max iteration budget is reached
        stable_steps = 0 # Total number of stable steps so far
        falling = 0 #time steps since ???
        converged = False #bool has the solver converged?
    
        # Stretches -> fails -> forces -> accs -> vels -> disps -> stretches
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
            rep = '\nCompleted iteration {}    {}'.format(it, version)
            t[7] = datetime.now() - gst

            
            print(rep)
#            self.t_rec = [ tt/(it+1) + tc for tt, tc in zip(t, self.t_rec) ]
            self.t_rec = t
            self.print_times()
            self.prune_mats()
            
            # Save results and check for convergence every save_every steps
            if it % save_every == 0:
                
                # Define local variables
                bonds_remaining = int(self.geom.conn.nnz/2) 
                total_bonds = np.floor(self.geom.num_bonds/2)
                perc_damage = 1.0 - bonds_remaining/total_bonds
                change_damage = perc_damage - prev_perc_damage
                prev_perc_damage =  perc_damage
                # Calculate the average node displacement
                av_node_zdisplacement = self.calc_loaded_disp()
                
                
                # Save results
                print('SAVING RESULTS \n{} out of {} bonds remaining ({:4.2f}%)'.format(bonds_remaining, total_bonds, 100.0*(1.0-perc_damage)))
                csvParse.save_results(args.out_dir, it, self.geom, load_tag, n_num_nodes)
                
                
                # Check if has converged or spent iteration budget
                print('AVERAGE NODE DISPLACEMENT:{}'.format(av_node_zdisplacement)) # Average displacement in the z direction
                
                if abs(change_damage) < tol:
                    stable_steps += save_every
                else:
                    stable_steps = 0
                
                if stable_steps >= 2* self.settings.build_up:
                    if av_node_zdisplacement * np.sum(self.geom.body_forces, axis=0)[2] < 0: # average node displacement in z direction * sum of body forces in z direction
                        falling = 0
                        if self.calc_abs_disp() <= tol:
                            converged = True
                            print('SIMULATION CONVERGED IN {} TIMESTEPS!'.format(it))
                            print('AVERAGE NODE DISPLACEMENT: {}'.format(av_node_zdisplacement))
                            print('BOND DAMAGE: {4.2f}%'.format(perc_damage*100.0))
                            break
                    else:
                        falling += save_every
                    if falling >= 2* self.settings.build_up:
                        converged = False
                        print('BEAM BROKE, FALLING')
                        print('BOND DAMAGE: {4.2f}%'.format(perc_damage*100.0))
                        # break
                            
                
                # If over budget and not converged beam has failed
                if it >= max_its:
                    converged = False
                    print('SIMULARTION DID NOT CONVERGE, REACHED MAX ITERATION BUDGET OF {} ITERATIONS'.format(max_its))
                    print('BOND DAMAGE: {4.2f}%'.format(perc_damage*100.0))
                    break
                
        return converged
                
            
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
    
    def calc_loaded_disp(self):
        #TODO check this function behaves correctly and is what Preben wanted
        
        # Disp of all nodes since start
        tot_disps = self.geom.coors - self.geom.coors_0
        x_disp, y_disp, z_disp = zip(*tot_disps)
        
        # Calculate average displacement of the nodes in each direction
        #x_av_disp = np.sum(x_disp)/self.geom.num_loaded
        #y_av_disp = np.sum(y_disp)/self.geom.num_loaded
        z_av_disp = np.sum(z_disp)/self.geom.num_loaded
        
        return z_av_disp # Only returns Z disp at the moment
    
    def calc_abs_disp(self):
        #TODO check this function behaves correctly and is what Preben wanted
        
        # Disp of all nodes since start
        tot_disps = self.geom.coors - self.geom.coors_0
        x_disp, y_disp, z_disp= zip(*tot_disps)
        
        # Pythag abs distance of each of the N nodes, gives (N * 1) vector of abs disp
        abs_disp = pow(x_disp, 2) + pow(y_disp, 2) + pow(z_disp, 2)
        abs_disp = pow(abs_disp, 0.5)
        
        # Return scalar value for the average absolute displacement
        av_abs_disp = np.sum(abs_disp)/self.geom.num_loaded
        
        return av_abs_disp
    
    def print_times(self):
        tot = sum(tt.microseconds for tt in self.t_rec)
        tot = self.t_rec[-1].microseconds
        rep = 'Timing report\n'
        for f, tt in zip(self.func_order[:-1], self.t_rec[:-1]):
            rep += 'Func: {:20}    Time[ms]: {:8d}    Of total: {:4.2f}%\n'.format(f, tt.microseconds, 100*tt.microseconds/tot)
#        rep += 'Total: {:>40}'.format(tot)
        rep += '\n'
        
        print(rep)
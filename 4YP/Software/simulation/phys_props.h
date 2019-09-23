#ifndef PHYS_PROPS_H
#define PHYS_PROPS_H

#include "classes.h"			// defines geometry class
//#include <string>

// FUNCTIONS
// Materials
void init_mat_flags(Geometry &geom, Material mat);							// set all as steel initially
void add_steel_bar(Geometry &geom, Material Steel, std::string test_tag);	// add steel bars
void init_mat_props(Material &mat, int flag);								// initialise material properties

// Simulation settings
void init_simsets(Sim_settings &config, Geometry geom,\
				  Material &Conc, Material &Steel);							// initialise and store Sim_settings object

// Forces
void set_body_forces(Geometry &geom, Sim_settings set, float scale, std::string test_tag);

// Constraints
void set_constraints(Geometry &geom, int constr_type);

// Corrections and other fudge factors
void calc_vol_corrs(Geometry &geom, Sim_settings set);
void set_bond_types_stiff_facs(Geometry &geom, Sim_settings set, \
								Material Conc, Material Steel);

#endif
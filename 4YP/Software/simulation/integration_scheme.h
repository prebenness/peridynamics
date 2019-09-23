#ifndef INTEGRATION_SCHEME
#define INTEGRATION_SCHEME

#include "classes.h"
#include <Eigen/Eigen>
#include <string>

using namespace Eigen;

// MASTER TIME STEPPING LOOP
bool run_time_integration(Geometry &geom, Sim_settings set, Material Conc, Material Steel, std::string out, std::string test_tag);

// COORDINATES -> STRETCH
void calc_bond_stretch(Geometry &geom);
void calc_bond_stretch_2(Geometry &geom);

// STRETCHES -> UPDATE FAILS
void calc_fails(Geometry &geom);

// STRETCHES -> FORCES
void calc_node_forces(Geometry &geom, Sim_settings set);

// FORCES AND VELS -> ACCS -> VELS -> DISPS
void calc_accs(Geometry &geom, Sim_settings set);
void calc_vels(Geometry &geom, Sim_settings set);
void calc_disps(Geometry &geom, Sim_settings set);

// GET DISP OF LOADED NODES AT END
float calc_loaded_disp(Geometry geom);



#endif
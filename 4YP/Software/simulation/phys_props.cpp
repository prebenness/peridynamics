// These functions set material flags, calculate various physical properties

#include <cmath>
#include <string>
#include <iostream>
#include <Eigen/Eigen>
#include <omp.h>
#include <cstdlib>

#include "phys_props.h"	// main header
#include "classes.h"	// class definitions, enumerations
#include "main.h"		// for N_NUM_NODES and other constants

using namespace Eigen;

// SET ALL FLAGS AS CONCRETE INITIALLY
void init_mat_flags(Geometry &geom, Material mat) {

	MatrixXi mat_flags;	mat_flags.resize(1, N_NUM_NODES);
	mat_flags.setOnes(1, N_NUM_NODES);							// use col vecs, here only one entry per col
	MatrixXf dens; dens.resize(1, N_NUM_NODES);
	dens.setOnes(1, N_NUM_NODES);								// density values for all nodes

	// use this rather than += MatrixXi::Constant, (as N_NUM_NODES no longer 32 bit int)
	mat_flags *= concrete;
	dens *= mat.dens;
	//mat_flags += MatrixXi::Constant(1, N_NUM_NODES, concrete);	// init all as concrete
	//dens += MatrixXf::Constant(1, N_NUM_NODES, mat.dens);		// concrete density
	geom.mat_flags = mat_flags; geom.dens = dens;				// store

	//std::cout << "init_mat_flags: dens.cols(): " << dens.cols() << std::endl;
}

// ADD STEEL BARS ALONG LONGITUDINAL DIRECTION
void add_steel_bar(Geometry &geom, Material Steel, std::string test_tag) {
	// As cannot switch on string, so use this hack
	int switch_case = -1;
	if (test_tag == "TEST_0") { switch_case = 0; }
	if (test_tag == "TEST_1A") { switch_case = 1; }
	if (test_tag == "TEST_1B") { switch_case = 2; }


	switch (switch_case) {
	case 0:
	{
		// Define center positions and radius of bars
		MatrixXf pos_0(2, 2);
		pos_0 << 0.025F, 0.075F, // [y_1, y_2,
				 0.025F, 0.025F; //  z_1, z_2]
		float radius_0 = 0.005F;							// meters, radius of bars
		MatrixXi tmp_mat_0 = geom.mat_flags;				// concrete initially
		MatrixXf tmp_dens_0 = geom.dens;

		// add steel bars
		double cost_1_0, cost_2_0;
		for (int i = 0; i < tmp_mat_0.cols(); i++) {
			// bars define two cylindrical regions, if point within this region, cost will be negative

			cost_1_0 = pow(geom.coors(1, i) - pos_0(0, 0), 2) + pow(geom.coors(2, i) - pos_0(1, 0), 2) - pow(radius_0, 2);
			cost_2_0 = pow(geom.coors(1, i) - pos_0(0, 1), 2) + pow(geom.coors(2, i) - pos_0(1, 1), 2) - pow(radius_0, 2);

			if (cost_1_0 <= 0 || cost_2_0 <= 0) {
				tmp_mat_0(0, i) = steel;
				tmp_dens_0(0, i) = Steel.dens;
			}
		}
		geom.mat_flags = tmp_mat_0; geom.dens = tmp_dens_0;			// store in Geometry object
		break;
	}
	case 1:
	{
		// Define center positions and radius of bars
		MatrixXf pos_1A(2, 4);
		pos_1A <<	0.03100F, 0.03825F, 0.21900F, 0.21175F,			// [y_1, y_2, y_3, y_4
					0.03100F, 0.56900F, 0.03100F, 0.56900F;			//  z_1, z_2, z_3, z_4]
		float radius_c_1A = 0.00600F; float radius_t_1A = 0.01325F;	// meters, radius of bars
		MatrixXi tmp_mat_1A = geom.mat_flags;						// concrete initially
		MatrixXf tmp_dens_1A = geom.dens;

		// add steel bars
		double cost_1_1A, cost_2_1A, cost_3_1A, cost_4_1A;
		for (int i = 0; i < tmp_mat_1A.cols(); i++) {
			// bars define two cylindrical regions, if point within this region, cost will be negative

			cost_1_1A = pow(geom.coors(1, i) - pos_1A(0, 0), 2) + pow(geom.coors(2, i) - pos_1A(1, 0), 2) - pow(radius_c_1A, 2);
			cost_2_1A = pow(geom.coors(1, i) - pos_1A(0, 1), 2) + pow(geom.coors(2, i) - pos_1A(1, 1), 2) - pow(radius_t_1A, 2);
			cost_3_1A = pow(geom.coors(1, i) - pos_1A(0, 2), 2) + pow(geom.coors(2, i) - pos_1A(1, 2), 2) - pow(radius_c_1A, 2);
			cost_4_1A = pow(geom.coors(1, i) - pos_1A(0, 3), 2) + pow(geom.coors(2, i) - pos_1A(1, 3), 2) - pow(radius_t_1A, 2);

			if (cost_1_1A <= 0 || cost_2_1A <= 0 || cost_3_1A <= 0 || cost_4_1A <= 0) {
				tmp_mat_1A(0, i) = steel;
				tmp_dens_1A(0, i) = Steel.dens;
			}
		}
		geom.mat_flags = tmp_mat_1A; geom.dens = tmp_dens_1A;	// store in Geometry object
		break;
	}
	case 2:
	{
		// Define center positions and radius of bars
		MatrixXf pos_1B(2, 4);
		pos_1B <<	0.03100F, 0.03825F, 0.21900F, 0.21175F,			// [y_1, y_2, y_3, y_4
					0.03100F, 0.56900F, 0.03100F, 0.56900F;			//  z_1, z_2, z_3, z_4]

		float radius_c_1B = 0.00600F; float radius_t_1B = 0.01325F;	// meters, radius of bars
		MatrixXi tmp_mat_1B = geom.mat_flags;						// concrete initially
		MatrixXf tmp_dens_1B = geom.dens;

		// add steel bars
		double cost_1_1B, cost_2_1B, cost_3_1B, cost_4_1B;
		for (int i = 0; i < tmp_mat_1B.cols(); i++) {
			// bars define two cylindrical regions, if point within this region, cost will be negative

			cost_1_1B = pow(geom.coors(1, i) - pos_1B(0, 0), 2) + pow(geom.coors(2, i) - pos_1B(1, 0), 2) - pow(radius_c_1B, 2);
			cost_2_1B = pow(geom.coors(1, i) - pos_1B(0, 1), 2) + pow(geom.coors(2, i) - pos_1B(1, 1), 2) - pow(radius_t_1B, 2);
			cost_3_1B = pow(geom.coors(1, i) - pos_1B(0, 2), 2) + pow(geom.coors(2, i) - pos_1B(1, 2), 2) - pow(radius_c_1B, 2);
			cost_4_1B = pow(geom.coors(1, i) - pos_1B(0, 3), 2) + pow(geom.coors(2, i) - pos_1B(1, 3), 2) - pow(radius_t_1B, 2);

			if (cost_1_1B <= 0 || cost_2_1B <= 0 || cost_3_1B <= 0 || cost_4_1B <= 0) {
				tmp_mat_1B(0, i) = steel;
				tmp_dens_1B(0, i) = Steel.dens;
			}
		}
		geom.mat_flags = tmp_mat_1B; geom.dens = tmp_dens_1B;	// store in Geometry object
		break;
	}
	default:
		using namespace std;
		cout << "Does not recognise case for adding steel bars" << endl;
		throw runtime_error("steel bar switch fail");
	}

}

// INITIALISE MATERIAL PROPERTIES
void init_mat_props(Material &mat, int flag) {
	// currently handles concrete and steel
	if (flag == concrete) {
		mat.E_mod = 22e9F;
		mat.v_ratio = 0.2F;
		mat.G_mod = 8.8e9F;
		mat.calc_eff_mod();
		mat.dens = 2400.0F;
		mat.crit_ts = 0.000533F;	// ORIGINALLY 0.000533F

	}
	else if (flag == steel) {
		mat.E_mod = 210e9F;
		mat.v_ratio = 0.3F;
		mat.G_mod = 78e9F;
		mat.calc_eff_mod();
		mat.dens = 8000.0F;
		mat.crit_ts = 0.01F;		// ORIGINALLY 0.01F

	}
	else {
		using namespace std;
		throw invalid_argument("material flag does not correspond to known material");
	}
}

// SET AND STORE BODY FORCE MATRIX IN GEOMTERY OBJECT
void set_body_forces(Geometry &geom, Sim_settings set, float scale, std::string test_tag) {

	// As cannot switch on string, so use this hack
	int switch_case = -1;
	if (test_tag == "TEST_0") { switch_case = 0; }
	if (test_tag == "TEST_1A") { switch_case = 1; }
	if (test_tag == "TEST_1B") { switch_case = 2; }

	int N_loaded = 0;

	switch (switch_case) {
	case 0:
	{
		// UNIFORM LOADING THROUGHOUT VOLUME

		float z_force = float(-MAX_REACTION * scale / (double(N_NUM_NODES)*set.volume));	// negative as pointing down
		MatrixXf body_force_0(3, 1);
		body_force_0 << 0, 0, z_force;														// load all nodes in neg z dir
		MatrixXf body_forces_0 = body_force_0.replicate<1, N_NUM_NODES>();					// repeat for all nodes
		N_loaded = body_forces_0.cols();													// all nodes loaded

		geom.body_forces = body_forces_0;													// store in Geometry object
		geom.nodes_loaded = N_loaded;
		break;
	}
	case 1:
	{
		// LOADING SPECIFIED IN TEST 1A DOCUMENTATION (~DISTRIBUTED)
		float x_min_1_1A = 0.20625F - 0.040F; float x_max_1_1A = 0.20625F + 0.040F;
		float x_min_2_1A = x_min_1_1A + 0.4125F; float x_max_2_1A = x_max_1_1A + 0.4125F;
		float x_min_3_1A = x_min_2_1A + 0.4125F; float x_max_3_1A = x_max_2_1A + 0.4125F;
		float x_min_4_1A = x_min_3_1A + 0.4125F; float x_max_4_1A = x_max_3_1A + 0.4125F;
		float x_min_5_1A = x_min_4_1A + 0.4125F; float x_max_5_1A = x_max_4_1A + 0.4125F;
		float x_min_6_1A = x_min_5_1A + 0.4125F; float x_max_6_1A = x_max_5_1A + 0.4125F;
		float x_min_7_1A = x_min_6_1A + 0.4125F; float x_max_7_1A = x_max_6_1A + 0.4125F;
		float x_min_8_1A = x_min_7_1A + 0.4125F; float x_max_8_1A = x_max_7_1A + 0.4125F;

		float y_min_1A = 0.125F - 0.080F; float y_max_1A = 0.125F + 0.080F;
		float z_min_1A = 0.600F - 2.0F*set.delta; float z_max_1A = 0.600F;					// applying load 2 deltas into beam

		bool x_bounds_1A, y_bounds_1A, z_bounds_1A;

		MatrixXf body_forces_1A(3, N_NUM_NODES); body_forces_1A.setZero();
		for (int i = 0; i < geom.coors.cols(); i++) {
			x_bounds_1A = (geom.coors(0, i) > x_min_1_1A) && (geom.coors(0, i) < x_max_1_1A);
			x_bounds_1A = (geom.coors(0, i) > x_min_2_1A) && (geom.coors(0, i) < x_max_2_1A) || (x_bounds_1A);
			x_bounds_1A = (geom.coors(0, i) > x_min_3_1A) && (geom.coors(0, i) < x_max_3_1A) || (x_bounds_1A);
			x_bounds_1A = (geom.coors(0, i) > x_min_4_1A) && (geom.coors(0, i) < x_max_4_1A) || (x_bounds_1A);
			x_bounds_1A = (geom.coors(0, i) > x_min_5_1A) && (geom.coors(0, i) < x_max_5_1A) || (x_bounds_1A);
			x_bounds_1A = (geom.coors(0, i) > x_min_6_1A) && (geom.coors(0, i) < x_max_6_1A) || (x_bounds_1A);
			x_bounds_1A = (geom.coors(0, i) > x_min_7_1A) && (geom.coors(0, i) < x_max_7_1A) || (x_bounds_1A);
			x_bounds_1A = (geom.coors(0, i) > x_min_8_1A) && (geom.coors(0, i) < x_max_8_1A) || (x_bounds_1A);
			y_bounds_1A = (geom.coors(1, i) > y_min_1A) && (geom.coors(1, i) < y_max_1A);
			z_bounds_1A = (geom.coors(2, i) > z_min_1A) && (geom.coors(2, i) < z_max_1A);

			if ((x_bounds_1A) && (y_bounds_1A) && (z_bounds_1A)) {
				body_forces_1A.col(i) << 0.00F,
								 		 0.00F,
									    -1.00F;
				N_loaded += 1;
			}
		}
		body_forces_1A *= MAX_REACTION * scale / (float(N_loaded)*set.volume);				// distribute load over loaded volume
		geom.body_forces = body_forces_1A;													// store in geometry object
		geom.nodes_loaded = N_loaded;
		break;
	}
	case 2:
	{
		// LOADING SPECIFIED IN TEST 1B DOCUMENTATION (~TIP LOADING)
		float x_min_1B = 1.550F; float x_max_1B = 1.750F;									// define a "box" within which nodes are loaded
		float y_min_1B = 0.000F; float y_max_1B = 0.250F;
		float z_min_1B = 0.600F - 2.0F*set.delta; float z_max_1B = 0.600F;					// now loading 2 horizons into beam
		bool x_bounds_1B, y_bounds_1B, z_bounds_1B;

		MatrixXf body_forces_1B(3, N_NUM_NODES); body_forces_1B.setZero();
		for (int i = 0; i < geom.coors.cols(); i++) {
			x_bounds_1B = (geom.coors(0, i) > x_min_1B) && (geom.coors(0, i) < x_max_1B);
			y_bounds_1B = (geom.coors(1, i) > y_min_1B) && (geom.coors(1, i) < y_max_1B);
			z_bounds_1B = (geom.coors(2, i) > z_min_1B) && (geom.coors(2, i) < z_max_1B);

			if ((x_bounds_1B) && (y_bounds_1B) && (z_bounds_1B)) {
				body_forces_1B.col(i) << 0.00F,
										 0.00F,
					                    -1.00F;
				N_loaded += 1;
			}
		}
		body_forces_1B *= MAX_REACTION * scale / (float(N_loaded)*set.volume);				// distribute load over loaded volume
		geom.body_forces = body_forces_1B;													// store in geometry object
		geom.nodes_loaded = N_loaded;
		break;
	}
	default:
		using namespace std;
		cout << "Failed on switch statement when setting body forces" << endl;
		throw runtime_error("Forces switch statement fail");
	}

	if (N_loaded == 0) {
		using namespace std;
		cout << "No nodes were loaded" << endl;
		throw runtime_error("Error loading nodes");

	}
	else
	{
		std::cout << "Loaded " << N_loaded << " out of " << N_NUM_NODES << " nodes\n";
	}

}

// INITIALISE SIMULATION SETTINGS
void init_simsets(Sim_settings &config, Geometry geom, \
	Material &Conc, Material &Steel) {
	// Stores values for all simulation settings, as well as calculating stiffness values for steel and conc bonds

	// loading parameters
//	config.damping = 2.5e6F			// ORIGINAL
	config.build_up = 500;			// build up load over N timesteps
	config.damping = 2.5e6F;		// damping coefficient

	// peridynamic constants
//	congig.delta = M_PI*geom.d_xyz[0];												// ORIGINAL
//	config.volume = pow(geom.d_xyz[0], 3);											// ORIGINAL
//	config.radij = geom.d_xyz[0] / 2.0F;											// ORIGINAL

	config.delta = M_PI * (geom.d_xyz[0] + geom.d_xyz[1] + geom.d_xyz[2]) / 3.0F;	// mean spacing times pi
	config.neigh_vol = 4.0F * M_PI*pow(config.delta, 3) / 3.0F;						// neighbourhood volume for node
	config.volume = geom.d_xyz[0] * geom.d_xyz[1] * geom.d_xyz[2];					// cell volume
	config.radij = (geom.d_xyz[0] + geom.d_xyz[1] + geom.d_xyz[2]) / (3.0F*2.0F);	// material point radius
	Conc.c_bond = (12 * Conc.E_mod) / (M_PI*pow(config.delta, 4));					// bond stiffness
	Steel.c_bond = (12 * Steel.E_mod) / (M_PI*pow(config.delta, 4));				// bond stiffness

	// time step
	config.saf_fac = 1.5F;															// 2 in mark's code
	config.calc_dt(Conc.dens, geom.d_xyz[0], Conc.c_bond);							// calc crit timestep
	config.nt = 20000;																// 10000 iterations originally

}

// SET CONSTRAINTS
void set_constraints(Geometry &geom, int constr_type) {
	// Set all points as unconstrained initially
	MatrixXf constrs(3, N_NUM_NODES); constrs.setOnes(3, N_NUM_NODES);

	switch (constr_type) {

		case cantilever :
			// fully clamp one end. Nodes within CLAMP_DISTANCE of wall in x_dir are immobile

			int i;
			for (i = 0; i < geom.coors.cols(); i++) {					// loop through all points, not speed critical:
				if (geom.coors.col(i)(0, 0) <= CLAMP_DIST) {			// done once at start of simulation.
					constrs.col(i) *= 0;								// set all constraints to 0 (fully clamped)
				}
			}
			geom.constraints = constrs;
			break;

		default :
			using namespace std;
			cout << "Does not recognise contraint type, exiting." << endl;
			throw runtime_error("Not implemented");
	}
}

// CALCULATE VOLUME CORRECTIONS
void calc_vol_corrs(Geometry &geom, Sim_settings set) {
	geom.vol_corrs.resize(N_NUM_NODES, N_NUM_NODES);									// Initialise storage
	geom.vol_corrs.reserve(VectorXi::Constant(N_NUM_NODES, AV_BOND_NUM));				// RESERVE SPACE BEFORE STARTING TO FILL!

	Eigen::initParallel();
	#pragma omp parallel if(RUN_IN_PARALLEL)
	{
		//std::cout << "PARALLEL" << std::endl;											// will output N times if threading is correct
		std::ptrdiff_t node, neigh; float len_ij, temp_fac;								// give each thread its own copy of these
		#pragma omp for schedule(dynamic, 256)
		for (int k = 0; k < geom.connect.outerSize(); ++k) {							// This nested loop iterates over the non-zero...
			for (SparseMatrix<bool>::InnerIterator it(geom.connect, k); it; ++it) {		// elements of the connectivity matrix. Giving access...
				node = it.row(); neigh = it.col();										// to bonds, and neighbour pairs.
				len_ij = (geom.coors.col(node) - geom.coors.col(neigh)).norm();			// len_ij is euclidean distance from node to neigh

				// correction factor depends on distance between neighbours len_ij
				if (len_ij <= set.delta - set.radij)
				{																		// cor factor = 1.0
					geom.vol_corrs.insert(node, neigh) = 1.0F;
				}
				else if ((len_ij > set.delta - set.radij) && (len_ij <= set.delta))
				{																		// 0.0 <= cor factor <= 1.0
					temp_fac = (set.delta + set.radij - len_ij) / (2.0F * set.radij);
					geom.vol_corrs.insert(node, neigh) = temp_fac;
				}
				else
				{																		// cor factor = 0.0
					geom.vol_corrs.insert(node, neigh) = 0.0F;
				}
			}
		}
	}
}

// SET BOND TYPES AND STIFFNESS FACTORS FOR ALL BONDS
// STORE VALUE OF MAX ALLOWED STRETCH FOR EACH BOND
void set_bond_types_stiff_facs(Geometry &geom, Sim_settings set, Material Conc, Material Steel) {
	// REWRITE TO ONLY EDIT GEOM OBJECT DIRECTLY! NEED TO SAVE MEMORY

	geom.bond_types.resize(N_NUM_NODES, N_NUM_NODES);
	geom.bond_types.reserve(VectorXi::Constant(N_NUM_NODES, AV_BOND_NUM));			// RESERVE APPROXIMATE SIZE!
																					// MEMORY USAGE CRITICAL HERE
	geom.stiff_facs.resize(N_NUM_NODES, N_NUM_NODES);
	geom.stiff_facs.reserve(VectorXi::Constant(N_NUM_NODES, AV_BOND_NUM));

	SparseMatrix<bool> connect_full = geom.connect.transpose();
	connect_full = connect_full || geom.connect;									// get full matrix for counting non-zeros
																					// want to ONLY use column vectors in loop

	geom.fail_stretches.resize(N_NUM_NODES, N_NUM_NODES);							// set size of sparse bond matrix
	geom.fail_stretches.reserve(VectorXi::Constant(N_NUM_NODES, AV_BOND_NUM));		// RESERVE APPROXIMATE SIZE!

	Eigen::initParallel();
	#pragma omp parallel if(RUN_IN_PARALLEL)
	{
		float neigh_vol_node, neigh_vol_neigh;										// neigh vols for pair (neigh called cnode in MATLAB)
		float stiff_fac, c_tmp;														// tmp values for calculations
		SparseVector<bool> col;
		std::ptrdiff_t node, neigh;												// Eigen::index is of type std::ptrdiff_t


		#pragma omp for schedule(static)
		for (int k = 0; k < geom.connect.outerSize(); ++k) {						// loop over all bonds (non-zero elements of connect)
			for (SparseMatrix<bool>::InnerIterator it(geom.connect, k); it; ++it) {	// iterate over bonds only once
				node = it.row(); neigh = it.col();									// get node pairs

				if ((geom.mat_flags(0, node) == concrete) && (geom.mat_flags(0, neigh) == concrete))	// concrete to concrete bond
				{
					c_tmp = Conc.c_bond;
					geom.bond_types.insert(node, neigh) = ctc;
					geom.fail_stretches.insert(node, neigh) = Conc.crit_ts;
				}
				else if((geom.mat_flags(0, node) == concrete) && (geom.mat_flags(0, neigh) == steel))	// concrete to steel bond
				{
					c_tmp = Conc.c_bond;
					geom.bond_types.insert(node, neigh) = cts;
					geom.fail_stretches.insert(node, neigh) = Conc.crit_ts*3.0F;						// mult force and crit_ts by 3
				}
				else if ((geom.mat_flags(0, node) == steel) && (geom.mat_flags(0, neigh) == concrete))	// steel to concrete bond
				{
					c_tmp = Conc.c_bond;
					geom.bond_types.insert(node, neigh) = stc;
					geom.fail_stretches.insert(node, neigh) = Conc.crit_ts*3.0F;
				}
				else if ((geom.mat_flags(0, node) == steel) && (geom.mat_flags(0, neigh) == steel))		// steel to steel bond
				{
					c_tmp = Steel.c_bond;
					geom.bond_types.insert(node, neigh) = sts;
					geom.fail_stretches.insert(node, neigh) = Steel.crit_ts;
				}
				else																// unrecognised, throw error
				{
					using namespace std;
					cout << "Pair of nodes tagged incorrectly, exiting." << endl;
					throw runtime_error("Did not recognise bond connection type");
				}

				// FASTEST OPTION FOUND: THINK BECAUSE ACCESSING ROW VECTORS OF COLUMN MAJOR MATRIX VERY COSTLY
				// neigh_vol = num neighs * vol of cell
				col = connect_full.col(node);
				neigh_vol_node = col.nonZeros() * set.volume;
				col = connect_full.col(neigh);
				neigh_vol_neigh = col.nonZeros() * set.volume;

				stiff_fac = (2.0F*set.neigh_vol)*c_tmp / (neigh_vol_node + neigh_vol_neigh);
				geom.stiff_facs.insert(node, neigh) = stiff_fac;
			}
		}
	}
}

// Functions relating to the actual step by step simulation

#include "integration_scheme.h"
#include "classes.h"
#include "main.h"
#include "csv_parse.h"
#include <iostream>
#include <iomanip>					// make output look nice
#include <chrono>					// time script
#include <cmath>					// for pow(a, b)

// MAIN FUNCTION FOR TIME STEP INTEGRATION
bool run_time_integration(Geometry &geom, Sim_settings set, Material Conc, Material Steel, std::string out, std::string load_tag) {
	using namespace std;
	int w = 70;																	// width of cout text column
	bool verbose = false;														// if true: cout timing info
	chrono::steady_clock::time_point start = chrono::steady_clock::now();		// begin timer
	chrono::steady_clock::time_point end = chrono::steady_clock::now();			// initialise now counter

	int max_its = set.nt; string time_tag;										// maximum number of iterations
	int k = 1; int save_every = 50;												// it count and save every N its
	bool converged = false;	float tolerance = 1e-5F; float av_abs_disp, av_disp;// if average node disp smaller: has converged (over loaded nodes)
//	float prev_av_disp = 0.0F;													// compare current
	float perc_damage, prev_perc_damage;										// use to store percentage of bonds broken for current step and prev step
	prev_perc_damage = 0.0F;
	int stable_steps = 0; int falling = 0;										// iterations bond damage unchanged, iterations disp in same dir as force


	while (true) {

		// USE POSITION TO WORK OUT STRETCH OF ALL BONDS USING CURRENT AND STARTING COORDINATES
		start = chrono::steady_clock::now();
		calc_bond_stretch(geom);
		//calc_bond_stretch_2(geom);
		end = chrono::steady_clock::now();

		if (verbose) {
			std::cout << left << setw(w) << "CALCULATED BOND STRETCHES";
			std::cout << "TIME TAKEN: (SECS): " << chrono::duration_cast<chrono::microseconds>(end - start).count() / pow(10, 6) << std::endl;
		}


		// USE STRETCHES TO UPDATE FAILS
		start = chrono::steady_clock::now();
		calc_fails(geom);
		end = chrono::steady_clock::now();

		if (verbose) {
			std::cout << left << setw(w) << "UPDATED FAILS";
			std::cout << "TIME TAKEN: (SECS): " << chrono::duration_cast<chrono::microseconds>(end - start).count() / pow(10, 6) << std::endl;
		}

		// USE STRETCH TO CALC FORCES ON ALL NODES
		start = chrono::steady_clock::now();
		calc_node_forces(geom, set);
		end = chrono::steady_clock::now();

		// IF IN BUILD-UP PHASE, SCALE FORCES ACCORDINGLY
		if ((k <= set.build_up) && set.build_up > 0) {
			//cout << "ATTEMPTING TO SCALE FORCES: current force sum: " << geom.forces.sum() << endl;
			geom.forces *= pow(float(k) / float(set.build_up), 2.0F);									// squared rather than linear build up
			//cout << "SCALED FORCE SUM: " << geom.forces.sum() << endl;
		}


		if (verbose) {
			std::cout << left << setw(w) << "CALCULATED FORCES ON NODES";
			std::cout << "TIME TAKEN: (SECS): " << chrono::duration_cast<chrono::microseconds>(end - start).count() / pow(10, 6) << std::endl;
		}

		// FROM FORCES GET ACCS, VELS, DISPS FOR ALL NODES
		start = chrono::steady_clock::now();
		calc_accs(geom, set);
		end = chrono::steady_clock::now();

		if (verbose) {
			std::cout << left << setw(w) << "CALCULATED NODE ACCELERATIONS";
			std::cout << "TIME TAKEN: (SECS): " << chrono::duration_cast<chrono::microseconds>(end - start).count() / pow(10, 6) << std::endl;
		}

		start = chrono::steady_clock::now();
		calc_vels(geom, set);
		end = chrono::steady_clock::now();

		if (verbose) {
			std::cout << left << setw(w) << "CALCULATED NODE VELOCITIES";
			std::cout << "TIME TAKEN: (SECS): " << chrono::duration_cast<chrono::microseconds>(end - start).count() / pow(10, 6) << std::endl;
		}

		start = chrono::steady_clock::now();
		calc_disps(geom, set);
		end = chrono::steady_clock::now();

		if (verbose) {
			std::cout << left << setw(w) << "CALCULATED NODE DISPS";
			std::cout << "TIME TAKEN: (SECS): " << chrono::duration_cast<chrono::microseconds>(end - start).count() / pow(10, 6) << std::endl;
		}

		// USE DISPS AND OLD POSITIONS TO GET NEW POSITIONS
		start = chrono::steady_clock::now();
		geom.coors += geom.disps;

		end = chrono::steady_clock::now();
		if (verbose) {
			std::cout << left << setw(w) << "UPDATED POSITIONS";
			std::cout << "TIME TAKEN: (SECS): " << chrono::duration_cast<chrono::microseconds>(end - start).count() / pow(10, 6) << std::endl;
		}

		// SAVE RESULTS AND CHECK FOR CONVERGENCE EVERY N STEPS
		if (k % save_every == 0) {

			// find percentage of broken bonds and update record
 			perc_damage = 1.0F - float(geom.connect.cast<int>().sum()) / float(geom.connect_0.cast<float>().sum());
			float change_damage = perc_damage - prev_perc_damage;
			prev_perc_damage = perc_damage;

			// save results
			cout << "COMPLETED TIMESTEP " << k << " SAVING RESULTS " << endl;
			cout << "PERCENTAGE OF BONDS INTACT: " << 100.0F*(1.00F - perc_damage) << "%" <<endl;
			time_tag = to_string(k);
			save_results(out, time_tag, geom, load_tag);

			// check if has converged or spent iteration budget
			// average absolute displacements over loaded nodes:
			{
				Matrix<bool, 3, N_NUM_NODES> tmp_bool = geom.body_forces.cast<bool>();			// = 1 if node loaded (regardless of sign)
				Matrix<float, 3, N_NUM_NODES> tmp_float = tmp_bool.cast<float>();				// cast back to float
				tmp_float = (tmp_float.array() * geom.disps.array()).matrix();					// now only disps of loaded nodes (in loading dir)
				av_abs_disp = tmp_float.cwiseAbs().sum() / float(geom.nodes_loaded);			// average absolute displacement
				av_disp = tmp_float.sum() / float(geom.nodes_loaded);							// average disp
			}

			cout << "AVERAGE ABSOLUTE NODE DISPLACEMENT: " << av_abs_disp << endl;
			cout << "AVERAGE NODE DISPLACEMENT: " << av_disp << endl;


			if (abs(change_damage) < 1e-5) {													// damage percentage stabilised, falling or vibrating?
				stable_steps += save_every;														// store how many timesteps damage has remained stable for
			}
			else
			{
				stable_steps = 0;
			}

			if (stable_steps >= 2*set.build_up)
			{																					// now damage is stable, but is beam vibrating or falling?
				if (av_disp * geom.body_forces.sum() <= 0.0F) {									// if displacement not in same direction as applied force
					falling = 0;																// not falling
					if (av_abs_disp <= tolerance)												// if movements also small, we have convergence
					{
						converged = true;
						cout << "SIMULATION CONVERGED IN " << k << " TIMESTEPS!" << endl;
						cout << "AVERAGE ABSOLUTE DIPLACEMENT OVER LOADED NODES: " << av_abs_disp << endl;
						cout << "BOND DAMAGE: " << perc_damage * 100.0F << "%" << endl;
						break;
					}
				}
				else
				{
					falling += save_every;														// could be falling, need to see if persistent
				}

				if (falling >= 2 * set.build_up) {
					converged = false;
					cout << "BEAM BROKE, FALLING" << endl;
					cout << "BOND DAMAGE: " << perc_damage * 100.0F << "%" << endl;
				//	break;																		// don't break here, so can save animation
				}
			}

			// if over budget and not converged beam has failed
			if (k >= max_its)
			{
				converged = false;
				cout << "SIMULATION DID NOT CONVERGE, REACHED MAX ITERATION BUDGET OF " << max_its << " ITERATIONS" << endl;
				cout << "BOND DAMAGE: " << perc_damage * 100.0F << "%" << endl;
				break;
			}

		}

		// INCREMENT ITERATION COUNTER
		k += 1;

	}
	return converged;

}

// COORDINATES -> STRETCHES
void calc_bond_stretch(Geometry &geom) {
	geom.bond_stretches.resize(N_NUM_NODES, N_NUM_NODES);								// set size of matrices
	//geom.bond_stretches.reserve(VectorXi::Constant(N_NUM_NODES, AV_BOND_NUM));
	//geom.abs_stretches.resize(N_NUM_NODES, N_NUM_NODES);								// NOTE: look into using reserve properly
	//geom.abs_stretches.reserve(VectorXi::Constant(N_NUM_NODES, AV_BOND_NUM));

	float cur_len, init_len, stretch;
	std::ptrdiff_t node, neigh;
	for (int k = 0; k < geom.connect.outerSize(); ++k) {
		for (SparseMatrix<bool>::InnerIterator it(geom.connect, k); it; ++it) {
			node = it.row(); neigh = it.col();											// i,j of non-zero element
			cur_len = (geom.coors.col(node) - geom.coors.col(neigh)).norm();			// current bond length
			//init_len = geom.init_bond_lens.coeff(node, neigh);						// inital bond length
			init_len = (geom.coors_0.col(node) - geom.coors_0.col(neigh)).norm();

			// THIS IS HOW STRETCH IS CALCULATED
			stretch = (cur_len - init_len) / init_len;

			geom.bond_stretches.insert(node, neigh) = stretch;							// insert v coeffRef?

		}
	}

}

void calc_bond_stretch_2(Geometry &geom) {
	// get full connectivity matrix as sparse integer matrix
	SparseMatrix<float> con = geom.connect.transpose().cast<float>();
	con += geom.connect.cast<float>(); con.makeCompressed();

	MatrixXf Xs_row = geom.coors.row(0); MatrixXf Xs_col = Xs_row.transpose();
	//ArrayXXf Ys_row = geom.coors.row(1).array(); ArrayXXf Ys_col = Ys_row.transpose();
	//ArrayXXf Zs_row = geom.coors.row(2).array(); ArrayXXf Zs_col = Zs_row.transpose();
	SparseMatrix<float> X = con.cwiseProduct(Xs_row.replicate<N_NUM_NODES, 1>());
	X -= con.cwiseProduct(Xs_col.replicate<1, N_NUM_NODES>());

	// map to array to allow elementwise operations
	//Map<ArrayXf> connect(con.valuePtr(), con.nonZeros());
	//SparseMatrix<float> X = (connect * (Xs_col.replicate<N_NUM_NODES, 1>()) - connect * (Xs_row.replicate<1, N_NUM_NODES>())).matrix().sparseView();
	//SparseMatrix<float> Y = (connect * Ys_col.replicate<N_NUM_NODES, 1>() - connect * Xs_row.replicate<1, N_NUM_NODES>()).matrix().sparseView();
	//SparseMatrix<float> Z = (connect * Zs_col.replicate<N_NUM_NODES, 1>() - connect * Xs_row.replicate<1, N_NUM_NODES>()).matrix().sparseView();
	//SparseMatrix<float> X = connect_full

}

// STRETCHES -> FORCES
void calc_node_forces(Geometry &geom, Sim_settings set) {

	float bf_multi,	tmp_force;															// empirical fudge factor for cts and stc bonds
	VectorXf force_vec(3, 1);															// force vector for bond
	VectorXf vec_node_to_neigh(3, 1);													// direction vector from node to neigh
	geom.forces.resize(3, N_NUM_NODES); geom.forces = geom.body_forces;					// NOTE: SET BODY FORCES FIRST
	int num_bonds_loaded = 0;

	std::ptrdiff_t node, neigh;
	for (int k = 0; k < geom.connect.outerSize(); ++k) {
		for (SparseMatrix<bool>::InnerIterator it(geom.connect, k); it; ++it) {
			node = it.row(); neigh = it.col();

			if (geom.connect.coeffRef(node, neigh) == false) {
				continue;
			}

			if (geom.bond_types.coeffRef(node, neigh) == cts || geom.bond_types.coeffRef(node, neigh) == stc)
			{
				bf_multi = 3.0;
			}
			else
			{
				bf_multi = 1.0;
			}

			// force magnitude
			tmp_force = bf_multi * geom.stiff_facs.coeffRef(node, neigh) *\
						geom.bond_stretches.coeffRef(node, neigh) *\
						set.volume * geom.vol_corrs.coeffRef(node, neigh);				// now have magnitude of bond force

			// force direction
			vec_node_to_neigh = geom.coors.col(neigh) - geom.coors.col(node);			// get direction vector from node to neigh
			vec_node_to_neigh /= vec_node_to_neigh.norm();								// normalise

			force_vec = vec_node_to_neigh * tmp_force;									// with pos strain, points from node to neigh

			// now add to force sum on node, and add negative to neighbour
			geom.forces.col(node)  += force_vec;
			geom.forces.col(neigh) -= force_vec;
			num_bonds_loaded += 1;
		}
	}

	int bonds_connect_matrix = geom.connect.cast<int>().sum();
	if (num_bonds_loaded != bonds_connect_matrix) {
		using namespace std;
		cout << "NUMBER OF LOADED BONDS: " << num_bonds_loaded << endl;
		cout << "NUMBER OF BONDS IN CONNECTIVITY MATRIX: " << bonds_connect_matrix << endl;
		throw runtime_error("INCORRECT NUMBER OF BONDS LOADED");
	}
}

// FORCES AND VELOCITIES TO ACCELERATIONS
void calc_accs(Geometry &geom, Sim_settings set) {
	geom.accs.resize(3, N_NUM_NODES); geom.accs.setZero();
	geom.accs -= geom.vels*set.damping;									// damping term
	geom.accs += geom.forces;											// add force term
	geom.accs = (geom.accs.array() / geom.dens.array()).matrix();		// divide by density to get acceleration
}

// ACCELERATIONS TO VELOCITIES
void calc_vels(Geometry &geom, Sim_settings set) {
	MatrixXf old_vels = geom.vels;										// v_t
	geom.vels = old_vels + geom.accs*set.dt;							// v_t+1 = v_t + a_t*dt
}

// VELOCITIES TO DISPLACEMENTS
void calc_disps(Geometry &geom, Sim_settings set) {
	geom.disps = ((geom.vels*set.dt).array()*geom.constraints.array()).matrix();

	//for (int i = 900; i < 1050; i++) {
	//	std::cout << geom.disps.col(i) << std::endl;
	//}
}

// USE VALUES OF STRETCH TO COMPUTE BONDS THAT HAVE FAILED
void calc_fails(Geometry &geom) {
	SparseMatrix<float> bond_healths(N_NUM_NODES, N_NUM_NODES);			// max stretch - current stretch
	bond_healths = geom.fail_stretches - geom.bond_stretches.cwiseAbs();// if <0 bond has failed
																		// bit awkward as no .abs() function for sparse matrices
	// loop through all health values, check if < 0
	std::ptrdiff_t node, neigh;
	for (int k = 0; k < bond_healths.outerSize(); ++k) {
		for (SparseMatrix<float>::InnerIterator it(bond_healths, k); it; ++it) {
			node = it.row(); neigh = it.col();
			if (it.value() < 0.0F) {

				// BOND HAS FAILED
				geom.connect.coeffRef(node, neigh) = false;				// remove from connect matrix
				//geom.fail_stretches.coeffRef(node, neigh) = 0.0F;		// set max stretch to zero to ensure bond never reset

				// use .prune() maybe? explicitly removes zero values and makes compressed

 			}

		}
	}

}

// CALCULATE DEFLECTION OF LOADED NODES
float calc_loaded_disp(Geometry geom) {

	// DISP OF ALL NODES SINCE START
	Matrix<float, 3, N_NUM_NODES> tot_disps = geom.coors - geom.coors_0;

	// AVERAGE OVER LOADED NODES
	MatrixXf vec_filter = geom.body_forces.colwise().any();							// find loaded nodes
	ArrayXXf filter = vec_filter.replicate<3, 1>().array();							// col is vec of 1s if node loaded in any dir

	ArrayXXf filtered_disps = tot_disps.array() * filter;							// only non-zero disps for loaded nodes

	Vector3f av_disp = filtered_disps.rowwise().sum() / float(geom.nodes_loaded);

//	tmp_float = (tmp_float.array() * geom.disps.array()).matrix();					// now only disps of loaded nodes (in loading dir)
//	av_abs_disp = tmp_float.cwiseAbs().sum() / float(geom.nodes_loaded);			// average absolute displacement
//	av_disp = tmp_float.sum() / float(geom.nodes_loaded);							// average disp

	// ONLY RETURNS Z DISP AT MOMENT
	return av_disp(2,0);
}

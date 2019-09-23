// Peridynamics Eigen implementation
// started Jan 2019 -- Preben Monteiro Ness

#include "main.h"					// import constants
#include "classes.h"				// class definitions,
#include "csv_parse.h"				// for reading and writing .csv data
#include "phys_props.h"				// for setting materials, correction factors
#include "integration_scheme.h"		// for actual simulation

#include <iostream>					// input/output to terminal
#include <fstream>					// read/write files
#include <iomanip>					// make output look nice
#include <chrono>					// time script
#include <cmath>					// for pow(a, b)
#include <boost/filesystem.hpp>		// for organising directories

using namespace std; using namespace Eigen;

int main(int argc, char *argv[]) {
	chrono::steady_clock::time_point begin = chrono::steady_clock::now();		// begin timer
	chrono::steady_clock::time_point now = chrono::steady_clock::now();			// initialise now counter

	//TESTING AND DEBUGGING AREA	###########################################
	{
	//	MatrixXf b(3, 2);
	//	b << 0.0F, -1.0F,
	//		 0.0F, -3.6F,
	//		 0.0F, -0.0F;
	//	cout << "test matrix b" << b << endl;
	//	MatrixXf b_ar = b.colwise().any();
	//	cout << "colwise any" << b_ar << "\n rows cols of above: " << b_ar.rows() << ", " << b_ar.cols() << endl;
	//	ArrayXXf b_arrr = b_ar.replicate<3, 1>().array();
	//	cout << "replicated back to correct dims" << b_arrr << endl;
	//
	//	SparseMatrix<bool> test; test.resize(4, 4);
	//	test.insert(0, 0) = true; test.insert(1, 1) = false;
	//	SparseMatrix<int> test_int = test.cast<int>();
	//	cout << "this is test as bool:\n" << test << "\npaused\n";
	//	cout << "num non zeros in bool:\n" << test.nonZeros() << endl;
	//	cout << "test as int:\n" << test_int << endl;
	//	cout << "num non zeros in int:\n" << test_int.nonZeros() << endl;
	//	cout << "sum of int:\n" << test.cast<int>().sum() << endl;
	//	char pause = cin.get();
	//	MatrixXf b(3, 2);
	//	b << 0.0F, -1.0F,
	//		 0.1F, -3.6F,
	//		 0.0F, -0.0F;
	//
	//	cout << "this is the float matrix b:\n" << b << endl;
	//	Matrix<bool, 3, 2> bb;
	//	bb = b.cast<bool>();
	//	bb = b.cast<bool>();
	//	cout << "this is the bool matrix b:\n" << bb << endl;
	//	b = bb.cast<float>();
	//	cout << "this is b back to float:\n" << b << endl;
	//	char a = cin.get();
	}
	// END OF TESTING AND DEBUGGING AREA ######################################

	int w = 70;									// sets width of print field for cout

	if (RUN_IN_PARALLEL) {
		#pragma omp parallel
		{
			if (omp_get_thread_num() == 0) {
				cout << "Running some sections in parallel, NUM_THREADS: " << omp_get_num_threads() << endl;
			}
		}
	}

	// PARSE COMMAND LINE ARGUMENTS
	string csv_path, out_path;					// directory of .csv files, and output directory
	if (argc != 3) {
		cout << "Incorrect number of command line arguments specified! Exiting\n" << endl;
		cout << "argc: " << argc << endl;
		throw runtime_error("argc_error");
	}
	else
	{
		// Specify location of csv files as command line argument
		csv_path = argv[1];						// directory with .csv files
		out_path = argv[2];						// directory where output is stored
		cout << "Loading .csv files from this directory: " << csv_path << endl;
		cout << "Storing output in this directory: " << out_path << endl;
	}

	// CREATE DIRECTORY STRUCTURE TO STORE RESULTS
	//boost::filesystem::path out_csv_path = out_path + "/" + TEST_TAG + "/" + to_string(N_NUM_NODES) + "-nodes";
	string out_csv_path = out_path + "/" + TEST_TAG + "/" + to_string(N_NUM_NODES) + "-nodes";
	boost::filesystem::create_directories(out_csv_path);

//	char pause = cin.get();

	// INITIALISE MATERIAL PROPERTIES
	Material Concrete, Steel;
	init_mat_props(Concrete, concrete); init_mat_props(Steel, steel);

	now = chrono::steady_clock::now();
	cout << left << setw(w) << "Initialised materials";
	cout << "Time elapsed (secs): " << chrono::duration_cast<chrono::microseconds>(now - begin).count() / pow(10, 6) << endl;

	// LOAD GEOMETRY FROM .CSV FILES INTO GEOMETRY CLASS
	Geometry geom;								// geom stores current state, coors_0 and connect_0 for t = 0

	now = chrono::steady_clock::now();
	cout << left << setw(w) << "Attempting to load .csv files";
	cout << "Time elapsed (secs): " << chrono::duration_cast<chrono::microseconds>(now - begin).count() / pow(10, 6) << endl;

	load_geom(csv_path, geom, TEST_TAG);
	now = chrono::steady_clock::now();
	cout << left << setw(w) << "Successfully loaded .csv files and set initial bond lengths";
	cout << "Time elapsed (secs): " << chrono::duration_cast<chrono::microseconds>(now - begin).count() / pow(10, 6) << endl;

	// SET SIMULATION SETTINGS AND PARAMETERS
	Sim_settings peri_config;
	init_simsets(peri_config, geom, Concrete, Steel);

	now = chrono::steady_clock::now();
	cout << left << setw(w) << "Set simulation settings";
	cout << "Time elapsed (secs): " << chrono::duration_cast<chrono::microseconds>(now - begin).count() / pow(10, 6) << endl;

	// SET MATERIAL FLAGS AND DENSITY
	// NOTE: EXPORT THIS TO PREPROCESSING STAGE AND STORE AS .CSV FILE IN FUTURE
	init_mat_flags(geom, Concrete);				// initialises all nodes as concrete
	add_steel_bar(geom, Steel, TEST_TAG);		// adds steel bars to beam

	MatrixXf tmp = geom.dens.replicate<3, 1>();
	geom.dens = tmp;							// repeat each row to get a dens value for each dimension
	tmp.resize(0, 0);							// memory is precious

	now = chrono::steady_clock::now();
	cout << left << setw(w) << "Set material flags";
	cout << "Time elapsed (secs): " << chrono::duration_cast<chrono::microseconds>(now - begin).count() / pow(10, 6) << endl;

	// SET FORCES
	set_body_forces(geom, peri_config, 1.0F, TEST_TAG);	// store body force matrix in Geometry object, scale by 1.0

	now = chrono::steady_clock::now();
	cout << left << setw(w) << "Set applied body forces";
	cout << "Time elapsed (secs): " << chrono::duration_cast<chrono::microseconds>(now - begin).count() / pow(10, 6) << endl;

	// SET CONSTRAINTS AND NO FAIL REGION
	set_constraints(geom, cantilever);			// fully clamped cantilever
	geom.set_no_fails();						// no fail bond iff at least one node in a bond is constrained

	now = chrono::steady_clock::now();
	cout << left << setw(w) << "Set constraints and no fail flags";
	cout << "Time elapsed (secs): " << chrono::duration_cast<chrono::microseconds>(now - begin).count() / pow(10, 6) << endl;

	// CALCULATE VOLUME CORRECTIONS
	calc_vol_corrs(geom, peri_config);

	now = chrono::steady_clock::now();
	cout << left << setw(w) << "Calculated volume corrections for all bonds";
	cout << "Time elapsed (secs): " << chrono::duration_cast<chrono::microseconds>(now - begin).count() / pow(10, 6) << endl;

	// SET BOND TYPES, CALCULATE AND STORE STIFFNESS FACTORS
	// STORE VALUE OF MAX STRETCH FOR EACH BOND
	set_bond_types_stiff_facs(geom, peri_config, Concrete, Steel);

	now = chrono::steady_clock::now();
	cout << left << setw(w) << "Stored all bond types, calculated and stored all stiffness factors";
	cout << "Time elapsed (secs): " << chrono::duration_cast<chrono::microseconds>(now - begin).count() / pow(10, 6) << endl;

	// START SIMULATION
	geom.vels.resize(3, N_NUM_NODES); geom.vels.setZero();		// initialise velocities to zero for t = 0
	calc_bond_stretch(geom);									// initialise stretches to zero for t = 0


	// BINARY SEARCH FOR FAILURE
	bool bi_conv = false; float tol = 0.01F;					// binary search converged 0.01F for +- 2.5% on end value
	float hi = 5.0F; float lo = -5.0F;							// scale forces by pow(10,x)
	float scale; bool safe;	int num_tested = 0;					// scale forces, check if load breaks beam
	string load_tag, tmp_path, loaded_node_disp;				// for specifying output directory and files
	while (!bi_conv) {
		scale = (hi + lo) / 2.0F;											// set mean as test value for force
		load_tag = "L" + to_string(MAX_REACTION * pow(10, scale)) + "N";	// store info about loading

		// STORE RESULTS IN NEW DIRECTORY
		tmp_path = out_csv_path + "/" + load_tag;
		boost::filesystem::create_directories(tmp_path);
		save_results(out_path, "0", geom, load_tag);			// make sure simulation saved for t = 0

		// CHECK MEAN VALUE
		geom.coors = geom.coors_0;								// set coordinates to initial values
		geom.vels.resize(3, N_NUM_NODES); geom.vels.setZero();	// set velocities to zero
		geom.connect = geom.connect_0;							// set connectivity matrix to inital value
		set_body_forces(geom, peri_config, pow(10.0F, scale), TEST_TAG);
		safe = run_time_integration(geom, peri_config, Concrete, Steel, out_path, load_tag);
		num_tested += 1;

		// ADJUST BOUNDS
		if (safe) {
			lo = scale;											// if survived load, this is new lower bound
			loaded_node_disp = to_string(calc_loaded_disp(geom));
		}
		else {
			hi = scale;											// if broken, this is new upper bound
			loaded_node_disp = "BEAM FAILED";
		}

		// STORE VALUES OF LOADED NODE DISPLACEMENTS
		ofstream disps(out_csv_path + "/" + load_tag + "/loaded_node_disp.txt");
		if (!disps.is_open()) {
			cout << "Could not open disp file to store results, exiting\n";
			throw runtime_error("failed on results save");
		}
		else
		{
			disps << "AVERAGE LOADED NODE DISPLACEMENT: " << loaded_node_disp << "\n";
		}
		disps.close();

		// DISPLAY PROGRESS AND CHECK CONVERGENCE
		if (hi - lo <= tol) { bi_conv = true; }					// if bounds close, terminate loop

		cout << "\nNUMBER TESTED: " << num_tested << endl;
		cout << " BOUNDS: " << "\nUB: " << MAX_REACTION * pow(10, hi) << "\nLB: " << MAX_REACTION * pow(10, lo) << endl << endl;

	}
	scale = (hi + lo) / 2.0F;									// save final value
	float fail_load = MAX_REACTION * pow(10.0F, scale);			// store fail load this is distributed over all nodes
	float sum_body_forces = geom.body_forces.sum() * geom.nodes_loaded * peri_config.volume / geom.nodes_loaded;

	cout << "scale = (hi-lo)/2: " << scale << endl;
	cout << "\nFailure load from summed body forces: " << sum_body_forces << endl;
	cout << "\n\nFailure load from applied loading: " << fail_load << " N" << endl;

	// STORE RESULTS ON FAILURE LOAD AND BOUNDS
	ofstream summary(out_csv_path + "summary.txt");				// attempt to open output file
	if (!summary.is_open()) {
		cout << "Could not open summary file to store results, exiting\n";
		throw runtime_error("failed on results save");
	}
	else
	{
		summary << "BEAM TEST COMPLETED FOR: " << TEST_TAG << "\n";
		summary << "NUMBER TESTED: " << num_tested << "\n";
		summary << "BOUNDS: " << "\nUB: " << MAX_REACTION * pow(10, hi) << "\nLB: " << MAX_REACTION * pow(10, lo) << "\n";
	}
	summary.close();

	char a = cin.get();											// because visual studio is being horrible on windows...
	return 0;
}

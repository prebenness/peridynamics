#ifndef MAIN_H
#define MAIN_H

// LIBRARIES
#include <Eigen/Eigen>
#include <vector>
#include <string>

// CONSTANTS
// MARK'S BEAM:
//const std::string TEST_TAG = "TEST_0";
//const int N_NUM_NODES = 2500;				// debug
//const int N_NUM_NODES = 67500;			// full geometry

// TEST 1A
const std::string TEST_TAG = "TEST_1A";
//const int N_NUM_NODES = 952;
const int N_NUM_NODES = 2970;
//const int N_NUM_NODES = 4392;
//const int N_NUM_NODES = 6048;
//const int N_NUM_NODES = 12376;
//const int N_NUM_NODES = 25080;

// TEST 1B
//const std::string TEST_TAG = "TEST_1B";
//const int N_NUM_NODES = 792;
//const int N_NUM_NODES = 2652;
//const int N_NUM_NODES = 4095;
//const int N_NUM_NODES = 6256;
//const int N_NUM_NODES = 12760;
//const int N_NUM_NODES = 23328;

const int AV_BOND_NUM = N_NUM_NODES / 15;	// each node connected to ~3.6% of other nodes (for N = 2500), OVERESTIMATE!

const bool RUN_IN_PARALLEL = false;			// if set to true, will make use of multithreading where possible

// LOADING PARAMETERS
const double MAX_REACTION = 4.0e5;			// sum of all forces [Newton] 400kN a good fail estimate
const float CLAMP_DIST = 0.05F;				// x-distance from wall of clamped region in cantilever

const float M_PI = 3.14159265359F;			// "math.h" seems to fail to set this



#endif

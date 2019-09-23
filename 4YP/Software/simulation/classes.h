#ifndef CLASSES_H
#define CLASSES_H

#include <Eigen/Eigen>
#include <cstdlib>

using namespace Eigen;

class Geometry {
	// The Geometry class stores node coordinates, bonds, material flags, volume corrections, stiffness corrections
public:
	// BONDS: ALL MATRICES SPARSE AND UPPER TRIANGULAR
	// LOOK INTO USING <int> or even <bool> connectivity matrix to save memory
	alignas(32) SparseMatrix<bool> connect;			// [NxN] connectivity matrix: [i,j] = TRUE if bond from node i to j
	alignas(32) SparseMatrix<bool> connect_0;		// connect at t = 0
	alignas(32) SparseMatrix<float> vol_corrs;		// [NxN] volume corrections: [i,j] = correction factor i to j if bonded
	alignas(32) SparseMatrix<int> bond_types;		// [NXN] stores bond types for all bonds
	alignas(32) SparseMatrix<float> stiff_facs;		// [NxN] stiffening factor for all bonds
	alignas(32) SparseMatrix<float> init_bond_lens;	// [NxN] absolute value of bond lens at t=0, used for stretch calc
	alignas(32) SparseMatrix<bool> no_fails;		// [NxN] (i,j) set to TRUE iff i and j bonded and at least one of i,j constrained
	alignas(32) SparseMatrix<float> bond_stretches;	// [NxN] scalar value of stretch for bond i,j if bonded
	alignas(32) SparseMatrix<float> fail_stretches;	// [NxN] store value of stretch at which bond fails

	// NODES: MATRICES DENSE
	// SCALARS
	MatrixXi mat_flags;					// [1xN] material flags for all nodes
	MatrixXf dens;						// [3xN] density values for all nodes, for each dimension
	// VECTORS
	MatrixXf coors;						// [3xN] coordinates matrix: coors.cols(i) = {X_i, Y_i, Z_i}'
	MatrixXf coors_0;					// [3xN] coordinates matrix at time t = 0
	MatrixXf constraints;				// [3xN] matrix of constraints: {1,0,0}' means free to move in x, not in y or z
	MatrixXf body_forces;				// [3xN] matrix of applied body forces on each node
	MatrixXf forces;					// [3xN] matrix of total node forces: forces.cols(i) = {Fx_i, Fy_i, Fz_i}'
	MatrixXf accs;						// [3xN] matrix of total node accelerations
	MatrixXf vels;						// [3xN] matrix of total node velocities
	MatrixXf disps;						// [3xN] matrix of total node displacements

	// SINGLE VALUES
	//int N_num_nodes;					// number of nodes in simulation
	std::vector<float> d_xyz;			// {dx, dy, dz}
	int nodes_loaded;					// store number of nodes loaded

	// FUNCTIONS
	void set_no_fails();				// Fill no fail matrix
	void calc_init_bond_lens();			// fill init_bond_lens matrix
};

class Material {
public:
	float E_mod;						// youngs modulus
	float v_ratio;						// poissons ratio
	float G_mod;						// shear modulus
	float eff_mod;						// effective modulus
	float dens;							// density
	float crit_ts;						// critical tensile strain
	float c_bond;						// bond stiffness
	void calc_eff_mod();				// func to calc eff mod
};

class Sim_settings {
public:
	float delta;						// horizon
	float neigh_vol;					// neighbourhood volume for node
	float volume;						// cell volume
	float radij;						// material point radius
	float saf_fac;						// safety factor in crit timestep
	float dt;							// time step
	int nt;							// number of timesteps
	float max_force;					// max force per node
	int build_up;						// build load over N timesteps
	float damping;						// damping coefficient
	void calc_dt(float dens, float dx, float c_bond);
};

// ENUMERATIONS
// Materials
enum material {			// material flags
	concrete = 0,
	steel = 1
};

// Connection types
enum mat_con {			// bond type flags
	ctc = 1,			// concrete to concrete
	cts = 2,			// concrete to steel
	stc = 3,			// steel to concrete
	sts = 4				// steel to steel
};

// Constraints
enum constraint_types {
	cantilever = 0
};

#endif

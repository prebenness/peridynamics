// Classes, enumerations used in peridynamics simulation

#include "classes.h"
#include "main.h"

#include <cmath>

// Geometry
void Geometry::set_no_fails() {
	// TODO:	STOP ACCESSING ROW VECTORS FROM CONENCT MATRIX
	// Set bond i,j to 1 (no fail) iff i,j bonded and at least one of i,j constrained

	no_fails.resize(N_NUM_NODES, N_NUM_NODES);
	no_fails.reserve(VectorXi::Constant(N_NUM_NODES, AV_BOND_NUM));

	SparseMatrix<bool> connect_full = connect.transpose();
	connect_full = connect_full || connect;										// connect_full += connect (just with booleans)

	// Loop through all nodes and if constrained set ALL connected nodes to 1 in tmp_no_fails
	int i;
	for (i = 0; i < N_NUM_NODES; i++) {
		if (!constraints.col(i).all()) {										// if not all elements set to 1, node is constrained
			no_fails.col(i) = connect_full.col(i);								// connect up. tri., so this finds all neighbours
		}
	}
	no_fails = no_fails.triangularView<Upper>();								// symmetric matrix, so need store only up. tria. half
}

// Fill the init_bond_len matrix
void Geometry::calc_init_bond_lens() {
	// resize and reserve memory
	//init_bond_lens.resize(N_NUM_NODES, N_NUM_NODES);
	//init_bond_lens.reserve(VectorXi::Constant(N_NUM_NODES, AV_BOND_NUM));

	// loop over all bonds once
	std::ptrdiff_t node, neigh;
	for (int k = 0; k < connect.outerSize(); ++k) {
		for (SparseMatrix<bool>::InnerIterator it(connect, k); it; ++it) {
			node = it.row(); neigh = it.col();
			//init_bond_lens.insert(node, neigh) = (coors_0.col(node) - coors_0.col(neigh)).norm();
		}

	}
}

// Material
void Material::calc_eff_mod() {
	eff_mod = E_mod / ((1 - 2 * v_ratio)*(1 + v_ratio));
}

// Simulation settings
void Sim_settings::calc_dt(float dens, float dx, float c_bond) {
	dt = (0.8F*pow(2.0F * dens*dx / (M_PI*pow(delta, 2.0F)*dx*c_bond), 0.5F)) / saf_fac;
}

// SUBFUNCTIONS FOR GENERATING DISCRETISATIONS AND BONDS

#include "node_bonds.h"
#include "main.h"

#include <iostream>
#include <Eigen/Eigen>

using namespace Eigen;

// GENERATE UNIFORM DISCRETISAION
Eigen::MatrixXf gen_uniform_disc(Vector3f lens, Vector3i divs, Vector3f &spacings) {

	// WORK OUT SPACINGS IN ALL THREE DIRECTIONS
	//Vector3f spacings(3, 1);
	spacings << lens(0, 0) / divs(0, 0),
				lens(1, 0) / divs(1, 0),
				lens(2, 0) / divs(2, 0);
	
	//delta = spacings.sum() * M_PI / 3.0F;								// horizon is average spacing * M_PI
	
	int N = (divs(0, 0) + 1)*(divs(1, 0) + 1)*(divs(2, 0) + 1);			// number of nodes
	MatrixXf coors(3, N); coors.setZero();								// coordinates matrix to fill

	// USE NESTED FOR LOOPS TO GET POINT COORDINATES
	int node_count = 0; Vector3f tmp_point;
	for (int i = 0; i <= divs(0, 0); i++) {								// NOTE: for D divisions in some direction, get D+1 nodes
		for (int j = 0; j <= divs(1, 0); j++) {
			for (int k = 0; k <= divs(2, 0); k++) {
				
				tmp_point << spacings(0, 0)*float(i),
							 spacings(1, 0)*float(j),
							 spacings(2, 0)*float(k);
				
				coors.col(node_count) = tmp_point;
				node_count += 1;
			}
		}
	}
	return coors;
}

Eigen::SparseMatrix<bool> find_neighbours(MatrixXf coors, float delta) {
	int N = coors.cols();												// number of nodes
	SparseMatrix<bool> connect;											// declare and reserve memory for connect matrix
	connect.resize(N, N); connect.reserve(VectorXi::Constant(N, N / 15));
	float delta_2 = delta * delta;

	// NAIVE EXHAUSTIVE SEARCH FOR NEIGHBOURS
	for (int i = 0; i < coors.cols(); i++) {
		for (int j = 0; j < coors.cols(); j++) {
			if (i < j) {												// save some time by only checking each bonds once
				if ((coors.col(i) - coors.col(j)).squaredNorm() <= delta_2) {
					connect.insert(i, j) = true;
				}
			}
		}
	}
	return connect;
}
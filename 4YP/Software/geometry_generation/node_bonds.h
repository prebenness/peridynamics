#ifndef NODE_BONDS_H
#define NODE_BONDS_H

#include <Eigen/Eigen>

// FUNCTIONS
// DISCRETISATIONS
Eigen::MatrixXf gen_uniform_disc(Eigen::Vector3f lens, Eigen::Vector3i divs, Eigen::Vector3f &spacings);

// BOND CREATION
Eigen::SparseMatrix<bool> find_neighbours(Eigen::MatrixXf coors, float delta);

#endif
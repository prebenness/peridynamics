#ifndef MAIN_H
#define MAIN_H

#include <string>
#include <Eigen/Eigen>

const double M_PI = 3.14159265359;		// for convenience as cannot seem to get this from math.h lib

// FUNCTIONS
// STORE TO .CSV FILES
void store_coors(std::string out_dir, Eigen::MatrixXf coors, std::string tag);
void store_connect(std::string out_dir, Eigen::SparseMatrix<bool> connect, std::string tag);
void store_spacings(std::string out_dir, Eigen::Vector3f spacings, std::string tag, int N);

#endif
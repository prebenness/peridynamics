// FUNCTIONS FOR READING AND WRITING CSV FILES

#include <iostream>
#include <fstream>
#include <sstream>
#include <windows.h>		// for managing directories

#include <Eigen/Eigen>
#include <vector>

#include "classes.h"
#include "main.h"

using namespace Eigen; using namespace std;

// LOAD GEOMETRY FILES
void load_geom(std::string csv_path, Geometry &geom, std::string test_tag) {
	/* This function reads the .csv files specified in the directory csv_path, 
	fills bond connectivity matrix, node coordinates matrix, d_xyz values.
	*/
	
	// LOAD COORDINATES
	// Open file to read
	string coors_path = csv_path + "/coors_" + test_tag + "_" + to_string(N_NUM_NODES) + "-nodes.csv";
	ifstream coor_infile(coors_path);					// attempt to open coors .csv file
	if (!coor_infile.is_open()) {
		throw runtime_error("could not open file to load coordinates");
	}

	// Create temporary and permanent storage
	vector< vector<float> > point_coors;				// used to store coordinates for all points
	vector<float> tmp_point;							// used to store coordinates for a single point (ie from one .csv line)
	float coor_float;									// used to store value of string converted to float using stof

	// Load into std::vector, then move to Eigen storage at end
	// Read from line 0: one point per line, three coordinates per point.
	int lin_num = 0;									//count line number, indexed from 0
	while (coor_infile)
	{
		string coor_line = "";							// reset coor_line here before calling getline
		if (!getline(coor_infile, coor_line)) break;	// if at end of file, break

		istringstream ss(coor_line);					// change string line to stringstream so can use getline() on it
		point_coors.push_back(vector<float>());			// "initialise" entry i as empty vector
		tmp_point.clear();								// clear previous values
		while (ss)
		{
			string coor = "";
			if (!getline(ss, coor, ',')) break;			// split on "," if at end of line break inner loop
			coor_float = stof(coor);					// convert string to float
			tmp_point.push_back(coor_float);			// store single coordinate
		}
		point_coors[lin_num] = tmp_point;				// store this point
		lin_num += 1;
	}
	
	coor_infile.close();

	MatrixXf coors(3, lin_num);							// convert to Eigen storage scheme
	int i,j;											// using for loops as speed is not critical here
	for (i = 0; i < coors.cols(); i++) {
		for (j = 0; j < point_coors[i].size(); j++) {
			coors(j, i) = point_coors[i][j];			// NOTE: coors [3xN], point_coors [Nx3]
		}
	}
	geom.coors = coors;									// store coordinates in Geometry object
	geom.coors_0 = coors;

	coors.resize(0, 0);									// memory is precious

	// LOAD CONNECTIVIY MATRIX
	// Load into std::vector, then copy to Eigen storage
	string neighs_path = csv_path + "/neighbours_" + test_tag + "_" + to_string(N_NUM_NODES) + "-nodes.csv";
	ifstream neigh_infile(neighs_path);
	if (!neigh_infile.is_open()) {
		throw runtime_error("could not open file to load neighbours");
	}

	vector<int> tmp_point_neighs;						// store neighbour list of single node
	vector< vector<int> > neighbours;					// store all neighbour lists

	while (neigh_infile)
	{
		string neigh_line;								// important to reset value of neigh line here before calling getline
		if (!getline(neigh_infile, neigh_line)) break;	// gets next line, if at end of file, break

		istringstream ss(neigh_line);					// change string line to stringstream so can use getline()
		tmp_point_neighs.clear();
		while (ss)
		{
			string neigh;
			if (!getline(ss, neigh, ',')) break;
			tmp_point_neighs.push_back(stoi(neigh));
		}
		neighbours.push_back(tmp_point_neighs);
	}
	neigh_infile.close();

	// copy into sparse Eigen matrix
	geom.connect.resize(N_NUM_NODES, N_NUM_NODES);
	geom.connect.reserve(VectorXi::Constant(N_NUM_NODES, AV_BOND_NUM));	// IMPORTANT: reserve size for approx number of non-zeros
	int node_a, node_b;													// NOTE: neighbours[i] = {node_id1, node_id2, ...}			
	for (i = 0; i < neighbours.size(); i++) {
		node_a = i;
		for (j = 0; j < neighbours[i].size(); j++) {
			node_b = neighbours[i][j];
			if (node_b < node_a) {										// connectivity matrix upper triangular for efficiency
				geom.connect.insert(node_a, node_b) = true;				// if a connected to b, b connected to a
				//std::cout << "inserted element in connectivity matrix!" << std::endl;
			}
			//geom.connect.insert(node_a, node_b);
		}
	}

	//vector<bool> ones(N_NUM_NODES, 1); ones.setOne();
	//if (geom.connect.row(0).nonZeros() == 0) {									// to deal with different formattings
	//	SparseMatrix<bool> tmp_transpose = geom.connect.transpose();
	//	geom.connect = tmp_transpose();
	//
	//}
	//geom.connect = geom.connect.
	geom.connect_0 = geom.connect;

	// LOAD VALUES OF dx, dy, dz
	string d_xyz_path = csv_path + "/d_xyz_" + test_tag + "_" + to_string(N_NUM_NODES) + "-nodes.csv";
	ifstream d_xyz_infile(d_xyz_path);
	if (!d_xyz_infile.is_open()) {
		throw runtime_error("could not open file to load d_xyz");
	}

	vector<float> tmp_d_xyz;
	string d_xyz_line;

	getline(d_xyz_infile, d_xyz_line);
	istringstream ss(d_xyz_line);
	while (ss) {
		string d = "";
		if (!getline(ss, d, ',')) break;
		tmp_d_xyz.push_back(stof(d));
	}
	d_xyz_infile.close();
	geom.d_xyz = tmp_d_xyz;								// store in Geometry class

}

// STORE RESULTS FOR PLOTTING
void save_results(std::string out_dir, std::string time_tag, Geometry geom, std::string load_tag) {
	// out_dir is root output directory, test_tag = TEST_1A_LOAD400000 for example
	// time_tag is eg "150"
	// STORE ALL NODE COORDINATES
	float X, Y, Z;
	// eg out_dir\TEST_1A\coors_TEST_1A_LOAD400000_2970-nodes_500.csv
	//string out_path = out_dir + "/" + TEST_TAG + "/coors_" + test_tag + "_" + to_string(N_NUM_NODES) + "-nodes.csv." + time_tag;
	string out_path = out_dir + "/" + TEST_TAG + "/" + to_string(N_NUM_NODES) + "-nodes/" + load_tag + "/" + "coors.csv." + time_tag;

	ofstream out(out_path);								// attempt to open output file
	if (!out.is_open()) {
		cout << "this is the load tag: " << load_tag << endl;
		cout << "attempted to create file: " << out_path << endl;
		cout << "Could not open output file to store results, exiting\n";
		throw runtime_error("failed on results save");
	} 
	else
	{
		out << "X,Y,Z,\n";
	}

	for (int i = 0; i < geom.coors.cols(); i++) {
		X = geom.coors(0, i);
		Y = geom.coors(1, i);
		Z = geom.coors(2, i);
		out << X << ", " << Y << ", " << Z << ",\n";
	}
	out.close();

	// STORE A "POINT" FOR EACH BOND AND AN ASSOSIATED STRAIN VALUE
	//out_path = out_dir + "/" + TEST_TAG + "/bonds_" + test_tag + "_" + to_string(N_NUM_NODES) + "-nodes.csv." + time_tag;
	out_path = out_dir + "/" + TEST_TAG + "/" + to_string(N_NUM_NODES) + "-nodes/" + load_tag + "/" + "bonds.csv." + time_tag;

	ofstream out_bonds(out_path);
	if (!out_bonds.is_open()) {
		cout << "attempted to create file: " << out_path << endl;
		cout << "Could not open output file to store results, exiting\n";
		throw runtime_error("failed on results save");
	}
	else
	{
		out_bonds << "X,Y,Z,strain,health,\n";
	}

	using namespace Eigen;
	float strain, health;
	std::ptrdiff_t node, neigh;													// Eigen::index is of type std::ptrdiff_t
	for (int k = 0; k < geom.connect.outerSize(); ++k) {						// loop over all bonds (non-zero elements of connect)
		for (SparseMatrix<bool>::InnerIterator it(geom.connect, k); it; ++it) {	// iterate over bonds only once
			node = it.row(); neigh = it.col();
			X = geom.coors(0, node) + geom.coors(0, neigh); X = X * 0.50F;		// take average coordinate as bond location
			Y = geom.coors(1, node) + geom.coors(1, neigh); Y = Y * 0.50F;
			Z = geom.coors(2, node) + geom.coors(2, neigh); Z = Z * 0.50F;

			strain = geom.bond_stretches.coeffRef(node, neigh);
			health = 1.0F - abs(strain) / geom.fail_stretches.coeffRef(node, neigh);

			if (health <= 0.0F) {
				health = 0.0F;
				strain = geom.fail_stretches.coeffRef(node, neigh);
			}

			// ONLY STORE BONDS THAT ARE NOT BROKEN
			if (it.value() == true) {
				out_bonds << X << ", " << Y << ", " << Z << ", " << strain << ", " << health << ",\n";
			}
			
		}
	}
}


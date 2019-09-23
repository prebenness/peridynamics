// MAIN SCRIPT FOR GENERATING INPUT GEOMETRIES

#include "main.h"
#include "node_bonds.h"

#include <iostream>
#include <cmath>
#include <fstream>


int main(int argc, char *argv[]) {

	// SET APPROX NUMBER OF NODES WANTED FOR DISCRETISATION	
	int approx_N = 500;

	// SET COMMON NAMESPACES
	using namespace std; 
	using namespace Eigen;
	
	int expected_argc = 2;

	// PARSE COMMAND LINE ARGUMENTS
	string out_path;							// output directory for .csv files
	if (argc != expected_argc) 
	{
		cout << "Incorrect number of command line arguments specified! Expected: " << expected_argc << "\nExiting\n" << endl;
		cout << "argc: " << argc << endl;
		throw runtime_error("argc_error");
	}
	else
	{
		// Specify location of csv files as command line argument
		out_path = argv[1];						// directory where output is stored
		cout << "Storing output in this directory: " << out_path << endl;
	}

	// SET DIMENSIONS OF BEAM
	Vector3f beam_lens(3, 1);

	// ## TEST 1A
	beam_lens << 3.300F,		// X length (length)
				 0.250F,		// Y length	(width)
				 0.600F;		// Z length	(height)
	string tag = "TEST_1A";

	// ## TEST 1B
//	beam_lens << 1.650F,		// X length (length)
//				 0.250F,		// Y length	(width)
//				 0.600F;		// Z length	(height)
//	string tag = "TEST_1B";

	// SET SPACING IN X, Y, Z DIRECTION (FOR ROUGHLY UNIFORM DISCRETISATION)
	// spacing ~ (X*Y*Z/N)^(1/3) = (volume per node)^(1/3)
	float approx_S = pow(beam_lens(0, 0)*beam_lens(1, 0)*beam_lens(2, 0) / float(approx_N), 1.0F/3.0F);
	
	VectorXi num_divs(3, 1);
	num_divs << round(beam_lens(0, 0) / approx_S),
				round(beam_lens(1, 0) / approx_S),
				round(beam_lens(2, 0) / approx_S);

	// DISCRETISE
	Vector3f spacings(3, 1);					// store value of horizon -> average dx, dy, dz * M_PI
	MatrixXf coors = gen_uniform_disc(beam_lens, num_divs, spacings);
	int N_nodes = coors.cols();
	float delta = spacings.sum()*M_PI / 3.0F;
	cout << "Finished generating coordinates, used " << N_nodes << " nodes" << endl;
	cout << "Got this value for the horizon delta: " << delta << endl;

	// GENERATE NEIGHBOURS - SLOW
	cout << "Started generating connectivity matrix\n";
	SparseMatrix<bool> connect = find_neighbours(coors, delta);
	cout << "Finished generating connectivity matrix\n";

	cout << "Nodes, bonds, pos bonds: " << coors.cols() << ", " << connect.nonZeros() << ", " << coors.cols()*coors.cols() / 2 << endl;

	// STORE RESULTS
	cout << "Storing matrices as .csv files\n";
	store_coors(out_path, coors, tag); 
	cout << "Done storing coordinates\n";
	store_connect(out_path, connect, tag);
	cout << "Done storing connectivity matrix\n";
	store_spacings(out_path, spacings, tag, coors.cols());
	cout << "Done storing all, exiting\n";

	char pause = cin.get();						// because Visual Studio is being difficult
}

// STORE GEOMETRIES AS .CSV FILES
// COORDINATES
void store_coors(std::string out_dir, Eigen::MatrixXf coors, std::string tag) {
	using namespace std;
	int N = coors.cols();
	string out_path = out_dir + "\\coors_" + tag + "_" + to_string(N) + "-nodes.csv";
	
	ofstream out(out_path);								// attempt to open output file
	if (!out.is_open()) {
		throw runtime_error("Could not open output file to store coordinates, exiting");
	}

	for (int i = 0; i < coors.cols(); i++) {
		out << coors(0, i) << "," << coors(1, i) << "," << coors(2, i) << "\n";
	}
	out.close();
}

// BONDS
void store_connect(std::string out_dir, Eigen::SparseMatrix<bool> connect, std::string tag) {
	using namespace std;
	int N = connect.cols();
	string out_path = out_dir + "\\neighbours_" + tag + "_" + to_string(N) + "-nodes.csv";

	ofstream out(out_path);
	if (!out.is_open()) {
		throw runtime_error("Could not open output file to store neighbours, exiting");
	}

	for (int k = 0; k < connect.outerSize(); ++k) {
		for (Eigen::SparseMatrix<bool>::InnerIterator it(connect, k); it; ++it) {
			out << it.row() << ",";
		}
		out << "\n";
	}
	out.close();
}

// SPACINGS 
void store_spacings(std::string out_dir, Eigen::Vector3f spacings, std::string tag, int N) {
	using namespace std;
	string out_path = out_dir + "\\d_xyz_" + tag + "_" + to_string(N) + "-nodes.csv";

	ofstream out(out_path);
	if (!out.is_open()) {
		throw runtime_error("Could not open output file to store spacings, exiting\n");
	}

	out << spacings(0, 0) << "," << spacings(1, 0) << "," << spacings(2, 0) << "\n";
	out.close();
}
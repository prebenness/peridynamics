### Preben Monteiro Ness May 2019 | King's College ###
### Modelling the Failure Behaviour of Concrete Structures Using Peridynamics ###
### SUpervisor: Dr John Orr ###
### pmn26@cam.ac.uk or preben.monteiro@gmail.com ###

# OVERVIEW:
All scripts used to run simulations and produce input beam geometries are found in the directory "Software". A full code listing can also be found at https://github.com/prebenness/peridynamics

Software/simulation contains scripts for running the simulation,
Software/geometry_generation contains scripts for generating .csv input geometry files
SOftware/example_csvs shows some example input and output .csv files

# BEFORE YOU START:
The simulation requires two external C++ libraries in order to run:
1.) Eigen
	- linear algebra library, easy to install
	- download link and install instructions here:
	https://eigen.tuxfamily.org/dox/GettingStarted.html
2.) Boost
	- For dealing with creating directories and storing .csv results
	- download link and install instructions here:
	https://www.boost.org/doc/libs/1_67_0/doc/html/quickbook/install.html

# NOTE: when running scripts make sure the C++ compiler you are using is pointed correctly to the location of the Eigen and Boost libraries above. See the links for further instructions on how to do this

# GENERATING GEOMETRIES:
- Include all the script files in Software/geometry_generation in your project, or use the -I argument on the commnand line when compiling using eg g++. 
- The main script takes exactly one argument which is the location of the output directory where the geometry files are stored.
- Changing the variable approx_N sets the approximate number of nodes in output geometry.

# RUNNING SIMULATIONS:
- Include all the script files in Software/simulation in your project, or use the -I argument on the commnand line when compiling using eg g++.
- The main script takes exactly two arguments which specify the location of the input .csv files and the location where simulation .csv results are stored
	- eg specifying "\home\PD_out" will create the following directory structure when running simulations:
		- \home\PD_out\TEST_1A\6256-nodes
				      \792-nodes
				      \...
				          \L400000N
					  \L300000N
					  \...
						\bonds.csv.0
						\bonds.csv.50
						\...
						\coors.csv.0
						\...
						\loaded_node_disp.txt
						\node_disp.txt
						\bond_damage.txt
		- The bonds.csv.t stores the position of all nodes (average of two nodes connected) for timestep t, information on node health and strain is also stored.
		- The coors.csv.t stores all node positions at time t
		- node_disp.txt stores the change in \Delta u for all loaded nodes at time step for all time steps
		- bond_damage.txt stores the change in the proportion of bonds broken for all time steps
		- loaded_node_disp.txt stores the average displacement of loaded nodes if the simulation converged
		- coors.csv.t and bonds.csv.t files can be analysed using paraview to generate animations
			- paraview can be downloaded for free from here:
			https://www.paraview.org/download/

END OF FILE 
import argparse, os
from phys_props import Geometry, Material, SimulationSettings
from simulation import Simulation
from csv_parse import csvParse

models = ['1A', '1B', 'debug']
verbose = False

parser = argparse.ArgumentParser()
parser.add_argument('--in_coor_file', required=True, dest='in_coor_file')
parser.add_argument('--out_dir', required=True, dest='out_dir')
parser.add_argument('--v', action='store_true', required=False)
args = parser.parse_args()

TEST_TAG = os.path.split(args.in_coor_file)[-1].split('_')[2]
assert(TEST_TAG in models)
if args.v:
    verbose = True

with open(args.in_coor_file, 'r') as r:
    N = len(r.readlines()) - 1

def main():
    # Load and store geometry information
    geom = Geometry(N=N, verbose=verbose)
    geom.load_geom(args.in_coor_file)
    
    # Load csv parser
    csv = csvParse(geom)
    
    conc = Material('concrete')
    steel = Material('steel')
    
    # Initialise simulation
    peri_config = SimulationSettings(geom, bulk_mat=conc, rebar_mat=steel)
    
    geom.set_body_forces(settings=peri_config, tag=TEST_TAG)
    geom.set_constraints(settings=peri_config, const_tag='cantilever')
    
    geom.calc_vol_corrs(settings=peri_config)
    geom.set_bond_types_stiff_facs(settings=peri_config, bulk_mat=conc, rebar_mat=steel)
    
    geom.init_sim()
    
    # Start simulation
    sim = Simulation(geom=geom, settings=peri_config, verbose=verbose)
    sim.run_time_integration(bulk_mat=conc, rebar_mat=steel, args=args, load_tag=TEST_TAG, n_num_nodes=N)
    
    # Binary search for failure
    bi_conv = False # True is failure
    tol = 0.01 # Tolerance
    hi = 5.0
    lo = -5.0
    
    # Failure checks
    while not bi_conv:
        #scale = (hi + lo )/2.0  # set mean as test value for force
        #load_tag = 'L' + str(peri_config.MAX_REAC * pow(10, scale)) + 'N' # store info about loading
        
        # Store results in new directory
        tmp_path = args.out_dir
        csv.save_results(tmp_path, "0", geom, TEST_TAG, N) #make sure simulation saved for t=0
        bi_conv = True
    print('done')
    
if __name__ == '__main__':
    main()
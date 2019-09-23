import time, argparse, csv, os
import numpy as np
import itertools
from datetime import datetime

from scipy import spatial, constants, sparse

from plot import *
from geom_utils import *

models = ['1A', '1B', 'debug']

parser = argparse.ArgumentParser()
parser.add_argument('--out_dir', required=True, dest='out_dir')
parser.add_argument('--model', required=True, choices=models, dest='model')
parser.add_argument('--num_nodes', required=True, type=int, dest='num_nodes')
args = parser.parse_args()


start = datetime.now()
N_aprx = args.num_nodes

if args.model == '1A':
    dims = [3.3, 0.250, 0.6]
elif args.model == '1B':
    dims = [1.770, 0.250, 0.6]
elif args.model == 'debug':
    dims = [0.1, 0.2, 1.0]
else:
    raise NotImplementedError('Did not recognise model argument: {}\nModel must be one of: {}'.format(args.model, models))

s = np.power(np.prod(dims)/N_aprx, 1/3)
n_divs = [ int(round(d/s)) for d in dims ]

print('Nodes per dim: {}\nSpacing per dim: {}\nN: {}'.format(n_divs, [d/n for d, n in zip(dims, n_divs)], np.prod(n_divs)))

## Generate coordinates
# Evenly spaced points
N = np.prod(n_divs)
c = [ np.linspace(start=0, stop=d, num=n, dtype=np.float64) for d,n in zip(dims, n_divs) ]
coors = np.array( [ [x,y,z] for x in c[0] for y in c[1] for z in c[2] ] )

# Set material flags
mat_flags = set_rebar(coors, args.model)

## Find bonds
delta = constants.pi*s
tree = spatial.cKDTree(coors)
conn = sparse.csr_matrix((N,N))

# Find up to k nearest neighbours within horizon delta
# Store in sparse matrix conn
n = 0
for p in coors:
    _, res = tree.query(p, k=N, distance_upper_bound=delta, n_jobs=-1)
    row = res[res != N]
    col = np.full( (row.shape[0]), n)
    data = np.ones(row.shape[0])
    conn += sparse.csr_matrix( (data, (row, col)), shape=(N, N) )
    n += 1

print('Finished creating coors and connectivity matrix for N = {}\nRuntime: {}'.format(N, datetime.now()-start))

## Get data in csv friendly format
# Node coordinates
coordata = [['x', 'y', 'z']]
for p in coors:
    coordata.append([ c for c in p ])

# Bonds: loop over nonzero elements
neigh_list = [ [] for i in range(N) ]
cx = conn.tocoo()
for i,j in zip(cx.row, cx.col):
    if i != j:
        neigh_list[i].append(j)
bonddata = [['neighbours']] + neigh_list

# Material flags
matdata = [['material', 'material flag']] + [[materials[m]] for m in mat_flags]

## Write csv files
name = 'coors_TEST_{}_{}-nodes.csv'.format(args.model, N)
with open(os.path.join(args.out_dir, name), 'w') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerows(coordata)

name = 'neighbours_TEST_{}_{}-nodes.csv'.format(args.model, N)
with open(os.path.join(args.out_dir, name), 'w') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerows(bonddata)

name = 'd_xyz_TEST_{}_{}-nodes.csv'.format(args.model, N)
with open(os.path.join(args.out_dir, name), 'w') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerows([['dx', 'dy', 'dz'], [d/n for d, n in zip(dims, n_divs)]])

name = 'mat_flags_TEST_{}_{}-nodes.csv'.format(args.model, N)
with open(os.path.join(args.out_dir, name), 'w') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerows(matdata)

import numpy as np

materials = {
    1 : 'steel',
    2 : 'concrete',
    }

inv_materials = {
    'steel' : 1,
    'concrete' : 2,    
    }

mat_cons = {
	'btb' : 1,
	'btr' : 3,
	'rtb' : 3,
	'rtr' : 1,
    }

#inv_mat_cons = {
#    1 : 'btb',
#    2 : 'btr',
#    3 : 'rtb',
#    4 : 'rtr',
#    }

def set_rebar(X, model):
    checks = {
        '1A' : in_rebar_1A,
        '1B' : in_rebar_1B,
        }
    check = checks[model]

    mat_flags = np.zeros((X.shape[0]))
    n = 0
    for p in X:
        if check(p):
            mat_flags[n] = 1
        else:
            mat_flags[n] = 2
        n += 1
    return mat_flags

'''
Rebar funcitons expect p to be a row vector, typically [1x3]
'''
def in_rebar_1A(p):
    p = p[1:]

    # [y, z]
    bar_centers = [
        # Compressive bars
        np.array((0.031, 0.031)),
        np.array((0.219, 0.031)),

        # Tensile bars
        np.array((0.03825, 0.569)),
        np.array((0.21175, 0.569)),
        ]

    rad_c = 0.006
    rad_t = 0.01325

    radii = [
        rad_c,
        rad_c,
        rad_t,
        rad_t,
        ]

    costs = [ np.sum(np.square(cent - p) - (np.square(rad))) for cent, rad in zip(bar_centers, radii) ]
    if any( c <= 0 for c in costs ):
        return True
    else:
        return False

def in_rebar_1B(p):
    p = p[1:]

    # [y, z]
    bar_centers = [
        # Compressive bars
        np.array((0.031, 0.031)),
        np.array((0.219, 0.031)),

        # Tensile bars
        np.array((0.03825, 0.569)),
        np.array((0.21175, 0.569)),
        ]

    rad_c = 0.006
    rad_t = 0.01325

    radii = [
        rad_c,
        rad_c,
        rad_t,
        rad_t,
        ]

    costs = [ np.sum(np.square(cent - p) - (np.square(rad))) for cent, rad in zip(bar_centers, radii) ]
    if any(c <= 0 for c in costs):
        return True
    else:
        return False

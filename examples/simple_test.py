"""
Created on Thu Mar  5 20:44:28 2026

Simple test case for OT-BGS solver using LC1 data and interpolation of masses
using a Gaussian grid in latitude. Adjust the integer parameter 'skip' to 
change resolution of the target measure. Bigger skip means lower resolution.
Zero skip means full resolution.
"""

import numpy as np
import atmosphere_bgs

# Load input data
data_type = 'lc1low'
date = '200101'
dtime = '0000'

filepath = '../data/'
filename = 'bs_'+data_type+date+dtime
p00 = 1e+5
pmin = 10
input_data = atmosphere_bgs.DataLoader(filepath+filename,p00=p00,load_all=True,pmin=pmin)

# solve ot problem
solv = atmosphere_bgs.OTSolver(input_data)
solv.ot_tol = 1e-4
solv.get_bgs(verbose=True)
ld = solv.ld
solv.ld.sim.pmax=140000

# plot resulting Laguerre diagrams coloured by zonal angular momentum and potential temperature
th = ld.ys[:,1].copy()
idx_top = th == np.max(th)
th[idx_top] = np.nan

z = ld.ys[:,0].copy()
z[idx_top] = np.nan

_, ax = atmosphere_bgs.plot_lag_tess(ld,val=z,invert_yaxis=True)
ax.set_title('Potential temperature')
ax.set_xlabel('s')
ax.set_ylabel('p')

_, ax = atmosphere_bgs.plot_lag_tess(ld,val=th,invert_yaxis=True)
ax.set_title('Zonal angular momentum')
ax.set_xlabel('s')
ax.set_ylabel('p')

u = atmosphere_bgs.get_u(ld)
u[idx_top] = np.nan
_, ax = atmosphere_bgs.plot_lag_tess(ld,val=u,invert_yaxis=True)
ax.set_title('Zonal wind')
ax.set_xlabel('s')
ax.set_ylabel('p')

# rasterize theta
H = 6500

s2lat = lambda s : np.rad2deg(np.arcsin(s))
p2h = lambda p : -H*np.log(p/solv.pp.p00)

lat2s = lambda lat : np.sin(np.deg2rad(lat))
h2p = lambda h : solv.pp.p00*np.exp(-h/H)

transform = lambda x : [s2lat(x[0]), p2h(x[1])]
rast = ld.get_rasterizer(transform)

res = [1000,1000]
rv = rast.rasterize(solv.ys[:,1],res,calc_inverse=True)

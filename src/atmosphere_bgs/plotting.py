import numpy as np
import shapely
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.ndimage import gaussian_filter

def plot_Lag_tess(ld,val,res=[1000,1000], bw = 0, contour_levels = 0, plot_lines = True, title = None, background = False):
    '''ld - Laguerre diagram class
    val - value to assign to each cell
    res - resolution of rasterized image
    bw - bandwidth of Gaussian blur kernel
    '''    
    rast = ld.get_rasterizer(0)
            
    val = np.concatenate([val, [np.nan]])
    rv = rast.rasterize(val, res)
    
    if bw > 0:
        tmp0 = rv.copy()
        tmp0[np.isnan(tmp0)] = 0
        tmp1 = 1. - np.isnan(rv)
        rvc = gaussian_filter(tmp0, bw) / gaussian_filter(tmp1, bw)
        rvc[np.isnan(rv)] = np.nan
        rv = rvc

    fig, ax = plt.subplots(1, 1, figsize=(10,5), dpi=500, tight_layout = True)

    sp = ld.sim
    im = plt.imshow(rv.T, origin="lower", aspect="auto",
               extent=[sp.spmin[0], sp.spmax[0], sp.spmin[1], sp.spmax[1]])

    if plot_lines:
        lc = LineCollection([[p.x for p in e.ls.points] for e in ld.edglist], 
                            alpha=0.8, linewidth=0.2, colors="white")
        ax.add_collection(lc)
        
    if contour_levels > 0:
        
        cs = plt.contour(rv.T, 
                        levels=contour_levels, 
                        colors = ["w"] * (contour_levels//2 + 1) + ["k"] * contour_levels,
                        extent=[sp.spmin[0], sp.spmax[0], sp.spmin[1], sp.spmax[1]])
        plt.clabel(cs, inline=True, fontsize=10)


    plt.colorbar(im)
    ax.invert_yaxis()
    
    if background:
        fig.patch.set_facecolor('xkcd:white')
    
    if title is not None:
        plt.title(title)
    
    return rv

def plot_pressure_surface(ld,title='Pressure surface'):
    # Get points that are on the pressure surface and delete duplicates
    # By convention, points on the pressure surface have index pj greater than or equal to N+6 because
    # the first six 'edges' are the sides of a 3d cube in the lifted space
    N = ld.ys.shape[0]
    top_edges = [e for e in ld.edglist if e.pj >= N+6]
    edge_pts = np.unique(np.array([p.x for te in top_edges for p in te.ls.points]),axis=0)
    
    # Sort points by s-coordinate 
    idx = np.argsort(edge_pts[:,0])
    pressure_surface = edge_pts[idx]
    
    fig, ax = plt.subplots(1, 1, figsize=(10,5), dpi=500)
    ax.plot(pressure_surface[:,0],pressure_surface[:,1])
    ax.set_xlim([ld.sim.spmin[0],ld.sim.spmax[0]])
    ax.set_ylim([ld.sim.spmin[1],ld.sim.spmax[1]])
    ax.set_title(title)
    
def get_u(ld):

    # compute zonal wind at cell centroids
    pp = ld.phys
    y = ld.ys
    n = ld.ys.shape[0]

    polys = [ld.get_poly(i) for i in range(n)]
    centroids = [shapely.centroid(shapely.Polygon(p)) for p in polys]
    svals = np.array([p.x for p in centroids])
    cosvals = np.sqrt(1-svals**2)
    u_vals = y[:,0]/(pp.a*cosvals) - pp.Omega*pp.a*cosvals # zonal wind at centroids

    return u_vals

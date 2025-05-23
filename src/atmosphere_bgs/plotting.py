import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.ndimage import gaussian_filter

def plot_lag_tess(ld, val = None, res=[1000,1000], bw = 0, contour_levels = 0, transform = None, invert_yaxis=False, plot_lines = True, line_clr = 'white', title = None, background = False):
    '''
    ld - Laguerre diagram class
    val - value to assign to each cell
    res - resolution of rasterized image
    bw - bandwidth of Gaussian blur kernel
    '''    

    rast = None

    if val is not None:
        
        if transform is not None:
            rast = ld.get_rasterizer(transform)
        else:
            rast = ld.get_rasterizer()
          
        rv = rast.rasterize(val, res).copy()
        fill = rast.fill
        rv = np.where(fill > 0.5, np.divide(rv, fill, where=(fill > 0.5)), np.nan)
        
        if bw > 0:
            tmp0 = rv.copy()
            tmp0[np.isnan(tmp0)] = 0
            tmp1 = 1. - np.isnan(rv)
            rvc = gaussian_filter(tmp0, bw) / gaussian_filter(tmp1, bw)
            rvc[np.isnan(rv)] = np.nan
            rv = rvc

    fig, ax = plt.subplots(1, 1, figsize=(10,5), dpi=500, tight_layout = True)

    sp = ld.sim

    if val is not None:
        im = plt.imshow(rv.T, origin="lower", aspect="auto",
                   extent=rast.bounds) #[sp.spmin[0], sp.spmax[0], sp.spmin[1], sp.spmax[1]])

        if contour_levels > 0:
        
            cs = plt.contour(rv.T, 
                            levels=contour_levels, 
                            colors = ["w"] * (contour_levels//2 + 1) + ["k"] * contour_levels,
                            extent=rast.bounds)#extent=[sp.spmin[0], sp.spmax[0], sp.spmin[1], sp.spmax[1]])
            plt.clabel(cs, inline=True, fontsize=10)

        plt.colorbar(im)

    if plot_lines:
        lc = LineCollection([[p for p in e.ls.x] for e in ld.edglist], 
                            alpha=0.8, linewidth=0.2, colors=line_clr)
        ax.add_collection(lc)

    if rast is not None:
        ax.set_xlim(rast.bounds[:2])
        ax.set_ylim(rast.bounds[2:])

    if invert_yaxis:
        ax.invert_yaxis()
    
    if background:
        fig.patch.set_facecolor('xkcd:white')
    
    if title is not None:
        plt.title(title)
    
    if val is not None:
        return rv, ax
    else:
        return ax

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

def centroid(p):
    """
    Compute the centroid of a polygon.

    Parameters:
    p - numpy array of shape [n, 2] representing polygon vertices

    Returns:
    c - numpy array of shape [2] representing centroid of p
    """
    p = np.concatenate((p, p[0][None,:]), axis=0)
    x0 = p[:-1,0]
    x1 = p[1:,0]
    y0 = p[:-1,1]
    y1 = p[1:,1]
    A = 0.5 * np.sum(x0 * y1 - x1 * y0)
    C = 1. / (6 * A) * np.sum((p[1:] + p[:-1]) * (x0 * y1 - x1 * y0)[:,None], axis=0)
    return C
    
def get_u(ld):
    """
    Compute zonal wind at cell centroids.
    """

    pp = ld.phys
    y = ld.ys
    n = ld.ys.shape[0]

    polys = [ld.get_poly(i) for i in range(n)]
    centroids = [centroid(p) for p in polys]
    svals = np.array([p[0] for p in centroids])
    cosvals = np.sqrt(1-svals**2)
    u_vals = y[:,0]/(pp.a*cosvals) - pp.Omega*pp.a*cosvals # zonal wind at centroids

    return u_vals

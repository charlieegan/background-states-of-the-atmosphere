import numpy as np
import matplotlib.pyplot as plt

from scipy import optimize
from scipy.special import softmax
from scipy import sparse
from scipy.interpolate import PchipInterpolator

class PostProcessor:
    def __init__(self,solv,bscirc,pvlevs,thlevs,pvmaxth,\
                 res=[200,200],eps=5e-3,nb_threshold=None,nbs=None,nb_type='Lag_cell_8',n_nb_th_levs=None,\
                 make_surface_adjustments=False,adjust_weights=False,verbose=False):
        
        self.ld = solv.ld
        self.pp = solv.pp
        self.sp = solv.sp
        
        self.eps = eps
        self.nb_type = nb_type
        if nb_threshold is not None:
            self.nb_threshold = nb_threshold
        self.res = res
        
        self.g = 9.81 # these should be in pp
        self.H = 6500 # these should be in pp
        
        # Adjust the surface zonal angular momentum to get rid of spurious
        # oscillations in the surface zonal wind
        if make_surface_adjustments:
            if verbose:
                print('Adjusting surface Z')
            Z_adj_surf, i_surf, s_mids = self.adjust_Z_surf(d=2,verbose=verbose)
            Z_adj = self.ld.ys[:,0].copy()
            Z_adj[i_surf] = Z_adj_surf
            self.surf_s_mids = s_mids
            self.i_surf = i_surf
            self.Z_adj = Z_adj
            
            # Adjust the surface potential temperature to get rid of steps
            if verbose:
                print('Adjusting surface theta')
            th_adj_surf, i_surf, s_mids = self.adjust_th_surf(d=2,verbose=verbose)
            th_adj = self.ld.ys[:,0].copy()
            th_adj[i_surf] = th_adj_surf
            self.surf_s_mids = s_mids
            self.i_surf = i_surf
            self.th_adj = th_adj
            
        # Smooth geopotential and obtain Z and theta from balance relations
        if verbose:
            print('Obtaining smooth fields')
        self.apply_smoothing(eps=eps,nb_threshold=nb_threshold,nbs=nbs,nb_type=nb_type,n_nb_th_levs=n_nb_th_levs,adjust_weights=adjust_weights,verbose=verbose);
        
        # Define smooth u from smooth Z
        self.u_eps = self.get_u(self.Z_eps,self.sg)
        
        # Get smooth surface pressure
        if verbose:
            print('Obtaining smooth surface pressure')
        self.get_smooth_lower_surface_pressure()
        
        # Get smooth surface variables
        if verbose:
            print('Obtaining smooth lower surface variables')
        self.get_smooth_lower_surface_vars(eps=eps)
        
        if make_surface_adjustments:
            self.interpolate_adjusted_lower_surface_vars()
        
        # Push to isentropic coordinates
        self.bscirc = bscirc
        self.pvlevs = pvlevs
        self.thlevs = thlevs
        self.pvmaxth = pvmaxth
        self.push_forward_to_insentropic_coords()
        self.get_q_isentropic()
        
        
    def get_u(self,Z,s):
        '''Compute zonal wind for zonal angular momentum Z and sin-of-lat s'''
        Omega = self.pp.Omega
        a = self.pp.a
        return Z/a/np.sqrt(1-s**2) - Omega*a*np.sqrt(1-s**2)
    
    def adjust_Z_surf(self,d=2,verbose=False):
        '''
        Function adjusts Laguerre-cell-based surface zonal angular momentum to
        remove spurious oscillations in surface wind. The (approximate) L2-norm of
        the d-th derivative of zonal wind is minimised subject to the corresponding
        surface Z values being constrained between midpoints of adjacent initial
        surface Z values.

        Parameters
        ----------
        solv : solution of the OT problem
        verbose : boolean
        plot : boolean

        Returns
        -------
        Z_adj : numpy array of adjusted surface Z values
        i_surf : integer numpy array of surface cell indices
        s_mids : midpoints in s of boundary-cell intersections
        '''
        # Get left endpoints of cell-surface intersections
        ld = self.ld
        n = len(ld.ys)
        s_end_pts_all = np.vstack([np.min(e.ls.x[:,0]) for e in ld.edglist if e.pj>n+6])
        i_surf = np.hstack([e.pi-6 for e in ld.edglist if e.pj>n+6])
        s_end_pts = np.hstack([np.min(s_end_pts_all[np.where(i_surf==i)[0]]) for i in np.unique(i_surf)])
        i_idx_sort = np.argsort(s_end_pts)
        i_surf = np.unique(i_surf)[i_idx_sort]
        s_l_end_pts = s_end_pts[i_idx_sort]
        s_r_end_pts = np.append(s_l_end_pts[1:],ld.sim.smax)
        s_mids = (s_l_end_pts + s_r_end_pts)/2
        
        # Get surface Z values and define bounds on u using modpoints of neighbouring surface values
        Z = ld.ys[:,0]
        Z_surf = Z[i_surf]

        Z_l_nb = np.hstack([np.roll(Z_surf,shift=-1)[:-1],np.min(Z)])
        Z_u_nb = np.hstack([np.max(Z),np.roll(Z_surf,shift=1)[1:]])

        Z_lb = .5*(Z_surf + Z_l_nb)
        Z_ub = .5*(Z_surf + Z_u_nb)
        
        u_lb = self.get_u(Z_lb,s_mids)
        u_ub = self.get_u(Z_ub,s_mids)
        
        bds = optimize.Bounds(lb=u_lb,ub=u_ub)
            
        # Define loss as the integral of third derivative of u squared
        K = len(s_mids)
        if d == 1:
            D = -np.eye(K-1,K,k=0) + np.eye(K-1,K,k=1)
        elif d == 2:
            D = np.eye(K-2,K,k=0) - 2*np.eye(K-2,K,k=1) + np.eye(K-2,K,k=2)
        elif d == 3:
            D = -np.eye(K-3,K,k=0) + 3*np.eye(K-3,K,k=1) - 3*np.eye(K-3,K,k=2) + np.eye(K-3,K,k=3)
        Q = D.T @ D
        
        loss = lambda u_adj : 0.5 * u_adj @ Q @ u_adj
        jac = lambda u_adj : Q @ u_adj
        hess = lambda u_adj : Q
        
        # Initialise using 10-point rolling average u adjusted to satisfy the bounds
        u_orig = self.get_u(Z_surf,s_mids)
        max_shift = 5
        u_init = np.mean([np.roll(u_orig,shift=shift) for shift in np.arange(-max_shift,max_shift+1)],axis=0)
        u_init[:max_shift] = u_orig[:max_shift]
        u_init[-max_shift:] = u_orig[-max_shift:]
        
        i_ub_ex = np.where(u_init>u_ub)[0]
        i_lb_ex = np.where(u_init<u_lb)[0]
        u_init[i_ub_ex] = u_ub[i_ub_ex]
        u_init[i_lb_ex] = u_lb[i_lb_ex]
        
        tol = 1e-3
        if verbose:
            print('Optimising for surface zonal wind')
        opt_res = optimize.minimize(loss,x0=u_init.copy(),bounds=bds,jac=jac,hess=hess,method='trust-constr',tol=tol)
        
        if verbose:
            print(f"Adjustment of surface zonal angular momentum successful: {opt_res.success}")
        
        # Define adjusted surface Z
        u_adj = opt_res.x
        Omega = self.pp.Omega
        a = self.pp.a
        Z_adj = Omega*a**2*(1-s_mids**2) + u_adj*a*np.sqrt(1-s_mids**2)
        
        # fig, ax = plt.subplots(1,1,dpi=200)
        # ax.plot(s_mids,u_orig,lw=1,label='Original u',alpha=0.4)
        # ax.plot(s_mids,u_adj,label='Adjusted u',lw=1,alpha=0.9)
        # ax.plot(s_mids,u_ub,label='Max u',lw=1,alpha=0.4,ls='--',c='k')
        # ax.plot(s_mids,u_lb,label='Min u',lw=1,alpha=0.4,ls='--',c='k')
        # ax.set_ylabel('$u_s$')
        # ax.set_xlabel('$s$')
        # ax.legend()    
        # ax.set_title('Minimising 2nd derivative')
                
        self.surf_adjustments_u = {'orig' : u_orig,
                                   'adj' : u_adj,
                                   'lb' : u_lb,
                                   'ub' : u_ub}
            
        return Z_adj, i_surf, s_mids
    
    def adjust_th_surf(self,d=1,verbose=False):
        '''
        Function adjusts Laguerre-cell-based surface theta to
        remove big steps. The (approximate) L2-norm of
        the d-th derivative of surface theta is minimised subject to each theta 
        value being constrained between neighbouring theta levels
        
        Parameters
        ----------
        solv : solution of the OT problem
        verbose : boolean
        plot : boolean

        Returns
        -------
        th_adj : numpy array of adjusted surface theta values
        i_surf : integer numpy array of surface cell indices
        s_mids : midpoints in s of boundary-cell intersections
        '''
        # Get left endpoints of cell-surface intersections
        ld = self.ld
        n = len(ld.ys)
        s_end_pts_all = np.vstack([np.min(e.ls.x[:,0]) for e in ld.edglist if e.pj>n+6])
        i_surf = np.hstack([e.pi-6 for e in ld.edglist if e.pj>n+6])
        s_end_pts = np.hstack([np.min(s_end_pts_all[np.where(i_surf==i)[0]]) for i in np.unique(i_surf)])
        i_idx_sort = np.argsort(s_end_pts)
        i_surf = np.unique(i_surf)[i_idx_sort]
        s_l_end_pts = s_end_pts[i_idx_sort]
        s_r_end_pts = np.append(s_l_end_pts[1:],ld.sim.smax)
        s_mids = (s_l_end_pts + s_r_end_pts)/2
        
        # Get surface theta values and define bounds on theta using midpoints to neighbouring theta levels
        th = ld.ys[:,1]
        th_surf_orig = th[i_surf]
        
        th_levs = np.unique(th)
        th_diffs = np.diff(th_levs)
        th_diffs = np.hstack([th_diffs[0],th_diffs,th_diffs[-1]])
        idx = np.where(th_surf_orig[:,None]==th_levs)[1]
        th_ub = th_surf_orig + th_diffs[idx+1]/2
        th_lb = th_surf_orig - th_diffs[idx]/2
        
        bds = optimize.Bounds(lb=th_lb,ub=th_ub)
            
        # Define loss as the integral of third derivative of u squared
        K = len(s_mids)
        if d == 1:
            D = -np.eye(K-1,K,k=0) + np.eye(K-1,K,k=1)
        elif d == 2:
            D = np.eye(K-2,K,k=0) - 2*np.eye(K-2,K,k=1) + np.eye(K-2,K,k=2)
        elif d == 3:
            D = -np.eye(K-3,K,k=0) + 3*np.eye(K-3,K,k=1) - 3*np.eye(K-3,K,k=2) + np.eye(K-3,K,k=3)
        Q = D.T @ D
        
        loss = lambda th_adj : 0.5 * th_adj @ Q @ th_adj
        jac = lambda th_adj : Q @ th_adj
        hess = lambda th_adj : Q
        
        # Initialise using 10-point rolling average u adjusted to satisfy the bounds
        max_shift = 5
        th_init = np.mean([np.roll(th_surf_orig,shift=shift) for shift in np.arange(-max_shift,max_shift+1)],axis=0)
        th_init[:max_shift] = th_surf_orig[:max_shift]
        th_init[-max_shift:] = th_surf_orig[-max_shift:]
        
        i_ub_ex = np.where(th_init>th_ub)[0]
        i_lb_ex = np.where(th_init<th_lb)[0]
        th_init[i_ub_ex] = th_ub[i_ub_ex]
        th_init[i_lb_ex] = th_lb[i_lb_ex]
        
        tol = 1e-3
        opt_res = optimize.minimize(loss,x0=th_init.copy(),bounds=bds,jac=jac,hess=hess,method='trust-constr',tol=tol)
        
        if verbose:
            print(f"Adjustment of surface theta successful: {opt_res.success}")
        
        # Define adjusted surface Z
        th_adj = opt_res.x
        
        # fig, ax = plt.subplots(1,1,dpi=200)
        # ax.plot(s_mids,th_surf_orig,lw=1,label='Original $\\theta$',alpha=0.4)
        # ax.plot(s_mids,th_adj,label='Adjusted $\\theta$',lw=1,alpha=0.9)
        # ax.plot(s_mids,th_ub,label='Max $\\theta$',lw=1,alpha=0.4,ls='--',c='k')
        # ax.plot(s_mids,th_lb,label='Min $\\theta$',lw=1,alpha=0.4,ls='--',c='k')
        # ax.set_ylabel('$\\theta_s$')
        # ax.set_xlabel('$s$')
        # ax.legend()    
        # ax.set_title('Minimising 2nd derivative')
        
        self.surf_adjustments_th = {'orig' : th_surf_orig,
                                    'adj' : th_adj,
                                    'lb' : th_lb,
                                    'ub' : th_ub}
        
        return th_adj, i_surf, s_mids
    
    def get_cost(self,s,p,Z,th):
        '''Evaluate cost function at ((s,p),(Z,th))'''
        cosines = np.sqrt(1-s**2)
        pp = self.pp
        return .5*(Z[None,:]/pp.a/cosines[:,None] - pp.Omega*pp.a*cosines[:,None])**2 + pp.cp*th[None,:]*(p[:,None]/pp.p00)**pp.kappa
    
    def centroid(self,p):
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

    def get_smoothing_neighbourhoods(self,threshold=1e-5,eps=1e-3):
        '''Define neighbourhoods for which to compute the SoftMax weights for each Laguerre cell
        For each Laguerre cell index i, the geopotential kernels 
        
            k_ij = (psi_j - c((s_i,p_i),(Z_j,th_j)))/eps/g/H
        
        are computed for each Laguerre cell index j. The i-th neighbourhood is defined to
        contain all indices j for which
        
            exp(k_ij)/exp(k_ii) > threshold.
        
        In other words, index j then has a relatively large contribution to the SoftMax
        at the centroid of cell i.
        '''
        # Get centroids of Laguerre cells in (s,p)
        ld = self.ld
        lag_polys = [ld.get_poly(i) for i in np.arange(len(ld.ys))]
        centroids = np.array([self.centroid(p) for p in lag_polys])
        
        s = centroids[:,0]
        p = centroids[:,1]
        
        # Evaluate all geopotential kernels at all centroids
        Z = ld.ys[:,0]
        th = ld.ys[:,1]
        
        psi = ld.duals[:-1]
        
        costs = self.get_cost(np.ravel(s),np.ravel(p),Z,th)
        kernels = (psi - costs)/eps/self.g/self.H
        
        # For each index i, find j for which the ratio exp(kernels[i,i])/exp(kernels[i,j])
        # is greater than the given threshold
        nbs = [np.where(kernels[i]-kernels[i,i] > np.log(threshold))[0] for i in np.arange(len(ld.ys))]
        
        return nbs

    def get_Lag_cell_nbs(self):
        '''Given a Laguerre diagram, this returns a list of Laguerre cell neighourhoods,
        where the neighbourhood of cell i is the set of indices j distinct from i
        such that the intersection of the i-th and j-th cells is non-empty'''
        ld = self.ld
        
        # Set up list to store Laguerre cell neihbourhoods
        imax = np.max([e.pj for e in ld.edglist])
        nbs_Lag = [[] for _ in np.arange(imax+1)]
        
        # For each edge, record the adjacent cells as being neighbours
        for e in ld.edglist:
            nbs_Lag[e.pi].append(e.pj-6) # there are six boundary indices to be ignored
            nbs_Lag[e.pj].append(e.pi-6)
        
        # Ignore the six boundary indices
        n = len(self.ld.ys)
        nbs_Lag = [np.array(nb) for nb in nbs_Lag[6:(n+6)]]
        nbs_Lag = [nb[nb>=0] for nb in nbs_Lag]
        
        return nbs_Lag
    
    def get_neighbours_up_to_degree(self,nbs, d):
        """
        Finds all neighbours of degree d and lower using powers of adjacency matrix.

        Parameters:
        nbs (list of lists): Adjacency list where nbs[i] contains first-degree neighbours of i.
        d (int): Maximum degree of separation.

        Returns:
        list of lists: nbs_d, where nbs_d[i] contains neighbours of degree <= d.
        """
        n = len(nbs)

        coords = np.hstack([np.array([i*np.ones(len(nb)),nb]) for i, nb in enumerate(nbs)])
        A = sparse.coo_array((np.ones(coords.shape[1],dtype=int),coords),[n,n]).tocsr()

        I = sparse.eye(n,dtype=int,format='csr')
        A_weighted = A + I
        
        d_count = 1
        A_d = A_weighted
        while d_count<d: #  could use sparse.linalg.matrix_power for this
            d_count = d_count+1
            A_d = A_d@A_weighted

        nbs_d = []
        for i in np.arange(n):
            row_start = A_d.indptr[i]
            row_end = A_d.indptr[i+1]
            row_neighbours = A_d.indices[row_start:row_end]

            nbs_d.append(row_neighbours)

        return nbs_d
    
    
    # # For checking against finite difference
    # i0 = sr==sr[0]
    # i1 = np.argsort(pr[i0])
    # prr=pr[i0][i1]
    # thr = th_eps[i0][i1]
    # dth_dp_fd = np.diff(thr)/np.diff(prr)
    # dth_dp_r = dth_dp[i0][i1]
    
    # smr = sm[i0][i1]
    # dsm_dp_fd = (np.diff(smr,axis=0).T/np.diff(prr)).T
    # dsm_dp_r = dsm_dp[i0][i1][:-1]
    
    # print((dsm_dp_fd[0]-dsm_dp_r[0])/dsm_dp_fd[0])
    
    # i0 = pr==pr[2]
    # i1 = np.argsort(pr[i0])
    # srr=sr[i0][i1]
    # Zr = Z_eps[i0][i1]**2
    # dZ_ds_fd = np.diff(Zr)/np.diff(srr)
    # dZ_ds_r = dZ2_ds[i0][i1]

    def apply_smoothing_local(self,s,p,isub=None,eps=5e-3,adjust_weights=False):
        '''Function to apply smoothing on a given set of points (s,p) using the
        Laguerre cell neighbourhood isub to compute softmax with smoothing parameter
        eps. Returns evaluation of analytic formulas for smooth variables at (s,p).
        These variables are: 
            
            Phi_eps  - geopotential
            Z_eps    - zonal angular momentum
            th_eps   - potential temperature
            r_eps    - isentropic density, r = -(1/g)/(dtheta/dp)
            
        The smooth zonal angular momentum and potential temperature are defined to
        satisfy gradient wind and hydrostatic balance with the smooth geopotential.
            
        If adjust_weights is set to True then adjusted values of surface Z are used
        to compute softmax weights. Otherwise, the original values of Z are used.
        '''
        # Extract variables and parameters from the Laguerre diagram
        ld = self.ld
        if adjust_weights and hasattr(self,'Z_adj'):
            Z = self.Z_adj
        else:
            Z = self.ld.ys[:,0]
        th = self.ld.ys[:,1]
        Omega = self.pp.Omega
        a = self.pp.a
        kappa = self.pp.kappa
        cp = self.pp.cp
        p00 = self.pp.p00
        g = self.g
        H = self.H
        
        # Due to smoothing Z^2 could be negative. Define minimum positive value
        # of Z to be planetary angular momentum at midpoint between smax and 1
        # and use wherever this is the case
        Z_min = Omega*a**2*(1-.5*(1+self.sp.smax**2))
        
        # Default to using all target points to compute softmax
        if isub is None:
            isub = np.arange(len(th))
            
        # Extract relevant elements from Z, theta and duals 
        psi = ld.duals[isub]
        Z = Z[isub]
        th = th[isub]
        
        # Compute kernel (psi_i - c((s,p),(Z_i,th_i)))/eps/g/H for each cell index
        # i at all (s,p)
        costs = self.get_cost(s, p, Z, th)#np.random.rand(npts,nnbs)
        ker = (psi - costs)/eps/g/H
    
        # compute the softmax of the kernels at each grid point
        sm = softmax(ker,axis=1)
    
        Z_eval = Z
        ker_eval = ker
        Phi_eps = eps*g*H*np.sum(ker_eval*sm,axis=1)
    
        # Ravel the s and p coordinate arrays
        sr = np.ravel(s)
        pr = np.ravel(p)
    
        # Define kernel derivatives
        #dk_ds = -(1/eps/g/H)*(Z[None,:]**2/a**2/(1-sr[:,None]**2)**2 - Omega**2*a**2)*sr[:,None]
        #dkeval_ds = -(1/eps/g/H)*(Z_eval[None,:]**2/a**2/(1-sr[:,None]**2)**2 - Omega**2*a**2)*sr[:,None]
        dk_dp = -(1/eps/g/H)*(kappa*cp/p00**kappa)*pr[:,None]**(kappa-1)*th[None,:]
    
        # Define softmax derivatives
        #dsm_ds = dk_ds*sm - (sm.T*np.sum(dk_ds*sm,axis=1)).T
        dsm_dp = dk_dp*sm - (sm.T*np.sum(dk_dp*sm,axis=1)).T
    
        # Define smooth theta
        fth0 = sm@th
        fth1 = 1-np.sum(ker_eval*sm,axis=1)
        th_eps = fth0*fth1 + np.sum(sm*ker_eval*th,axis=1)
    
        # Define smooth dth/dp
        dth_dp = fth1*(dsm_dp@th) - fth0*(np.sum(dk_dp*sm,axis=1) + np.sum(ker_eval*dsm_dp,axis=1)) \
                + np.sum(dsm_dp*ker_eval*th,axis=1) + np.sum(sm*dk_dp*th,axis=1)
                
        # Define smooth dth/ds (used to compute zeta)
        # dth_ds = fth1*(dsm_ds@th) - fth0*(np.sum(dk_ds*sm,axis=1) + np.sum(ker_eval*dsm_ds,axis=1)) \
        #         + np.sum(dsm_ds*ker_eval*th,axis=1) + np.sum(sm*dk_ds*th,axis=1)
    
        # Define smooth isentropic density
        r_eps = -(1/g)/dth_dp
    
        # Define smooth Z
        Z2_eval = Z_eval**2
        Z2 = Z**2
        Zp2 = Omega**2*a**4*(1-sr**2)**2
    
        fZ0 = sm@Z2
        fZ1 = np.sum(ker_eval*sm,axis=1)
        fZ2 = np.sum(sm.T*Zp2,axis=0)
        
        Z2_eps = sm@Z2_eval + (ker_eval*sm)@Z2 \
                - np.sum((ker_eval*sm).T*Zp2,axis=0) + fZ1*(fZ2 - fZ0)
    
        Z_eps = np.max([np.sqrt(Z2_eps),Z_min*np.ones(len(Z2_eps))],axis=0)
    
        # Define smooth dZ/ds. First get d(Z^2)/ds and use identity dZ/ds = .5*(d(Z^2)/ds)/Z
        # fZ3 = np.sum(dkeval_ds*sm,axis=1)
        # fZ4 = np.sum(ker_eval*dsm_ds,axis=1)
        # dZp2_ds = -4*Omega**2*a**4*(1-sr**2)*sr
    
        # dZ2_ds = dsm_ds@Z2_eval \
        #         + np.sum(dkeval_ds*sm*Z2,axis=1) + np.sum(ker_eval*dsm_ds*Z2,axis=1) \
        #         - np.sum((ker_eval*sm).T*dZp2_ds,axis=0) - np.sum((dkeval_ds*sm).T*Zp2,axis=0)- np.sum((ker_eval*dsm_ds).T*Zp2,axis=0) \
        #         - fZ1*(dsm_ds@Z2) - fZ0*(fZ3 + fZ4) \
        #         + fZ1*(np.sum(sm.T*dZp2_ds,axis=0) + np.sum(dsm_ds.T*Zp2,axis=0)) + fZ2*(fZ3 + fZ4)
    
        # dZ_ds = .5*dZ2_ds/Z_eps
        
        # The above optimised version of the code can be checked against the following simpler
        # but slower implementation using numpy.einsum
        #
        ## Compute theta
        # wth = th*(1+ker_eval[:,None,:]-ker_eval[:,:,None])
        # th_eps_einsum = np.einsum(wth,[...,1,2],sm,[...,1],sm,[...,2])
        #
        ## Compute r
        # dk_dp = -(1/eps/g/H)*(kappa*cp/p00**kappa)*pr[:,None]**(kappa-1)*th[None,:]
        # dSMSM_dp = dk_dp[:,:,None,None] + dk_dp[:,None,:,None] - 2*dk_dp[:,None,None,:] # terms in derivative of SoftMax*SoftMax
        # dth_dp_einsum = np.einsum(th,[1],dk_dp[:,:,None]-dk_dp[:,None,:],[...,1,2],sm,[...,1],sm,[...,2]) + \
        #             np.einsum(wth,[...,2,1],sm,[...,1],sm,[...,2],sm,[...,3],dSMSM_dp,[...,1,2,3])
        # r_eps_einsum = -(1/g)/dth_dp_einsum
        #
        ## Compute Z
        # wZ = Z_eval[None,None,:]**2 \
        #     + (Z[None,None,:]**2 - Omega**2*a**4*(1-sr[:,None,None]**2)**2)*(ker_eval[:,None,:]-ker_eval[:,:,None])
        # Z2_eps_einsum = np.einsum(wZ,[...,1,2],sm,[...,1],sm,[...,2])
        # Z_eps_einsum = np.max([np.sqrt(Z2_eps_einsum),Z_min*np.ones(len(Z2_eps_einsum))],axis=0)
        #
        ## Compute zeta
        # dkeval_ds = -(1/eps/g/H)*(Z_eval[None,:]**2/a**2/(1-sr[:,None]**2)**2 - Omega**2*a**2)*sr[:,None]
        # dk_ds = -(1/eps/g/H)*(Z[None,:]**2/a**2/(1-sr[:,None]**2)**2 - Omega**2*a**2)*sr[:,None]
        # dSMSM_ds = dk_ds[:,:,None,None] + dk_ds[:,None,:,None] - 2*dk_ds[:,None,None,:]
        # dZ2_ds_einsum = np.einsum(wZ,[...,2,1],sm,[...,1],sm,[...,2],sm,[...,3],dSMSM_ds,[...,1,2,3]) \
        #           - 4*Omega**2*a**4*(sr**3-sr)*np.einsum(dkeval_ds[:,:,None]-dkeval_ds[:,None,:],[...,1,2],sm,[...,1],sm,[...,2]) \
        #           + np.einsum(Z[None,:]**2 - Omega**2*a**4*(1-sr[:,None]**2)**2,[...,1],dkeval_ds[:,:,None]-dkeval_ds[:,None,:],[...,1,2],sm,[...,1],sm,[...,2])
        # dZ_ds_einsum = .5*dZ2_ds_einsum/Z_eps_einsum
        # zeta_eps_einsum = -dZ_ds_einsum/a**2
        #
        ## Compute q
        # q_eps_einsum = zeta_eps_einsum/r_eps_einsum
        
        if len(r_eps)==1:
            return Phi_eps[0], Z_eps[0], th_eps[0], r_eps[0]
        
        return Phi_eps, Z_eps, th_eps, r_eps


    def apply_smoothing(self,eps=5e-3,nb_threshold=1e-5,nbs=None,nb_type=None,n_nb_th_levs=None,adjust_weights=True,verbose=False):

        ld = self.ld
        th = ld.ys[:,1]
        pp = self.pp
        res = self.res
        n = len(th)

        # Get smoothing neighbourhoods
        if verbose:
            print('Getting smoothing neighbourhoods')
        if nbs is not None:
            nb_type = 'fixed, specified'
            nbs = n*[nbs]
        elif nb_type == 'threshold':
            # Threshold by the kernel values at the Laguerre cell centroids
            nbs = self.get_smoothing_neighbourhoods(ld,threshold=nb_threshold,eps=eps)
        elif 'Lag_cell' in nb_type:
            # Use Laguerre cells neighbourhoods with specifed degrees of separation
            nbs_Lag = self.get_Lag_cell_nbs()
            nbs = [nb[nb<n] for nb in nbs_Lag]
            nbs = self.get_neighbours_up_to_degree(nbs, float(nb_type[-1]))
            
            # keep only n_nb_th_levs-many theta levels either side of the given cell if specified 
            if n_nb_th_levs is not None:
                nbs_anisotropic = []
                thlevs = np.unique(th)

                for i, nb in enumerate(nbs):
                    i_lev = np.where(thlevs==th[i])[0][0]
                    thmin = thlevs[np.max([i_lev-n_nb_th_levs,0])]
                    thmax = thlevs[np.min([i_lev+n_nb_th_levs,len(thlevs)-1])]
                    idx = np.where((th[nb]>=thmin)*(th[nb]<=thmax))[0]
                    nbs_anisotropic.append(nb[idx])
                
                nbs = nbs_anisotropic
            
        # Get indices of Laguerre cells that are on the lower boundary (these
        # are the cells whose Laguerre cell neighbourhood contains an index bigger 
        # than the number of Laguerre cells)
        i_lower = np.where([np.max(nb)>=n for nb in nbs_Lag])[0]
        
        # Define bounds of pressure such that all cells sit within them, excluding 
        # fictitious cells and redefine limits of ld
        if verbose:
            print('Setting up grid for smoothing')
            
        i_genuine = np.where(th<np.max(th))[0]
        verts = np.vstack([ld.get_poly(i) for i in i_genuine])
        pmin = 0.9*np.min(verts[:,1])
        pmax = 1.1*np.max(verts[:,1])
        pmin_orig = ld.sim.pmin
        pmax_orig = ld.sim.pmax
        ld.sim.pmin = pmin
        ld.sim.pmax = pmax

        # Define regular rectangular grid in (lat,h) and find which Laguerre cell
        # each grid point belongs to
        s2lat = lambda s : np.rad2deg(np.arcsin(s))
        p2h = lambda p : -self.H*np.log(p/pp.p00)

        lat2s = lambda lat : np.sin(np.deg2rad(lat))
        h2p = lambda h : pp.p00*np.exp(-h/self.H)

        smin = ld.sim.smin
        smax = ld.sim.smax
        
        lg = np.linspace(s2lat(smin),s2lat(smax),res[0]+1)
        hg = np.linspace(p2h(pmin),p2h(pmax),res[1]+1)
        lg = lg[:-1] + np.diff(lg)/2 # midpoints
        hg = hg[:-1] + np.diff(hg)/2 # midpoints
        sg = lat2s(lg)
        pg = h2p(hg)
        sgu = sg.copy()
        pgu = pg.copy()
        
        lg, hg = np.meshgrid(lg,hg)
        sg, pg = np.meshgrid(sg,pg)
        sg = np.ravel(sg)
        pg = np.ravel(pg)
        
        transform = lambda x : [s2lat(x[0]), p2h(x[1])]
        rast = ld.get_rasterizer(transform)
        pts = np.vstack([np.ravel(lg),np.ravel(hg)]).T
        cell_idx = rast.indices(pts)

        ld.sim.pmin = pmin_orig
        ld.sim.pmax = pmax_orig
        
        # Set up matrices to store smooth fields
        Phi_eps = np.nan*np.zeros(res[0]*res[1])
        Z_eps = np.nan*np.zeros(res[0]*res[1])
        th_eps = np.nan*np.zeros(res[0]*res[1])
        r_eps = np.nan*np.zeros(res[0]*res[1])
        zeta_eps = np.nan*np.zeros(res[0]*res[1])
        q_eps = np.nan*np.zeros(res[0]*res[1])

        # Define smooth geopotential using softmax of cost minus Kantorovich 
        # potential over cell neighbour indices, and then define smooth zonal
        # angular momentum and potential temperature substituting exact derivatives
        # of this smooth geopotential into balance equations        
        if verbose:
            print('Smoothing geopotential')
            
        # Set up list for recording Laguerre cell index as a function of latitude
        # on the lower surface
        s_lower = []
            
        for i in i_genuine: # exclude fictitious cells
                
            # Get the points that are in the ith cell
            ig = np.where(cell_idx==i)[0]
            si = sg[ig]
            pi = pg[ig]
            
            # If a cell on the lower boundary, extend p to pmax on each selected s-level.
            # This is so that the smooth fields are defined beyond the surface pressure
            # that is given by the Laguerre cells and then a smooth surface pressure
            # can be defined by finding the 0-level set of the smooth geopotential
            if i in i_lower:
                su = np.unique(si)
                igextra = np.ones(0,dtype=np.int64)
                poly = ld.get_poly(i)
                smin, smax = np.min(poly[:,0]), np.max(poly[:,0])
                pmin, pmax = np.min(poly[:,1]), np.max(poly[:,1])
                
                if len(su)>0:
                    for sii in su:
                        imax = np.max(ig[np.where(si==sii)[0]])
                        jmax = np.max(np.where(sg==sii)[0])
                        igextra = np.concatenate([igextra,np.arange(imax+res[1],jmax+1,res[1],dtype=np.int64)])
                        
                    ig = np.concatenate([ig,igextra])
                    si = np.concatenate([si,sg[igextra]])
                    pi = np.concatenate([pi,pg[igextra]])
                    
                    # It can be that a boundary cell contains some points but 
                    # is so shallow at the ends that some columns are not 
                    # extended beyond the boundary. This next block of code deals
                    # with this case
                    j=np.where(np.isin(pgu,pg[igextra]))
                    
                    igleft = np.where((smin<sgu)*(sgu<np.min(si)))[0]
                    igright = np.where((np.max(si)<sgu)*(sgu<smax))[0]
                    
                    if any(igleft):
                        i_s, i_p = np.meshgrid(igleft,j)
                        igextra = np.ravel(i_p*200+i_s)
                        ig = np.concatenate([ig,igextra])
                        si = np.concatenate([si,sg[igextra]])
                        pi = np.concatenate([pi,pg[igextra]])
                    if any(igright):
                        i_s, i_p = np.meshgrid(igright,j)
                        igextra = np.ravel(i_p*200+i_s)
                        ig = np.concatenate([ig,igextra])
                        si = np.concatenate([si,sg[igextra]])
                        pi = np.concatenate([pi,pg[igextra]])    
                
                    # Record Laguerre cell index as a function of latitude
                    # on the lower surface 
                    s_lower.append([np.unique(si),i,nbs[i]])
                    
                # if no grid point is in the cell, find latitude lines that intersect the cell
                # and use cell nbhd to define smooth field below lower boundary
                else:
                    poly = ld.get_poly(i)
                    smin, smax = np.min(poly[:,0]), np.max(poly[:,0])
                    pmin, pmax = np.min(poly[:,1]), np.max(poly[:,1])
                    ig = np.where((sg>smin)*(sg<smax)*pg>pmax)[0]
                    si = sg[ig]
                    pi = pg[ig]
                    
                    # Record Laguerre cell index as a function of latitude
                    # on the lower surface 
                    s_lower.append([np.unique(si),i,nbs[i]])
                
            # compute the softmax using only the neighbour cells
            if len(si)>0:
                Phi_eps_i, Z_eps_i, th_eps_i, r_eps_i = self.apply_smoothing_local(si,pi,isub=nbs[i],eps=eps,adjust_weights=adjust_weights)
                
                # assign smooth values to matrices
                Phi_eps[ig] = Phi_eps_i
                Z_eps[ig] = Z_eps_i
                th_eps[ig] = th_eps_i
                r_eps[ig] = r_eps_i
                
        # reshape the smoothed fields and translate geopotential
        Phi_eps = np.reshape(Phi_eps,res) - ld.duals[-1]
        Z_eps = np.reshape(Z_eps,res)
        th_eps = np.reshape(th_eps,res)
        r_eps = np.reshape(r_eps,res)
            
        # reshape the grid point coordinate vectors
        sg = np.reshape(sg,res)
        pg = np.reshape(pg,res)
        lg = np.reshape(lg,res)
        hg = np.reshape(hg,res)
        
        # Force Z_eps to decrease with latitude and th_eps to decrease with p
        i_Z = np.where(np.isnan(Z_eps))
        Z_eps[i_Z] = -np.inf
        Z_eps = np.maximum.accumulate(Z_eps[:,::-1],axis=1)[:,::-1]
        Z_eps[i_Z] = np.nan
        i_th = np.where(np.isnan(th_eps))
        th_eps[i_th] = -np.inf
        th_eps = np.maximum.accumulate(th_eps[::-1],axis=0)[::-1]
        th_eps[i_th] = np.nan
        
        # assign variables to the class
        self.Phi_eps = Phi_eps
        self.Z_eps = Z_eps
        self.th_eps = th_eps
        self.r_eps = r_eps
        self.sg = sg
        self.pg = pg
        self.lg = lg
        self.hg = hg
        self.s_lower = s_lower
        
        return Phi_eps, Z_eps, th_eps, r_eps, sg, pg, lg, hg, s_lower
            
    
    def get_smooth_lower_surface_pressure(self):
        '''
        Get smooth pressure on the lower boundary. First, translate the smoothed 
        geopotenital by the final Kantorovich dual. Then define the smooth lower
        pressure surface as the zero level set of the translated smooth geopotential.
        At each latitude, the zero of geopotential is found by monotone interpolation of
        pressure as a function of geopotential.
        '''
        pg = self.pg
        Phi_eps = self.Phi_eps
        
        # Translate geopotential by final Kantorovich dual so that zero corresponds
        # to ground level
        n_lats = pg.shape[1]
        p_lower_eps = np.nan*np.zeros(n_lats)
        surf_mask = np.ones(pg.shape)
        
        # For each latitude level, find where smooth geopotential is either side of
        # zero and interpolate using isothermal relation to get pressure at zero geopotential
        pu = np.unique(pg)

        for i in np.arange(n_lats):
            
            i0 = np.argmax(Phi_eps[:,i]<0)
            i1 = i0 - 1
            
            p0 = pu[i0]
            p1 = pu[i1]
            
            phi0 = Phi_eps[i0,i]
            phi1 = Phi_eps[i1,i]
            
            Hi = (phi0-phi1)/self.g/np.log(p1/p0)
            p_surf_i = (p0 + p1)/(np.exp(-phi0/self.g/Hi) + np.exp(-phi1/self.g/Hi))
            p_lower_eps[i] = p_surf_i
            
            # define mask to mask out grid cells where smooth geopotential is negative
            surf_mask[i0:,i] = np.nan
            
        self.p_lower_eps = p_lower_eps
        self.surf_mask = surf_mask
        
        return p_lower_eps, surf_mask

    def get_smooth_lower_surface_vars(self,eps=5e-3):
        '''
        Get smooth geopotential, zonal angular momentum and potential temperature
        on the lower boudary by evaluating smooth fields at surface pressure
        '''
        
        sg = self.sg
        s_lower = self.s_lower
        p_lower_eps = self.p_lower_eps
        ld = self.ld
    
        n_lats = sg.shape[1]
        Phi_lower_eps = np.nan*np.zeros(n_lats)
        th_lower_eps = np.nan*np.zeros(n_lats)
        Z_lower_eps = np.nan*np.zeros(n_lats)
        u_lower_eps = np.nan*np.zeros(n_lats)
        r_lower_eps = np.nan*np.zeros(n_lats)
        
        s_surf = [sl[0] for sl in s_lower if len(sl[0])>0]
        nbs_surf = [sl[2] for sl in s_lower if len(sl[0])>0]
        
        s_surf.reverse()
        nbs_surf.reverse()
        
        # Compute smooth surface Z and theta cell-wise using SoftMax over appropriate indices
        su = np.unique(sg)
        idx1 = np.array([])    
        
        for si, nbsi in zip(s_surf,nbs_surf):
            # Find indices of latitudes (the same index can appear in several boundary cells
            # because of how s_lower is defined - ignore duplicates. Lists of latitudes
            # and neighbours are reversed for consistency in how variables are overwritten)
            si = np.sort(si)
            idx0 = np.where(np.isin(su,si))[0]
            not_duplicate = ~np.isin(idx0,idx1)
            idx0 = idx0[not_duplicate]
            si = su[idx0]
            pi = p_lower_eps[idx0]
            Phi_eps_i, Z_eps_i, th_eps_i, r_eps_i = self.apply_smoothing_local(si,pi,isub=nbsi,eps=eps)
            
            Phi_lower_eps[idx0] = Phi_eps_i - ld.duals[-1]
            Z_lower_eps[idx0] = Z_eps_i
            th_lower_eps[idx0] = th_eps_i
            u_lower_eps[idx0] = self.get_u(Z_eps_i,si)
            r_lower_eps[idx0] = r_eps_i
            
            idx1 = idx0
            
        self.Phi_lower_eps = Phi_lower_eps
        self.Z_lower_eps = Z_lower_eps
        self.th_lower_eps = th_lower_eps
        self.u_lower_eps = u_lower_eps
        self.r_lower_eps = r_lower_eps
                        
        return Phi_lower_eps, Z_lower_eps, th_lower_eps, u_lower_eps, r_lower_eps
    
    def interpolate_adjusted_lower_surface_vars(self):
        '''Function interpolates adjusted lower boundary variables
        Z, theta and u onto grid and sets as class attributes'''
        pp = self.pp
        sgu = np.unique(self.sg) # grid points in s
        
        i_surf = self.i_surf # surface Laguerre cell indices
        Z_adj = self.Z_adj[i_surf] # surface zonal angular momentum at Lag-cell midpoints
        th_adj = self.th_adj[i_surf] # surface theta at Lag-cell midpoints
        smids = self.surf_s_mids
        
        # interpolate Z using pseudo-s coordinate (1-z/Omega/a^2)^{1/2}, which is approximately linear in s
        s_pseudo = np.sqrt(1-Z_adj/pp.Omega/pp.a**2)
        s_pseudo_interp = np.interp(sgu,smids,s_pseudo)
        Z_lower_adj = pp.Omega*pp.a**2*(1-s_pseudo_interp**2)
        
        # interpolate theta linearly in s
        th_lower_adj = np.interp(sgu,smids,th_adj)
        
        # define u from Z and s-coordinates
        u_lower_adj = self.get_u(Z_lower_adj, sgu)
        
        # assign variables to class
        self.Z_lower_adj = Z_lower_adj
        self.th_lower_adj = th_lower_adj
        self.u_lower_adj = u_lower_adj
        
        return
        
    def push_forward_to_insentropic_coords(self):

        # Get the unique (sine-of-)latitudes, pressure levels and theta levels
        su = self.sg[0,:]
        pu = self.pg[:,0]
        thlevs = self.thlevs
        
        # Define decreasing interpolation of smooth pressure onto theta levels
        # at each latitude individually. Note that theta is extended smoothly below
        # ground and these values are used in the interpolation
        p_interps = []

        for i in np.arange(self.res[1]):
            th_eps_i = self.th_eps[:,i]
            idx = np.where(~np.isnan(th_eps_i))[0]
            p_interp = PchipInterpolator(np.flip(th_eps_i[idx]),np.flip(self.pg[idx,i]),extrapolate=True)
            p_interps.append(p_interp)
        
        # Interpolate pressure onto theta levels at each latitude individually.
        # Derive interpolation coefficients such that interpolation of pressure
        # is linear in log(p) between successive sample points and use these to
        # interpolate Z and r. Then define u from Z
        p_isentropic = np.nan*np.ones([len(thlevs),len(su)])
        Z_isentropic = np.nan*np.ones([len(thlevs),len(su)])
        r_isentropic = np.nan*np.ones([len(thlevs),len(su)])
        u_isentropic = np.nan*np.ones([len(thlevs),len(su)])
        
        for k, p_interp in enumerate(p_interps):
            # interpolate pressure onto theta levels
            idx = thlevs>self.th_lower_eps[k]
            p_isentropic_k = p_interp(thlevs[idx])
            p_isentropic[idx,k] = p_isentropic_k
            
            # get interpolation coefficients for linear interpolation in log(p) 
            i_sort = np.searchsorted(pu,np.sort(p_isentropic_k))
            i_sort[i_sort==len(pu)] = len(pu)-1
            lmbda = (np.log(np.sort(p_isentropic_k))-np.log(pu[i_sort]))/(np.log(pu[i_sort-1])-np.log(pu[i_sort]))
            
            # linearly interpolate Z using interpolation coefficients
            Z_isentropic_k = lmbda*self.Z_eps[i_sort,k] + (1-lmbda)*self.Z_eps[i_sort,k]
            Z_isentropic[idx,k] = Z_isentropic_k
            
            # define zonal wind from Z and latitude
            u_isentropic[idx,k] = self.get_u(Z_isentropic_k,su[k])
            
            # linearly interpolate isentropic density using interpolation coefficients
            r_isentropic_k = lmbda*self.r_eps[i_sort,k] + (1-lmbda)*self.r_eps[i_sort,k]
            r_isentropic[idx,k] = r_isentropic_k
    
        self.p_isentropic = p_isentropic
        self.Z_isentropic = Z_isentropic
        self.u_isentropic = u_isentropic
        self.r_isentropic = r_isentropic
        
        return p_isentropic, Z_isentropic, u_isentropic, r_isentropic
    
    def get_q_isentropic(self):
        '''Function obtains PV on isentropic levels by linear interpolation of 
        PV as a function of circulation on each isentropic level. Known data 
        points come from the 3-d data. Points of evaluation represent zonal angular 
        momentum of the MLM, which is proportional to circulation by zonal symmetry.
        '''
        eartharea = 4*np.pi*self.pp.a**2
        zam_scale_factor = eartharea/2/np.pi
        
        pvlevs = np.flip(self.pvlevs)
        zmom = zam_scale_factor*np.flip(self.bscirc,axis=0)
        q_isentropic = np.nan*np.ones([len(self.thlevs),self.res[1]])
        
        # interpolate PV as function of circulation on each theta level and
        # evaluate interpolant at 2*pi*Z_isentropic
        for k, thlev in enumerate(self.thlevs):
            zmom_k = zmom[:,k]
            i0 = np.where(zmom_k>0)
            zmom_k = np.hstack([0,zmom_k[i0]])
            pvlevs_k = np.hstack([self.pvmaxth[k],pvlevs[i0]])
            q_isentropic[k] = np.interp(self.Z_isentropic[k],zmom_k,pvlevs_k)
            
        self.q_isentropic = q_isentropic
            
        return q_isentropic
    

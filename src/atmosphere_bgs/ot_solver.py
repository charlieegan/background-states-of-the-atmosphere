import numpy as np
import _atmosphere_bgs
import matplotlib.pyplot as plt
from scipy.sparse import coo_array, dia_array
from scipy.sparse.linalg import spsolve
import time

class OTSolver:
    def __init__(self, input_data, ot_tol=1e-4, boundary_res=10000, initial_weights=None):
        '''
        initial_weights -  specified how to initialise damped Newton method
                           can be None (use defualt depending on whether data
                           has been interpolated onto a grid) or 'Voro' (use 
                           translation and rescaling corresponding to Voronoi
                           tessellation).
        '''
        
        # assign point masses and locations to solver
        self.y = input_data.y
        self.tm = input_data.tm
        self.tmn = input_data.tmn
        
        # assign physical and simulation parameters to solver
        self.pp = input_data.pp
        self.sp = _atmosphere_bgs.SimulationParameters(area_tolerance=.2*ot_tol, 
                                                      line_tolerance=1e-4,
                                                      min_line_resolution=3,
                                                      boundary_res=boundary_res,
                                                      max_refine_steps=1000, 
                                                      smin=input_data.smin,
                                                      smax=input_data.smax,
                                                      pmin=input_data.pmin,
                                                      pmax=self.pp.p00 * 1e6)
        
        # if data has been interpolated onto grid then assign this to the solver
        # for use when initialising weights
        self.interpolate_onto_grid = input_data.interpolate_onto_grid
        if self.interpolate_onto_grid:
            self.zam_grid = input_data.zam_grid
            self.th_grid = input_data.th_grid
            self.mass_grid = input_data.mass_grid
        
        # assign solver tolerance and number of target masses to solver and
        # intialise runstats as None
        self.ot_tol = ot_tol
        self.n = self.y.shape[0]
        self.runstats = None
        self.initial_weights = initial_weights
        
    def initialise_weights_grid(self,top_scale=2,initial_weights_perturbation=1e-14):
        '''Initialises weights for target measure with support on a rectangular grid.
        Initialisation corresponds to mass density in (Z,theta) being uniform
        between each successive pair of Z-levels and being uniform between each
        successive pair of theta levels.
        Z values in the matrix should be constant along each row (except at 
        estimated intersection with ground) and increasing with column index.
        Theta values in matrix should be constant along each column and 
        decreasing with row index.
        '''
        # Get data in array-grid format
        zam_grid = self.zam_grid
        th_grid = self.th_grid
        mass_grid = self.mass_grid
        J, K = zam_grid.shape
        
        # Compute cumulative masses
        idx = np.isnan(zam_grid)
        mass_grid[idx] = 0
        mass_zam = np.cumsum(np.sum(mass_grid,axis=1))
        mass_th = np.cumsum(np.sum(mass_grid,axis=0))
        total_mass = np.max(mass_th)
        
        # Initialise weights psi as a matrix with rows corresponding to Z-levels
        # and columns corresponding to Theta-level
        psi = np.zeros(zam_grid.shape)
        
        for k in np.arange(K):
            if k>0:
                p = self.sp.pmin + (mass_th[k-1]/total_mass)*(self.pp.p00-self.sp.pmin)
                psi[0,k] = psi[0,k-1] + self.c_enth(p,th_grid[0,k]) - self.c_enth(p,th_grid[0,k-1])
            if k == K-1 and ~np.isnan(zam_grid[0,k]):
                psi_ext = psi[0,k] + self.c_KE(self.sp.smax,zam_grid[0,k]) - self.pp.cp*th_grid[0,k]*top_scale**self.pp.kappa
            for j in np.arange(J-1):
                s = self.sp.smax - (mass_zam[j]/total_mass)*(self.sp.smax - self.sp.smin)
                psi[j+1,k] = psi[j,k] + self.c_KE(s,zam_grid[j+1,k]) - self.c_KE(s,zam_grid[j,k])
                if k==K-1 and ~np.isnan(zam_grid[j+1,k]) and ~np.isnan(psi[j+1,k]):
                    psi_ext = np.max([psi_ext, psi[j+1,k] + self.c_KE(s,zam_grid[j+1,k]) - self.pp.cp*th_grid[0,k]*top_scale**self.pp.kappa])
        
        # filter out weights for points with no mass or zonal angular momentum
        idx = ~np.isnan(zam_grid)
        psi = psi[idx]
        psi = psi + .5*np.random.uniform(-1,1,len(psi))*np.max(psi - np.min(psi))*initial_weights_perturbation
        psi = np.append(psi,psi_ext)
        
        return psi
    
    def c_KE(self,s,zam):
        '''Pointwise kinetic energy as function of s=sin(lat) and zonal angular momentum'''
        a = self.pp.a
        Omega = self.pp.Omega
        return (1/(1-s**2))*(zam**2/2/a**2) - Omega*zam + (1-s**2)*Omega**2*a**2/2

    def c_enth(self,p,th):
        '''Pointwise enthalpy as function of pressure p and potential temperature theta'''
        cp = self.pp.cp
        p00 = self.pp.p00
        kappa = self.pp.kappa
        return cp*th*(p/p00)**kappa
    
    def initialise_weights_Voro(self, phi_0=-1e+10):
        '''
        Given an (N,2) array of seeds y in (Z,Theta) coordinates, this function returns a weight vector phi such that 
        all cells in the dot-product-Laguerre tessellation generated by phi and the transformed seeds E 
        (defined in this function) are non-empty.
        '''
        sp = self.sp
        pp=self.pp
        
        # define weights
        E = np.array([(0.5 / self.pp.a**2) * self.y[:,0]**2, self.pp.cp * self.y[:,1]]).T
        E_min, E_max = np.min(E, axis=0), np.max(E, axis=0)

        slim = np.array([sp.smin,sp.smax])
        plim = np.array([sp.pmin,pp.p00])
        x0lim = 1/(1-slim**2)
        x1lim=(plim/pp.p00)**pp.kappa

        t_0 = 0.5 * E_max
        t_1 = np.array([x0lim[0],x1lim[0]])

        lmbd = np.min(np.array([np.diff(x0lim)/(E_max[0]/2-E_min[0]/2),np.diff(x1lim)/(E_max[1]/2-E_min[1]/2)]))
        
        t = lmbd * t_0 + t_1

        phi = E @ t - lmbd * 0.25 * np.linalg.norm(E, axis=1)**2  + pp.Omega*self.y[:,0]

        phi = np.append(phi, phi_0)
        
        # Laguerre tesselation generated by phi is equal to Voronoi 
        # tessellation generated by seeds G = -lmbd*E/2 + t, which are defined to
        # lie in the source domain. Here is a plot to check that they lie in the source domain
        plot = False
        
        if plot:
            fig, ax = plt.subplots(1,1,dpi=300)
            
            s=0.1
            
            ax.scatter(-E[:,0]/2,-E[:,1]/2,s=s,alpha=0.3)
            ax.scatter(-E[:,0]/2+t_0[0],-E[:,1]/2+t_0[1],s=s,alpha=0.3)
            ax.scatter(lmbd*(-E[:,0]/2+t_0[0]),lmbd*(-E[:,1]/2+t_0[1]),s=s,alpha=0.3)
            ax.scatter(lmbd*(-E[:,0]/2+t_0[0])+t_1[0],lmbd*(-E[:,1]/2+t_0[1])+t_1[1],s=s,alpha=0.3)
            
            
            fig, ax = plt.subplots(1,1,dpi=300)
            
            ax.scatter(lmbd*(-E[:,0]/2+t_0[0])+t_1[0],lmbd*(-E[:,1]/2+t_0[1])+t_1[1],s=s,alpha=0.3)
            ax.plot(x0lim,2*[x1lim[0]],c='k')
            ax.plot(x0lim,2*[x1lim[1]],c='k')
            ax.plot(2*[x0lim[0]],x1lim,c='k')
            ax.plot(2*[x0lim[1]],x1lim,c='k')
        
        return phi

    def get_bgs(self,
                use_long_double=False, verbose=False,
                max_its=1000, descent_accept_thresh=0.01,
                min_area=1e-10, max_lost_areas=None,
                lr_up=2.0, lr_down=0.5, lr_max=1.0, lr_min=1e-20, lr_init=1e-5,
                max_loss_fraction=1e-4, initial_weights_perturbation=1e-14, top_scale=2):
        """
        run the modified damped Newton solver
        
        Parameters:
        use_long_double : bool -- whether to use long double precision in critical parts of the calculation (default False)
        verbose : bool -- output information about the procedure during execution (default False)
        max_its : int -- maximal number of iterations to run the solver for (default 1000)
        descent_accept_thresh : float < 1 -- threshold of relative descent to accept step (default 0.01)
        min_area : float >= 0.0 -- minimal area to treat as non-empty (default 0)
        max_lost_areas : int >= 0 -- maximal number of cells that may be lost during an accepted step (default None, chosen based on problem size)
        lr_up : float > 1.0 -- factor by which to increase the learning rate (aka. step size) at the start of each iteration (default 2.0)
        lr_down : float < 1.0 -- factor by which to decrease the learning rate after a failed step (default 0.5)
        lr_max : float <= 1.0 -- upper bound for the learning rate (default 1.0)
        lr_min : float > 0 -- lower bound for the learning rate, when it falls below this, the solver is terminated (default 1e-20)
        lr_init : 0 < float <= 1 -- initial learning rate (default 1e-5)
        """
        
        err_goal = self.ot_tol/2
        
        # define functions returning mean error
        get_err = lambda a0, a1 : np.mean(np.abs(a0-a1)/a0)
        
        # maximum number of areas to allow to be lost in an accepted step
        if max_lost_areas is None:
            max_lost_areas = int(np.ceil(self.n))*max_loss_fraction
            
        lr = lr_init
        
        # initialise weights
        if self.initial_weights is None:
            if self.interpolate_onto_grid:
                psi = self.initialise_weights_grid(top_scale=top_scale,initial_weights_perturbation=initial_weights_perturbation)
            else:
                psi = self.initialise_weights_Voro()
        elif type(self.initial_weights) == str:
            psi = self.initialise_weights_Voro()
        else:
            psi = self.initial_weights

        if use_long_double:
            psi = psi.astype(np.float128)
            LD = _atmosphere_bgs.LaguerreDiagram_float128
        else:
            LD = _atmosphere_bgs.LaguerreDiagram_float64

        self.runstats = {"na" : [],  "lr" : [], "maxerr" : [], "l2err" : [], "meanerr" : [], "good" : [], "bad" : [], "tries" : []}
        self.timer = _atmosphere_bgs.Timer([])

        if verbose:
            print(f'initial duals calculated', flush=True)
        
        ld = LD(self.y, psi, self.pp, self.sp)
        self.ld = ld

        if verbose:
            print(f'try to fix {np.sum(ld.areas < min_area)} bad areas', flush=True)
            
        psi = ld.touching_dual(randomize=True)
        self.timer += ld.time
        self.timer += ld.hs.time
        
        psi -= np.mean(psi)
        
        # initialise list to keep track of weights
        psis = [psi]
        
        ld = LD(ld, psi)
        self.ld = ld
        err = get_err(self.tmn,ld.areas)#np.abs(self.tmn - ld.areas)
        good_areas = (ld.areas > min_area)
        if verbose:
            print(f'it={-1}, lr={lr:.2e}, good_areas={np.sum(good_areas)}/{good_areas.shape[0]}', flush=True)

        t00 = time.time()
        for it in range(max_its):

            jac = coo_array(ld.jac_coo(), shape=(self.n+1, self.n+1)).tocsr()
            if jac[self.n, self.n] == 0:
                raise RuntimeError("ERROR: top cell became empty")
            jac = jac[:-1,:-1]
            if self.sp.negative_area_scaling <= 0:
                jac += dia_array((1.0 - good_areas[None,:], [0]), shape=(self.n, self.n))
            d = spsolve(jac, (self.tmn - ld.areas))
            if self.sp.negative_area_scaling <= 0:
                d[~good_areas] = np.mean(d[good_areas]) - 1e9

            cnt_lost_areas_prev = self.n
            lr = np.minimum(lr_max, lr_up * lr)
            cnt_tries = 0
            while lr > lr_min:
                cnt_tries += 1
                # calculate step given step size lr
                psi2 = psi + lr * np.concatenate([d, [0]])

                if cnt_tries > 1:
                    self.timer += ld2.time
                    self.timer += ld2.hs.time
                
                # calculate ld after step
                ld2 = LD(ld, psi2)
                #err2 = get_err(self.tmn,ld2.areas) #np.abs(self.tmn - ld2.areas)
                good_areas2 = (ld2.areas > min_area)

                # if self.sp.negative_area_scaling <= 0:
                #     aerr = np.sum((err/self.tmn)[good_areas])
                #     aerr2 = np.sum((err2/self.tmn)[good_areas])
                # else:
                #     aerr = np.sum((err/self.tmn))
                #     aerr2 = np.sum((err2/self.tmn))

                cnt_lost_areas = np.sum(good_areas & ~good_areas2)
                
                if self.sp.negative_area_scaling <= 0:
                    areas_good = cnt_lost_areas <= max_lost_areas # max_loss_fraction * self.n or (cnt_lost_areas >= cnt_lost_areas_prev and cnt_lost_areas <= max_lost_areas)
                else:
                    areas_good = True
                cnt_lost_areas_prev = cnt_lost_areas

                # areas' ~= areas + lr * jac @ d
                #         = areas + lr * jac @ (jac-1 (tmn - areas))
                #         = areas + lr * (tmn - areas)
                #         = (1 - lr) * areas + lr * tmn
                # => ||areas' - tmn||_1 ~= (1 - lr) ||areas - tmn||_1
                if self.sp.negative_area_scaling <= 0:
                    err = get_err(self.tmn[good_areas],ld.areas[good_areas])
                    err2 = get_err(self.tmn[good_areas2],ld2.areas[good_areas2])
                else:
                    err = get_err(self.tmn,ld.areas)
                    err2 = get_err(self.tmn,ld2.areas)
                    
                descent_good = err2 < (1 - descent_accept_thresh * min(0.1, lr)) * err

                if areas_good and descent_good:
                    self.timer += ld.time
                    self.timer += ld.hs.time
                    # accept step
                    psi = psi2 - np.mean(psi2)
                    ld = ld2
                    ld.detach()
                    err = err2
                    good_areas = good_areas2
                    
                    psis.append(psi)
                    
                    # if some areas are lost and step is accepted, reset learning rate to its initial value
                    if cnt_lost_areas > 0:
                        lr = lr_init
                        
                    break

                else:
                    # reject step
                    if verbose:
                        print(f"failed step at it {it}, lr {lr:.2e}, error {err:.10e} -> {err2:.10e}, areas_good: {areas_good} ({np.sum(good_areas)}, {np.sum(good_areas2)}, {np.sum(good_areas & ~good_areas2)}), descent_good: {descent_good}", flush=True)
                    lr *= lr_down
                    continue

            if self.sp.negative_area_scaling <= 0 and not np.all(good_areas):
                good_areas_prev = good_areas
                if verbose:
                    print(f"try to fix {np.sum(~good_areas_prev)} bad area(s)", flush=True)
                    
                psi = ld.touching_dual(randomize=True)
                self.timer += ld.time
                self.timer += ld.hs.time
                ld = LD(ld, psi, False)
                good_areas = (ld.areas > min_area)

                if verbose:
                    print(f"managed to fix {np.sum(good_areas & ~good_areas_prev)} and broke {np.sum(good_areas_prev & ~good_areas)}", flush=True)

            if lr <= lr_min:
                raise RuntimeError("can not find sufficiently small learning rate")

            self.ld = ld
            
            abs_errs = np.abs(self.tmn - ld.areas)
            self.runstats["lr"] += [lr]
            self.runstats["maxerr"] += [np.max(abs_errs / self.tmn)]
            self.runstats["l2err"] += [np.sum((abs_errs / self.tmn)**2)**.5]
            self.runstats["meanerr"] += [np.mean(abs_errs / self.tmn)]
            self.runstats["good"] += [np.sum(good_areas)]
            self.runstats["bad"] += [good_areas.shape[0] - self.runstats["good"][-1]]
            self.runstats["tries"] += [cnt_tries]

            if verbose:
                print(f'it={it}, lr={lr:.2e}, good_areas={self.runstats["good"][-1]}/{good_areas.shape[0]}, meanerr={self.runstats["meanerr"][-1]:.6e}, l2err = {self.runstats["l2err"][-1]:.6e}, max_err={self.runstats["maxerr"][-1]:.6e}', flush=True)

            if err < err_goal:
                break

        t11 = time.time()
        if verbose:
            print(f"finished in {t11 - t00:.2f}s", flush=True)

        self.timer += ld.time
        self.timer += ld.hs.time
            
        # assign variables to class
        ld.detach()
        self.ld = ld
        return ld, psis
    
    def plot_runstats(self):
        """
        plot the behavior of the last run
        """
        
        if self.runstats is None:
            return
        
        # plot convergence behavior
        fig, axs = plt.subplots(3, 1, figsize=(10, 2*3), dpi=100, sharex=True, tight_layout=True)

        pi = 0
        axs[pi].set_yscale("log")
        axs[pi].set_ylabel("lr")
        axs[pi].plot(self.runstats["lr"], label="step size")
        axs[pi].legend(loc="upper right")

        pi += 1
        axs[pi].set_ylabel("err")
        axs[pi].plot(self.runstats["maxerr"], label="max")
        axs[pi].plot(self.runstats["l2err"], label="l2")
        axs[pi].plot(self.runstats["meanerr"], label="mean")
        axs[pi].set_yscale("log")
        axs[pi].legend(loc="upper right")

        pi += 1
        axs[pi].set_ylabel("count")
        axs[pi].plot(self.runstats["bad"], label="bad areas")
        axs[pi].plot(self.runstats["tries"], label="step tries")
        axs[pi].legend(loc="upper right")
        

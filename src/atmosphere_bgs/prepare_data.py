import numpy as np
import _atmosphere_bgs
import re

class DataLoader:
    
    def __init__(self,path,pmin=None,p00=None,load_all=True,interpolate_onto_grid=True,skip=0,split_param=None):
        
        # parse the input data text file and make a dictionary of data arrays
        with open(path) as f:
            self.s = f.read()
        
        self.float_pattern = r'[-+]?(\d+([.,]\d*)?|[.,]\d+)([eE][-+]?\d+)?'
        self.path = path
        self.data_dict = {}

        try:
            self.latitudes_cnt = int(next(re.compile(r"(\d+)\s+LATITUDES ON GAUSSIAN GRID").finditer(self.s)).group(1))
        except StopIteration:
            raise RuntimeError("did not find number of \"LATITUDES ON GAUSSIAN GRID\"")
        
        try:
            self.tracer_mixing_ratio_contour_cnt = int(next(re.compile(r"(\d+)\s+TRACER MIXING RATIO CONTOURS").finditer(self.s)).group(1))
        except StopIteration:
            raise RuntimeError("did not find number of \"TRACER MIXING RATIO CONTOURS\"")

        try:
            self.isentropic_level_cnt = int(next(re.compile(r"(\d+)\s+ISENTROPIC LEVELS").finditer(self.s)).group(1))
        except StopIteration:
            raise RuntimeError("did not find number of \"ISENTROPIC LEVELS\"")

        self.parse_block(self.latitudes_cnt, "LATITUDES ON GAUSSIAN GRID")
        self.parse_block(self.isentropic_level_cnt, "ISENTROPIC LEVELS")
        self.parse_block(self.tracer_mixing_ratio_contour_cnt, "TRACER MIXING RATIO CONTOURS")
        self.parse_block(self.isentropic_level_cnt, "FACTOR TO CONVERT FROM LAIT TO ERTEL PV")
        self.parse_block(self.isentropic_level_cnt * self.tracer_mixing_ratio_contour_cnt, "MASS INTEGRALS IN PV-THETA COORDINATES")
        self.parse_block(self.isentropic_level_cnt * self.tracer_mixing_ratio_contour_cnt, "AREA INTEGRALS IN PV-THETA COORDINATES")
        self.parse_block(self.isentropic_level_cnt * self.tracer_mixing_ratio_contour_cnt, "CIRCULATION INTEGRALS IN PV-THETA COORDINATES")
        self.parse_block(self.latitudes_cnt, "BACKGROUND PRESSURE ON TOP BOUNDARY")

        if load_all:
            try:
                self.data_dict["DATA BASE TIME"] = int(next(re.compile(r"DATA BASE TIME IS\s+(\d+)").finditer(self.s)).group(1))
            except StopIteration:
                print("WARNING: did not find \"DATA BASE TIME\"")
                
            try:
                self.data_dict["TOP BOUNDARY IN ISENTROPIC COORDS"] = float(next(re.compile(r"(" + self.float_pattern + r")\s+IS TOP BOUNDARY IN ISENTROPIC COORDS").finditer(self.s)).group(1))
            except StopIteration:
                print("WARNING: did not find \"TOP BOUNDARY IN ISENTROPIC COORDS\"")
            
            try:
                match = next(re.compile(r"(" + self.float_pattern + r")\s+(" + self.float_pattern + r")\s+MAX AND MIN VALUES OF SURFACE THETA").finditer(self.s))
                self.data_dict["MAX VALUE OF SURFACE THETA"] = float(match.group(1))
                self.data_dict["MIN VALUE OF SURFACE THETA"] = float(match.group(2))
            except StopIteration:
                print("WARNING: did not find \"MAX AND MIN VALUES OF SURFACE THETA\"")

            try:
                self.data_dict["TOTAL PVS ENCLOSED BY LOWEST VALUE TRACER CONTOUR"] = float(next(re.compile(r"TOTAL PVS ENCLOSED BY LOWEST VALUE TRACER CONTOUR\s+(" + self.float_pattern + r")").finditer(self.s)).group(1))
            except StopIteration:
                print("WARNING: did not find \"TOTAL PVS ENCLOSED BY LOWEST VALUE TRACER CONTOUR\"")
            
            try:
                self.data_dict["TOTAL ATM MASS ENCLOSED BY LOWEST VALUE TRACER CONTOUR"] = float(next(re.compile(r"TOTAL ATM MASS ENCLOSED BY LOWEST VALUE TRACER CONTOUR\s+(" + self.float_pattern + r")").finditer(self.s)).group(1))
            except StopIteration:
                print("WARNING: did not find \"TOTAL ATM MASS ENCLOSED BY LOWEST VALUE TRACER CONTOUR\"")
        
            self.parse_block(self.isentropic_level_cnt, "MAX PV ON THETA LEVELS")
            self.parse_block(self.isentropic_level_cnt, "MIN PV ON THETA LEVELS")
            self.parse_block(self.isentropic_level_cnt, "AREA INTEGRAL OVER POLAR SHELLS THETA COORDINATES")
            self.parse_block(self.isentropic_level_cnt, "MASS INTEGRALS OVER POLAR SHELLS IN THETA COORDINATES")
            self.parse_block(self.isentropic_level_cnt, "CIRCULATION INTEGRALS OVER POLAR SHELLS IN THETA COORDINATES")
            self.parse_block(self.latitudes_cnt, "BACKGROUND SURFACE GEOPOTENTIAL")
            self.parse_block(self.isentropic_level_cnt, r"BACKGROUND u\*cos\(phi\) AT EQUATOR ON THETA LEVELS")
        
        # get physical parameters
        self.pp = _atmosphere_bgs.PhysicalParameters()

        # assign minimum pressure level and reference pressure to class
        if pmin is None:
            # Define minimum pressure to be the mean value of the zonal average
            # pressure on the top isentropic level, copmuted as an integral over
            # [0,1] using the trapezium rule
            ptop = self.data_dict['BACKGROUND PRESSURE ON TOP BOUNDARY']
            ptop = np.flip(ptop) # flip to be ordered accending with latitude
            sgrid = np.sin(np.deg2rad(self.data_dict['LATITUDES ON GAUSSIAN GRID']))
            sgrid = np.flip(sgrid) # order to be acsending
            
            area_top = np.diff(sgrid)*((ptop[:-1]-self.pmin)+(ptop[1:] - self.pmin))/2
            pmin = np.sum(area_top)
        self.pmin = pmin
        
        if p00 is not None:
            self.pp.p00 = p00
        
        # get target measure
        self.get_target_measure(interpolate_onto_grid=interpolate_onto_grid,skip=skip,split_param=split_param)
        
        # record whether or not interpolation has been used
        self.interpolate_onto_grid = interpolate_onto_grid
        
    def _multifloat_pattern(self, mincnt, maxcnt=None):
        if maxcnt is None:
            maxcnt = mincnt
        return "(" + self.float_pattern + r"\s+){" + str(mincnt-1) + "," + str(maxcnt-1) + "}" + self.float_pattern
    
    def parse_block(self, cnt, title):
        """
        find block of `cnt` floats starting with title `title`
        """
        
        try:
            r = re.compile(title)
            match = next(r.finditer(self.s))
            
            try:
                r = re.compile(self._multifloat_pattern(1, cnt))
                vals = [float(x) for x in next(r.finditer(self.s, pos=match.end())).group(0).split()]
                self.data_dict[title] = np.array(vals)
            except StopIteration:
                print(f"WARNING: could not find values for block {title}")
                
            if self.data_dict[title].shape[0] != cnt:
                print(f"WARNING: could not find all {cnt} values for block {title}, only found {self.data_dict[title].shape[0]}")

        except StopIteration:
            print(f"WARNING: could not find block {title}")
            
    def get_target_measure(self,interpolate_onto_grid=True,skip=0,split_param=None):
        
        '''Function that defines target measure to be used in the OT problem. 
        
        Default is to interpolate mass on each theta level linearly in the 'pseudo-s'
        coordinate
        
            spseudo = (1-Z/Omega/a**2)**(1/2).
            
        This is chosen because it is equal to s if u=0, and mass is linear in (s,p) 
        because it is proportional to area. The interpolation nodes in 
        pseudo-s are from a Gaussian grid in latitude. If the minimum speudo on
        a given theta level is greater than the minimum grid point then it is 
        used as an interpolation node and smaller interpolation nodes are discarded.
        This effectively deals with theta levels that intersect the ground.
        smin and smax that define the lateral extent of the (s,p) domain that is
        to be tessellated are also defined here.
        smax is defined to to avoid large negative pressure gradients at the pole.
        Specifically, it is chosen as the turning point of the pressure 'kernel'
        for the minimum Z value on the lowest theta level. This excludes a 
        spherical cap around the pole. smin is defined to exclude a band from the 
        equator of the same mass as this cap.
        
        Data is typically truncated at some theta level. A row of masses is 
        added on a theta level that is above the top theta level of the data 
        to represent the missing data. The input data contains the zonal average 
        pressure on the top theta level from the 3-d state as a function of latitude. 
        This is used is used to approximate the mass above the top theta level
        and between sucessive grid points in s.
        For each grid point s, a corresponding mass point is added at Z satisfying 
        
            s=(1-Z/Omega/a**2)**(1/2).
        
        If mass is not interpolated then point masses are defined directly from
        input data by differencing cumulative mass integrals to get masses and
        by differencing and rescaling circulation integrals to get zonal angular
        momentum values. Theta levels are retained.
        
        To do:
            - interpolate PV and laitPV and assign these variables to the class instance
            - add masses for top boundary to non-interpolant version
        
        Inputs:
            
            interpolate_onto_grid - boolean; determines whether or not input data is interpolated onto a grid
            skip - integer; number of nodes to skip in interpolation in order to reduce the resolution of the interpolation
            split_param - float or None; number between 0 and 1 used to force mass 
                          and circulation integrals to be strictly decreasing with 
                          PV on every theta level; bigger value means bigger enforced increase
                          
        Outputs (assigned to class instance):
            
            if interpolate_onto_grid is True:
                zam_grid - (nnodes,(nthlev+1)) numpy array; gridded zonal angular momentum values; constant along each row (except at estimated intersection with ground) and increasing with column index
                th_grid - (nnodes,(nthlev+1)) numpy array; gridded theta values; constant along each column and decreasing with row index
                mass_grid - (nnodes,(nthlev+1)) numpy array; gridded masses
            in any case:
                smin - float; minimum s
                smax - float; maximum s
                y - (N,2) numpy array; 2-d locations of point masses in (Z,theta)
                tmn - (N,) numpy array; target masses associated to points y normalised to sum to area of (s,p) domain [smin,smax]x[pmin,pmax]
        '''
        
        # extract input data from the dictionary
        latitudes = self.data_dict['LATITUDES ON GAUSSIAN GRID']
        pvlev = self.data_dict['TRACER MIXING RATIO CONTOURS']
        thlev = self.data_dict['ISENTROPIC LEVELS']
        lait2pv = self.data_dict['FACTOR TO CONVERT FROM LAIT TO ERTEL PV']
        bscirc = self.data_dict['CIRCULATION INTEGRALS IN PV-THETA COORDINATES']
        bsmass = self.data_dict['MASS INTEGRALS IN PV-THETA COORDINATES']
        
        # get physical and simulation parameters
        pp = self.pp
        a = pp.a
        eartharea = 4*np.pi*pp.a**2
        Omega = pp.Omega
        
        # reshape data to have rows as pvlevels and columns as theta levels
        npvlev = pvlev.shape[0]
        nthlev = thlev.shape[0]
        bscirc = np.reshape(bscirc,(nthlev, npvlev)).T
        bsmass = np.reshape(bsmass,(nthlev,npvlev)).T
        
        # keep only theta levels that have some mass and circulation
        idx = (bsmass[0] > 0)*(bscirc[0] > 0)
        bscirc = bscirc[:,idx]
        bsmass = bsmass[:,idx]
        thlev = thlev[idx]
        nthlev = thlev.shape[0]
        
        # Define mass scale factors on each level (final theta level depth is repeated)
        dth = np.absolute(np.append(np.diff(thlev),thlev[-1]-thlev[-2]))
        mass_scale_factor = eartharea*dth
        
        # Define angular momentum scale factor for converting circulation
        zam_scale_factor = eartharea/2/np.pi
        
        # if interpolating, make circulation integrals strictly decreasing on every theta level
        # otherwise, ignore flat parts by setting split_param to 0
        if interpolate_onto_grid:
            if split_param is None or split_param == 0:
                split_param = 1e-2
        elif split_param is None:
            split_param = 0
            
        ## force difference in cumulative masses on each theta level to be positive
        for k, mass_k in enumerate(bsmass.T):
            if np.any(np.diff(mass_k)>0):
                idx = np.where(np.diff(mass_k)>0)[0]
                for i in np.flip(idx):
                    mass_k[idx] = mass_k[idx+1]
            bsmass[:,k] = mass_k
            
        ## force difference in circulations on each theta level to be positive
        for k, circ_k in enumerate(bscirc.T):
            if np.any(np.diff(circ_k)>0):
                idx = np.where(np.diff(circ_k)>0)[0]
                for i in np.flip(idx):
                    circ_k[idx] = circ_k[idx+1]
            bscirc[:,k] = circ_k
            
        ## force both mass and circulation to be strictly decreasing on each theta level
        for k, (circ_k, mass_k) in enumerate(zip(bscirc.T,bsmass.T)):
            # find where cumulative mass or circulation are flat as function of PV
            zero_mass = np.hstack([False,np.diff(mass_k) == 0])
            zero_circ = np.hstack([False,np.diff(circ_k) == 0])
            zero_mass_or_circ = np.where(zero_mass+zero_circ)[0]
            
            if len(zero_mass_or_circ)>0:
                # Could have several places on same theta level where this happens
                # Split indices into blocks of consecutive numbers
                split = np.where((zero_mass_or_circ - np.roll(zero_mass_or_circ,shift=1))>1)[0]
                zms = np.split(zero_mass_or_circ,split)
            
                # for each block of consecutive indices, share mass and circulation across them
                for zms_i in zms:
                    # only split mass where it exists
                    if zms_i[-1]<len(mass_k)-2:
                        mass_to_split = split_param*(mass_k[zms_i[0]]-mass_k[zms_i[-1]+1])
                        circ_to_split = split_param*(circ_k[zms_i[0]]-circ_k[zms_i[-1]+1])
                        
                        split_mass = mass_to_split*np.arange(1,len(zms_i)+1)/len(zms_i)
                        split_circ = circ_to_split*np.arange(1,len(zms_i)+1)/len(zms_i)
                        
                        mass_k[zms_i] = mass_k[zms_i] - split_mass
                        circ_k[zms_i] = circ_k[zms_i] - split_circ
                        
            bscirc[:,k] = circ_k
            bsmass[:,k] = mass_k
            
        self.bscirc = bscirc
        self.pvlevs = pvlev
        self.thlevs = thlev
        
        if not interpolate_onto_grid:
            # define matrix of theta values
            th = np.vstack(npvlev*[thlev])
            # define masses by differencing cumulative mass and rescaling
            mass = np.vstack([-np.diff(bsmass,axis=0),bsmass[-1][None,:]])*mass_scale_factor
            # define zonal angular momentum by rescaling circulation
            bscirc_diff = np.diff(bscirc,axis=0)
            bscirc_diff = np.vstack([bscirc_diff,bscirc_diff[-1]])
            zam = zam_scale_factor*(bscirc + .5*bscirc_diff)
            # define smin and smax
            smin = np.min(np.sin(np.deg2rad(latitudes)))
            smax = np.max(np.sin(np.deg2rad(latitudes)))
            # define vectorised seeds, masses and normalised masses to be used in OT routine
            y = np.vstack([np.ravel(zam),np.ravel(th)]).T
            tm = np.ravel(mass)
            source_mass = (self.pp.p00 - self.pmin)*(smax - smin)
            tmn = source_mass*tm/np.sum(tm)
            # Retain only points where mass and zonal angular momentum are positive
            idx = np.where((tmn>0)*(y[:,0]>0))
            y = y[idx]
            tm = tm[idx]
            tmn = tmn[idx]
        else:
            
            # Interpolate mass linearly in pseudo-s coordinate onto Gaussian grid of latitudes
            ## Get maximum of Omega*a**2 and maximum zonal angular momentum from the data
            bscirc_max = np.max([np.max(bscirc),Omega*a**2/zam_scale_factor])
                
            ## Define pseudo-s coordinate
            spseudo = np.sqrt(1-(bscirc/bscirc_max))
            
            ## Define interpolation nodes to be sine of the given latitudes in ascending order
            snodes = np.flip(np.sin(np.deg2rad(latitudes)))
            if skip>0:
                snodes = snodes[::skip]
                
            ## set up matrices to store angmom and mass with one extra level for theta above the data
            nnodes = len(snodes)+1
            zam_grid = np.nan*np.zeros([nnodes,nthlev+1])
            mass_grid = np.nan*np.zeros([nnodes,nthlev+1])
            
            ## do interpolation on each theta level
            for k, cum_mass_k in enumerate(bsmass.T):
                # Get points in pseudo-s on k-th theta level, select interpolation
                # nodes that are above the minimum and append the minimium.
                # If the gap to the minimum node is small, this can create 
                # prohibitively small masses. In case, discard the minimum node.
                spseudo_k = spseudo[:,k]
                sk_min = spseudo_k[0]
                idx = snodes>sk_min
                snodes_k = snodes[idx]
                if (snodes_k[0]-sk_min)/(snodes_k[1]-snodes_k[0])>.1:
                    snodes_k = np.hstack([sk_min,snodes_k])
                else:
                    snodes_k = snodes_k[1:] # discard minimum node
                    snodes_k = np.hstack([sk_min,snodes_k])
                    
                # append zero mass at pole for interpolation and differencing
                spseudo_k = np.append(spseudo_k,1)
                cum_mass_k = np.append(cum_mass_k,0)
                snodes_k = np.append(snodes_k,1)
                
                # interpolate cumulative mass linearly in s onto interpolation nodes and
                # difference cumulative mass and rescale to get mass density
                cm_interp_k = np.interp(snodes_k,spseudo_k,cum_mass_k)
                mass_k = -mass_scale_factor[k]*np.diff(cm_interp_k)
                
                # get midpoints in s and define corresponding angular momentum values
                smids_k = .5*(snodes_k + np.roll(snodes_k,shift=1))[1:]
                zam_k = zam_scale_factor*bscirc_max*(1-smids_k**2)
                
                # flip angmom and mass vectors to be ordered increasing in zonal angular
                # momentum (effectively from pole to equator) and assign to matrices
                nk = len(mass_k)
                zam_grid[:nk,k+1] = np.flip(zam_k)
                mass_grid[:nk,k+1] = np.flip(mass_k)
                
                # check that no mass is lost
                rel_mass_loss_k = (np.sum(mass_k)/mass_scale_factor[k] - cum_mass_k[0])/cum_mass_k[0]
                if np.abs(rel_mass_loss_k) > 1e-11:
                    raise ValueError('Some mass is lost during interpolation on level '+str(k))
                
                # add masses to represent the missing theta levels using the same
                # interpolation nodes as for the top theta level
                if thlev[k] == np.max(thlev):
                    # estimate masses and zonal angular momentum using zonal average pressure from 3-d data
                    ptop = self.data_dict['BACKGROUND PRESSURE ON TOP BOUNDARY']
                    ptop = np.flip(ptop) # flip to be ordered accending with latitude
                    sgrid = np.sin(np.deg2rad(self.data_dict['LATITUDES ON GAUSSIAN GRID']))
                    sgrid = np.flip(sgrid) # order to be acsending
                    
                    # interpolate zonal average pressure onto interpolation nodes
                    ptop = np.interp(snodes_k,sgrid,ptop)
                    
                    # Compute areas using trapezium rule
                    mass_top = np.diff(snodes_k)*((ptop[:-1]-self.pmin)+(ptop[1:] - self.pmin))/2
                    
                    ## define zonal angular momentum values for extra theta level
                    s_top = smids_k
                    zam_top = zam_scale_factor*bscirc_max*(1-s_top**2)
    
                    ## define theta value at top to be the maximum theta value plus a level depth
                    th_top = thlev[0] + (thlev[0]-thlev[1])
                    
                    # flip angmom and mass vectors to be ordered increasing in zonal angular momentum (effectively from pole to equator)
                    nt = len(mass_top)
                    #zam_grid[:nt,k] = np.flip(zam_top)
                    #mass_grid[:nt,k] = np.flip(mass_top)
                    
            # define theta matrix such that each column corresponds to a theta level and
            # each row corresponds (roughly) to a Z-level
            th_grid = np.vstack(nnodes*[np.hstack([th_top,thlev])])
                    
            # check total mass is conserved
            total_mass = np.sum(bsmass[0]*mass_scale_factor)
            masses = mass_grid[:,1:][~np.isnan(mass_grid[:,1:])]
            rel_mass_loss = np.abs((total_mass - np.sum(masses))/total_mass)
            
            if rel_mass_loss > 1e-11:
                raise ValueError('Some mass is lost during interpolation')
            
            # Define smin and smax such that smax is at turning point of p-surface
            # corresponding to minimal zonal angular momentum on lowest theta level
            idx = ~np.isnan(zam_grid[:,-1])
            zmin = np.min(zam_grid[idx,-1])
            smax = np.sqrt(1-zmin/Omega/a**2)
            
            # Define smin so that areas in (s,p) excluded at pole and equator are the same
            smin = 1-smax
            
            # rescale mass and add top row
            source_mass = (self.pp.p00 - self.pmin)*(smax - smin)
            frac_top = np.sum(mass_top)/source_mass
            mass_grid = mass_grid*source_mass*(1-frac_top)/np.sum(mass_grid[~np.isnan(mass_grid)])
            mass_grid[:nt,0] = np.flip(mass_top)
            zam_grid[:nt,0] = np.flip(zam_top)
            
            # Define vectorised seeds, masses and normalised masses to be used in OT routine
            idx = ~np.isnan(zam_grid)
            zam = zam_grid[idx]
            th = th_grid[idx]
            
            y = np.vstack([zam,th]).T
            tmn = mass_grid[idx]
            
            # Assign variables to the class
            self.zam_grid = zam_grid
            self.th_grid = th_grid
            self.mass_grid = mass_grid
        
        # Assign variables to the class
        self.smin = smin
        self.smax = smax
        self.y = y
        self.tmn = tmn
        
        return y, tmn
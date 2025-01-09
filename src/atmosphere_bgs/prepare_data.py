import numpy as np
from numpy import matlib
import _atmosphere_bgs
import re

class DataLoader:
    
    def __init__(self, path, pmin=10, nextra=0, load_all=False):
        
        # parse the input data text file and make a dictionary of data arrays
        with open(path) as f:
            self.s = f.read()
        
        self.float_pattern = '[-+]?(\d+([.,]\d*)?|[.,]\d+)([eE][-+]?\d+)?'
        self.path = path
        self.data_dict = {}

        try:
            self.latitudes_cnt = int(next(re.compile("(\d+)\s+LATITUDES ON GAUSSIAN GRID").finditer(self.s)).group(1))
        except StopIteration:
            raise RuntimeError("did not find number of \"LATITUDES ON GAUSSIAN GRID\"")
        
        try:
            self.tracer_mixing_ratio_contour_cnt = int(next(re.compile("(\d+)\s+TRACER MIXING RATIO CONTOURS").finditer(self.s)).group(1))
        except StopIteration:
            raise RuntimeError("did not find number of \"TRACER MIXING RATIO CONTOURS\"")

        try:
            self.isentropic_level_cnt = int(next(re.compile("(\d+)\s+ISENTROPIC LEVELS").finditer(self.s)).group(1))
        except StopIteration:
            raise RuntimeError("did not find number of \"ISENTROPIC LEVELS\"")

        self.parse_block(self.latitudes_cnt, "LATITUDES ON GAUSSIAN GRID")
        self.parse_block(self.isentropic_level_cnt, "ISENTROPIC LEVELS")
        self.parse_block(self.tracer_mixing_ratio_contour_cnt, "TRACER MIXING RATIO CONTOURS")
        self.parse_block(self.isentropic_level_cnt, "FACTOR TO CONVERT FROM LAIT TO ERTEL PV")
        self.parse_block(self.isentropic_level_cnt * self.tracer_mixing_ratio_contour_cnt, "MASS INTEGRALS IN PV-THETA COORDINATES")
        self.parse_block(self.isentropic_level_cnt * self.tracer_mixing_ratio_contour_cnt, "CIRCULATION INTEGRALS IN PV-THETA COORDINATES")

        if load_all:
            try:
                self.data_dict["DATA BASE TIME"] = int(next(re.compile("DATA BASE TIME IS\s+(\d+)").finditer(self.s)).group(1))
            except StopIteration:
                print("WARNING: did not find \"DATA BASE TIME\"")
                
            try:
                self.data_dict["TOP BOUNDARY IN ISENTROPIC COORDS"] = float(next(re.compile("(" + self.float_pattern + ")\s+IS TOP BOUNDARY IN ISENTROPIC COORDS").finditer(self.s)).group(1))
            except StopIteration:
                print("WARNING: did not find \"TOP BOUNDARY IN ISENTROPIC COORDS\"")
            
            try:
                match = next(re.compile("(" + self.float_pattern + ")\s+(" + self.float_pattern + ")\s+MAX AND MIN VALUES OF SURFACE THETA").finditer(self.s))
                self.data_dict["MAX VALUE OF SURFACE THETA"] = float(match.group(1))
                self.data_dict["MIN VALUE OF SURFACE THETA"] = float(match.group(2))
            except StopIteration:
                print("WARNING: did not find \"MAX AND MIN VALUES OF SURFACE THETA\"")

            try:
                self.data_dict["TOTAL PVS ENCLOSED BY LOWEST VALUE TRACER CONTOUR"] = float(next(re.compile("TOTAL PVS ENCLOSED BY LOWEST VALUE TRACER CONTOUR\s+(" + self.float_pattern + ")").finditer(self.s)).group(1))
            except StopIteration:
                print("WARNING: did not find \"TOTAL PVS ENCLOSED BY LOWEST VALUE TRACER CONTOUR\"")
            
            try:
                self.data_dict["TOTAL ATM MASS ENCLOSED BY LOWEST VALUE TRACER CONTOUR"] = float(next(re.compile("TOTAL ATM MASS ENCLOSED BY LOWEST VALUE TRACER CONTOUR\s+(" + self.float_pattern + ")").finditer(self.s)).group(1))
            except StopIteration:
                print("WARNING: did not find \"TOTAL ATM MASS ENCLOSED BY LOWEST VALUE TRACER CONTOUR\"")
        
            self.parse_block(self.isentropic_level_cnt, "MAX PV ON THETA LEVELS")
            self.parse_block(self.isentropic_level_cnt, "MIN PV ON THETA LEVELS")
            self.parse_block(self.isentropic_level_cnt * self.tracer_mixing_ratio_contour_cnt, "AREA INTEGRALS IN PV-THETA COORDINATES")
            self.parse_block(self.isentropic_level_cnt, "AREA INTEGRAL OVER POLAR SHELLS THETA COORDINATES")
            self.parse_block(self.isentropic_level_cnt, "MASS INTEGRALS OVER POLAR SHELLS IN THETA COORDINATES")
            self.parse_block(self.isentropic_level_cnt, "CIRCULATION INTEGRALS OVER POLAR SHELLS IN THETA COORDINATES")
            self.parse_block(self.latitudes_cnt, "BACKGROUND PRESSURE ON TOP BOUNDARY")
            self.parse_block(self.latitudes_cnt, "BACKGROUND SURFACE GEOPOTENTIAL")
            self.parse_block(self.isentropic_level_cnt, "BACKGROUND u\*cos\(phi\) AT EQUATOR ON THETA LEVELS")
        
        # get physical parameters
        self.pp = _atmosphere_bgs.PhysicalParameters()
        
        # get max and min s values
        self.smax = np.sin(2*np.pi*np.max(self.data_dict['LATITUDES ON GAUSSIAN GRID'])/360)
        self.smin = np.sin(2*np.pi*np.min(self.data_dict['LATITUDES ON GAUSSIAN GRID'])/360)
        self.pmin = pmin
        self.nextra = nextra
        
        # get seeds and masses
        self.get_bgs_target_measure(nextra=nextra)
        
    def _multifloat_pattern(self, mincnt, maxcnt=None):
        if maxcnt is None:
            maxcnt = mincnt
        return "(" + self.float_pattern + "\s+){" + str(mincnt-1) + "," + str(maxcnt-1) + "}" + self.float_pattern
    
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
            
    def get_bgs_target_measure(self,nextra=0):
        '''
        This function returns seed locations and target masses given input data, physical parameters, and simulation parameters.
        '''
        # extract input data from the dictionary
        latitudes = self.data_dict['LATITUDES ON GAUSSIAN GRID']
        pv_lev = self.data_dict['TRACER MIXING RATIO CONTOURS']
        th_lev = self.data_dict['ISENTROPIC LEVELS']
        lait_to_pv = self.data_dict['FACTOR TO CONVERT FROM LAIT TO ERTEL PV']
        bs_circ = self.data_dict['CIRCULATION INTEGRALS IN PV-THETA COORDINATES']
        bs_mass = self.data_dict['MASS INTEGRALS IN PV-THETA COORDINATES']

        # get physical and simulation parameters
        pp = self.pp
        
        earth_radius = pp.a
        earth_area = 4*np.pi*pp.a**2
        Omega = pp.Omega

        # parameter used to place extra mass points near the North Pole in target space
        ang_mom_min = Omega*earth_radius**2*(1-self.smax**2) # planetary zonal ang mom at polar cap latitude

        # number of pv and theta levels
        num_pv_lev = pv_lev.shape[0]; num_th_lev = th_lev.shape[0]

        # reshape and rescale circulation and areas vectors
        bs_circ = np.reshape(bs_circ,(num_th_lev, num_pv_lev)).T
        bs_mass = np.reshape(bs_mass,(num_th_lev,num_pv_lev)).T
        
        d_mass_all = np.diff(bs_mass,axis=0)
        if (d_mass_all > 0).any():
            None
            #raise ValueError('Mass data is inconsistent')

        # find theta half-levels and layer depth in terms of theta
        d_th = -np.diff(th_lev) # this should be positive
        d_th = np.ravel(np.append(d_th,d_th[-1])) # dth is the weighting to obtain mass of each isentropic layer
        th_lev_h = th_lev + 0.5*d_th
        th_lev_h = np.ravel(np.append(th_lev_h,th_lev[-1]-0.5*d_th[-1]))

        # get target masses and angular momentum (horizontal coordinates of seeds)
        epv_field = np.outer(pv_lev,lait_to_pv) # (num_pv_lev,num_th_lev) numpy array of Ertel PV values

        ang_mom = np.zeros([num_pv_lev + nextra -1,num_th_lev])
        mass_weight = np.zeros([num_pv_lev + nextra -1,num_th_lev])
        pv_mid = np.zeros([num_pv_lev + nextra - 1,num_th_lev])
        
        neg_mass_1 = np.zeros([num_pv_lev + nextra -1,num_th_lev])
        neg_mass_2 = np.zeros([num_pv_lev + nextra -1,num_th_lev])
        
        for m in np.arange(num_th_lev):

            d_pv = np.diff(epv_field[:,m])
            d_mass = np.diff(bs_mass[:,m])
            d_circ = np.diff(bs_circ[:,m])

            pv_mid[np.arange(num_pv_lev-1),m] = epv_field[np.arange(num_pv_lev-1),m] + 0.5*d_pv
            c_mid = bs_circ[np.arange(num_pv_lev-1),m] + 0.5*d_circ

            ang_mom[np.arange(num_pv_lev-1),m] = earth_area*c_mid/(2*np.pi)
            mass_weight[np.arange(num_pv_lev-1),m] = -earth_area*d_mass*d_th[m] # d_mass should be negative
            
            if mass_weight[np.arange(num_pv_lev-1),m].any() < 0:
                None
                #raise ValueError('negative mass at location 0')

            # Find lowest full PV level for which bs_circ=0, if it exists (i.e., at North Pole).
            # Introduce a new mass point at position Z_min to make sure
            # that there is a column of points next to pole.
            # Shift the position of the closest existing point so that the mass is 
            # re-distributed across these two points (total mass is unchanged).
            elz = np.where(bs_circ[:,m] == 0)
            if elz[0].shape[0]>0:
                kp = int(elz[0][0])
            else:
                kp = 0

            if kp > 0:
                zkpm1 = earth_area*bs_circ[kp-1,m]/(2*np.pi)
                if zkpm1 > ang_mom_min:
                    dmom = zkpm1/(nextra+1.)
                    pvtop = epv_field[kp,m]
                    massscal = bs_mass[kp-1,m]/zkpm1

                    qscal = (epv_field[kp-1,m]-pvtop)/zkpm1
                    zextra = np.zeros([nextra+2])
                    mextra = np.zeros([nextra+2])

                    # Find angular momentum at new introduced points
                    # between Z(kp-1) and Z_min.

                    for i in np.arange(nextra+2):
                        zextra[i] = zkpm1 - i*dmom

                    #m,thlev(m),zextra
                    mextra = massscal*zextra

                    # Find mass weights, and midpoint Z and Q corresponding to each
                    # introduced interval.

                    for i in np.arange(nextra+1):
                        zmid = zextra[i]-0.5*dmom
                        ang_mom[kp-1+i,m] = zmid
                        pv_mid[kp-1+i,m] = pvtop + qscal*zmid
                        mass_weight[kp-1+i,m] = earth_area*(mextra[i]-mextra[i+1])*d_th[m]
                        if np.min(mass_weight) < 0:
                            None
                            #raise ValueError('negative mass location 1')

            # Find repeated points in (M, theta) space.
            # Occurs where same mass and circ assigned to a range of PV values.
            elz = np.argwhere((d_circ == 0) & (c_mid != 0))

            num_elz = elz.shape
            if num_elz[0] > 1:
                # Where this occurs, sum the mass
                sum_mass_weight = np.sum(mass_weight[elz,m])            
                mass_weight[elz,m] = 0
                mass_weight[elz[0],m] = sum_mass_weight
                
                if np.min(mass_weight) < 0:
                    None
                    #raise ValueError('negative mass location 2')

        mass_check = np.sum(mass_weight)
        #print('Total mass of point masses: ', mass_check)

        # potential temperature
        th = np.matlib.repmat(th_lev,num_pv_lev+nextra-1,1) # ((num_pv_lev+nextra-1),num_th_lev) numpy array of potential temperature values

        # ravel the data
        pv_mid = np.ravel(pv_mid)
        ang_mom = np.ravel(ang_mom)
        th = np.ravel(th)
        mass_weight = np.ravel(mass_weight)

        # find indices of seeds with zero angular momentum
        idx_mom = ang_mom == 0

        # find indices of seeds corresponding to zero target masses
        idx_mass = mass_weight == 0

        # delete seeds with zero momentum or zero mass
        idx = idx_mom | idx_mass
        pv_mid = np.delete(pv_mid,idx)
        ang_mom = np.delete(ang_mom,idx)
        th = np.delete(th,idx)
        mass_weight = np.delete(mass_weight,idx)

        # create seed and mass arrays
        y = np.append(ang_mom[:,None],th[:,None],1)
        tm = mass_weight

        # eliminate duplicate seeds
        _, i = np.unique(y,axis = 0,return_index = True)
        i = np.sort(i)
        y = y[i]
        tm = tm[i]
        
        # normalise the masses
        tmn = tm / np.sum(tm) * (pp.p00 - self.pmin) * (self.smax - self.smin)
        
        # assign seeds and masses to the class instance
        self.epv = pv_mid
        self.y = y
        self.tm = tm
        self.tmn = tmn
        
        return y, tm, tmn

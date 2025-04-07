import pytest
import numpy as np
import pathlib
import atmosphere_bgs

@pytest.mark.parametrize("ot_tol", [1e-1, 1e-2, 1e-3, 1e-4, 1e-5])
def test_precision_n20(ot_tol):
    data = atmosphere_bgs.DataLoader(pathlib.Path(__file__).parent.resolve() / "../data/bs_lc1low2001010000", pmin=1000, nextra=17)

    # restrict to subsample
    np.random.seed(42)
    use = np.random.choice(np.arange(data.y.shape[0]), 20, replace=False)
    data.y = data.y[use,:]
    data.tm = data.tm[use]
    data.tmn = data.tmn[use] * np.sum(data.tmn) / np.sum(data.tmn[use])

    solv = atmosphere_bgs.OTSolver(data, ot_tol=ot_tol, boundary_res=100)

    # set max_line_resolution very high so it is not hit
    solv.sp.max_line_resolution = 1000000

    solv.get_bgs()

    # relative area error for discretization
    assert(np.max(solv.ld.areaerrs / solv.ld.areas) <= solv.sp.area_tolerance)

    # areas calculated with high precision
    A = np.array([1279.5947346549979, 18613.640026409474 ,  1340.4692988341799,
                  1340.288465098682 ,  2977.7059503186906, 18666.43084603698  ,
                  2543.338176998129 ,  1827.9896267266265,   824.1960985971008,
                  1238.117775420699 , 10868.314017967821 ,  5277.591663615872 ,
                  1257.34701529699  ,  5054.471261028026 ,  1669.126053508578 ,
                  2238.0370784916854, 16237.558821434208 ,  1582.7569462733984,
                  4487.674091413537 ,   691.7534464234041])
    
    # areas are within absolute error
    assert(np.max((solv.ld.areas - A) / solv.tmn) < ot_tol)


@pytest.mark.parametrize("ot_tol, n", [(1e-1, 20), (1e-2, 20), (1e-3, 20), (1e-4, 20), (1e-5, 20),
                                       (1e-1, 50), (1e-2, 30), (1e-3, 20), (1e-4, 10), (1e-5, 5)])
def test_precision_bootstrap(ot_tol, n):
    data = atmosphere_bgs.DataLoader(pathlib.Path(__file__).parent.resolve() / "../data/bs_lc1low2001010000", pmin=1000, nextra=17)

    # restrict to subsample
    np.random.seed(42)
    use = np.random.choice(np.arange(data.y.shape[0]), 20, replace=False)
    data.y = data.y[use,:]
    data.tm = data.tm[use]
    data.tmn = data.tmn[use] * np.sum(data.tmn) / np.sum(data.tmn[use])

    # solve with increased precision
    solv_prec = atmosphere_bgs.OTSolver(data, ot_tol = ot_tol * 1e-3, boundary_res=100)
    solv_prec.sp.max_line_resolution = 1000000
    solv_prec.get_bgs()
    assert(np.max(solv_prec.ld.areaerrs / solv_prec.ld.areas) <= solv_prec.sp.area_tolerance) # relative area error for discretization

    solv = atmosphere_bgs.OTSolver(data, ot_tol=ot_tol, boundary_res=100)
    solv.sp.max_line_resolution = 1000000
    solv.get_bgs()
    assert(np.max(solv.ld.areaerrs / solv.ld.areas) <= solv.sp.area_tolerance)

    # areas are within absolute error
    assert(np.max((solv.ld.areas - solv_prec.ld.areas) / solv.tmn) < ot_tol * (1 + 1e-3))

import pytest
import numpy as np
import atmosphere_bgs

@pytest.mark.parametrize("ot_tol", [1e-1, 1e-2, 1e-3, 1e-4, 1e-5])
def test_precision_n20(ot_tol):
    data = atmosphere_bgs.DataLoader("data/bs_lc1low2001010000", pmin=1000, nextra=17)

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
    A = np.array([3.5025063175828626e+02, 1.6230217936286797e+03, 4.0725270099405101e+01, 3.1512841352633779e+02,
                  8.2335361016268166e+01, 4.4443558239278502e+01, 9.7776027784522594e+03, 2.3411595834304197e+01,
                  1.9104743250731888e+03, 5.5638682120010890e+02, 5.0317779183617622e+03, 2.8781630836723963e+01,
                  3.8418502961235863e+03, 1.8627818396485368e+03, 6.8757831016173775e+04, 4.6147320764476053e+03,
                  1.9486908652979946e+02, 8.0612078753658295e+02, 5.7484628186238766e+01, 9.6391568210418370e+01])

    # areas are within absolute error
    assert(np.max((solv.ld.areas - A) / solv.tmn) < ot_tol)


@pytest.mark.parametrize("ot_tol, n", [(1e-1, 20), (1e-2, 20), (1e-3, 20), (1e-4, 20), (1e-5, 20),
                                       (1e-1, 50), (1e-2, 30), (1e-3, 20), (1e-4, 10), (1e-5, 5)])
def test_precision_bootstrap(ot_tol, n):
    data = atmosphere_bgs.DataLoader("data/bs_lc1low2001010000", pmin=1000, nextra=17)

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

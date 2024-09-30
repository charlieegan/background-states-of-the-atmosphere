import pytest
import numpy as np
import itertools
from _atmosphere_bgs import Rasterizer

def meanpool(a, ps):
    res = np.zeros((a.shape[0] // ps[0], a.shape[1] // ps[1]))
    J, I = np.meshgrid(np.repeat(np.arange(res.shape[1]), ps[1]), np.repeat(np.arange(res.shape[0]), ps[0]))
    print(np.repeat(res.shape[1], ps[1]))
    np.add.at(res, (I.ravel(), J.ravel()), a.ravel())
    res /= ps[0] * ps[1]
    return res

def get_start(seg, x0):
    return [i for i in range(len(seg)) if seg[i][0][0] == x0 or seg[i][1][0] == x0]

@pytest.fixture
def geometry_simple():
    res = {}
    res["val"] = np.array([1, 2, 3, 0])
    res["seg"] = [((0,0), (4,0), -1), 
                  ((4,0), (10,0), -1), 
                  ((4,0), (6,2), 1), 
                  ((6,2), (8,1), 1), 
                  ((8,1), (10,3), 1), 
                  ((4,6), (6,2), 0), 
                  ((0,7), (4,6), 0), 
                  ((4,6), (8,7), 2), 
                  ((8,7), (10,7), 2), 
                  ((0,10), (10,10), 3)]
    res["bounds"] = np.array([0, 10, 0, 10])
    return res

@pytest.fixture
def geometry_simple_10x10_truth():
    return np.array([[1.  , 1.  , 1.  , 1.  , 1.  , 1.  , 0.875, 0.  , 0.  , 0.  ],
                     [1.  , 1.  , 1.  , 1.  , 1.  , 1.  , 0.625, 0.  , 0.  , 0.  ],
                     [1.  , 1.  , 1.  , 1.  , 1.  , 1.  , 0.375, 0.  , 0.  , 0.  ],
                     [1.  , 1.  , 1.  , 1.  , 1.  , 1.  , 0.125, 0.  , 0.  , 0.  ],
                     [1.5 , 1.  , 1.  , 1.  , 1.5 , 2.5 , 0.375, 0.  , 0.  , 0.  ],
                     [2.  , 1.5 , 1.5 , 2.5 , 3.  , 3.  , 1.125, 0.  , 0.  , 0.  ],
                     [2.  , 2.25, 3.  , 3.  , 3.  , 3.  , 1.875, 0.  , 0.  , 0.  ],
	             [2.  , 2.75, 3.  , 3.  , 3.  , 3.  , 2.625, 0.  , 0.  , 0.  ],
                     [2.  , 2.5 , 3.  , 3.  , 3.  , 3.  , 3.   , 0.  , 0.  , 0.  ],
                     [2.  , 2.  , 2.5 , 3.  , 3.  , 3.  , 3.   , 0.  , 0.  , 0.  ]])

def test_10x10_simple(geometry_simple, geometry_simple_10x10_truth):
    res = [10,10]
    rast = Rasterizer(seg = geometry_simple["seg"],
                      bounds = geometry_simple["bounds"]);
    v = rast.rasterize(val = geometry_simple["val"], res = res)
    assert np.max(np.abs(v - geometry_simple_10x10_truth)) < 1e-9

def test_consecutive_simple(geometry_simple, geometry_simple_10x10_truth):
    res = [10,10]
    rast = Rasterizer(seg = geometry_simple["seg"],
                      bounds = geometry_simple["bounds"]);

    v = rast.rasterize(val = geometry_simple["val"], res = res)
    assert np.max(np.abs(v - geometry_simple_10x10_truth)) < 1e-9
    
    v = rast.rasterize(val = geometry_simple["val"], res = res)
    assert np.max(np.abs(v - geometry_simple_10x10_truth)) < 1e-9

@pytest.mark.parametrize("seed", np.arange(50))
def test_segment_direction_simple(geometry_simple, geometry_simple_10x10_truth, seed):
    np.random.seed(seed)
    seg = [(s[1], s[0], s[2]) if np.random.random() > 0.5 else s for s in geometry_simple["seg"]]
    res = [10, 10]
    rast = Rasterizer(seg = seg,
                      bounds = geometry_simple["bounds"]);
    v = rast.rasterize(val = geometry_simple["val"], res = res)
    assert np.max(np.abs(v - geometry_simple_10x10_truth)) < 1e-9

@pytest.mark.parametrize("rep", list(itertools.product(range(1,5), range(1,5))))
def test_downsample_simple(geometry_simple, geometry_simple_10x10_truth, rep):
    res = [rep[0] * 10, rep[1] * 10]
    rast = Rasterizer(seg = geometry_simple["seg"],
                      bounds = geometry_simple["bounds"]);
    v = rast.rasterize(val = geometry_simple["val"], res = res)
    assert np.max(np.abs(meanpool(v, rep) - geometry_simple_10x10_truth)) < 1e-9

def test_consecutive_different_values_simple(geometry_simple, geometry_simple_10x10_truth):
    res = [10,10]
    rast = Rasterizer(seg = geometry_simple["seg"],
                      bounds = geometry_simple["bounds"]);

    v = rast.rasterize(val = geometry_simple["val"], res = res)
    assert np.max(np.abs(v - geometry_simple_10x10_truth)) < 1e-9
    
    v = rast.rasterize(val = 2 * geometry_simple["val"], res = res)
    print(v)
    assert np.max(np.abs(v - 2 * geometry_simple_10x10_truth)) < 1e-9

def test_consecutive_different_sizes_simple(geometry_simple, geometry_simple_10x10_truth):
    res = [10,10]
    rast = Rasterizer(seg = geometry_simple["seg"],
                      bounds = geometry_simple["bounds"]);

    v = rast.rasterize(val = geometry_simple["val"], res = res)
    assert np.max(np.abs(v - geometry_simple_10x10_truth)) < 1e-9
    
    v = rast.rasterize(val = geometry_simple["val"], res = 2 * np.array(res))
    assert np.max(np.abs(meanpool(v, (2,2)) - geometry_simple_10x10_truth)) < 1e-9
    
def test_values_one_by_one_simple(geometry_simple, geometry_simple_10x10_truth):
    res = [10,10]
    rast = Rasterizer(seg = geometry_simple["seg"],
                      bounds = geometry_simple["bounds"]);
    v = np.zeros(res);
    
    for i in range(len(geometry_simple["val"])):
        val = np.zeros(len(geometry_simple["val"]))
        val[i] = geometry_simple["val"][i]
        v += rast.rasterize(val = val, res = res)
        
    assert np.max(np.abs(v - geometry_simple_10x10_truth)) < 1e-9

def test_ending_2():
    seg = [((0,0), (2,0), -1),
           ((0,0), (1,1), 1),
           ((0,2), (1,1), 0),
           ((0,2), (2,2), 1)]
    bounds = [0, 2, 0, 2]

    val = [0, 1]
    res = [2, 2]
    rast = Rasterizer(seg = seg, bounds = bounds)
    v = rast.rasterize(val = val, res = res)

    tv = np.array([[0.5, 0.5],
                   [1.0, 1.0]])
    assert np.max(np.abs(v - tv)) < 1e-9

def test_ending_3():
    seg = [((0,0), (2,0), -1),
           ((0,0), (1,1), 2),
           ((0,1), (1,1), 0),
           ((0,2), (1,1), 1),
           ((0,2), (2,2), 2)]
    bounds = [0, 2, 0, 2]

    val = [0, 1, 2]
    res = [2, 2]
    rast = Rasterizer(seg = seg, bounds = bounds)
    v = rast.rasterize(val = val, res = res)

    tv = np.array([[1.0, 1.5],
                   [2.0, 2.0]])
    assert np.max(np.abs(v - tv)) < 1e-9

def test_double_vertical():
    seg = [((0,0), (1,0), -1),
           ((1,0), (2,0), -1),
           ((1,0), (1,1), 1),
           ((1,2), (1,1), 1),
           ((0,2), (1,2), 0),
           ((1,2), (2,2), 1)]
    bounds = [0, 2, 0, 2]

    val = [0, 1]
    res = [2, 2]
    rast = Rasterizer(seg = seg, bounds = bounds)
    v = rast.rasterize(val = val, res = res)

    tv = np.array([[0.0, 0.0],
                   [1.0, 1.0]])
    assert np.max(np.abs(v - tv)) < 1e-9

def test_starting_2():
    seg = [((0,0), (2,0), -1),
           ((2,0), (1,1), 0),
           ((2,2), (1,1), 1),
           ((0,2), (2,2), 0)]
    bounds = [0, 2, 0, 2]

    val = [0, 1]
    res = [2, 2]
    rast = Rasterizer(seg = seg, bounds = bounds)
    v = rast.rasterize(val = val, res = res)

    tv = np.array([[0.0, 0.0],
                   [0.5, 0.5]])
    assert np.max(np.abs(v - tv)) < 1e-9

def test_simple_intersect():
    seg = [((0,0), (1,0), -1),
           ((0,0), (1,1), 0),
           ((0,1), (1,0), 1),
           ((0,1), (1,1), 2)]
    bounds = [0, 1, 0, 1]
    val = [0, 1, 2]
    res = [4, 4]
    rast = Rasterizer(seg = seg, bounds = bounds)
    v = rast.rasterize(val = val, res = res)

    tv = np.array([[0.5, 1. , 1. , 1.5],
                   [0. , 0.5, 1.5, 2. ],
                   [1. , 0.5, 1. , 2. ],
                   [0.5, 0. , 0. , 1. ]])
    assert np.max(np.abs(v - tv)) < 1e-9
    

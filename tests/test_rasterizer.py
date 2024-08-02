import pytest
import numpy as np
import itertools
from _atmosphere_bgs import Rasterizer

@pytest.fixture
def geometry_simple():
    res = {}
    vals = [1, 2, 3, 0]
    res["seg"] = [((0,0), (4,0), np.nan), 
                  ((4,0), (10,0), np.nan), 
                  ((4,0), (6,2), vals[1]), 
                  ((6,2), (8,1), vals[1]), 
                  ((8,1), (10,3), vals[1]), 
                  ((4,6), (6,2), vals[0]), 
                  ((0,7), (4,6), vals[0]), 
                  ((4,6), (8,7), vals[2]), 
                  ((8,7), (10,7), vals[2]), 
                  ((0,10), (10,10), vals[3])]
    res["con"] = [[0,1,2], [2,3,5], [3,4], [6,5,7], [7,8]]
    res["start"] = [0, 6, 9]
    res["bounds"] = [0, 10, 0, 10]
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
    rast = Rasterizer(**geometry_simple, res=res);
    assert np.max(np.abs(rast.out - geometry_simple_10x10_truth)) < 1e-9

@pytest.mark.parametrize("seed", np.arange(50))
def test_segment_direction_simple(geometry_simple, geometry_simple_10x10_truth, seed):
    np.random.seed(seed)
    seg = [(s[1], s[0], s[2]) if np.random.random() > 0.5 else s for s in geometry_simple["seg"]]
    con = [np.array(c)[np.random.permutation(len(c))] for c in geometry_simple["con"]]
    res = [10, 10]
    rast = Rasterizer(seg = seg, con = con, start = geometry_simple["start"], bounds = geometry_simple["bounds"], res=res)
    assert np.max(np.abs(rast.out - geometry_simple_10x10_truth)) < 1e-9

def meanpool(a, ps):
    res = np.zeros((a.shape[0] // ps[0], a.shape[1] // ps[1]))
    J, I = np.meshgrid(np.repeat(np.arange(res.shape[1]), ps[1]), np.repeat(np.arange(res.shape[0]), ps[0]))
    print(np.repeat(res.shape[1], ps[1]))
    np.add.at(res, (I.ravel(), J.ravel()), a.ravel())
    res /= ps[0] * ps[1]
    return res

@pytest.mark.parametrize("rep", list(itertools.product(range(1,5), range(1,5))))
def test_downsample_simple(geometry_simple, geometry_simple_10x10_truth, rep):
    res = [rep[0] * 10, rep[1] * 10]
    rast = Rasterizer(**geometry_simple, res=res)
    assert np.max(np.abs(meanpool(rast.out, rep) - geometry_simple_10x10_truth)) < 1e-9



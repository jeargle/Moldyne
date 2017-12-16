# John Eargle
# 2017


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits import axes_grid1
from mpl_toolkits.mplot3d import Axes3D
import scipy as sp
from scipy import misc
from scipy.integrate import simps, quad, nquad, tplquad
from scipy.special import sph_harm, eval_genlaguerre, expi

import random


__all__ = ["markov_pi"]

def markov_pi(N, delta):
    """
    Reject the move if it falls outside the unit square.
    Count a hit if it falls within the unit circle.
    """
    x, y = 1.0, 1.0
    n_hits = 0
    n_accepts = 0
    for i in range(N):
        del_x, del_y = random.uniform(-delta, delta), random.uniform(-delta, delta)
        if abs(x + del_x) < 1.0 and abs(y + del_y) < 1.0:
            x, y = x + del_x, y + del_y
            n_accepts += 1
        if x**2 + y**2 < 1.0:
            n_hits += 1
    return n_hits, n_accepts

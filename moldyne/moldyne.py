# John Eargle
# 2017

import math
import random

import numpy as np
import pylab
import scipy as sp
from scipy.special import gamma


__all__ = ["markov_pi", "markov_pi_all_data", "markov_pi_all_data2",
           "direct_disks_box"]

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

def markov_pi_all_data(N, delta):
    """
    """
    x, y = 1.0, 1.0
    data_sum = 0.0
    data_sum_sq = 0.0
    for i in range(N):
        del_x, del_y = random.uniform(-delta, delta), random.uniform(-delta, delta)
        if abs(x + del_x) < 1.0 and abs(y + del_y) < 1.0:
            x, y = x + del_x, y + del_y
        if x ** 2 + y ** 2 < 1.0:
            data_sum += 4.0
            data_sum_sq += 4.0 ** 2
    return data_sum / float(N), data_sum_sq / float(N)

def markov_pi_all_data2(N, delta):
    """
    """
    x, y = 1.0, 1.0
    data = []
    for i in range(N):
        del_x, del_y = random.uniform(-delta, delta), random.uniform(-delta, delta)
        if abs(x + del_x) < 1.0 and abs(y + del_y) < 1.0:
            x, y = x + del_x, y + del_y
        if x**2 + y**2 < 1.0:
            data.append(4.0)
        else:
            data.append(0.0)
    return data

def direct_disks_box(N, sigma):
    """
    Place N non-overlapping disks in a box.
    N: number of disks to place
    sigma: disk radius
    """
    overlap = True
    while overlap == True:
        L = [(random.uniform(sigma, 1.0 - sigma), random.uniform(sigma, 1.0 - sigma))]
        for k in range(1, N):
            a = (random.uniform(sigma, 1.0 - sigma), random.uniform(sigma, 1.0 - sigma))
            min_dist_sq = min(((a[0] - b[0])**2 + (a[1] - b[1])**2) for b in L)
            if min_dist_sq < 4.0 * sigma**2:
                overlap = True
                break
            else:
                overlap = False
                L.append(a)
    return L

def direct_disks_box2(N, sigma):
    condition = False
    att = 0
    acc = 0
    while condition == False:
        L = [(random.uniform(sigma, 1.0 - sigma), random.uniform(sigma, 1.0 - sigma))]
        for k in range(1, N):
            att += 1
            a = (random.uniform(sigma, 1.0 - sigma), random.uniform(sigma, 1.0 - sigma))
            min_dist = min(np.sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2) for b in L)
            if min_dist < 2.0 * sigma:
                condition = False
                break
            else:
                L.append(a)
                acc += 1
                condition = True
    return att, acc, L

def markov_disks_box(L, sigma, delta):
    """
    Place one disk in a box with existing disks.  The placement
    is a jump from one of the existing disks.
    L: array of existing disk locations
    sigma: disk radius
    delta: max jump distance
    """
    sigma_sq = sigma**2
    accepted = False
    a = random.choice(L)
    b = [a[0] + random.uniform(-delta, delta), a[1] + random.uniform(-delta, delta)]
    min_dist = min((b[0] - c[0])**2 + (b[1] - c[1])**2 for c in L if c != a)
    box_cond = min(b[0], b[1]) < sigma or max(b[0], b[1]) > 1.0 - sigma

    if not (box_cond or min_dist < 4.0 * sigma_sq):
        a[:] = b
        accepted = True

    return accepted

def wall_time(pos_a, vel_a, sigma):
    """
    Determine amount of time before disk motion in 1D hits a wall.
    """
    if vel_a > 0.0:
        del_t = (1.0 - sigma - pos_a) / vel_a
    elif vel_a < 0.0:
        del_t = (pos_a - sigma) / abs(vel_a)
    else:
        del_t = float('inf')
    return del_t

def pair_time(pos_a, vel_a, pos_b, vel_b, sigma, sigma_sq):
    """
    Determine amount of time disk motions in 2D cause a pair of disks
    to hit each other.
    """
    del_x = [pos_b[0] - pos_a[0], pos_b[1] - pos_a[1]]
    del_x_sq = del_x[0]**2 + del_x[1]**2
    del_v = [vel_b[0] - vel_a[0], vel_b[1] - vel_a[1]]
    del_v_sq = del_v[0]**2 + del_v[1]**2
    scal = del_v[0] * del_x[0] + del_v[1] * del_x[1]
    upsilon = scal**2 - del_v_sq * (del_x_sq - 4.0 * sigma_sq)
    if upsilon > 0.0 and scal < 0.0:
        del_t = - (scal + np.sqrt(upsilon)) / del_v_sq
    else:
        del_t = float('inf')
    return del_t

def disk_dist(x, y):
    """
    Distance between two disks in a 1x1 square with periodic boundary
    conditions.
    """
    d_x = abs(x[0] - y[0]) % 1.0
    d_x = min(d_x, 1.0 - d_x)
    d_y = abs(x[1] - y[1]) % 1.0
    d_y = min(d_y, 1.0 - d_y)

    return  np.sqrt(d_x**2 + d_y**2)

def phi6(phi):
    """
    """
    return (sum([np.exp(6.0j * ( phi + (x * np.pi/3.0)))
                 for x in range(6)])
            / 6.0)

def V(x, cubic, quartic):
    """
    """
    return (x**2 / 2) + (cubic * x**3) + (quartic * x**4)

def Energy(n, cubic, quartic):
    """
    """
    return n + 0.5 - (15.0 / 4.0 * cubic**2 * (n**2 + n + 11.0/30.0)) \
         + (3.0 / 2.0 * quartic * (n**2 + n + 1.0/2.0))

def Z(cubic, quartic, beta, n_max):
    """
    """
    return sum(np.exp(-beta * Energy(n, cubic, quartic))
               for n in range(n_max + 1))

def Z2(beta):
    """
    """
    return 1.0/(2*math.sinh(beta/2.0))



# ====================
# Sphere Volume
# ====================

def direct_pi(N, dimension):
    """
    Direct randomized calculation of PI.
    """
    n_hits = 0
    for i in range(N):
        x_len = 0.0
        for j in range(dimension):
            x = random.uniform(-1.0, 1.0)
            x_len += x**2
            if x_len > 1.0:
                break
        if x_len <= 1.0:
            n_hits += 1
    return n_hits

def sphere_volume(dimension):
    """
    Volume of N-dimensional sphere.
    """
    return np.pi**(dimension / 2.0) / gamma(dimension / 2.0 + 1.0)

def sample_sphere(n_trials, d):
    """
    Markov chain sampling of sphere volume.
    Each trial adjusts value of one dimension chosen at random.
    n_trials : number of trials
    d: sphere dimension
    """
    x = [0] * d
    r_sqs = []
    for i in range(n_trials):
        r_sq = sum([j**2 for j in x])
        k = random.randint(0, d - 1)
        x_new_k = x[k] + random.uniform(-1.0, 1.0)
        # Sum of squares including x_new_k
        r_sq_new = r_sq - x[k]**2 + x_new_k**2
        # Accept move if point is within unit sphere.
        # sqrt not needed since radius bound is 1.0
        if r_sq_new < 1.0:
            x[k] = x_new_k
            r_sqs.append(r_sq_new)
            # n_hits += 1
        else:
            r_sqs.append(r_sq)
        # points.append(x[:])

    return r_sqs

def sample_cylinder_old(n_trials, d):
    """
    Markov sampling of sphere within a cylinder.  Moves are accepted
    if they lie within the cylinder and counted towards the sphere
    volume calculation if they also lie within the sphere.  Dimension
    d represents cylinder height.
    n_trials : number of Markov trials
    d : dimension of hypercylinder
    """
    x = [0] * (d+1)
    n_hits = 0
    r_sqs = []
    for i in range(n_trials):
        r_sq = sum([j**2 for j in x])
        k = random.randint(0, d - 1)
        x_new_k = x[k] + random.uniform(-1.0, 1.0)
        x_supp = random.uniform(-1.0, 1.0)
        # Sum of squares including x_new_k but not x[d], the height
        # dimension.
        r_sq_new = r_sq - x[k]**2 + x_new_k**2 - x[d]**2
        # Check that position is within cylinder.
        if np.sqrt(r_sq_new) < 1.0:
            x[k] = x_new_k
            x[d] = x_supp
            r_sq = r_sq_new + x[d]**2
        # Check that cylinder position is within sphere.
        if r_sq < 1.0:
            n_hits += 1
        # r_sqs.append(r_sq)
        # points.append(x[:])

    return n_hits

def sample_cylinder(n_trials, d):
    """
    Markov sampling of sphere within a cylinder.  Moves are accepted
    if they lie within the cylinder and counted towards the sphere
    volume calculation if they also lie within the sphere.  Dimension
    d represents cylinder height.
    n_trials : number of Markov trials
    d : dimension of hypercylinder
    """
    x = [0] * (d+1)
    n_hits = 0
    r_sq = 0.0
    r_sq_new = 0.0
    r_sqC = 0.0
    for i in range(n_trials):
        # if i%100 == 0:
        #     r_sq2 = sum([j**2 for j in x[:-1]])
        #     print('|%f - %f| = %f' % (r_sq, r_sq2, abs(r_sq-r_sq2)))
        k = random.randint(0, d - 1)
        x_new_k = x[k] + random.uniform(-1.0, 1.0)
        x_supp = random.uniform(-1.0, 1.0)

        # Sum of squares including x_new_k
        r_sq_new = r_sq - x[k]**2 + x_new_k**2
        if np.sqrt(r_sq_new) < 1.0:
            x[k] = x_new_k
            x[d] = x_supp
            r_sq = r_sq_new

        # Sum of squares including the cylinder x_supp term
        r_sqC = r_sq + x[d]**2
        if r_sqC < 1.0:
            n_hits += 1

    return n_hits



# ====================
# Plotting
# ====================

def show_conf(L, sigma, title, fname=None):
    """
    """
    pylab.axes()
    for [x, y] in L:
        for ix in range(-1, 2):
            for iy in range(-1, 2):
                cir = pylab.Circle((x + ix, y + iy), radius = sigma,  fc = 'r')
                pylab.gca().add_patch(cir)
    pylab.axis('scaled')
    pylab.title(title)
    pylab.axis([0.0, 1.0, 0.0, 1.0])
    if fname is not None:
        pylab.savefig(fname)
    pylab.show()

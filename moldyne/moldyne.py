# John Eargle
# 2017

import random

import numpy as np


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
    min_dist = min((b[0] - c[0]) ** 2 + (b[1] - c[1]) ** 2 for c in L if c != a)
    box_cond = min(b[0], b[1]) < sigma or max(b[0], b[1]) > 1.0 - sigma

    if not (box_cond or min_dist < 4.0 * sigma_sq):
        a[:] = b
        accepted = True

    return accepted

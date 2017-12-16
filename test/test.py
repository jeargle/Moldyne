# John Eargle
# 2017


from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import pylab
from scipy.integrate import simps, quad, nquad
from scipy.optimize import minimize

import moldyne.moldyne as md


def markov_pi_test1():
    """
    """
    # n_runs = 1000
    # n_trials = 1000
    n_runs = 100
    n_trials = 100
    deltas = [(x+1)/10.0 for x in range(50)]
    acceptance_ratios = []

    for delta in deltas:
        ars = 0.0
        for run in range(n_runs):
            n_hits, n_accepts = md.markov_pi(n_trials, delta)
            ars +=  n_accepts / float(n_trials)
        acceptance_ratios.append(ars / n_runs)
        print('%.1f, %f' % (delta, ars/n_runs))

    pylab.plot(deltas, acceptance_ratios, 'o')
    # pylab.gca().set_xscale('log')
    # pylab.gca().set_yscale('log')
    pylab.xlabel('$\delta$')
    pylab.ylabel('acceptance ratio')
    pylab.title('1/2 Rule: Acceptance ratio as a function of $\delta$')
    # pylab.savefig('acceptance_ratio.png')
    pylab.show()


def markov_pi_test2():
    """
    """
    # n_runs = 1000
    # n_trials = 1000
    n_runs = 100
    n_trials = 100
    deltas = [(x+1)/10.0 for x in range(50)]
    sigmas = []

    for delta in deltas:
        sigma = 0.0
        for run in range(n_runs):
            n_hits, n_accepts = md.markov_pi(n_trials, delta)
            pi_est = 4.0 * n_hits / float(n_trials)
            sigma += (pi_est - np.pi)**2
        sigmas.append(np.sqrt(sigma / n_runs))
        print('%.1f, %f' % (delta, np.sqrt(sigma/n_runs)))


    pylab.plot(deltas, sigmas, 'o')
    # pylab.gca().set_xscale('log')
    # pylab.gca().set_yscale('log')
    pylab.xlabel('$\delta$')
    pylab.ylabel('$\sigma$')
    pylab.title('Performance: Standard deviation $\sigma$ as a function of $\delta$')
    # pylab.title('1/2 Rule: Acceptance ratio as a function of $\delta$')
    # pylab.savefig('performance.png')
    pylab.show()



if __name__=='__main__':

    print('*********************')
    print('*** MOLDYNE TESTS ***')
    print('*********************')

    # ====================
    # Markov Pi tests
    # ====================

    markov_pi_test1()
    markov_pi_test2()

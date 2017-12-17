# John Eargle
# 2017


from __future__ import print_function

import numpy as np
import pylab

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


def markov_pi_test3():
    """
    """
    n_trials = 2 ** 14
    delta = 0.1
    n_parties = 100
    inside_error_bar = 0
    for iteration in range(n_parties):
        mean, mean_square = md.markov_pi_all_data(n_trials, delta)
        naive_error = np.sqrt(mean_square  - mean ** 2) / np.sqrt(n_trials)
        error =  abs(mean - np.pi)
        if error < naive_error: inside_error_bar += 1
        print(mean, error, naive_error)
    print(inside_error_bar / float(n_parties), 'fraction: error bar including pi')



if __name__=='__main__':

    print('*********************')
    print('*** MOLDYNE TESTS ***')
    print('*********************')

    # ====================
    # Markov Pi tests
    # ====================

    # markov_pi_test1()
    # markov_pi_test2()
    markov_pi_test3()

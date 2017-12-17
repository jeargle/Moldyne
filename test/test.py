# John Eargle
# 2017


from __future__ import print_function

import numpy as np
import pylab

import moldyne.moldyne as md


def markov_pi_test1():
    """
    Trial acceptance ratio
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
    Accuracy
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
    Error
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


def markov_pi_test4():
    """
    Bunching
    """
    power = 14
    n_trials = 2 ** power
    delta = 0.1
    data = md.markov_pi_all_data2(n_trials, delta)
    errors  = []
    bunches = []

    for i in range(power):
        new_data = []
        mean = 0.0
        mean_sq = 0.0
        N = len(data)
        while data != []:
            x = data.pop()
            y = data.pop()
            mean += x + y
            mean_sq += x ** 2 + y ** 2
            new_data.append((x + y) / 2.0 )
        errors.append(np.sqrt(mean_sq / N - (mean / N) ** 2) / np.sqrt(N))
        bunches.append(i)
        data = new_data[:]

    pylab.plot(bunches, errors, 'o')
    pylab.xlabel('iteration')
    pylab.ylabel('naive error')
    pylab.title('Bunching: naive error vs iteration number')
    pylab.savefig('apparent_error_bunching.png', format='PNG')
    pylab.show()



if __name__=='__main__':

    print('*********************')
    print('*** MOLDYNE TESTS ***')
    print('*********************')

    # ====================
    # Markov Pi tests
    # ====================

    # markov_pi_test1()
    # markov_pi_test2()
    # markov_pi_test3()
    markov_pi_test4()

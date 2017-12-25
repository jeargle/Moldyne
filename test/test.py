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


def disk_test1():
    """
    Test direct placement of 4 disks in a box.
    """
    N = 4
    sigma = 0.1196
    # n_runs = 1000000
    n_runs = 100000
    histo_data = []
    for run in range(n_runs):
        pos = md.direct_disks_box(N, sigma)
        for k in range(N): histo_data.append(pos[k][0])

    pylab.hist(histo_data, bins=100, normed=True)
    pylab.xlabel('x')
    pylab.ylabel('frequency')
    pylab.title('x-coordinates for 1e6 runs of direct_disks_box\nwith 4 disks of radius 0.1196')
    pylab.grid()
    #pylab.savefig('direct_disks_histo.png')
    pylab.show()


def disk_test2():
    """
    Test markov placement of 1 disk in a box with 4 existing disks.
    """
    L = [[0.25, 0.25], [0.75, 0.25], [0.25, 0.75], [0.75, 0.75]]
    sigma = 0.1196
    delta = 0.18   # 0.5 acceptance ratio

    # n_steps = 2000000
    n_steps = 20000
    accept_count = 0
    histo_data = []
    for steps in range(n_steps):
        acc = md.markov_disks_box(L, sigma, delta)
        if acc:
            accept_count += 1
        for k in L: histo_data.append(k[0])

    print('  acceptance ratio: %f' % (1.0*accept_count/n_steps))

    pylab.hist(histo_data, bins=100, normed=True)
    pylab.xlabel('x')
    pylab.ylabel('frequency')
    pylab.title('x-coordinates for 2e6 runs of markov_disks_box\nwith 4 disks of radius 0.1196 and $\delta$=0.18')
    pylab.grid()
    # pylab.savefig('markov_disks_histo.png')
    pylab.show()


def disk_test3():
    """
    Molecular dynamics of four disks in a box.
    """
    pos = [[0.25, 0.25],
           [0.75, 0.25],
           [0.25, 0.75],
           [0.75, 0.75]]
    vel = [[0.21, 0.12],
           [0.71, 0.18],
           [-0.23, -0.79],
           [0.78, 0.1177]]
    singles = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (3, 0), (3, 1)]
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    sigma = 0.1196
    sigma_sq = sigma**2
    t = 0.0
    n_events = 200000
    histo_data = []

    for event in range(n_events):
        wall_times = [md.wall_time(pos[k][l], vel[k][l], sigma)
                      for k, l  in singles]
        pair_times = [md.pair_time(pos[k], vel[k], pos[l], vel[l], sigma, sigma_sq)
                      for k, l in pairs]
        next_event = min(wall_times + pair_times)
        t += next_event
        for k, l in singles:
            pos[k][l] += vel[k][l] * next_event

        if min(wall_times) < min(pair_times):
            collision_disk, direction = singles[wall_times.index(next_event)]
            vel[collision_disk][direction] *= -1.0
        else:
            a, b = pairs[pair_times.index(next_event)]
            del_x = [pos[b][0] - pos[a][0], pos[b][1] - pos[a][1]]
            abs_x = np.sqrt(del_x[0] ** 2 + del_x[1] ** 2)
            e_perp = [c / abs_x for c in del_x]
            del_v = [vel[b][0] - vel[a][0], vel[b][1] - vel[a][1]]
            scal = del_v[0] * e_perp[0] + del_v[1] * e_perp[1]
            for k in range(2):
                vel[a][k] += e_perp[k] * scal
                vel[b][k] -= e_perp[k] * scal

        print('event', event)
        print('time', t)
        # print('pos', pos)
        # print('vel', vel)
        for k in pos: histo_data.append(k[0])

    pylab.hist(histo_data, bins=100, normed=True)
    pylab.xlabel('x')
    pylab.ylabel('frequency')
    pylab.title('x-coordinates for 2e5 events of event_disks_box\nwith 4 disks of radius 0.1196')
    pylab.grid()
    # pylab.savefig('event_disks_histo.png')
    pylab.show()


def disk_test4():
    """
    Molecular dynamics of four disks in a box.
    """
    pos = [[0.25, 0.25],
           [0.75, 0.25],
           [0.25, 0.75],
           [0.75, 0.75]]
    vel = [[0.21, 0.12],
           [0.71, 0.18],
           [-0.23, -0.79],
           [0.78, 0.1177]]
    singles = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (3, 0), (3, 1)]
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    sigma = 0.1196
    sigma_sq = sigma**2
    t = 0.0
    n_events = 500000
    # n_events = 100000
    histo_data = []

    for event in range(n_events):
        wall_times = [md.wall_time(pos[k][l], vel[k][l], sigma)
                      for k, l  in singles]
        pair_times = [md.pair_time(pos[k], vel[k], pos[l], vel[l], sigma, sigma_sq)
                      for k, l in pairs]
        next_event = min(wall_times + pair_times)
        t += next_event

        for k, l in singles: pos[k][l] += vel[k][l] * next_event

        if min(wall_times) < min(pair_times):
            collision_disk, direction = singles[wall_times.index(next_event)]
            vel[collision_disk][direction] *= -1.0
        else:
            a, b = pairs[pair_times.index(next_event)]
            del_x = [pos[b][0] - pos[a][0], pos[b][1] - pos[a][1]]
            abs_x = np.sqrt(del_x[0]**2 + del_x[1]**2)
            e_perp = [c / abs_x for c in del_x]
            del_v = [vel[b][0] - vel[a][0], vel[b][1] - vel[a][1]]
            scal = del_v[0] * e_perp[0] + del_v[1] * e_perp[1]
            for k in range(2):
                vel[a][k] += e_perp[k] * scal
                vel[b][k] -= e_perp[k] * scal
            # only recording pairwise disk collisions
            for k in pos: histo_data.append(k[0])
        print('event', event)
        print('time', t)
        # print('pos', pos)
        # print('vel', vel)

    pylab.hist(histo_data, bins=100, normed=True)
    pylab.xlabel('x')
    pylab.ylabel('frequency')
    pylab.title('x-coordinates for 5e5 events of event_disks_box\nwith 4 disks of radius 0.1196 (only non-wall collisions plotted)')
    pylab.grid()
    # pylab.savefig('event_disks_histo2.png')
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
    # markov_pi_test4()

    # ====================
    # Disk placement tests
    # ====================

    # disk_test1()
    # disk_test2()
    # disk_test3()
    disk_test4()

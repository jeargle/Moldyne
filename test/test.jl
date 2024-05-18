# John Eargle (mailto: jeargle at gmail.com)
# Moldyne test

using moldyne

using Plots
using Printf


function print_test_header(test_name)
    border = repeat("*", length(test_name) + 4)
    println(border)
    println("* ", test_name, " *")
    println(border)
end


function test_structure()
    print_test_header("Structure")

    pos1 = read_pdb_file("test2.pdb", 2)
    structure1 = Structure("structure1", 2, pos1)

    pos2 = read_pdb_file("test2.pdb", 3)
    structure2 = Structure("structure2", 3, pos2)

    println(structure1)
    println(structure2)
    println()
end


# Trial acceptance ratio
function test_markov_pi1()
    # n_runs = 100
    # n_trials = 100
    n_runs = 1000
    n_trials = 1000
    deltas = [(x+1)/10.0 for x in 1:50]
    acceptance_ratios = []

    for delta in deltas
        ars = 0.0
        for run in 1:n_runs
            n_hits, n_accepts = markov_pi(n_trials, delta)
            ars +=  n_accepts / float(n_trials)
        end
        push!(acceptance_ratios, ars / n_runs)
        @printf "%.1f, %f" delta (ars/n_runs)
    end

    p = plot(deltas,
             acceptance_ratios,
             title="1/2 Rule: Acceptance ratio as a function of Δ",
             xlabel="Δ",
             ylabel="acceptance ratio",
             legend=false,
             marker=:circle)
    savefig(p, "test_markov_pi1.svg")

end


# Accuracy
function test_markov_pi2()
    # n_runs = 100
    # n_trials = 100
    n_runs = 1000
    n_trials = 1000
    deltas = [(x+1)/10.0 for x in 1:50]
    sigmas = []

    for delta in deltas
        sigma = 0.0

        for run in 1:n_runs
            n_hits, n_accepts = markov_pi(n_trials, delta)
            pi_est = 4.0 * n_hits / float(n_trials)
            sigma += (pi_est - pi)^2
        end

        push!(sigmas, sqrt(sigma / n_runs))
        @printf "%.1f, %f" delta sqrt(sigma/n_runs)
    end

    p = plot(deltas,
             sigmas,
             title="Performance: Standard deviation σ as a function of Δ",
             xlabel="Δ",
             ylabel="σ",
             legend=false,
             marker=:circle)
    savefig(p, "test_markov_pi2.svg")
end


# Error
function test_markov_pi3()
    n_trials = 2^14
    delta = 0.1
    n_parties = 100
    inside_error_bar = 0

    for iteration in 1:n_parties
        mean, mean_square = markov_pi_all_data(n_trials, delta)
        naive_error = sqrt(mean_square  - mean^2) / sqrt(n_trials)
        error =  abs(mean - pi)

        if error < naive_error
            inside_error_bar += 1
        end

        println(mean, " ,", error, " ,", naive_error)
    end

    println(inside_error_bar / float(n_parties), " fraction: error bar including pi")
end


# Bunching
function test_markov_pi4()
    power = 14
    n_trials = 2^power
    delta = 0.1
    data = markov_pi_all_data2(n_trials, delta)
    errors  = []
    bunches = []

    for i in 1:power
        new_data = []
        mean = 0.0
        mean_sq = 0.0
        N = length(data)

        while data != []
            x = pop!(data)
            y = pop!(data)
            mean += x + y
            mean_sq += x^2 + y^2
            push!(new_data, (x + y) / 2.0)
        end

        push!(errors, sqrt(mean_sq/N - (mean/N)^2) / sqrt(N))
        push!(bunches, i)
        # data = new_data[:]
        data = new_data
    end

    p = plot(bunches,
             errors,
             title="Bunching: naive error vs iteration number",
             xlabel="iteration",
             ylabel="naive error",
             legend=false,
             marker=:circle)
    savefig(p, "test_markov_pi4.svg")
end


# Test direct placement of 4 disks in a box.
function test_disk1()
    N = 4
    sigma = 0.1196
    # n_runs = 100000
    n_runs = 1000000
    histo_data = []

    for run in 1:n_runs
        pos = direct_disks_box(N, sigma)
        for k in 1:N
            push!(histo_data, pos[k][1])
        end
    end

    p = histogram(histo_data,
                  bins=100,
                  title="x-coordinates for 1e6 runs of direct_disks_box\nwith 4 disks of radius 0.1196",
                  xlabel="x",
                  ylabel="frequency",
                  legend=false)
    savefig(p, "test_disk1.svg")
end


# Test markov placement of 1 disk in a box with 4 existing disks.
function test_disk2()
    L = [[0.25, 0.25], [0.75, 0.25], [0.25, 0.75], [0.75, 0.75]]
    sigma = 0.1196
    delta = 0.18   # 0.5 acceptance ratio

    n_steps = 2000000
    # n_steps = 20000
    accept_count = 0
    histo_data = []

    for steps in 1:n_steps
        acc, pos = markov_disks_box(L, sigma, delta)

        if acc
            accept_count += 1
            push!(histo_data, pos[1])
        end
    end

    @printf "  acceptance ratio: %f\n" accept_count/n_steps

    p = histogram(histo_data,
                  bins=100,
                  title="x-coordinates for 2e6 runs of markov_disks_box\nwith 4 disks of radius 0.1196 and Δ=0.18",
                  xlabel="x",
                  ylabel="frequency",
                  legend=false)
    savefig(p, "test_disk2.svg")
end


function test_disk3()
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
    singles = [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2), (4, 1), (4, 2)]
    pairs = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    sigma = 0.1196
    sigma_sq = sigma^2
    t = 0.0
    n_events = 200000
    histo_data = []

    for event in 1:n_events
        wall_times = [wall_time(pos[k][l], vel[k][l], sigma)
                      for (k, l) in singles]
        pair_times = [pair_time(pos[k], vel[k], pos[l], vel[l], sigma, sigma_sq)
                      for (k, l) in pairs]
        next_event = minimum([wall_times; pair_times])
        t += next_event
        for (k, l) in singles
            pos[k][l] += vel[k][l] * next_event
        end

        if minimum(wall_times) < minimum(pair_times)
            collision_disk, direction = singles[findfirst(==(next_event), wall_times)]
            vel[collision_disk][direction] *= -1.0
        else
            a, b = pairs[findfirst(==(next_event), pair_times)]
            del_x = [pos[b][1] - pos[a][1], pos[b][2] - pos[a][2]]
            abs_x = sqrt(del_x[1] ^ 2 + del_x[2] ^ 2)
            e_perp = [c / abs_x for c in del_x]
            del_v = [vel[b][1] - vel[a][1], vel[b][2] - vel[a][2]]
            scal = del_v[1] * e_perp[1] + del_v[2] * e_perp[2]
            for k in 1:2
                vel[a][k] += e_perp[k] * scal
                vel[b][k] -= e_perp[k] * scal
            end
        end

        println("event ", event)
        println("time ", t)
        # println("pos ", pos)
        # println("vel ", vel)
        for k in pos
            push!(histo_data, k[1])
        end
    end

    p = histogram(histo_data,
                  bins=100,
                  title="x-coordinates for 2e5 events of event_disks_box\nwith 4 disks of radius 0.1196",
                  xlabel="x",
                  ylabel="frequency",
                  legend=false)
    savefig(p, "test_disk3.svg")
end


# Molecular dynamics of four disks in a box.
function test_disk4()
    pos = [[0.25, 0.25],
           [0.75, 0.25],
           [0.25, 0.75],
           [0.75, 0.75]]
    vel = [[0.21, 0.12],
           [0.71, 0.18],
           [-0.23, -0.79],
           [0.78, 0.1177]]
    singles = [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2), (4, 1), (4, 2)]
    pairs = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    sigma = 0.1196
    sigma_sq = sigma^2
    t = 0.0
    n_events = 500000
    # n_events = 100000
    histo_data = []

    for event in 1:n_events
        wall_times = [wall_time(pos[k][l], vel[k][l], sigma)
                      for (k, l) in singles]
        pair_times = [pair_time(pos[k], vel[k], pos[l], vel[l], sigma, sigma_sq)
                      for (k, l) in pairs]
        next_event = minimum([wall_times; pair_times])
        t += next_event

        for (k, l) in singles
            pos[k][l] += vel[k][l] * next_event
        end

        if minimum(wall_times) < minimum(pair_times)
            collision_disk, direction = singles[findfirst(==(next_event), wall_times)]
            vel[collision_disk][direction] *= -1.0
        else
            a, b = pairs[findfirst(==(next_event), pair_times)]
            del_x = [pos[b][1] - pos[a][1], pos[b][2] - pos[a][2]]
            abs_x = sqrt(del_x[1]^2 + del_x[2]^2)
            e_perp = [c / abs_x for c in del_x]
            del_v = [vel[b][1] - vel[a][1], vel[b][2] - vel[a][2]]
            scal = del_v[1] * e_perp[1] + del_v[2] * e_perp[2]
            for k in 1:2
                vel[a][k] += e_perp[k] * scal
                vel[b][k] -= e_perp[k] * scal
            end
            # only recording pairwise disk collisions
            for k in pos
                push!(histo_data, k[1])
            end
        end
        println("event ", event)
        println("time ", t)
        # println("pos ", pos)
        # println("vel ", vel)
    end

    p = histogram(histo_data,
                  bins=100,
                  title="x-coordinates for 5e5 events of event_disks_box\nwith 4 disks of radius 0.1196 (only non-wall collisions plotted)",
                  xlabel="x",
                  ylabel="frequency",
                  legend=false)
    savefig(p, "test_disk4.svg")
end


# Molecular dynamics of four disks in a box.
function test_disk5()
    pos = [[0.25, 0.25],
           [0.75, 0.25],
           [0.25, 0.75],
           [0.75, 0.75]]
    # vel = [[0.21, 0.12],
    #        [0.71, 0.18],
    #        [-0.23, -0.79],
    #        [0.78, 0.1177]]
    vel = [[0.021, 0.012],
           [0.071, 0.018],
           [-0.023, -0.079],
           [0.078, 0.01177]]
    singles = [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2), (4, 1), (4, 2)]
    pairs = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    sigma = 0.1196
    sigma_sq = sigma^2
    t = 0.0
    n_events = 1000000
    # n_events = 100000
    histo_data = []

    for event in 1:n_events
        wall_times = [wall_time(pos[k][l], vel[k][l], sigma)
                      for (k, l) in singles]
        pair_times = [pair_time(pos[k], vel[k], pos[l], vel[l], sigma, sigma_sq)
                      for (k, l) in pairs]
        next_event = minimum([wall_times; pair_times])
        println("next_event ", next_event)

        t_previous = t
        for inter_times in round(Int, t + 1):round(Int, t + next_event + 1)
            del_t = inter_times - t_previous
            for (k, l) in singles
                pos[k][l] += vel[k][l] * del_t
            end
            t_previous = inter_times
            for k in 1:4
                push!(histo_data, pos[k][1])
            end
            println("time ", inter_times)
        end

        t += next_event
        for (k, l) in singles
            pos[k][l] += vel[k][l] * (t - t_previous)
        end

        if minimum(wall_times) < minimum(pair_times)
            collision_disk, direction = singles[findfirst(==(next_event), wall_times)]
            vel[collision_disk][direction] *= -1.0
        else
            a, b = pairs[findfirst(==(next_event), pair_times)]
            del_x = [pos[b][1] - pos[a][1], pos[b][2] - pos[a][2]]
            abs_x = sqrt(del_x[1]^2 + del_x[2]^2)
            e_perp = [c / abs_x for c in del_x]
            del_v = [vel[b][1] - vel[a][1], vel[b][2] - vel[a][2]]
            scal = del_v[1] * e_perp[1] + del_v[2] * e_perp[2]
            for k in 1:2
                vel[a][k] += e_perp[k] * scal
                vel[b][k] -= e_perp[k] * scal
            end
        end
        # println("event ", event)
        # println("time ", t)
        # println("pos ", pos)
        # println("vel ", vel)
    end

    # pylab.hist(histo_data, bins=100, density=True)
    # pylab.xlabel("x")
    # pylab.ylabel("frequency")
    # pylab.title("x-coordinates at regular timesteps for 1e6 events\nof event_disks_box with 4 disks of radius 0.1196")
    # pylab.grid()
    # # pylab.savefig("event_disks_histo3.png")
    # pylab.show()

    p = histogram(histo_data,
                  bins=100,
                  title="x-coordinates at regular timesteps for 1e6 events\nof event_disks_box with 4 disks of radius 0.1196",
                  xlabel="x",
                  ylabel="frequency",
                  legend=false)
    savefig(p, "test_disk5.svg")
end


# Markov simulation of four disks in a box.
function test_disk6()
    # detailed balance (global balance)
    # irreducible
    # aperiodic
    accepts = 0
    attempts = 0
    N = 4
    sigma = 0.1
    # n_runs = 1000000
    n_runs = 100000
    conf_a = [(0.25, 0.25),
              (0.25, 0.75),
              (0.75, 0.25),
              (0.75,0.75)]
    conf_b = [(0.20, 0.20),
              (0.20, 0.80),
              (0.75, 0.25),
              (0.75,0.75)]
    conf_c = [(0.30, 0.20),
              (0.30, 0.80),
              (0.70, 0.20),
              (0.70,0.70)]
    hits = [0, 0, 0]
    # Total = 0
    del_xy = 0.1
    configuration = [conf_a, conf_b, conf_c]

    for run in 1:n_runs
        if run % 1000 == 0
            println("run ", run)
        end
        at, ac, x_vec = direct_disks_box2(N, sigma)
        attempts += at
        accepts += ac
        for c in 1:3
            cond = true
            for b in configuration[c]
                cond_b = minimum(maximum( [abs(a[1] - b[1]), abs(a[2] - b[2])] ) for a in x_vec) < del_xy
                cond *= cond_b
            end
            if cond
                hits[c] += 1
            end
        end
    end

    for c in 1:3
        println(hits[c] / float(n_runs), " proportion of confs in eight-dimensional volume element.")
    end

    println("acceptance ratio: ", 1.0*accepts/attempts)
end


# Markov simulation of N disks in a box.
# First time, place N disks in a lattice.
# Afterwards, start from a restart file with coords for all disks.
function test_disk7()
    filename = "disk_conf.txt"
    # eta = 0.72
    # N = 256
    eta = 0.42
    N = 64

    k = int(sqrt(N) + 0.5)
    k_offset = 1.0/k
    sigma = sqrt(eta/pi)/k  # radius
    sigma_sq = sigma^2
    # delta = 0.5 * sigma
    delta = 0.1 * sigma
    n_steps = 1000
    accept = 0
    reject = 0
    L = []

    # Set locations
    if os.path.isfile(filename)
        # from input file
        open(filename, "r") do f
            for line in eachline(f)
                a, b = line.split()
                # push!(L, [float(a), float(b)])
                push!(L, [parse(Float64, a), parse(Float64, b)])

            end
        end
        println("starting from file ", filename)
    else
        # place on lattice
        for x in 1:k
            for y in 1:k
                push!(L, [k_offset/2.0 + k_offset*x, k_offset/2.0 + k_offset*y])
            end
        end
        println("starting from scratch")
    end

    uniform_dist = Uniform(-delta, delta)
    for step in 1:n_steps
        println("step ",step)
        a = random.choice(L)
        b = [(a[1] + rand(uniform_dist)) % 1.0,
             (a[2] + rand(uniform_dist)) % 1.0]
        min_dist = minimum(disk_dist(b, c) for c in L if c != a)
        println(" ", min_dist)
        if !(min_dist < 2.0 * sigma)
            a[:] = b
            accept += 1
            println("  accept")
            # println(L)
        else
            reject += 1
            # println("  reject")
        end
    end

    println("Acceptance ratio: ", float(accept)/n_steps)

    f = open(filename, "w")
    for a in L
        f.write(str(a[1]) + " " + str(a[2]) + "\n")
    end
    f.close()

    println("sigma: ", sigma)

    # md.show_conf(L, sigma, "test graph", "four_disks_b2.png")
    # md.show_conf(L, sigma, "test graph")
    # md.show_conf(L, sigma, "N=%d, $\eta$=%.2f" % (N, eta))
    show_conf(L, sigma, @sprintf "N=%d, η=%.2f" N eta, "test_disk7.png")
end


# Draw a disk.
function test_plot1()
    L = [[0.9, 0.9]]
    sigma = 0.4
    # md.show_conf(L, sigma, "test graph", "one_disk.png")
    show_conf(L, sigma, "test graph", "test_plot1.png")
end


# Read and print 2D coordinates from an existing file or write 3
# random 2D coordinates to a new file.
function test_file_io1()
    filename = "disk_configuration.txt"
    L = []

    if os.path.isfile(filename)
        f = open(filename, "r")

        for line in f
            a, b = line.split()
            push!(L, [float(a), float(b)])
        end

        f.close()
        println("starting from file ", filename)
    else
        uniform_dist = Uniform(0.0, 1.0)
        for k in 1:3
            push!(L, [rand(uniform_dist), rand(uniform_dist)])
        end

        println("starting from scratch")
    end

    L[1][1] = 3.3
    f = open(filename, "w")

    for a in L
        f.write(str(a[1]) + " " + str(a[2]) + "\n")
    end

    f.close()
end


function test_phi1()
    a = phi6(0.0)
    b = phi6(pi/6.0)
    c = phi6(pi/12.0)

    println(a)
    println(b)
    println(c)
    println(a+b)
end


function test_volume1()
    for dimension in 1:20
        println(dimension, ": ", sphere_volume(dimension))
    end
end


function test_volume2()
    dimensions = 1:200
    volumes = []

    for dimension in dimensions
        vol = sphere_volume(dimension)
        push!(volumes, vol)
        println(dimension, " ", vol)
    end

    p = plot(dimensions,
             volumes,
             yaxis=:log,
             title="Hypersphere volume (unit radius) in d dimensions",
             xlabel="d (dimensions)",
             ylabel="volume",
             legend=false)
    savefig(p, "test_volume2.svg")
end


function test_volume3()
    n_trials = 1000000

    println("---------------------------------------------------------")
    println("n_trials=", n_trials, " used for all")
    println("d | estimation of volume(d) | volume(d) (exact) | n_hits")
    println("---------------------------------------------------------")

    for d in 1:12
        n_hits = direct_pi(n_trials, d)
        est_vol = 2.0^d * n_hits / float(n_trials)
        exact_vol = sphere_volume(d)
        @printf "%2d  %1.4f                    %1.4f              %7d\n" d est_vol exact_vol n_hits
    end
end


# Markov chain sampling of sphere volume for different
# dimensionalities.
function test_volume4()
    n_trials = 10000

    # print 2.0**d * n_hits / float(n_trials)
    println("Average distance from origin to sample")

    for dim in 1:10
        # print sum(r_sqs)
        # print n_trials
        r_sqs = sample_sphere(n_trials, dim)
        println("  ", dim, " ", sum(r_sqs)/n_trials)
    end

    # pylab.plot([i[0] for i in points], [i[1] for i in points], ".")
    # pylab.axis([-1.5, 1.5, -1.5, 1.5])
    # pylab.xlabel("x")
    # pylab.ylabel("y")
    # pylab.title("Markov chain")
    # # # pylab.savefig("a2.png")
    # pylab.show()
end


function main()
    # ====================
    # Structure
    # ====================

    # test_structure()

    # ====================
    # Markov Pi
    # ====================

    # test_markov_pi1()
    # test_markov_pi2()
    # test_markov_pi3()
    # test_markov_pi4()

    # ====================
    # Disk Placement
    # ====================

    # test_disk1()
    # test_disk2()
    # test_disk3()
    # test_disk4()
    # test_disk5()
    # test_disk6()
    # test_disk7()

    # ====================
    # Plotting
    # ====================

    # test_plot1()

    # ====================
    # File IO
    # ====================

    # test_file_io1()

    # ====================
    # Phi
    # ====================

    # test_phi1()

    # ====================
    # Volume
    # ====================

    # test_volume1()
    # test_volume2()
    # test_volume3()
    test_volume4()
    # test_volume5()
    # test_volume6()
    # test_volume7()
    # test_volume8()

end

main()

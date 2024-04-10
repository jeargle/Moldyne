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
                      for (k, l)  in singles]
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
        println("event", event)
        println("time", t)
        # println("pos", pos)
        # println("vel", vel)
    end

    p = histogram(histo_data,
                  bins=100,
                  title="x-coordinates for 5e5 events of event_disks_box\nwith 4 disks of radius 0.1196 (only non-wall collisions plotted)",
                  xlabel="x",
                  ylabel="frequency",
                  legend=false)
    savefig(p, "test_disk4.svg")
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
    test_disk4()
    # test_disk5()
    # test_disk6()
    # test_disk7()

end

main()

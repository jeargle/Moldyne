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


function main()
    test_structure()
    test_markov_pi1()
    test_markov_pi2()
end

main()

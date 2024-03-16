# John Eargle (mailto: jeargle at gmail.com)
# moldyne

module moldyne

using Distributions
using Plots
using Printf
using Random

export Structure, Trajectory, set_positions, read_pdb_file, read_xyz_file
export markov_pi, markov_pi_all_data, markov_pi_all_data2

# Molecular structure including atomic positions and pairwise bonds.
struct Structure
    name::AbstractString
    dimension::Int64
    positions::Array{Array{Float64, 1}, 1}

    Structure(name::AbstractString, dimension::Int64, positions) = new(name, dimension, positions)
end


#
struct Trajectory
    structure::Structure
    trajectory_filename::AbstractString
    dimension::Int64

    function Trajectory(structure::Structure, trajectory_filename::AbstractString, dimension::Int64)
        new(structure, trajectory_filename, dimension)
    end

end


#
struct MdSystem
    structure::Structure
    temperature::Float64
    dimension::Int64
    positions::Array{Array{Float64, 1}, 1}
    new_positions::Array{Array{Float64, 1}, 1}
    velocities::Array{Array{Float64, 1}, 1}
    forces::Array{Array{Float64, 1}, 1}
    energy::Float64
    kinetic_energy::Float64
    timestep::Float64

    function MdSystem(structure::Structure, temperature, dimension)
        positions = Array{Float64, 1}[]
        new_positions = Array{Float64, 1}[]
        velocities = Array{Float64, 1}[]
        forces = Array{Float64, 1}[]

        new(structure, temperature, dimension, positions, new_positions, velocities, forces, 0.0, 0.0, 0.01)
    end

end


function set_positions(structure::Structure)
    return 0
end


# Read trajectory file and extract atomic positions for each frame.
function read_pdb_file(filename, dimension)
    positions = Array{Float64, 1}[]
    open(filename, "r") do f
        for line in eachline(f)
            println(line)
            println(line[31:38])
            println(line[39:46])
            c1 = parse(Float64, line[31:38])
            c2 = parse(Float64, line[39:46])

            if dimension == 2
                push!(positions, [c1, c2])
            elseif dimension == 3
                c3 = parse(Float64, line[47:54])
                push!(positions, [c1, c2, c3])
            end
        end
    end

    return positions
end


# Read trajectory file and extract atomic positions for each frame.
function read_xyz_file(trajectory::Trajectory)
    trajectory_file = open(trajectory.trajectory_filename, "r")
    num_atoms = Int64(readline(trajectory_file))
    if num_atoms != length(trajectory.structure.potitions)
        println("Error: atom count mismatch")
    end

    for line in eachline(trajectory_file)
        print(line)
    end
    close(trajectory_file)
end


# Reject the move if it falls outside the unit square.
# Count a hit if it falls within the unit circle.
function markov_pi(N, delta)
    x, y = 1.0, 1.0
    n_hits = 0
    n_accepts = 0
    uniform_dist = Uniform(-delta, delta)

    for i in 1:N
        del_x = rand(uniform_dist)
        del_y = rand(uniform_dist)

        if abs(x + del_x) < 1.0 && abs(y + del_y) < 1.0
            x += del_x
            y += del_y
            n_accepts += 1
        end

        if x^2 + y^2 < 1.0
            n_hits += 1
        end
    end

    return n_hits, n_accepts
end

function markov_pi_all_data(N, delta)
    x, y = 1.0, 1.0
    data_sum = 0.0
    data_sum_sq = 0.0
    uniform_dist = Uniform(-delta, delta)

    for i in 1:N
        del_x = rand(uniform_dist)
        del_y = rand(uniform_dist)

        if abs(x + del_x) < 1.0 && abs(y + del_y) < 1.0
            x = x + del_x
            y = y + del_y
        end

        if x^2 + y^2 < 1.0
            data_sum += 4.0
            data_sum_sq += 16.0  # 4.0^2
        end
    end

    return data_sum / float(N), data_sum_sq / float(N)
end

function markov_pi_all_data2(N, delta)
    x, y = 1.0, 1.0
    data = []
    uniform_dist = Uniform(-delta, delta)

    for i in 1:N
        del_x = rand(uniform_dist)
        del_y = rand(uniform_dist)

        if abs(x + del_x) < 1.0 && abs(y + del_y) < 1.0
            x = x + del_x
            y = y + del_y
        end

        if x^2 + y^2 < 1.0
            push!(data, 4.0)
        else
            push!(data, 0.0)
        end
    end

    return data
end

end

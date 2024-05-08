# John Eargle (mailto: jeargle at gmail.com)
# moldyne

module moldyne

using Distributions
using Plots
using Printf
using Random
using SpecialFunctions

export Structure, Trajectory, set_positions, read_pdb_file, read_xyz_file
export markov_pi, markov_pi_all_data, markov_pi_all_data2
export direct_disks_box, direct_disks_box2, markov_disks_box
export wall_time, pair_time, disk_dist, phi6
export direct_pi, sphere_volume
export show_conf

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

# Place N non-overlapping disks in a box.
# N: number of disks to place
# sigma: disk radius
function direct_disks_box(N, sigma)
    overlap = true
    uniform_dist = Uniform(sigma, 1.0 - sigma)
    L = []

    while overlap
        L = [(rand(uniform_dist), rand(uniform_dist))]
        for k in 2:N
            a = (rand(uniform_dist), rand(uniform_dist))
            min_dist_sq = minimum( [((a[1] - b[1])^2 + (a[2] - b[2])^2) for b in L] )

            if min_dist_sq < 4.0 * sigma^2
                overlap = true
                break
            else
                overlap = false
                push!(L, a)
            end
        end
    end

    return L
end

function direct_disks_box2(N, sigma)
    condition = false
    att = 0
    acc = 0
    uniform_dist = Uniform(sigma, 1.0 - sigma)
    L = []

    while condition == false
        # L = [(random.uniform(sigma, 1.0 - sigma), random.uniform(sigma, 1.0 - sigma))]
        L = [(rand(uniform_dist), rand(uniform_dist))]

        for k in 1:N
            att += 1
            # a = (random.uniform(sigma, 1.0 - sigma), random.uniform(sigma, 1.0 - sigma))
            a = (rand(uniform_dist), rand(uniform_dist))
            min_dist = minimum( sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2) for b in L )
            if min_dist < 2.0 * sigma
                condition = false
                break
            else
                push!(L, a)
                acc += 1
                condition = true
            end
        end
    end

    return att, acc, L
end

# Place one disk in a box with existing disks.  The placement
# is a jump from one of the existing disks.
# L: array of existing disk locations
# sigma: disk radius
# delta: max jump distance
function markov_disks_box(L, sigma, delta)
    sigma_sq = sigma^2
    accepted = false
    uniform_dist = Uniform(-delta, delta)
    a = rand(L)
    b = [a[1] + rand(uniform_dist), a[2] + rand(uniform_dist)]
    min_dist = minimum( (b[1] - c[1])^2 + (b[2] - c[2])^2 for c in L if c != a )
    box_cond = min(b[1], b[2]) < sigma || max(b[1], b[2]) > (1.0 - sigma)

    if !(box_cond || min_dist < 4.0 * sigma_sq)
        accepted = true
    end

    return accepted, b
end

# Determine amount of time before disk motion in 1D hits a wall.
function wall_time(pos_a, vel_a, sigma)
    if vel_a > 0.0
        del_t = (1.0 - sigma - pos_a) / vel_a
    elseif vel_a < 0.0
        del_t = (pos_a - sigma) / abs(vel_a)
    else
        del_t = Inf
    end

    return del_t
end

# Determine amount of time disk motions in 2D cause a pair of disks
# to hit each other.
function pair_time(pos_a, vel_a, pos_b, vel_b, sigma, sigma_sq)
    del_x = [pos_b[1] - pos_a[1], pos_b[2] - pos_a[2]]
    del_x_sq = del_x[1]^2 + del_x[2]^2
    del_v = [vel_b[1] - vel_a[1], vel_b[2] - vel_a[2]]
    del_v_sq = del_v[1]^2 + del_v[2]^2
    scal = del_v[1] * del_x[1] + del_v[2] * del_x[2]
    upsilon = scal^2 - del_v_sq * (del_x_sq - 4.0 * sigma_sq)

    if upsilon > 0.0 && scal < 0.0
        del_t = -(scal + sqrt(upsilon)) / del_v_sq
    else
        del_t = Inf
    end

    return del_t
end


# Distance between two disks in a 1x1 square with periodic boundary
# conditions.
function disk_dist(x, y)
    d_x = abs(x[0] - y[0]) % 1.0
    d_x = min(d_x, 1.0 - d_x)
    d_y = abs(x[1] - y[1]) % 1.0
    d_y = min(d_y, 1.0 - d_y)

    return  sqrt(d_x^2 + d_y^2)
end


function phi6(phi)
    return sum([exp(6.0im * ( phi + (x * pi/3.0))) for x in 1:6]) / 6.0
end


# ====================
# Sphere Volume
# ====================

function direct_pi(N, dimension)
    """
    Direct randomized calculation of PI.
    """
    n_hits = 0
    uniform_dist = Uniform(-1.0, 1.0)

    for i in 1:N
        x_len = 0.0
        for j in 1:dimension
            x = rand(uniform_dist)
            x_len += x^2
            if x_len > 1.0
                break
            end
        end

        if x_len <= 1.0
            n_hits += 1
        end
    end

    return n_hits
end


# Volume of N-dimensional sphere.
function sphere_volume(dimension)
    return pi^(dimension / 2.0) / gamma(dimension / 2.0 + 1.0)
end


# ====================
# Plotting
# ====================

function show_conf(L, sigma, title, fname=None)
    # pylab.axes()
    # for [x, y] in L
    #     for ix in -1:1
    #         for iy in -1:1
    #             cir = pylab.Circle((x + ix, y + iy), radius = sigma,  fc = 'r')
    #             pylab.gca().add_patch(cir)
    #         end
    #     end
    # end

    # pylab.axis('scaled')
    # pylab.title(title)
    # pylab.axis([0.0, 1.0, 0.0, 1.0])

    # if fname is not None
    #     pylab.savefig(fname)
    # end

    # pylab.show()
end

end

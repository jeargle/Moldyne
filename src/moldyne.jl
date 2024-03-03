# John Eargle (mailto: jeargle at gmail.com)
# moldyne

module moldyne

using Plots
using Printf
using Random

export Structure, Trajectory, set_positions, read_pdb_file, read_xyz_file

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
    trajectoryFilename::AbstractString
    dimension::Int64

    function Trajectory(structure::Structure, trajectoryFilename::AbstractString, dimension::Int64)
        new(structure, trajectoryFilename, dimension)
    end

end


#
struct MdSystem
    structure::Structure
    temperature::Float64
    dimension::Int64
    positions::Array{Array{Float64, 1}, 1}
    newPositions::Array{Array{Float64, 1}, 1}
    velocities::Array{Array{Float64, 1}, 1}
    forces::Array{Array{Float64, 1}, 1}
    energy::Float64
    kineticEnergy::Float64
    timestep::Float64

    function MdSystem(structure::Structure, temperature, dimension)
        positions = Array{Float64, 1}[]
        newPositions = Array{Float64, 1}[]
        velocities = Array{Float64, 1}[]
        forces = Array{Float64, 1}[]

        new(structure, temperature, dimension, positions, newPositions, velocities, forces, 0.0, 0.0, 0.01)
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
    trajectoryFile = open(trajectory.trajectoryFilename, "r")
    numAtoms = Int64(readline(trajectoryFile))
    if numAtoms != length(trajectory.structure.potitions)
        println("Error: atom count mismatch")
    end

    for line in eachline(trajectoryFile)
        print(line)
    end
    close(trajectoryFile)
end

end

# John Eargle (mailto: jeargle at gmail.com)
# 2017-2018
# moldyne

module moldyne

export Structure, setPositions

# Molecular structure including atomic positions and pairwise bonds.
type Structure
    name::AbstractString
    dimension::Int64
    positions::Array{Array{Float64, 1}, 1}
    function Structure(name::AbstractString, dimension, filename)
        positions = Array{Float64, 1}[]
        open(filename, "r") do f
            for line in eachline(f)
                println(line)
                println(line[31:38])
                println(line[39:46])
                c1 = float(line[31:38])
                c2 = float(line[39:46])

                if dimension == 2
                    push!(positions, [c1, c2])
                elseif dimension == 3
                    c3 = float(line[47:54])
                    push!(positions, [c1, c2, c3])
                end
            end
        end
        new(name, dimension, positions)
    end
end


end

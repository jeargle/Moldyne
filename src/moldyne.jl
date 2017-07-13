# John Eargle (mailto: jeargle at gmail.com)
# 2017
# Moldyne

module Moldyne

# Molecular structure including atomic positions and pairwise bonds.
type Structure
    name::AbstractString
    dimension::Int64
    positions::Array{Float64, 1}
end

# Read PDB file and extract atomic positions.
# filename: input PDB file
function setPositions(structure, filename)
    numLines = 0

    open(filename, "r") do f
        for line in eachline(f)
            print(line)
            c1 = float(line[31:38])
            c2 = float(line[39:46])

            if structure.dimension == 2
                structure.positions[numLines] = [c1, c2]
            elseif structure.dimension == 3
                c3 = float(line[47:54])
                structure.positions[numLines] = [c1, c2, c3]
            end

            numLines += 1
        end
    end
end




class Structure

  attr_reader :positions

  # 
  def initialize(structureFilename, dimension)
    @structureFilename = structureFilename
    @dimension = dimension
    numPositions = getNumPositions
    @positions = Array.new(numPositions)
    setPositions
  end

  # @return Number of atoms in the simulation
  def numAtoms
    return @positions.length
  end

  # Get number of atomic positions defined in @structureFile
  # @return Number of lines in the file
  def getNumPositions
    numLines = 0

    File.readlines(@structureFilename).each do |line|
      numLines += 1
    end

    return numLines
  end

  # Read PDB file and extract atomic positions.
  def setPositions
    numLines = 0

    File.readlines(@structureFilename).each do |line|
      c1 = line[31,7].to_f
      c2 = line[39,7].to_f

      if @dimension == 2
        @positions[numLines] = Coord2d.new(c1, c2)
      elsif @dimension == 3
        c3 = line[47,7].to_f
        @positions[numLines] = Coord3d.new(c1, c2, c3)
      end

      numLines += 1
    end

  end

end

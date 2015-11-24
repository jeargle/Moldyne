# Author:: John Eargle (mailto: jeargle at gmail.com)
# 2015
# :title: Structure

require_relative 'Coord2d'
require_relative 'Coord3d'

# Molecular structure including atomic positions and pairwise bonds.
class Structure

  # 
  def initialize(structureFilename, dimension)
    @structureFilename = structureFilename
    @dimension = dimension
    numPositions = getNumPositions
    @positions = Array.new(numPositions)
    setPositions

  # @return Number of atoms in the simulation
  def numAtoms
    return @positions.length
  end

  # Get number of atomic positions defined in @structureFile
  # @return Number of lines in the file
  def getNumPositions
    structureFile = File.open(@structureFilename)
    numLines = 0

    structureFile.each_line do
      numLines += 1
    end

    structureFile.close
    return numLines
  end

  # Read PDB file and extract atomic positions.
  def setPositions
    structureFile = File.open(@structureFilename)
    numLines = 0

    structureFile.each_line do |line|
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

    structureFile.close
  end

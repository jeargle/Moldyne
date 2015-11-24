# Author:: John Eargle (mailto: jeargle at gmail.com)
# 2015
# :title: Trajectory

require_relative 'Coord2d'
require_relative 'Coord3d'
require_relative 'Structure'

#
class Trajectory

  # 
  def initialize(structure, trajectoryFilename, dimension)
    @structure = structure
    @trajectoryFilename = trajectoryFilename
    @dimension = dimension
  end

  # Read trajectory file and extract atomic positions for each frame.
  def readXyzFile
    numLines = 0
    
    trajectoryFile = File.open(@trajectoryFilename, "r")
    numAtoms = Integer(trajectoryFile.readline)
    if numAtoms != @structure.numAtoms then
      puts "Error: atom count mismatch"
    end
    
    trajectoryFile.each_line do |line|
      puts line
    end
    trajectoryFile.close
    
  end
      
end

# Author:: John Eargle (mailto: jeargle at gmail.com)
# 2015
# :title: Trajectory

require_relative 'Coord2d'
require_relative 'Coord3d'
require_relative 'Structure'

#
class Trajectory

  # 
  def initialize(structureFilename, trajectoryFilename, dimension)
    @structureFilename = structureFilename
    @trajectoryFilename = trajectoryFilename
    @dimension = dimension


    

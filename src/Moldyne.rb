#! /usr/bin/ruby

# Author:: John Eargle (mailto: jeargle at gmail.com)
# 2010-2015
# :title: Moldyne

require_relative 'MdSystem'

# 
class Moldyne

  def initialize()
    @structureFileName = ""
    @outFileName = ""
    @temperature = 1.0
    @dimension = 2
    @initialTimestep = 0
    @maxTimestep = 0
    parseArgs()
    @system = MDSystem.new(@structureFileName, @temperature, @dimension)
  end

  # Parse commandline arguments
  def parseArgs()
    if ARGV.length != 3
      printUsage
      exit
    end

    @maxTimestep = ARGV[0].to_i
    @structureFileName = ARGV[1]
    @outFileName = ARGV[2]
  end    

  # Print command-line usage information
  def printUsage
    print "usage -- Moldyne.rb <numTimesteps> <structureFile> <outFile>\n"
    print "  <numTimesteps>   - total number of timesteps to run\n"
    print "  <structureFile>  - the PDB file to read in\n"
    print "  <outFilePrefix>  - prefix for two output files (.out, .xyz)\n"
  end
      
  ########
  # MAIN #
  ########
  def main
    outFile = File.open(@outFileName + ".out","w")
    outFile.print "Setup\n"
    outFile.print "  pos: " + @system.getPositions()
    outFile.print "  vel: " + @system.getVelocities()
    
    xyzFile = File.open(@outFileName + ".xyz","w")
    xyzFile.print @system.numAtoms.to_s + "\n\n"
    xyzFile.print @system.getXyz()
    
    @initialTimestep.upto(@maxTimestep).each do |timestep|
      outFile.print "Timestep: #{timestep}\n"
      xyzFile.print "\n\n"
      print "Timestep: #{timestep}\n"
      @system.computeForces
      @system.integrate
      outFile.print @system.getSystemStats()
      outFile.print "  pos: " + @system.getPositions()
      outFile.print "  vel: " + @system.getVelocities()
      xyzFile.print @system.getXyz()

    end
    outFile.close
    xyzFile.close
  end
  
end



cmd = Moldyne.new()
cmd.main

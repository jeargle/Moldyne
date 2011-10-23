#! /usr/bin/ruby

# John Eargle
# October 2010, May 2011

require 'MdSystem'

# 
class CrimsonMd

  def initialize()
    @structureFileName = ""
    @outFileName = ""
    @temperature = 1.0
    @initialTimestep = 0
    @maxTimestep = 1000
    parseArgs()
    @system = MDSystem.new(@structureFileName,@temperature)
  end

  # Parse commandline arguments
  def parseArgs()
    if ARGV.length != 2
      printUsage
      exit
    end
    
    @structureFileName = ARGV[0]
    @outFileName = ARGV[1]
    
  end    

  # Print command-line usage information
  def printUsage
    
    print "usage -- main.rb <structureFile> <outFile>\n"
    print "  <structureFile>  - the PDB file to read in\n"
    print "  <outFile> - MD output\n"
    
  end
      
  ########
  # MAIN #
  ########
  def main
    #structureFile = "test1.pdb"
    
    outFile = File.open(@outFileName,"w")
    outFile.print "Setup\n"
    outFile.print "  pos: " + @system.getPositions()
    outFile.print "  vel: " + @system.getVelocities()
    
    @initialTimestep.upto(@maxTimestep).each do |timestep|
      outFile.print "Timestep: #{timestep}\n"
      print "Timestep: #{timestep}\n"
      @system.computeForces
      @system.integrate
      outFile.print @system.getSystemStats()
      outFile.print "  pos: " + @system.getPositions()
      outFile.print "  vel: " + @system.getVelocities()
    end
    outFile.close
  end
  
end



cmd = CrimsonMd.new()
cmd.main

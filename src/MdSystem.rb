# Author:: John Eargle (mailto: jeargle at gmail.com)
# October 2010, May 2011
# :title: MDSystem

require_relative 'Coord2d'


#
class MDSystem

  # 
  def initialize(structureFilename,temperature)
    @structureFilename = structureFilename
    numPositions = getNumPositions
    @positions = Array.new(numPositions)
    @newPositions = Array.new(numPositions)
    @velocities = Array.new(numPositions)
    @forces = Array.new(numPositions)
    @energy = 0.0
    @kineticEnergy = 0.0
    @timestep = 0.01
    @temperature = temperature
    @box = Coord2d.new(20.0,20.0)
    @cutoff = 14.0
    @cutoffEnergy = 4 * (1.0 / @cutoff**12 - 1.0 / @cutoff**6)
    setPositions
    setVelocities
  end

  # @return Number of atoms in the simulation
  def numAtoms
    return @positions.length
  end

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

  
  #
  def setPositions
    #0.upto(@newPositions.length-1).each do |i|
    #  @newPositions[i] = Coord2d.new(rand*@box.x,rand*@box.y)
    #end
    structureFile = File.open(@structureFilename)
    numLines = 0
    structureFile.each_line do |line|
      match = /(.+), (.+)/.match(line)
      @newPositions[numLines] = Coord2d.new(match[1].to_f,match[2].to_f)
      numLines += 1
    end
    structureFile.close
    
  end

  # 
  def setVelocities
    sumVelocities = Coord2d.new(0.0,0.0)
    sumVelocitiesSquared = 0
    0.upto(@velocities.length-1).each do |i|
      newVelocity = Coord2d.new(rand-0.5,rand-0.5)
      @velocities[i] = newVelocity
      sumVelocities = sumVelocities.plus(newVelocity)
      sumVelocitiesSquared += newVelocity.dot(newVelocity)
    end
    
    # Velocity center of mass
    sumVelocities = sumVelocities.times(1.0/@positions.length)
    # Kinetic energy
    sumVelocitiesSquared /= @positions.length
    scaleFactor = Math.sqrt((2*@temperature)/sumVelocitiesSquared)
    0.upto(@velocities.length-1).each do |i|
      @velocities[i] = @velocities[i].minus(sumVelocities).times(scaleFactor)
      @positions[i] = @newPositions[i].minus((@velocities[i]).times(@timestep))
    end
  end
  
  # 
  def computeForces
    #puts ">computeForces"
    @energy = 0.0
    0.upto(@forces.length-1) do |i|
      @forces[i] = Coord2d.new(0.0, 0.0)
    end
    # Loop over all pairs
    0.upto(@forces.length-2) do |i|
      (i+1).upto(@forces.length-1) do |j|	
	separation = @newPositions[i].minus(@newPositions[j])
	# Periodic boundary condition
	boxImage = separation.elementDivide(@box)
	boxImage = Coord2d.new(boxImage.x.round,boxImage.y.round)
	separation = separation.minus(@box.elementTimes(boxImage))
	sepLength = separation.length
	#puts "(#{i},#{j}) sepLength: #{sepLength}"
	sepSquared = separation.dot(separation)
	if sepLength < @cutoff then
	  sep2i = 1.0 / sepSquared
	  sep6i = sep2i ** 3.0
	  force = 48 * sep2i * sep6i * (sep6i-0.5)
	  #force = 48 * sep6i * (sep2i-0.5)
	  #puts "force: #{force}"
	  #@forces[i] = @forces[i].plus(separation.times(force))
	  #@forces[j] = @forces[j].minus(separation.times(force))
	  @forces[i] = @forces[i].plus(separation.times(force))
	  @forces[j] = @forces[j].minus(separation.times(force))
	  @energy += 4 * sep6i * (sep6i - 1.0) - @cutoffEnergy
	  #@energy += 4 * sep6i * (sep2i - 1.0) - @cutoffEnergy
	  #puts "energy: #{@energy}"
	end
      end
    end
    #puts "<computeForces"
  end


  # 
  def integrate
    sumVelocities = Coord2d.new(0.0,0.0)
    sumVelSquared = 0.0
    0.upto(@positions.length-1) do |i|
      newPos = @newPositions[i].times(2).minus(@positions[i]).plus(@forces[i].times(@timestep**2))
      @velocities[i] = newPos.minus(@positions[i]).times(1.0/(2*@timestep))
      newPos = newPos.elementMod(@box)
      sumVelocities = sumVelocities.plus(@velocities[i])
      sumVelSquared += @velocities[i].dot(@velocities[i])
      @positions[i] = @newPositions[i]
      @newPositions[i] = newPos
    end
    # K = 0.5 * m * v**2
    @kineticEnergy = 0.5 * sumVelSquared
    @temperature = @kineticEnergy / (2*@positions.length)
    
  end

  # 
  def getSystemStats
    stats = ""
    stats += "  PE: " + "%.3f" % @energy + ", "
    stats += "KE: " + "%.3f" % @kineticEnergy + ", "
    stats += "TOTAL: " + "%.3f" % (@energy + @kineticEnergy) + ", "
    particleEnergy = (@energy + @kineticEnergy) / @positions.length
    stats += "PARTICLE: " + "%.3f" % particleEnergy + ", "
    stats += "TEMP: " + "%.3f" % @temperature + "\n"
    return stats
  end

  #
  def getPositions
    posString = ""
    0.upto(@positions.length-2) do |i|
      posString += "#{i}: (" + "%.3f" % @newPositions[i].x + "," + "%.3f" % @newPositions[i].y + "), "
    end
    i = @positions.length-1
      posString += "#{i}: (" + "%.3f" % @newPositions[i].x + "," + "%.3f" % @newPositions[i].y + ")\n"
    return posString
  end

  # Return string for all atoms in xyz trajectory format.
  def getXyz
    xyzString = ""
    0.upto(@positions.length-1) do |i|
      xyzString += "#{i} " + "%.3f" % @newPositions[i].x + " " + "%.3f" % @newPositions[i].y + " 0.0\n"
    end
    return xyzString
  end

  #
  def getVelocities
    velString = ""
    0.upto(@velocities.length-2) do |i|
      velString += "#{i}: (" + "%.3f" % @velocities[i].x + "," + "%.3f" % @velocities[i].y + "), "
    end
    i = @velocities.length-1
    velString += "#{i}: (" + "%.3f" % @velocities[i].x + "," + "%.3f" % @velocities[i].y + ")\n"
    return velString
  end

end

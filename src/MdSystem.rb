# Author:: John Eargle (mailto: jeargle at gmail.com)
# 2010-2015
# :title: MDSystem

require_relative 'Coord2d'
require_relative 'Coord3d'
require_relative 'Structure'


#
class MDSystem

  # 
  def initialize(structure, temperature, dimension)
    # @structureFilename = structureFilename
    @structure = structure
    @temperature = temperature
    @dimension = dimension
    # @positions = @structure.positions
    @newPositions = @structure.positions
    numPositions = @newPositions.length
    @positions = Array.new(numPositions)
    # @newPositions = Array.new(numPositions)
    @velocities = Array.new(numPositions)
    @forces = Array.new(numPositions)
    @energy = 0.0
    @kineticEnergy = 0.0
    @timestep = 0.01

    if @dimension == 2
      @box = Coord2d.new(20.0, 20.0)
    elsif @dimension == 3
      @box = Coord3d.new(20.0, 20.0, 20.0)
    end

    @cutoff = 14.0
    @cutoffEnergy = 4 * (1.0 / @cutoff**12 - 1.0 / @cutoff**6)

    # setPositions
    setVelocities
  end

  # @return Number of atoms in the simulation
  def numAtoms
    return @positions.length
  end

  # Get number of atomic positions defined in @structureFile
  # @return Number of lines in the file
  # def getNumPositions
  #   numLines = 0

  #   File.readlines(@structureFilename).each do |line|
  #     numLines += 1
  #   end

  #   return numLines
  # end

  # Read PDB file and extract atomic positions.
  # def setPositions
  #   numLines = 0

  #   File.readlines(@structureFilename).each do |line|
  #     c1 = line[31,7].to_f
  #     c2 = line[39,7].to_f

  #     if @dimension == 2
  #       @newPositions[numLines] = Coord2d.new(c1, c2)
  #     elsif @dimension == 3
  #       c3 = line[47,7].to_f
  #       @newPositions[numLines] = Coord3d.new(c1, c2, c3)
  #     end

  #     numLines += 1
  #   end

  # end

  # Set initial velocities to random values.
  def setVelocities
    if @dimension == 2
      sumVelocities = Coord2d.new(0.0, 0.0)
    elsif @dimension == 3
      sumVelocities = Coord3d.new(0.0, 0.0, 0.0)
    end

    sumVelocitiesSquared = 0
    0.upto(@velocities.length-1).each do |i|
      if @dimension == 2
        newVelocity = Coord2d.new(rand-0.5, rand-0.5)
      elsif @dimension == 3
        newVelocity = Coord3d.new(rand-0.5, rand-0.5, rand-0.5)
      end
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
  
  # Compute all pairwise atomic forces for the current system state.
  def computeForces
    #puts ">computeForces"
    @energy = 0.0
    0.upto(@forces.length-1) do |i|
      if @dimension == 2
        @forces[i] = Coord2d.new(0.0, 0.0)
      elsif @dimension == 3
        @forces[i] = Coord3d.new(0.0, 0.0, 0.0)
      end
    end

    # Loop over all pairs
    0.upto(@forces.length-2) do |i|
      (i+1).upto(@forces.length-1) do |j|	
	separation = @newPositions[i].minus(@newPositions[j])
	# Periodic boundary condition
	boxImage = separation.elementDivide(@box)

        if @dimension == 2
	  boxImage = Coord2d.new(boxImage.x.round,
                                 boxImage.y.round)
        elsif @dimension == 3
	  boxImage = Coord3d.new(boxImage.x.round,
                                 boxImage.y.round,
                                 boxImage.z.round)
        end

	separation = separation.minus(@box.elementTimes(boxImage))
	sepLength = separation.length
	#puts "(#{i},#{j}) sepLength: #{sepLength}"

	if sepLength < @cutoff then
	  sepSquared = separation.dot(separation)
	  sep2i = 1.0 / sepSquared
	  sep6i = sep2i ** 3.0
	  force = 48 * sep2i * sep6i * (sep6i-0.5)
	  #force = 48 * sep6i * (sep2i-0.5)
	  #puts "force: #{force}"
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


  # Time integration of MD system.
  def integrate
    if @dimension == 2
      sumVelocities = Coord2d.new(0.0, 0.0)
    elsif @dimension == 3
      sumVelocities = Coord3d.new(0.0, 0.0, 0.0)
    end
    
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

  # Get statistics for the current state.
  def getSystemStats
    stats = "  "
    stats += "PE: " + "%.3f" % @energy + ", "
    stats += "KE: " + "%.3f" % @kineticEnergy + ", "
    stats += "TOTAL: " + "%.3f" % (@energy + @kineticEnergy) + ", "
    particleEnergy = (@energy + @kineticEnergy) / @positions.length
    stats += "PARTICLE: " + "%.3f" % particleEnergy + ", "
    stats += "TEMP: " + "%.3f" % @temperature + "\n"
    return stats
  end

  # Get positions for all atoms.
  def getPositions
    posString = ""

    if @dimension == 2
      0.upto(@positions.length-2) do |i|
        posString += "#{i}: (" + "%.3f" % @newPositions[i].x + "," +
                     "%.3f" % @newPositions[i].y + "), "
      end
      
      i = @positions.length-1
      posString += "#{i}: (" + "%.3f" % @newPositions[i].x + "," +
                   "%.3f" % @newPositions[i].y + ")\n"
    elsif @dimension == 3
      0.upto(@positions.length-2) do |i|
        posString += "#{i}: (" + "%.3f" % @newPositions[i].x + "," +
                     "%.3f" % @newPositions[i].y + "," +
                     "%.3f" % @newPositions[i].z + "), "        
      end
      
      i = @positions.length-1
      posString += "#{i}: (" + "%.3f" % @newPositions[i].x + "," +
                   "%.3f" % @newPositions[i].y + "," +
                   "%.3f" % @newPositions[i].z + ")\n"
    end
    
    return posString
  end

  # Return string for all atoms in xyz trajectory format.
  def getXyz
    xyzString = ""
    if @dimension == 2
      0.upto(@positions.length-1) do |i|
        xyzString += "#{i} " + "%.3f" % @newPositions[i].x + " " +
                     "%.3f" % @newPositions[i].y + " 0.0\n"
      end
    elsif @dimension == 3
      0.upto(@positions.length-1) do |i|
        xyzString += "#{i} " + "%.3f" % @newPositions[i].x + " " +
                     "%.3f" % @newPositions[i].y + " " +
                     "%.3f" % @newPositions[i].z + "\n"
      end
    end
    
    return xyzString
  end

  # Get velocities for all atoms.
  def getVelocities
    velString = ""
    if @dimension == 2
      0.upto(@velocities.length-2) do |i|
        velString += "#{i}: (" + "%.3f" % @velocities[i].x + "," +
                     "%.3f" % @velocities[i].y + "), "
      end
      
      i = @velocities.length-1
      velString += "#{i}: (" + "%.3f" % @velocities[i].x + "," +
                   "%.3f" % @velocities[i].y + ")\n"
    elsif @dimension == 3
      0.upto(@velocities.length-2) do |i|
        velString += "#{i}: (" + "%.3f" % @velocities[i].x + "," +
                     "%.3f" % @velocities[i].y + "," +
                     "%.3f" % @velocities[i].z + "), "
      end
      
      i = @velocities.length-1
      velString += "#{i}: (" + "%.3f" % @velocities[i].x + "," +
                   "%.3f" % @velocities[i].y + "," +
                   "%.3f" % @velocities[i].z + ")\n"
    end
    
    return velString
  end

end

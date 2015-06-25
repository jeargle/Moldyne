#! /usr/bin/ruby

# Author:: John Eargle (mailto: jeargle at gmail.com)
# 2011-2015
# :title: test


require_relative '../src/Coord2d'

# 
def computeForces(coord1, coord2, box)
  #puts ">computeForces"
  energy = 0.0
  forceV = Coord2d.new(0.0, 0.0)
  cutoff = 14.0
  
  # Loop over all pairs
  separation = coord1.minus(coord2)
  # Periodic boundary condition
  boxImage = separation.elementDivide(box)
  boxImage = Coord2d.new(boxImage.x.round,boxImage.y.round)
  separation = separation.minus(box.elementTimes(boxImage))
  sepLength = separation.length
  puts "sepLength: #{sepLength}"
  sepSquared = separation.dot(separation)
  if sepLength < cutoff then
    sep2i = 1.0 / sepSquared
    sep6i = sep2i ** 3.0
    # XXX - from book
    #force = 48 * sep2i * sep6i * (sep6i-0.5)
    #force = 48 * sep6i * (sep2i-0.5) * (1.0/sepLength)
    force = 48 * sep2i * sep6i * (sep6i-0.5)
    print "force: #{force}, "
    #@forces[i] = @forces[i].plus(separation.times(force))
    #@forces[j] = @forces[j].minus(separation.times(force))
    forceV = forceV.plus(separation.times(force))
    #@forces[j] = @forces[j].minus(separation.times(force))
    #@energy += 4 * sep6i * (sep6i - 1.0) - @cutoffEnergy
    energy = 4 * sep6i * (sep6i - 1.0)
    puts "energy: #{energy}"
  end
  #puts "<computeForces"
end

coord1 = Coord2d.new(0.0,0.0)
1.upto(40) do |i|
  coord2 = Coord2d.new(i*1.0/20,0.0)
  box = Coord2d.new(10.0,10.0)
  computeForces(coord1,coord2,box)
end


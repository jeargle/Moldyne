# Author:: John Eargle (mailto: jeargle at gmail.com)
# May 2011
# :title: Coord3d

# 
class Coord3d

  attr_reader :x, :y, :z

  def initialize(x,y,z)
    @x = x
    @y = y
  end

  # Scalar multiplication
  def times(scalar)
    return Coord2d.new(scalar*@x, scalar*@y)
  end

  # Element-wise multiplication
  def elementTimes(coord)
    return Coord2d.new(@x*coord.x, @y*coord.y)
  end

  # Element-wise division
  def elementDivide(coord)
    return Coord2d.new(@x/coord.x, @y/coord.y)
  end

  # Element-wise modulus
  def elementMod(coord)
    return Coord2d.new(@x%coord.x, @y%coord.y)
  end

  # Vector addition
  def plus(coord)
    return Coord2d.new(@x+coord.x, @y+coord.y)
  end

  # Vector subtraction
  def minus(coord)
    return Coord2d.new(@x-coord.x, @y-coord.y)
  end

  # Vector multiplication
  def dot(coord)
    return @x*coord.x + @y*coord.y
  end

  # Vector length
  def length()
    return Math.sqrt(@x*@x + @y*@y)
  end

end

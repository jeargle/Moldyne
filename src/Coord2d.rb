# DEPRECATED - classes not needed in julia

# Author:: John Eargle (mailto: jeargle at gmail.com)
# 2011-2015
# :title: Coord2d

# 2D coordinates
class Coord2d

  attr_reader :x, :y

  def initialize(x,y)
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

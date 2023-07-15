# DEPRECATED - classes not needed in julia

# Author:: John Eargle (mailto: jeargle at gmail.com)
# 2011-2015
# :title: Coord3d

# 3D coordinates
class Coord3d

  attr_reader :x, :y, :z

  def initialize(x,y,z)
    @x = x
    @y = y
    @z = z
  end

  # Scalar multiplication
  def times(scalar)
    return Coord3d.new(scalar*@x, scalar*@y, scalar*@z)
  end

  # Element-wise multiplication
  def elementTimes(coord)
    return Coord3d.new(@x*coord.x, @y*coord.y, @z*coord.z)
  end

  # Element-wise division
  def elementDivide(coord)
    return Coord3d.new(@x/coord.x, @y/coord.y, @z/coord.z)
  end

  # Element-wise modulus
  def elementMod(coord)
    return Coord3d.new(@x%coord.x, @y%coord.y, @z%coord.z)
  end

  # Vector addition
  def plus(coord)
    return Coord3d.new(@x+coord.x, @y+coord.y, @z+coord.z)
  end

  # Vector subtraction
  def minus(coord)
    return Coord3d.new(@x-coord.x, @y-coord.y, @z-coord.z)
  end

  # Vector multiplication
  def dot(coord)
    return @x*coord.x + @y*coord.y + @z*coord.z
  end

  # Vector length
  def length()
    return Math.sqrt(@x*@x + @y*@y + @z*@z)
  end

end

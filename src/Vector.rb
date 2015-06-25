# Author:: John Eargle (mailto: jeargle at gmail.com)
# May 2011
# :title: Vector


class Vector < Array

  def +(other)
    case other
    when Array
      raise "Incorrect Dimensions" unless self.size == other.size
      other = other.dup
      self.class.new(map{|i| i + other.shift})
    else
      super
    end
  end

  def -(other)
    case other
    when Array
      raise "Incorrect Dimensions" unless self.size == other.size
      other = other.dup
      self.class.new(map{|i| i - other.shift})
    else
      super
    end
  end

  def *(other)
    case other
    when Array
      raise "Incorrect Dimensions" unless self.size == other.size
      other = other.dup
      self.class.new(map{|i| i * other.shift})
    else
      super
    end
  end
  
  def /(other)
    case other
    when Array
      raise "Incorrect Dimensions" unless self.size == other.size
      other = other.dup
      self.class.new(map{|i| i / other.shift})
    else
      super
    end
  end

  def %(other)
    case other
    when Array
      raise "Incorrect Dimensions" unless self.size == other.size
      other = other.dup
      self.class.new(map{|i| i % other.shift})
    else
      super
    end
  end

  def dot(other)
    case other
    when Array
      raise "Incorrect Dimensions" unless self.size == other.size
      other = other.dup
      vec = self.class.new(map{|i| i * other.shift})
      sum = 0
      vec.each do |i|
	sum += i
      end
    else
      super
    end    
  end

end


class Array

  def to_vector
    Vector.new(self)
  end

end

# 
class Coord

  attr_reader :x, :y

  def initialize(x,y)
    @x = x
    @y = y
  end

  # Scalar multiplication
  def times(scalar)
    return Coord2d.new(scalar*@x, scalar*@y)
  end

  # Element-wise modulus
  def elementMod(coord)
    return Coord2d.new(@x%coord.x, @y%coord.y)
  end

end

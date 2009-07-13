class CoordinatesConverter
  attr_reader :x, :y, :zone, :northern_hemisphere, :lat, :lng

  MAJOR_AXIS = 6378137.0
  MINOR_AXIS = 6356752.3
  ECC = (MAJOR_AXIS * MAJOR_AXIS - MINOR_AXIS * MINOR_AXIS) / (MAJOR_AXIS * MAJOR_AXIS)
  ECC2 = ECC / (1.0 - ECC)
  K0 = 0.9996
  E4 = ECC * ECC
  E6 = ECC * E4

=begin 
  Y_ERROR = 0.0017442271873235882
	X_ERROR = 0.001317089164991181
=end
  # Initializes a new instance of CoordinatesConverter
  # == Parameters
  # utm
  # * :x => the utm_x coordinate, as a float. Default 0.
  # * :y => the utm_y coordinate, as a float. Default 0.
  # * :zone => the utm_zone coordinate, as an integer. Default 31 (Because I'm Catalan!).
  # * :northern_hemisphere => defines if the coordinate is at the northern hemisphere, as a boolean. Default true (Because I'm Catalan!).
  # latituge longitude
  # * :lat => the latitude coordinate, as a float. Default 0.
  # * :lng => the longitude coordinate, as a float. Default 0.
  def initialize(opts={})
    @x=opts[:x].to_f || 0
    @y=opts[:y].to_f || 0
    @zone=opts[:zone] || 31
    @northern_hemisphere=opts[:northern_hemisphere] || true
    @lat=(opts[:lat].to_f > -180 && opts[:lat].to_f < 180 ? opts[:lat].to_f : 0)
    @lng=(opts[:lng].to_f > -90 && opts[:lng].to_f < 90 ? opts[:lng].to_f : 0)
  end

  def x=(x)
    @x = x.to_f
  end

  def y=(y)
    @y = y.to_f 
  end

  def zone=(zone)
    @zone = zone.to_i
  end

  def northern_hemisphere=(hemisphere)
    @northern_hemisphere = (hemisphere == true)
  end

  def lat=(lat)
    @lat = lat.to_f if(lat.to_f > -180 && lat.to_f < 180)
  end

  def lng=(lng)
   @lng = lng.to_f if(lng.to_f > -90 && lng.to_f < 90)
  end


  def meridianDist(lat)
    c1 = MAJOR_AXIS * (1 - ECC / 4 - 3 * E4 / 64 - 5 * E6 / 256)
    c2 = -MAJOR_AXIS * (3 * ECC / 8 + 3 * E4 / 32 + 45 * E6 / 1024)
    c3 = MAJOR_AXIS * (15 * E4 / 256 + 45 * E6 / 1024)
    c4 = -MAJOR_AXIS * 35 * E6 / 3072

    (c1 * lat + c2 * Math.sin(lat * 2) + c3 * Math.sin(lat * 4) + c4 * Math.sin(lat * 6))
  end

  def latlng2utmxy
    centeralMeridian = Math::degree2radian -((30 - @zone) * 6 + 3)

    lat = Math::degree2radian(@lat + Y_ERROR)
    lon = Math::degree2radian(@lng + X_ERROR)

    latSin = Math.sin(lat)
    latCos = Math.cos(lat)
    latTan = latSin / latCos
    latTan2 = latTan * latTan
    latTan4 = latTan2 * latTan2

    n = MAJOR_AXIS / Math.sqrt(1 - ECC * (latSin*latSin))
    c = ECC2 * latCos*latCos
    a = latCos * (lon - centeralMeridian)
    m = meridianDist(lat)

    temp5 = 1.0 - latTan2 + c
    temp6 = 5.0 - 18.0 * latTan2 + latTan4 + 72.0 * c - 58.0 * ECC2
    temp11 = a**5#Math.pow(a, 5)

    @x = K0 * n * (a + (temp5 * a**3) / 6.0 + temp6 * temp11 / 120.0) + 500000

    temp7 = (5.0 - latTan2 + 9.0 * c + 4.0 * (c*c)) * a**4 / 24.0
    temp8 = 61.0 - 58.0 * latTan2 + latTan4 + 600.0 * c - 330.0 * ECC2
    temp9 = temp11 * a / 720.0

    @y = K0 * (m + n * latTan * ((a * a) / 2.0 + temp7 + temp8 * temp9))

    {:x => @x, :y => @y, :zone => @zone}
  end

  def utmxy2latlng
    centeralMeridian = Math::degree2radian -((30 - @zone) * 6 + 3)

    temp1 = Math.sqrt(1.0 - ECC)
    ecc1 = (1.0 - temp1) / (1.0 + temp1)
    ecc12 = ecc1 * ecc1
    ecc13 = ecc1 * ecc12
    ecc14 = ecc12 * ecc12

    @x -= 500000.0

    m = @y / K0
    um = m / (MAJOR_AXIS * (1.0 - (ECC / 4.0) - 3.0 * (E4 / 64.0) - 5.0 * (E6 / 256.0)))

    temp8 = (1.5 * ecc1) - (27.0 / 32.0) * ecc13
    temp9 = ((21.0 / 16.0) * ecc12) - ((55.0 / 32.0) * ecc14)

    latrad1 = um + temp8 * Math.sin(2 * um) + temp9 * Math.sin(4 * um) + (151.0 * ecc13 / 96.0) * Math.sin(6.0 * um)

    latsin1 = Math.sin(latrad1)
    latcos1 = Math.cos(latrad1)
    lattan1 = latsin1 / latcos1
    n1 = MAJOR_AXIS / Math.sqrt(1.0 - ECC * latsin1*latsin1)
    t2 = lattan1 * lattan1
    c1 = ECC2 * latcos1 * latcos1

    temp20 = (1.0 - ECC * latsin1 * latsin1)
    r1 = MAJOR_AXIS * (1.0 - ECC) / Math.sqrt(temp20 * temp20 * temp20)

    d1 = @x / (n1*K0)
    d2 = d1 * d1
    d3 = d1 * d2
    d4 = d2 * d2
    d5 = d1 * d4
    d6 = d3 * d3

    t12 = t2 * t2
    c12 = c1 * c1

    temp1 = n1 * lattan1 / r1
    temp2 = 5.0 + 3.0 * t2 + 10.0 * c1 - 4.0 * c12 - 9.0 * ECC2
    temp4 = 61.0 + 90.0 * t2 + 298.0 * c1 + 45.0 * t12 - 252.0 * ECC2 - 3.0 * c12
    temp5 = (1.0 + 2.0 * t2 + c1) * d3 / 6.0
    temp6 = 5.0 - 2.0 * c1 + 28.0 * t2 - 3.0 * c12 + 8.0 * ECC2 + 24.0 * t12

    @lat = ((latrad1 - temp1 * (d2 / 2.0 - temp2 * (d4 / 24.0) + temp4 * d6 / 720.0)) * 180 / Math::PI) - Y_ERROR
    @lng = ((centeralMeridian + (d1 - temp5 + temp6 * d5 / 120.0) / latcos1) * 180 / Math::PI) - X_ERROR

    @y += 500000.0

    {:lat => @lat, :lng => @lng, :zone => @zone}
  end
end
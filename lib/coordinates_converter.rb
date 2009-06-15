# CoordinatesConverter
require "math.rb"
class CoordinatesConverter

  attr_reader :x
  attr_reader :y
  attr_reader :zone
  attr_reader :northern_hemisphere
  attr_reader :lat
  attr_reader :lng
 
  PI = Math::PI
  UTM_SCALE_FACTOR = 0.9996
  # Ellipsoid model constants (actual values here are for WGS84) 
  SM_A = 6378137.0
  SM_B = 6356752.314
  SM_ECC_SQUARED = 6.69437999013e-03


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
    @x=opts[:x] || 0
    @y=opts[:y] || 0
    @zone=opts[:zone] || 31
    @northern_hemisphere=opts[:northern_hemisphere] || true
    @lat=opts[:lat] || 0
    @lng=opts[:lng] || 0
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

  def latlng2utmxy
    @zone = (((lng + 180) / 6) +1).to_f.floor
    
    xy=maplatlng2xy

    # Adjust easting and northing for UTM system. 
    @x = xy[0] * UTM_SCALE_FACTOR + 500000
    @y = xy[1] * UTM_SCALE_FACTOR
    @y += 10000000 if (@y < 0)

    {:x => @x, :y => @y, :zone => @zone}
  end

  def maplatlng2xy
    phi = Math::degree2radian(lat)
    lambda_ = Math::degree2radian(lng)
    lambda0 = utm_central_meridian

    xy=[]

    # Precalculate ep2 
    ep2 = ((SM_A ** 2) - (SM_B ** 2)) / (SM_B ** 2)

    # Precalculate nu2 
    nu2 = ep2 * (Math::cos(phi) ** 2)

    # Precalculate N 
    n = (SM_A ** 2) / (SM_B * Math::sqrt(1 + nu2))

    # Precalculate t 
    t = Math::tan(phi)
    t2 = t * t

    # Precalculate l 
    l = lambda_ - lambda0

    # Precalculate coefficients for l**n in the equations below so a normal human being can read the expressions for easting and northing -- l**1 and l**2 have coefficients of 1
    l3coef = 1 - t2 + nu2

    l4coef = 5 - t2 + 9 * nu2 + 4 * (nu2 * nu2)

    l5coef = 5 - 18 * t2 + (t2 * t2) + 14 * nu2 - 58 * t2 * nu2

    l6coef = 61 - 58 * t2 + (t2 * t2) + 270 * nu2 - 330 * t2 * nu2

    l7coef = 61 - 479 * t2 + 179 * (t2 * t2) - (t2 * t2 * t2)

    l8coef = 1385 - 3111 * t2 + 543 * (t2 * t2) - (t2 * t2 * t2)

    # Calculate easting (x) 
    xy << n * Math::cos(phi) * l + ((n / 6) * (Math::cos(phi) ** 3) * l3coef * (l ** 3)) + ((n / 120) * (Math::cos(phi) ** 5) * l5coef * (l ** 5)) + ((n / 5040) * (Math::cos(phi) ** 7) * l7coef * (l ** 7))
    # Calculate northing (y) 
    xy << arc_length_of_meridian + (t / 2 * n * (Math::cos(phi) ** 2) * (l ** 2)) + (t / 24 * n * (Math::cos(phi) ** 4) * l4coef * (l ** 4)) + (t / 720 * n * (Math::cos(phi) ** 6) * l6coef * (l ** 6)) + (t / 40320 * n * (Math::cos(phi) ** 8) * l8coef * (l ** 8))
    xy
  end

  def arc_length_of_meridian
    phi = Math::degree2radian(lat)
    # Precalculate n 
    n = (SM_A - SM_B) / (SM_A + SM_B)

    # Precalculate alpha 
    alpha = ((SM_A + SM_B) / 2) * (1 + ((n ** 2) / 4) + ((n ** 4) / 64))

    # Precalculate beta 
    beta = (-3 * n / 2) + (9 * (n ** 3) / 16) + (-3 * (n ** 5) / 32)

    # Precalculate gamma 
    gamma = (15 * (n ** 2) / 16) + (-15 * (n ** 4) / 32)

    # Precalculate delta 
    delta = (-35 * (n ** 3) / 48) + (105 * (n ** 5) / 256)

    # Precalculate epsilng 
    epsilng = (315 * (n ** 4) / 512)

    # Now calculate the sum of the series and return 
    (alpha * (phi + (beta * Math::sin(2 * phi)) + (gamma * Math::sin(4 * phi)) + (delta * Math::sin(6 * phi)) + (epsilng * Math::sin(8 * phi))))
  end

  def utmxy2latlng
    x = @x
    x -= 500000
    x /= UTM_SCALE_FACTOR
    
    y = @y
    y -= 10000000 unless @northern_hemisphere
    y /= UTM_SCALE_FACTOR
    mapxy2latlng(x,y)
    {:lat => @lat, :lng => @lng}
  end

  def utm_central_meridian
    Math::degree2radian(@zone*6 -183)
  end
 
  def mapxy2latlng(x,y)
    #philambda = []

    # Get the value of phif, the footpoint latitude. 
    phif = footpoint_latitude

    # Precalculate ep2 
    ep2 = ((SM_A ** 2) - (SM_B ** 2)) / (SM_B ** 2)

    # Precalculate cos (phif) 
    cf = Math::cos(phif)

    # Precalculate nuf2 
    nuf2 = ep2 * (cf ** 2)

    # Precalculate Nf and initialize Nfpow 
    nf = (SM_A ** 2) / (SM_B * Math::sqrt(1 + nuf2))
    nfpow = nf

    # Precalculate tf 
    tf = Math::tan(phif)
    tf2 = tf * tf
    tf4 = tf2 * tf2

    # Precalculate fractional coefficients for x**n in the equations below to simplify the expressions for latitude and lnggitude.
    x1frac = 1 / (nfpow * cf)

    nfpow *= nf # now equals nf**2) 
    x2frac = tf / (2 * nfpow)

    nfpow *= nf # now equals nf**3) 
    x3frac = 1 / (6 * nfpow * cf)

    nfpow *= nf # now equals nf**4) 
    x4frac = tf / (24 * nfpow)

    nfpow *= nf # now equals nf**5) 
    x5frac = 1 / (120 * nfpow * cf)

    nfpow *= nf # now equals nf**6) 
    x6frac = tf / (720 * nfpow)

    nfpow *= nf # now equals nf**7) 
    x7frac = 1 / (5040 * nfpow * cf)

    nfpow *= nf # now equals nf**8) 
    x8frac = tf / (40320 * nfpow)

    # Precalculate polynomial coefficients for x**n. -- x**1 does not have a polynomial coefficient.
    x2poly = -1 - nuf2

    x3poly = -1 - 2 * tf2 - nuf2

    x4poly = 5 + 3 * tf2 + 6 * nuf2 - 6 * tf2 * nuf2 - 3 * (nuf2 * nuf2) - 9 * tf2 * (nuf2 * nuf2)

    x5poly = 5 + 28 * tf2 + 24 * tf4 + 6 * nuf2 + 8 * tf2 * nuf2

    x6poly = -61 - 90 * tf2 - 45 * tf4 - 107 * nuf2 + 162 * tf2 * nuf2

    x7poly = -61 - 662 * tf2 - 1320 * tf4 - 720 * (tf4 * tf2)

    x8poly = 1385 + 3633 * tf2 + 4095 * tf4 + 1575 * (tf4 * tf2)

    # Calculate latitude 
    @lat = Math::radian2degree(phif + x2frac * x2poly * (x * x) + x4frac * x4poly * (x ** 4) + x6frac * x6poly * (x ** 6) + x8frac * x8poly * (x ** 8))

    # Calculate lnggitude 
    @lng = Math::radian2degree(utm_central_meridian + x1frac * x + x3frac * x3poly * (x ** 3) + x5frac * x5poly * (x ** 5) + x7frac * x7poly * (x ** 7))
  end 

  def footpoint_latitude
    # Precalculate n (Eq. 10.18) 
    n = (SM_A - SM_B) / (SM_A + SM_B)

    # Precalculate alpha_ (Eq. 10.22) 
    # (Same as alpha in Eq. 10.17) 
    alpha_ = ((SM_A + SM_B) / 2) * (1 + ((n ** 2) / 4) + ((n ** 4) / 64))

    # Precalculate y_ (Eq. 10.23) 
    y_ = y / alpha_

    # Precalculate beta_ (Eq. 10.22) 
    beta_ = (3 * n / 2) + (-27 * (n ** 3) / 32) + (269 * (n ** 5) / 512)

    # Precalculate gamma_ (Eq. 10.22) 
    gamma_ = (21 * (n ** 2) / 16) + (-55 * (n ** 4) / 32)

    # Precalculate delta_ (Eq. 10.22) 
    delta_ = (151 * (n ** 3) / 96) + (-417 * (n ** 5) / 128)

    # Precalculate epsilng_ (Eq. 10.22) 
    epsilng_ = (1097 * (n ** 4) / 512)

    # Now calculate the sum of the series (Eq. 10.21) 
    (y_ + (beta_ * Math::sin(2 * y_)) + (gamma_ * Math::sin(4 * y_)) + (delta_ * Math::sin(6 * y_)) + (epsilng_ * Math::sin(8 * y_)))
  end
end

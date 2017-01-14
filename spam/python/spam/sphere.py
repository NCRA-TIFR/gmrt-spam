###############################################################################

# import Python modules
from math import *

# import 3rd patry modules
from numpy import *

# import user modules
from acalc import *
try:
  import _sphere
except:
  __sphere = False
else:
  __sphere = True

###############################################################################
# be aware of the special cases -0 and .. 60 ..
# TODO: implement special cases in C-code
###############################################################################

def radec_to_name( radec ):
  rd = degdeg_to_hmsdms( radec )
  name = '%02d%02d' % ( rd[ 0 ], rd[ 1 ] )
  if ( radec[ 1 ] >= 0. ):
    name = name + '+%02d' % ( rd[ 3 ] )
  else:
    name = name + '%03d' % ( rd[ 3 ] )
  name = name + '%02d' % ( rd[ 4 ] )
  return name

###############################################################################

def convert_radec_from_j2000( radec, epoch, m = 3.075, n = 1.336 ):
  # TODO: make more accurate
  # TODO: special cases for dec close to poles
  [ ra, dec ] = radec
  depoch = epoch - 2000.
  mm = depoch * m * 15. / 3600.
  nn = depoch * n * 15. / 3600.
  dra = mm + nn * sin( radians( ra ) ) * tan( radians( dec ) )
  ddec = nn * cos( radians( ra ) )
  return [ amodulo( ra + dra, 360. ), dec + ddec ]

###############################################################################

def convert_b1950_to_j2000( radec ):
# based on Lieske (1976)
  [ ra1, dec1 ] = [ aradians( radec[ 0 ] ), aradians( radec[ 1 ] ) ]
  xyz1 = array( [ cos( ra1 ) * cos( dec1 ), sin( ra1 ) * cos( dec1 ), sin( dec1 ) ], dtype = float64 )
  rot = array( [ [   0.9999257079523629, - 0.0111789381377700, - 0.0048590038153592 ],
                 [   0.0111789381264276,   0.9999375133499888, - 0.0000271625947142 ],
                 [   0.0048590038414544, - 0.0000271579262585,   0.9999881946023742 ] ], dtype = float64 )
  xyz2 = dot( rot, xyz1 ).tolist()
  ra2 = amodulo( adegrees( atan2( xyz2[ 1 ], xyz2[ 0 ] ) ), 360. )
  dec2 = adegrees( max( - 1., min( 1., asin( xyz2[ 2 ] ) ) ) )
  return [ ra2, dec2 ]

###############################################################################

def convert_j2000_to_b1950( radec ):
# based on Lieske (1976)
  [ ra1, dec1 ] = [ aradians( radec[ 0 ] ), aradians( radec[ 1 ] ) ]
  xyz1 = array( [ cos( ra1 ) * cos( dec1 ), sin( ra1 ) * cos( dec1 ), sin( dec1 ) ], dtype = float64 )
  rot = array( [ [   0.9999257079523629,   0.0111789381264276,   0.0048590038414544 ],
                 [ - 0.0111789381377700,   0.9999375133499888, - 0.0000271579262585 ],
                 [ - 0.0048590038153592, - 0.0000271625947142,   0.9999881946023742 ] ], dtype = float64 )
  xyz2 = dot( rot, xyz1 ).tolist()
  ra2 = amodulo( adegrees( atan2( xyz2[ 1 ], xyz2[ 0 ] ) ), 360. )
  dec2 = adegrees( asin( max( - 1., min( 1., xyz2[ 2 ] ) ) ) )
  return [ ra2, dec2 ]

###############################################################################

def hmsdms_to_degdeg( hmsdms ):
#  if __sphere:
#    return _sphere.hmsdms_to_degdeg( [ float( x ) for x in hmsdms ] )
  ra_h = amodulo( hmsdms[ 0 ], 24. )
  ra = 15. * ( ra_h + ( hmsdms[ 1 ] / 60. ) + ( hmsdms[ 2 ] / 3600. ) )
  dec_d = asign( hmsdms[ 3 ] ) * degrees( asin( 
      max( - 1., min( 1., sin( radians( amodulo( fabs( hmsdms[ 3 ] ), 360. ) ) ) ) ) ) )
  dec = dec_d + asign( dec_d ) * ( ( hmsdms[ 4 ] / 60. ) + ( hmsdms[ 5 ] / 3600. ) )
  if ( dec > 90. ):
    dec = 90.
  elif ( dec < - 90. ):
    dec = - 90.
  return [ ra, dec ]

###############################################################################

def degdeg_to_hmsdms( degdeg, precision = None ):
#  if __sphere:
#    return _sphere.degdeg_to_hmsdms( [ float( x ) for x in degdeg ] )
  ra_deg = amodulo( degdeg[ 0 ], 360. )
  ra_h = floor( ra_deg / 15. )
  ra_m = floor( 60. * ( ( ra_deg / 15. ) - ra_h ) )
  ra_s = 3600. * ( ( ra_deg / 15. ) - ra_h - ( ra_m / 60. ) )
  dec_deg = asign( degdeg[ 1 ] ) * degrees( asin( 
      max( - 1., min( 1., sin( radians( amodulo( fabs( degdeg[ 1 ] ), 360. ) ) ) ) ) ) )
  dec_d = asign( dec_deg ) * floor( abs( dec_deg ) )
  dec_m = floor( 60. * abs( dec_deg - dec_d ) )
  dec_s = 3600. * ( abs( dec_deg - dec_d ) - ( dec_m / 60. ) )
  if ( not precision is None ):
    if ( len( shape( precision ) ) == 0 ):
      prec1 = int( precision )
      prec2 = int( precision )
    elif ( len( precision ) == 1 ):
      prec1 = int( precision[ 0 ] )
      prec2 = int( precision[ 0 ] )
    else:
      prec1 = int( precision[ 0 ] )
      prec2 = int( precision[ 1 ] )
    ra_s = around( ra_s, decimals = prec1 )
    dec_s = around( dec_s, decimals = prec2 )
  if ( ra_s >= 60. ):
    ra_s = ra_s - 60.
    ra_m = ra_m + 1.
  if ( ra_m >= 60. ):
    ra_m = ra_m - 60.
    ra_h = ra_h + 1.
  if ( ra_h >= 24. ):
    ra_h = ra_h - 24.
  if ( dec_s >= 60. ):
    dec_s = dec_s - 60.
    dec_m = dec_m + 1.
  if ( dec_m >= 60. ):
    dec_m = dec_m - 60.
    if ( asign( dec_deg ) > 0. ):
      dec_d = dec_d + 1.
      if ( dec_d == 90. ):
        dec_s = 0.
        dec_m = 0.
    else:
      dec_d = dec_d - 1.
      if ( dec_d == - 90. ):
        dec_s = 0.
        dec_m = 0.
  return [ ra_h, ra_m, ra_s, dec_d, dec_m, dec_s ]

###############################################################################

def dmsdms_to_degdeg( dmsdms ):
#  if __sphere:
#    return _sphere.dmsdms_to_degdeg( [ float( x ) for x in dmsdms ] )
  lon_d = amodulo( dmsdms[ 0 ], 360. )
  lon = lon_d + asign( dmsdms[ 0 ] ) * ( ( dmsdms[ 1 ] / 60. ) + ( dmsdms[ 2 ] / 3600. ) )
  lon = amodulo( lon + 180., 360. ) - 180.
  lat_d = asign( dmsdms[ 3 ] ) * degrees( asin( 
      max( - 1., min( 1., sin( radians( amodulo( fabs( dmsdms[ 3 ] ), 360. ) ) ) ) ) ) )
  lat = lat_d + asign( lat_d ) * ( ( dmsdms[ 4 ] / 60. ) + ( dmsdms[ 5 ] / 3600. ) )
  if ( lat > 90. ):
    lat = 90.
  elif ( lat < - 90. ):
    lat = - 90.
  return [ lon, lat ]

###############################################################################

def degdeg_to_dmsdms( degdeg, precision = None ):
#  if __sphere:
#    return _sphere.degdeg_to_dmsdms( [ float( x ) for x in degdeg ] )
  lon_deg = amodulo( degdeg[ 0 ] + 180., 360. ) - 180.
  lon_d = asign( lon_deg ) * floor( abs( lon_deg ) )
  lon_m = floor( 60. * abs( lon_deg - lon_d ) )
  lon_s = 3600. * ( abs( lon_deg - lon_d ) - ( lon_m / 60. ) )
  lat_deg = degrees( asin( 
      max( - 1., min( 1., sin( radians( amodulo( degdeg[ 1 ], 360. ) ) ) ) ) ) )
  lat_d = asign( lat_deg ) * floor( abs( lat_deg ) )
  lat_m = floor( 60. * abs( lat_deg - lat_d ) )
  lat_s = 3600. * ( abs( lat_deg - lat_d ) - ( lat_m / 60. ) )
  if ( not precision is None ):
    if ( len( shape( precision ) ) == 0 ):
      prec1 = int( precision )
      prec2 = int( precision )
    elif ( len( precision ) == 1 ):
      prec1 = int( precision[ 0 ] )
      prec2 = int( precision[ 0 ] )
    else:
      prec1 = int( precision[ 0 ] )
      prec2 = int( precision[ 1 ] )
    lon_s = around( lon_s, decimals = prec1 )
    lat_s = around( lat_s, decimals = prec2 )
  if ( lon_s >= 60. ):
    lon_s = lon_s - 60.
    lon_m = lon_m + 1.
  if ( lon_m >= 60. ):
    lon_m = lon_m - 60.
    lon_d = lon_d + 1.
  if ( lon_d >= 360. ):
    lon_d = lon_d - 360.
  if ( lat_s >= 60. ):
    lat_s = lat_s - 60.
    lat_m = lat_m + 1.
  if ( lat_m >= 60. ):
    lat_m = lat_m - 60.
    if ( asign( lat_deg ) > 0. ):
      lat_d = lat_d + 1.
      if ( lat_d == 90. ):
        lat_s = 0.
        lat_m = 0.
    else:
      lat_d = dec_d - 1.
      if ( lat_d == - 90. ):
        lat_s = 0.
        lat_m = 0.
  return( [ lon_d, lon_m, lon_s, lat_d, lat_m, lat_s ] )

###############################################################################

def calculate_angular_separation( degdeg1, degdeg2 ):
  if __sphere:
    return _sphere.calculate_angular_separation( [ float( x ) for x in degdeg1 ],
       [ float( x ) for x in degdeg2 ] )
  ra1 = amodulo( degdeg1[ 0 ], 360. )
  dec1 = degrees( asin( 
      max( - 1., min( 1., sin( radians( amodulo( degdeg1[ 1 ], 360. ) ) ) ) ) ) )
  ra2 = amodulo( degdeg2[ 0 ], 360. )
  dec2 = degrees( asin( 
      max( - 1., min( 1., sin( radians( amodulo( degdeg2[ 1 ], 360. ) ) ) ) ) ) )
  dra = amodulo( ( ra2 - ra1 ) + 180., 360. ) - 180.
  ddec = dec2 - dec1
  if ( ( dra == 0. ) or ( amodulo( dec1, 180. ) == 90. ) or ( amodulo( dec2, 180. ) == 90. ) ):
    radius = abs( ddec )
    if ( dec2 >= dec1 ):
      angle = 0.
    else:
      angle = 180.
#  elif ( ddec == 0. ):
#    radius = dra * cos( radians( dec1 ) )
#    if ( radius > 0. ):
#      angle = 90.
#    else:
#      angle = 270.
#    radius = abs( radius )
  else:
    # use Haversine formula (adapted from http://www.plutoproject.com/dist_pa2.cpp)
    hav_r = ( sin( radians( ddec / 2. ) )**2 + 
        cos( radians( dec1 ) ) * cos( radians( dec2 ) ) * sin( radians( dra / 2. ) )**2 )
    radius = degrees( 2. * asin( max( - 1., min( 1., sqrt( hav_r ) ) ) ) )
    t = cos( radians( dec1 ) ) * sin( radians( radius ) )
    hav_a = ( sin( radians( dec1 + radius ) ) - sin( radians( dec2 ) ) ) / ( 2. * t )
    angle = degrees( 2. * asin( min( 1., sqrt( hav_a ) ) ) )
    if ( sin( radians( dra ) ) < 0. ):
      angle = 360. - angle
  return [ radius, angle ]

###############################################################################

def calculate_offset_position( degdeg, radius, angle ):
# 0. <= radius <= 180.
  if __sphere:
    return _sphere.calculate_offset_position( [ float( x ) for x in degdeg ],
        float( radius ), float( angle ) )
  ra = degdeg[ 0 ]
  dec = degdeg[ 1 ]
  if ( radius <= 0. ):
    new_ra = ra
    new_dec = dec
  else:
    a = radians( radius )
    c = radians( 90. - dec )
    B = radians( - angle )
    b = acos( max( - 1., min( 1., sin( a ) * cos( B ) * sin( c ) + cos( a ) * cos( c ) ) ) )
    if ( b == 0. ):
      A = 0.
    else:
      A = asin( max( - 1., min( 1., sin( a ) * sin( B ) / sin( b ) ) ) )
      if ( ( ( cos( a ) * sin( c ) - sin( a ) * cos( B ) * cos( c ) ) / sin( b ) ) < 0. ):
        A = pi - A
    new_ra = amodulo( ra - degrees( A ), 360. )
    new_dec = 90. - degrees( b )
  return [ new_ra, new_dec ]

###############################################################################

def xyz_to_llr( xyz ):
  if __sphere:
    return _sphere.xyz_to_llr( [ float( x ) for x in xyz ] )
  x = xyz[ 0 ]
  y = xyz[ 1 ]
  z = xyz[ 2 ]
  lon = amodulo( degrees( atan2( y, x ) ) + 180., 360. ) - 180.
  lat = degrees( atan2( z, sqrt( x**2 + y**2 ) ) )
  rad = sqrt( x**2 + y**2 + z**2 )
  return [ lon, lat, rad ]

###############################################################################

def xyz_to_geo_llh( xyz, iterations = 4, a = 6378137., f = 1. / 298.257, e2 = 6.6943799013e-3 ):
# default Earth ellipticity definition (a,f) is WGS (1984)
# Note that longitude is defined as positive towards east, just like RA
  if __sphere:
    return _sphere.xyz_to_geo_llh( [ float( x ) for x in xyz ], int( iterations ),
        float( a ), float( f ), float( e2 ) )
  [ x, y, z ] = xyz
  glon = amodulo( degrees( atan2( y, x ) ) + 180., 360. ) - 180.
  glat = degrees( atan2( z, sqrt( x**2 + y**2 ) ) )
  gh = sqrt( x**2 + y**2 + z**2 ) - a * sqrt( 1. - f )
  if ( iterations > 0 ):
    phi = radians( glat )
    for i in range( iterations ):
      n = a / sqrt( 1. - e2 * ( sin( phi )**2 ) )
      gh = ( sqrt( x**2 + y**2 ) / cos( phi ) ) - n
      phi = atan( z / ( sqrt( x**2 + y**2 ) * ( 1. - e2 * ( n / ( n + gh ) ) ) ) )
    glat = degrees( phi )
  return [ glon, glat, gh ]

###############################################################################

def geo_llh_to_xyz( geo_llh, a = 6378137., f = 1. / 298.257, e2 = 6.6943799013e-3 ):
# default Earth ellipticity definition (a,f) is WGS (1984)
# Note that longitude is defined as positive towards east, just like RA
  if __sphere:
    return _sphere.geo_llh_to_xyz( [ float( x ) for x in geo_llh ], float( a ),
        float( f ), float( e2 ) )
  [ glon, glat, gh ] = geo_llh
  lamda = radians( glon )
  phi = radians( glat )
  n = a / sqrt( 1. - e2 * ( sin( phi )**2 ) )
  x = ( n + gh ) * cos( phi ) * cos( lamda )
  y = ( n + gh ) * cos( phi ) * sin( lamda )
  z = ( n * ( 1. - e2 ) + gh ) * sin( phi )
  return [ x, y, z ]

###############################################################################

def calculate_hour_angles_at_elevation_limit( lat, dec, elevation_limit = 0. ):
  if __sphere:
    return _sphere.calculate_hour_angles_at_elevation_limit( float( lat ), float( dec ),
        float( elevation_limit ) )
  if ( ( dec + lat >= 90. ) or ( dec + lat <= - 90. ) ): # check for circumpolar sources
    ha = 180.
  elif ( ( dec - lat >= 90. ) or ( dec - lat <= - 90. ) ): # check for non-visible sources
    ha = 0.
  else:
    a = radians( 90. - elevation_limit )
    b = radians( 90. - dec ) # 0 < b < 180
    c = radians( 90. - lat ) # 0 < c < 180
    A = acos( max( - 1., min( 1., ( cos( a ) - cos( b ) * cos( c ) ) / (  sin( b ) * sin( c ) ) ) ) )
    # 0 < A < 180 degrees
    ha = degrees( A )
  return [ - ha, ha ]

###############################################################################

def time_to_dhms( time ):
  if __sphere:
    return _sphere.time_to_dhms( float( time ) )
  res = abs( time )
  day = sign( time ) * floor( res )
  res = 24. * ( res - day )
  hour = floor( res )
  res = 60. * ( res - hour )
  mins = floor( res )
  sec = 60. * ( res - mins )
  return [ day, hour, mins, sec ]

###############################################################################

def round_dhms( dhms, decimals = 0 ):
  [ day, hour, mins, sec ] = dhms
  sec = round( sec, decimals )
  if ( sec == 60. ):
    mins = mins + 1.
    sec = 0.
  if ( mins == 60. ):
    hour = hour + 1.
    mins = 0.
  if ( hour == 24. ):
    day = sign( day ) * ( abs( day ) + 1. )
    hour = 0.
  return [ day, hour, mins, sec ]

###############################################################################

def dhms_to_time( dhms ):
#  if __sphere:
#    return _sphere.dhms_to_time( [ float( x ) for x in dhms ] )
  [ day, hour, mins, sec ] = dhms
  time = float( day ) + ( float( hour ) / 24. ) + ( float( mins ) / 1440. ) + ( float( sec ) / 86400. )
  return time

###############################################################################

def dhmsdhms_to_timetime( dhmsdhms ):
  timetime = [ dhms_to_time( dhmsdhms[ 0 : 4 ] ), dhms_to_time( dhmsdhms[ 4 : 8 ] ) ]
  return timetime

###############################################################################

def timetime_to_dhmsdhms( timetime ):
  dhmsdhms = time_to_dhms( timetime[ 0 ] ) + time_to_dhms( timetime[ 1 ] )
  return dhmsdhms

###############################################################################

def calculate_enu( ref_xyz, xyz ):
  rot_xyz = array( xyz, dtype = float64 )
  ref_geo_llh = xyz_to_geo_llh( ref_xyz )
  ref_lon = radians( ref_geo_llh[ 0 ] )
  ref_lat = radians( ref_geo_llh[ 1 ] )
  rot = array( [ [ - sin( ref_lon )                 ,   cos( ref_lon )                 ,             0. ], 
                 [ - cos( ref_lon ) * sin( ref_lat ), - sin( ref_lon ) * sin( ref_lat ), cos( ref_lat ) ],
                 [   cos( ref_lon ) * cos( ref_lat ),   sin( ref_lon ) * cos( ref_lat ), sin( ref_lat ) ] ],
                 dtype = float64 )
  rot_xyz = dot( rot, rot_xyz )
  return rot_xyz.tolist()

###############################################################################

def calculate_local_sky_position( geo_llh, radec, gst ):
# gst in degrees
  if __sphere:
    return _sphere.calculate_local_sky_position( [ float( x ) for x in geo_llh ],
        [ float( x ) for x in radec ], float( gst ) )
  [ ra, dec ] = radec
  [ lon, lat ] = geo_llh[ 0 : 2 ]
  if ( ( amodulo( lat, 180. ) == 90. ) or ( amodulo( dec, 180. ) == 90. ) ):
    za = amodulo( abs( dec - lat ), 180. )
    if ( ( dec == 90. ) or ( lat == - 90 ) ):
      az = 0.
    else:
      az = 180.
  else:
    ha = amodulo( gst - ra + lon, 360. )
    A = radians( -ha )
    b = radians( 90. - dec )
    c = radians( 90. - lat )
    a = acos( max( - 1., min( 1., sin( b ) * cos( A ) * sin( c ) + cos( b ) * cos( c ) ) ) )
    B = acos( max( - 1., min( 1., ( cos( b ) - cos( a ) * cos( c ) ) / ( sin( a ) * sin( c ) ) ) ) )
    za = degrees( a )
    az = degrees( B )
    if ( sin( A ) < 0. ):
      az = 360. - az
  return [ za, az ]

###############################################################################

def calculate_pierce_point( xyz, radec, gst, height = 400.e3, iterations = 4 ):
# height in meters
# radec at observing epoch
  if __sphere:
    [ pp_x, pp_y, pp_z, pp_za ] = _sphere.calculate_pierce_point( [ float( x ) for x in xyz ],
        [ float( x ) for x in radec ], float( gst ), float( height ), int( iterations ) )
    return [ [ pp_x, pp_y, pp_z ], pp_za ]

  # initialize some variables
  ant_xyz = array( xyz, dtype = float64 )
  ant_geo_llh = xyz_to_geo_llh( xyz, iterations = iterations )
  ant_lon = radians( ant_geo_llh[ 0 ] )
  ant_lat = radians( ant_geo_llh[ 1 ] )
  ant_lh = ant_geo_llh[ 2 ]
  if ( ant_lh > height ):
    raise error( 'specified location has a height larger than the pierce layer' )
  rot = array( [ [ - sin( ant_lon )                 ,   cos( ant_lon )                 ,             0. ], 
                 [ - cos( ant_lon ) * sin( ant_lat ), - sin( ant_lon ) * sin( ant_lat ), cos( ant_lat ) ],
                 [   cos( ant_lon ) * cos( ant_lat ),   sin( ant_lon ) * cos( ant_lat ), sin( ant_lat ) ] ],
                 dtype = float64 )
  ant_za_az = calculate_local_sky_position( ant_geo_llh, radec, gst )
  ant_za = radians( ant_za_az[ 0 ] )
  ant_az = radians( ant_za_az[ 1 ] )
  len2_ant_xyz = ( ant_xyz**2 ).sum()
  local_src_dxyz = array( [ sin( ant_za ) * sin( ant_az ), sin( ant_za ) * cos( ant_az ), cos( ant_za ) ],
      dtype = float64 )
  src_dxyz = dot( local_src_dxyz, rot )

  # determine xyz coordinates of pierce point through vector algebra
  B = 2. * ( ant_xyz * src_dxyz ).sum()
  len2_pp_xyz = ( sqrt( len2_ant_xyz ) + ( height - ant_lh ) )**2
  for i in range( iterations ):
    C = len2_ant_xyz - len2_pp_xyz # always < 0
    len_src_xyz = ( sqrt( B**2 - 4. * C ) - B ) / 2. # always > 0
    src_xyz = len_src_xyz * src_dxyz
    pp_xyz = ant_xyz + src_xyz
    len_pp_xyz = sqrt( ( pp_xyz**2 ).sum() )
    pp_geo_llh = xyz_to_geo_llh( pp_xyz.tolist(), iterations = iterations )
    dlen_pp_xyz = height - pp_geo_llh[ 2 ]
    len2_pp_xyz = ( len_pp_xyz + dlen_pp_xyz )**2
  C = len2_ant_xyz - len2_pp_xyz # always < 0
  len_src_xyz = ( sqrt( B**2 - 4. * C ) - B ) / 2. # always > 0
  src_xyz = len_src_xyz * src_dxyz
  pp_xyz = ant_xyz + src_xyz

  # determine zenith angle at pierce point
  pp_geo_llh = xyz_to_geo_llh( pp_xyz.tolist(), iterations = iterations )
  [ radius, angle ] = calculate_angular_separation( ant_geo_llh[ 0 : 2 ], pp_geo_llh[ 0 : 2 ] )
  pp_za = ant_za_az[ 0 ] - radius

  return [ pp_xyz.tolist(), float( pp_za ) ]

###############################################################################

def calculate_offset_position_dradec( degdeg, offset ):
# offset in degrees
  dr = sqrt( offset[ 0 ]**2 + offset[ 1 ]**2 )
  dp = degrees( atan2( offset[ 0 ], offset[ 1 ] ) )
  radec = calculate_offset_position( degdeg, dr, dp )
  return radec

###############################################################################

def calculate_offset_position( degdeg, radius, angle ):
# 0. <= radius <= 180.
  if __sphere:
    return _sphere.calculate_offset_position( [ float( x ) for x in degdeg ],
        float( radius ), float( angle ) )
  ra = degdeg[ 0 ]
  dec = degdeg[ 1 ]
  if ( radius <= 0. ):
    new_ra = ra
    new_dec = dec
  else:
    a = radians( radius )
    c = radians( 90. - dec )
    B = radians( - angle )
    b = acos( max( - 1., min( 1., sin( a ) * cos( B ) * sin( c ) + cos( a ) * cos( c ) ) ) )
    if ( b == 0. ):
      A = 0.
    else:
      A = asin( max( -1., min( 1., sin( a ) * sin( B ) / sin( b ) ) ) )
      if ( ( ( cos( a ) * sin( c ) - sin( a ) * cos( B ) * cos( c ) ) / sin( b ) ) < 0. ):
        A = pi - A
    new_ra = amodulo( ra - degrees( A ), 360. )
    new_dec = 90. - degrees( b )
  return [ new_ra, new_dec ]

###############################################################################

def equatorial_to_galactic( degdeg, epoch = 2000 ):
  a = radians( degdeg[ 0 ] )
  d = radians( max( min( degdeg[ 1 ], 90. ), -90. ) )
  if ( epoch == 2000 ):
    radec_G = [ 192.859508, 27.128336 ]
    radec_B = calculate_offset_position( radec_G, 90., 122.932 )
  elif ( epoch == 1950 ):
    radec_G = hmsdms_to_degdeg( [ 12,49,0, 27.4,0,0 ] )
    radec_B = hmsdms_to_degdeg( [ 17,42.4,0, -28.92,0,0 ] )
  else:
    raise error( 'unknown epoch: %s' % ( repr( epoch ) ) )
  aG = radians( radec_G[ 0 ] )
  dG = radians( radec_G[ 1 ] )
  aB = radians( radec_B[ 0 ] )
  dB = radians( radec_B[ 1 ] )
  BK = acos( max( -1., min( 1., 
      sin( dB ) * cos( dG ) - cos( dB ) * sin( dG ) * cos( aG - aB ) ) ) )
  b = asin( max( -1., min( 1., 
      sin( dG ) * sin( d ) + cos( dG ) * cos( d ) * cos( a - aG ) ) ) )
  if ( abs( degrees( b ) ) == 90. ):
    l = 0.
    b = asign( b ) * 90.
  else:
    BKml = asin( max( -1., min( 1., cos( d ) * sin( a - aG ) / cos( b ) ) ) )
    if ( ( cos( dG ) * sin( d ) - sin( dG ) * cos( d ) * cos( a - aG ) ) / cos( b ) < 0. ):
      BKml = pi - BKml
    l = BK - BKml
  l = amodulo( degrees( l ), 360. )
  b = degrees( b )
  return [ l,b ]

###############################################################################

def galactic_to_equatorial( degdeg, epoch = 2000 ):
  l = radians( degdeg[ 0 ] )
  b = radians( max( min( degdeg[ 1 ], 90. ), -90. ) )
  if ( epoch == 2000 ):
    radec_G = [ 192.859508, 27.128336 ]
    radec_B = calculate_offset_position( radec_G, 90., 122.932 )
  elif ( epoch == 1950 ):
    radec_G = hmsdms_to_degdeg( [ 12,49,0, 27.4,0,0 ] )
    radec_B = hmsdms_to_degdeg( [ 17,42.4,0, -28.92,0,0 ] )
  else:
    raise error( 'unknown epoch: %s' % ( repr( epoch ) ) )
  aG = radians( radec_G[ 0 ] )
  dG = radians( radec_G[ 1 ] )
  aB = radians( radec_B[ 0 ] )
  dB = radians( radec_B[ 1 ] )
  BK = acos( max( -1., min( 1., 
      sin( dB ) * cos( dG ) - cos( dB ) * sin( dG ) * cos( aG - aB ) ) ) )
  d = asin( max( -1., min( 1., 
      sin( dG ) * sin( b ) + cos( dG ) * cos( b ) * cos( BK - l ) ) ) )
  if ( abs( degrees( d ) ) >= 90. ):
    a = 0.
    d = asign( d ) * 90.
  else:
    amaG = asin( max( -1., min( 1., cos( b ) * sin( BK - l ) / cos( d ) ) ) )
    if ( ( cos( dG ) * sin( b ) - sin( dG ) * cos( b ) * cos( BK - l ) ) / cos( d ) < 0. ):
      amaG = pi - amaG
    a = amaG + aG
  a = amodulo( degrees( a ), 360. )
  d = degrees( d )
  return [ a,d ]

###############################################################################


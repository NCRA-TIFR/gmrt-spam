###############################################################################

# import Python modules
from math import pi, atan2, degrees, radians, atan
from pylab import randn

# import 3rd party modules
from numpy import *

# import user modules
from error import *
try:
  import _acalc
except:
  __acalc = False
else:
  __acalc = True

###############################################################################
# NOTE: there are several cases where one needs to preserve the sign of zero,
# e.g. +0 and -0. 
###############################################################################

def factorial( n ):
  nn = int( n )
  if ( nn < 0 ):
    fac = None
  else:
    if __acalc:
      fac = _acalc.factorial( nn )
    else:
      fac = long( 1 )
      while ( nn > 0 ):
        fac = fac * long( nn )
        nn = nn - 1
  return fac

###############################################################################

def binomial( n, k ):
  nn = int( n )
  kk = int( k )
  if ( kk < 0 ) or ( kk > nn ):
    binom = None
  else:
    if __acalc:
      binom = _acalc.binomial( nn, kk )
    else:
      binom = factorial( nn ) / ( factorial( kk ) * factorial( nn - kk ) )
  return binom

###############################################################################

def complex_to_r_phi( c ):
  if ( __acalc and ( not is_array( c ) ) ):
    [ r, phi ] = _acalc.complex_to_r_phi( [ c.real, c.imag ] )
  else:
    r = abs( c )
    sel = awhere( array( r ) > 0. )
#    phi = degrees( log( c ).imag )
    phi = degrees( atan2( c.imag, c.real ) )
  return [ r, phi ]

###############################################################################

def r_phi_to_complex( rp ):
  if ( __acalc and ( not is_array( rp ) ) ):
    [ cr, ci ] = _acalc.r_phi_to_complex( rp )
    c = complex( cr, ci )
  else:
    [ r, phi ] = rp
    c = r * exp( complex( 0., 1. ) * radians( phi ) )
  return c

###############################################################################

def is_array( a ):
  return isinstance( a, type( array( [ 1 ] ) ) )

###############################################################################

def azeros( x ):
  if ( len( shape( x ) ) == 0 ):
    zero = 0.
  else:
    zero = zeros( shape = x.shape, dtype = x.dtype )
  return zero

###############################################################################

def aones( x ):
  if ( len( x.shape ) == 0 ):
    one = 1.
  else:
    one = ones( shape = x.shape, dtype = x.dtype )
  return one

###############################################################################

def aatan2( y, x ):
  if ( shape( x ) != shape( y ) ):
    raise error( 'x and y have different shapes' )
  if ( len( shape( x ) ) == 0 ):
    z = atan2( y, x )
  else:
    xx = x.ravel()
    yy = y.ravel()
    zz = array( [ atan2( yy[ i ], xx[ i ] ) for i in range( len( xx ) ) ], dtype = x.dtype )
    z = zz.reshape( x.shape )
  return z

###############################################################################

def asign( x ):
# this function also separates between -0 and +0
  if ( not is_array( x ) ):
    s = ( - 2. * float( aatan2( x, x ) < 0. ) + 1. )
  else:
    s = ( - 2. * array( aatan2( x, x ) < 0., dtype = x.dtype ) + 1. )
  return s

###############################################################################

def amodulo( x, y ):
  if ( not is_array( x ) ):
    if ( not is_array( y ) ):
      m = x - y * floor( x / ( y + float( y == 0. ) ) )
    else:
      xx = x * aones( y )
      m = xx - y * floor( x / ( y + array( y == 0., dtype = y.dtype ) ) )
  else:
    if ( not is_array( y ) ):
      yy = y * aones( x )
      m = x - yy * floor( x / ( yy + array( yy == 0., dtype = yy.dtype ) ) )
    else:
      m = x - y * floor( x / ( y + array( y == 0., dtype = y.dtype ) ) )
  return m

###############################################################################

def aradians( x ):
  r = x * ( pi / 180. )
  return r

###############################################################################

def adegrees( x ):
  r = x * ( 180. / pi )
  return r

###############################################################################

def awhere( a ):
   return transpose( array( where( a ) ) )

###############################################################################

def aput( data, sel, sub ):
# data, sel and sub must be arrays

  # check input dimensions
  if ( len( sel ) == 0 ):
    return data.copy()
  if ( len( sel.ravel() ) == 0 ):
    return data.copy()
  if ( sel.shape[ 1 ] > len( data.shape ) ):
    raise error( 'number of dimensions of index array is higher than that of data array' )
  asub = array( sub )
  if ( len( asub.shape ) == len( data.shape ) - sel.shape[ 1 ] ):
    if ( asub.shape != data.shape[ sel.shape[ 1 ] : ] ):
      raise error( 'shape of subarray does not match selected data' )
    asub = resize( asub, [ sel.shape[ 0 ] ] + list( data.shape[ sel.shape[ 1 ] : ] ) )
  elif ( len( asub.shape ) == len( data.shape ) - sel.shape[ 1 ] + 1 ):
    if ( list( asub.shape ) != [ sel.shape[ 0 ] ] + list( data.shape[ sel.shape[ 1 ] : ] ) ):
      raise error( 'shape of subarray does not match selected data' )

  # collapse and write data
  coffset = [ int( product( data.shape[ i : sel.shape[ 1 ] ] ) ) for i in range( 1, 1 + sel.shape[ 1 ] ) ]
  coffset[ - 1 ] = 1
  csel = add.reduce( sel * array( coffset ), 1 )
  subsize = int( product( data.shape[ sel.shape[ 1 ] : ] ) )
  suboffset = resize( arange( subsize ), [ sel.shape[ 0 ], subsize ] )
  csel = ( transpose( resize( csel * subsize, [ subsize, sel.shape[ 0 ] ] ) ) + suboffset ).ravel()
  cdata = data.copy().ravel()
  csub = asub.ravel()
  put( cdata, csel, csub )

  return cdata.reshape( data.shape )

###############################################################################

def aget( data, sel ):
# data and sel must be arrays

  # check input dimensions
  if ( len( sel ) == 0 ):
    return array( [], dtype = data.dtype )
  if ( len( sel.ravel() ) == 0 ):
    return array( [], dtype = data.dtype )
  if ( sel.shape[ 1 ] > len( data.shape ) ):
    raise error( 'number of dimensions of index array is higher than that of data array' )
  
  # collapse data along sel axes
  cdata_len = int( product( data.shape[ 0 : sel.shape[ 1 ] ] ) )
  cdata_shape = [ cdata_len ] + list( data.shape[ sel.shape[ 1 ] : ] )
  cdata = data.reshape( cdata_shape )
  coffset = [ int( product( data.shape[ i : sel.shape[ 1 ] ] ) ) for i in range( 1, 1 + sel.shape[ 1 ] ) ]
  coffset[ - 1 ] = 1
  csel = add.reduce( sel * array( coffset ), 1 )

  return take( cdata, csel, axis = 0 ) 

###############################################################################

def acount( sel ):
# data and sel must be arrays
  
  # check input dimensions
  if ( len( sel ) == 0 ):
    return array( [], dtype = data.dtype )
  if ( len( sel.ravel() ) == 0 ):
    return array( [], dtype = data.dtype )
  if ( len( sel.shape ) != 2 ):
    raise error( 'number of dimensions of index array is not two' )
  
  # count occurences
  sel_code = array( [ 0 ] + list( ( sel + 1 ).max( axis = 0 ).cumprod()[ : -1 ] ) )
  sel_encode = dot( sel, sel_code )
  counts = aget( bincount( sel_encode ), sel_encode.reshape( ( len( sel_encode ), 1 ) ) )
  
  return counts

###############################################################################

def amean_phase_old( data ):
  # determine phase average (with trick to get around possible phase wrap problems)
  phases = array( [ data ], dtype = float64 ).ravel()
  offsets = [ 0., 90., 180., 270. ]
  mean_phases = [ ( amodulo( ( phases + offsets[ j ] ) + 180., 360. ) - 180. ).mean()
      for j in range( len( offsets ) ) ]
  var_phases = [ ( ( ( amodulo( ( phases + offsets[ j ] ) + 180., 360. ) - 180. ) - 
      mean_phases[ j ] )**2 ).mean() for j in range( len( offsets ) ) ]
  j = var_phases.index( min( var_phases ) )
  mean_phase = mean_phases[ j ] - offsets[ j ]
  return float( amodulo( mean_phase + 180, 360. ) - 180. )

###############################################################################

def amean_phase( data ):
  phases = array( [ data ], dtype = float64 ).ravel()
  # unwrap the phase (phase jumps greater than 180 are made smaller by adding n*360)
  if ( len( phases ) < 20 ):
    phase_array = adegrees( unwrap( aradians( data ) ) )
  else:
    phase_array = aunwrap_phase( phases, alpha = 0.01 )
    if sometrue( isnan( phase_array ) ):
      phase_array = aunwrap_phase( phases, alpha = 0.001 )
    if ( sometrue( isnan( phase_array ) ) or ( max( fabs( phase_array ) ) > 1.e4 ) ):
      phase_array = adegrees( unwrap( aradians( phases ) ) )
  if sometrue( isnan( phase_array ) ):
    return amean_phase_old( data )
  else:
    return float( amodulo( phase_array.mean() + 180, 360. ) - 180. )

###############################################################################

def amedian_phase_old( data ):
  # determine phase average (with trick to get around possible phase wrap problems)
  mean_phase = amean_phase( data )
  phases = array( [ data ], dtype = float64 ).ravel()
  median_phase = median( amodulo( ( phases - mean_phase ) + 180., 360. ) - 180. )
  median_phase = amodulo( ( median_phase + mean_phase ) + 180., 360. ) - 180.
  return float( median_phase )

###############################################################################

def amedian_phase( data ):
  phases = array( [ data ], dtype = float64 ).ravel()
  # unwrap the phase (phase jumps greater than 180 are made smaller by adding n*360)
  if ( len( phases ) < 20 ):
    phase_array = adegrees( unwrap( aradians( data ) ) )
  else:
    phase_array = aunwrap_phase( phases, alpha = 0.01 )
    if sometrue( isnan( phase_array ) ):
      phase_array = aunwrap_phase( phases, alpha = 0.001 )
    if ( sometrue( isnan( phase_array ) ) or ( max( fabs( phase_array ) ) > 1.e4 ) ):
      phase_array = adegrees( unwrap( aradians( phases ) ) )
  if sometrue( isnan( phase_array ) ):
    return amedian_phase_old( data )
  else:
    return float( amodulo( median( phase_array ) + 180, 360. ) - 180. )

###############################################################################

def aunwrap_phase( x, window = 10, alpha = 0.01, iterations = 3, 
    clip_range = [ 170., 180. ] ):
  xx = array( x, dtype = float64 )
#  a = zeros( ( window ), dtype = float64 )
#  a[ -1 ] = 1.
  if ( len( clip_range ) == 2 ):
    o = clip_range[ 0 ]
    s = ( clip_range[ 1 ] - clip_range[ 0 ] ) / 90.
  a = ones( ( window ), dtype = float64 ) / float( window )
  xs = xx[ 0 ]
  for j in range( 2 * iterations ):
    for k in range( window, len( x ) ):
      xi = xx[ k - window : k ]
      xp = dot( xi, a )
      e = xx[ k ] - xp
      e = amodulo( e + 180., 360. ) - 180.
      if ( len( clip_range ) == 2 ):
        if ( abs( e ) > o ):
          e = sign( e ) * ( s * degrees( atan( radians( ( abs( e ) - o ) / s ) ) ) + o )
      xx[ k ] = xp + e
#      a = a + xx[ k - window : k ] * alpha * e / 360.
      a = a + xi * alpha * e / ( dot( xi, xi ) + 1.e-6 )
    xx = xx[ : : -1 ].copy()
  xx = xx - xx[ 0 ] + xs
  return xx

###############################################################################

def gaussian_2d( amp, bmaj, bmin, bpa, r, phi ):
  theta = radians( phi - bpa )
  x = r * sin( theta )
  y = r * cos( theta )
  dx = bmin / sqrt( 8. * log( 2. ) )
  dy = bmaj / sqrt( 8. * log( 2. ) )
  alpha = ( x / dx )**2 + ( y / dy )**2
  f = amp * exp( - alpha / 2.)
  return f

###############################################################################

def draw_from_gaussian_distribution( mean, stddev, times = 10000 ):
  draw = float64( randn( times ) * stddev + mean )
  if ( times == 1 ):
    draw = draw[ 0 ]
  return draw

###############################################################################

def get_robust_mean_deviations( data, rejection = None, min_data = 10 ):
  sel_data = data.copy()
  while ( len( sel_data ) > min_data ):
#    print len( sel_data )
    data_median = median( sel_data )
    sel_plus = awhere( sel_data > data_median )
    data_plus = aget( sel_data, sel_plus ).tolist()
    data_plus.sort( cmp = lambda a, b: cmp( a, b ) )
    cut_plus = float( len( sel_plus ) ) * exp( -0.5 )
    cut_index = int( floor( cut_plus ) )
    cut_fraction = cut_plus - float( cut_index )
    sigma_plus = ( data_plus[ cut_index ] * ( 1. - cut_fraction ) + 
        data_plus[ cut_index + 1 ] * cut_fraction ) - data_median
    sel_minus = awhere( sel_data < data_median )
    data_minus = aget( sel_data, sel_minus ).tolist()
    data_minus.sort( cmp = lambda a, b: cmp( b, a ) )
    cut_minus = float( len( sel_minus ) ) * exp( -0.5 )
    cut_index = int( floor( cut_minus ) )
    cut_fraction = cut_minus - float( cut_index )
    sigma_minus = ( data_minus[ cut_index ] * ( 1. - cut_fraction ) + 
        data_minus[ cut_index + 1 ] * cut_fraction ) - data_median
    if ( rejection is None ):
      break
    sel = awhere( ( sel_data < data_median + rejection * sigma_plus ) &
        ( sel_data > data_median + rejection * sigma_minus ) )
    if ( len( sel ) == len( sel_data ) ):
      break
    sel_data = aget( sel_data, sel )
  if ( len( sel_data ) <= min_data ):
    raise error( 'data count too low for useful statistics' )
  return [ data_median, sigma_plus, sigma_minus ]

###############################################################################


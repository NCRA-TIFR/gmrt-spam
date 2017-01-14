###############################################################################

# import Python modules
from sys import *
from os import *
from datetime import *
from math import *
import pdb

# import 3rd party modules
from numpy import *
from pylab import *
#from scipy import *
#import scipy.signal as ss

# import user modules
from files import *
from aips import *
from acalc import *
from sphere import *
from parameter import *
from mpfit import *
from solutions import *
from plot import *
from error import *
#from instrumental import *
try:
  import _ionosphere
except:
  __ionosphere = False
else:
  __ionosphere = True

###############################################################################

def fit_phi_model( P, dojac = None, phi_model = None, X_table = None, 
    pza_table = None, ref_table = None, phase_table = None, error_table = None,
    normalize_weights = None, wrap_phases = None ):

  # handle input parameters
  if ( ( not dojac is None ) or ( phi_model is None ) or ( X_table is None ) or 
      ( pza_table is None ) or ( ref_table is None ) or ( phase_table is None ) or
      ( error_table is None ) or ( normalize_weights is None ) or
      ( wrap_phases is None ) ):
    return - 1, None, None

  # calculate phases at pierce points
  phi_table = phi_model( X_table, P ) / cos( aradians( pza_table ) )

  # split pierce points per source
  ref_list = []
  for ref in ref_table:
    if ( not ( ref in ref_list ) ):
      ref_list.append( ref )

  # calculate chi2 terms
  chi_list = []
  if normalize_weights:
    normalize_list = []
  for ref in ref_list:
    sel = awhere( ref_table == ref )
    dphase_table = aget( phase_table, sel )
    dphase_table = resize( dphase_table, ( len( sel ), len( sel ) ) )
    dphase_table = transpose( dphase_table ) - dphase_table
    dphi_table = aget( phi_table, sel )
    dphi_table = resize( dphi_table, ( len( sel ), len( sel ) ) )
    dphi_table = transpose( dphi_table ) - dphi_table
    dweight_table = aget( error_table, sel )
    dweight_table = resize( dweight_table**2, ( len( sel ), len( sel ) ) )
    dweight_table = transpose( dweight_table ) + dweight_table
    dweight_table = 1. / sqrt( dweight_table )
    if wrap_phases:
      chi_table = dweight_table * ( amodulo( ( dphase_table - dphi_table ) + 180.,
          360. ) - 180. )
    else:
      chi_table = dweight_table * ( dphase_table - dphi_table )
    dummy = [ [ chi_list.append( chi_table[ b, a ] )
        for a in range( b + 1, len( sel ) ) ] for b in range( len( sel ) ) ]
    if normalize_weights:
      dummy = [ [ normalize_list.append( dweight_table[ b, a ] )
          for a in range( b + 1, len( sel ) ) ] for b in range( len( sel ) ) ]

  # make (normalized) chi2 array
  chi_array = array( chi_list, dtype = chi_table.dtype )
  if normalize_weights:
    normalize_array = array( normalize_list, dtype = dweight_table.dtype )
    factor = 1. / sqrt( ( normalize_array**2 ).sum() )
    chi_array = factor * chi_array

  return 0, chi_array, None

###############################################################################

def fit_phi_mkl_model( P, dojac = None, U_table = None, ref_table = None,
    phase_table = None, error_table = None, normalize_weights = None ):
  
  # handle input parameters
  if ( ( not dojac is None ) or
       ( U_table is None ) or
       ( ref_table is None ) or
       ( phase_table is None ) or
       ( error_table is None ) or
       ( normalize_weights is None ) ):
    return - 1, None, None
  
  # calculate model phases
  phi_table = dot( U_table, P )

  # split pierce points per source
  ref_list = []
  for ref in ref_table:
    if ( not ( ref in ref_list ) ):
      ref_list.append( ref )
  
  # calculate chi2 terms
  chi_list = []
  if normalize_weights:
    normalize_list = []
  for ref in ref_list:
    index_list = []
    for j in range( len( ref_table ) ):
      if ( ref_table[ j ] == ref ):
        index_list.append( j )
    ref_count = len( index_list )
    dphase_table = array( [ phase_table[ index_list[ a ] ] 
        for a in range( ref_count ) ], dtype = phase_table.dtype )
    dphase_table = resize( dphase_table, ( ref_count, ref_count ) )
    dphase_table = transpose( dphase_table ) - dphase_table
    dphi_table = array( [ phi_table[ index_list[ a ] ]
        for a in range( ref_count ) ], dtype = phi_table.dtype )
    dphi_table = resize( dphi_table, ( ref_count, ref_count ) )
    dphi_table = transpose( dphi_table ) - dphi_table
    dweight_table = array( [ error_table[ index_list[ a ] ]
        for a in range( ref_count ) ], dtype = error_table.dtype )
    dweight_table = resize( dweight_table**2, ( ref_count, ref_count ) )
    dweight_table = transpose( dweight_table ) + dweight_table
    dweight_table = 1. / sqrt( dweight_table )
    chi_table = dweight_table * ( amodulo( ( dphase_table - dphi_table ) + 180.,
        360. ) - 180. )
    dummy = [ [ chi_list.append( chi_table[ b, a ] )
        for a in range( b + 1, ref_count ) ] for b in range( ref_count ) ]
    if normalize_weights:
      dummy = [ [ normalize_list.append( dweight_table[ b, a ] )
          for a in range( b + 1, ref_count ) ] for b in range( ref_count ) ]
  
  # make (normalized) chi2 array
  chi_array = array( chi_list, dtype = chi_table.dtype )
  if normalize_weights:
    normalize_array = array( normalize_list, dtype = dweight_table.dtype )
    factor = 1. / sqrt( ( normalize_array**2 ).sum() )
    chi_array = factor * chi_array
  
  return 0, chi_array, None

###############################################################################

def phi_gradient_model( X, p ):
  phi = dot( X, p )
  return phi

###############################################################################

def phi_poly_model( X, p ):
  Y = transpose( array( [ X[ : , 0 ], X[ : , 1 ],
      X[ : , 0 ]**2, X[ : , 0 ] * X[ : , 1 ], X[ : , 1 ]**2 ] ) )
  phi = dot( Y, p )
  return phi

###############################################################################

def phi_mkl_model( layer_weights, Xl_table, zal_table, Xpl_table, pzal_table, 
    Bl_table, F_table, beta = 5. / 3., r_0 = 1. ):
# B_table = (1/m)(1T)(C_table)(A_table)
# F_table = ( Ci_table )( U_table )( P_table )

  x_count = len( Xl_table[ 0 ] )
  p_count = len( Xpl_table[ 0 ] )

  for l in range( len( layer_weights ) ):

    # calculate structure matrix
    D_table = transpose( resize( Xl_table[ l ], ( p_count, x_count, 2 ) ),
        ( 1, 0, 2 ) )
    D_table = D_table - resize( Xpl_table[ l ], ( x_count, p_count, 2 ) )
    D_table = add.reduce( D_table**2, 2 )
    D_table = ( D_table / ( r_0**2 ) )**( beta / 2. )

    # calculate covariance matrix
    C_table = - D_table / 2.
    C_table = transpose( transpose( C_table ) - 
        ( add.reduce( C_table, 1 ) / float( p_count ) ) )
    C_table = C_table - Bl_table[ l ]

    # incorporate airmass functions and layer weights
    # add layer C to total C
    A_table = resize( 1. / cos( aradians( pzal_table[ l ] ) ), ( x_count, p_count ) )
    A_table = A_table * transpose( resize( 1. / cos( aradians( zal_table[ l ]) ),
        ( p_count, x_count ) ) )
    if ( l == 0 ):
      Cl_table = C_table * A_table * ( layer_weights[ l ]**2 )
    else:
      Cl_table = Cl_table + C_table * A_table * ( layer_weights[ l ]**2 )

  phi = dot( Cl_table, F_table )
  phi = reshape( phi, ( x_count ) )

  return phi

###############################################################################

def get_phase_calibration_data( uv, facets, facet_list = [], time_info = True,
    source_info = True, antenna_info = True, calibration_info = True,
    solution_version = 0, print_info = True ):

# TODO: possibly include delays into model fitting

  if ( source_info or calibration_info ):
    if ( len( facet_list ) > 0 ):
      source_count = len( facet_list )
      source_list = facet_list
    else:
      source_count = restore_parameter( facets, 'facet_count' )
      source_list = range( 1, 1 + source_count )

  # get time info
  time_list = get_time_list( uv )
  time_count = len( time_list )
  time_table = None
  if time_info:
    if print_info:
      print 'getting time info ...'
    gst_list = get_gst_list( uv, time_list = time_list )
    time_table = [ [ time_list[ n ], gst_list[ n ] ] for n in range( time_count ) ]

  # get source info
  center_table = None
  source_table = None
  if source_info:
    if print_info:
      print 'getting source info ...'
    obs_epoch = get_observing_epoch( uv )
    center_table = convert_radec_from_j2000( get_radec( uv ), obs_epoch )
    source_table = [ 
        convert_radec_from_j2000( get_radec( get_facet( facets, source ) ), obs_epoch )
        for source in source_list ]

  # get antenna info
  array_table = None
  antenna_table = None
  antenna_list = get_antenna_positions( uv )
  antenna_no_list = [ antenna[ 0 ] for antenna in antenna_list ]
  antenna_count = antenna_no_list[ - 1 ]
  if antenna_info:
    if print_info:
      print 'getting antenna info ...'
    array_xyz = get_mean_antenna_position( uv ) # get_array_position( uv )
    array_table = [ array_xyz, xyz_to_geo_llh( array_xyz ) ]
    antenna_table = []
    for antenna_no in range( 1, 1 + antenna_count ):
      if ( antenna_no in antenna_no_list ):
        antenna_index = antenna_no_list.index( antenna_no )
        antenna_xyz = antenna_list[ antenna_index ][ 1 ]
        antenna_table.append( [ antenna_xyz, xyz_to_geo_llh( antenna_xyz ) ] )
      else:
        antenna_table.append( [ [ 0., 0., 0. ], [ 0., 0., 0. ] ] )

  # read calibration info from peel facets
  reference_table = None
  phase_table = None
  error_table = None
  if calibration_info:
    if print_info:
      print 'getting calibration info ...'
    phase_array = 360. * ones( shape = ( time_count, source_count, antenna_count ),
        dtype = float64 )
    error_array = - 360. * ones( shape = ( time_count, source_count, antenna_count ),
        dtype = float64 )
    reference_array = zeros( shape = ( time_count, source_count ), dtype = int32 )
    k = 0
    for source in source_list:
      if print_info:
        print '... from source %d ...' % ( source )
      facet_k = get_facet( facets, source )
      if ( solution_version <= 0 ):
        sol_version = facet_k.table_highver( 'SN' ) + solution_version
      else:
        sol_version = solution_version
      solution_table = read_solution_table( facet_k, in_version = sol_version )
      n = 0
      w = 0
      for solution in solution_table:
        [ time, reference_antenna, dummy, dummy2 ] = solution[ 0 ]
        try:
          n = time_list.index( time )
        except:
          w = w + 1
        reference_array[ n ][ k ] = reference_antenna - 1
        [ ref_gain_real, ref_gain_imag, ref_delay, ref_weight ] = solution[ 
            reference_antenna ]
        ref_gain = complex( ref_gain_real, ref_gain_imag )
        if ( ( ref_weight == 0. ) or ( ref_gain == complex( 0., 0. ) ) ):
#          raise error( 'reference antenna has no solution for time instance' )
          # skip time instance altogether
          continue
        [ ref_amp, ref_phase ] = complex_to_r_phi( ref_gain )
        for i in range( antenna_count ):
          [ gain_real, gain_imag, delay, weight ] = solution[ i + 1 ]
          gain = complex( gain_real, gain_imag )
          if ( ( weight == 0. ) or 
              ( ( i != reference_array[ n ][ k ] ) and 
              ( gain == complex( 1., 0. ) ) ) or 
              ( gain == complex( 0., 0. ) ) ):
            phase = 360.
            phase_error = - 360.
          else:
            [ amp, phase ] = complex_to_r_phi( gain )
            phase_error = degrees( 1. / weight )
          phase_array[ n ][ k ][ i ] = phase
          error_array[ n ][ k ][ i ] = phase_error
      k = k + 1
      if ( print_info and ( w > 0 ) ):
        print ( 'WARNING: %d time(s) in solution table do(es) not ' % ( w ) + 
          'correspond to time(s) in observations' )
    reference_table = reference_array.tolist()
    phase_table = phase_array.tolist()
    error_table = error_array.tolist()

  return [ time_table, center_table, source_table, array_table, antenna_table, 
      reference_table, phase_table, error_table ]

###############################################################################

def calculate_pierce_coordinates( time, center_table, source_table, array_table,
    antenna_table, height = 400.e3, iterations = 4 ):
# source table radecs at observing epoch

  pierce_table = []

  # calculate Earth referenced coordinates of pierce points of array center
  # towards pointing center
  [ center_pxyz, center_pza ] = calculate_pierce_point( array_table[ 0 ], 
      center_table, time[ 1 ], height = height, iterations = iterations )
  center_p_geo_llh = xyz_to_geo_llh( center_pxyz, iterations = iterations )

  # loop over sources
  for k in range( len( source_table ) ):

    # loop over antennas
    for i in range( len( antenna_table ) ):

      # only handle active antennas
      if ( antenna_table[ i ] != [ [ 0., 0., 0. ], [ 0., 0., 0. ] ] ):

        # calculate Earth referenced coordinates of pierce points of antenna
        # towards peeled source
        [ pxyz, pza ] = calculate_pierce_point( antenna_table[ i ][ 0 ],
            source_table[ k ], time[ 1 ], height = height, iterations = iterations )
        p_geo_llh = xyz_to_geo_llh( pxyz, iterations = iterations )

        # calculate local angular coordinates of antenna pierce point 
        # ( x = East, y = North )
        [ radius, angle ] = calculate_angular_separation( center_p_geo_llh[ 0 : 2 ], 
            p_geo_llh[ 0 : 2 ] )
        X = [  radius * sin( radians( angle ) ), radius * cos( radians( angle ) ) ]

        # store model fit input data
        pierce_info = [ X, pza, [ k, i ] ]
        pierce_table.append( pierce_info )

  return pierce_table

###############################################################################

def calculate_phase_gradients( phase_table ):
#  if __ionosphere:
#    return _ionosphere.calculate_phase_gradients( 
#        [ phase[ 0 ].tolist() + phase[ 1 : 3 ] for phase in phase_table ] )

  phase_gradient_table = []
  phase_count = len( phase_table )
  for i in range( phase_count ):
    [ X_i, phase_i, error_i ] = phase_table[ i ]
    for j in range( i + 1, phase_count ):
      [ X_j, phase_j, error_j ] = phase_table[ j ]
      [ r_ij, angle_ij ] = calculate_angular_separation( X_i, X_j )
      phase_ij = amodulo( ( phase_j - phase_i ) + 180., 360. ) - 180.
      error_ij = sqrt( error_i**2 + error_j**2 )
      phase_gradient_table.append( [ r_ij, angle_ij, phase_ij, error_ij ] )
      [ r_ji, angle_ji ] = [ r_ij, amodulo( ( angle_ij + 180. ) + 180., 360. ) - 180. ]
      phase_ji = amodulo( ( phase_i - phase_j ) + 180., 360. ) - 180.
      error_ji = error_ij
      phase_gradient_table.append( [ r_ji, angle_ji, phase_ji, error_ji ] )
  return phase_gradient_table

###############################################################################

def estimate_phase_gradient_through_fft( uv, X_table, ref_table, phase_table,
    error_table, fft_cell_size = 200., fft_image_size = 1024, moment_radius = 128,
    peak_ratio = 0.75, warning = 0.05 ):
# this function returns the PROJECTED phase gradient
  
  # calculate phase differences
  gradient_table = []
  ref = ref_table[ 0 ]
  local_phase_table = []
  for j in range( len( X_table ) ):
    if ( ref_table[ j ] != ref ):
      if ( len( local_phase_table ) > 0 ):
        local_phase_gradient_table = calculate_phase_gradients( local_phase_table )
        for phase_gradient in local_phase_gradient_table:
          gradient_table.append( phase_gradient )
      ref = ref_table[ j ]
      local_phase_table = []
    local_phase_table.append( [ X_table[ j ], phase_table[ j ] , error_table[ j ] ] )
  if ( len( local_phase_table ) > 0 ):
    local_phase_gradient_table = calculate_phase_gradients( local_phase_table )
    for phase_gradient in local_phase_gradient_table:
      gradient_table.append( phase_gradient )

  # estimate gradient angle
  gradient_array = array( gradient_table, dtype = phase_table.dtype )
  sel = awhere( gradient_array[ : , 0 ] > 0. )
  w = abs( 1. / aget( gradient_array[ : , 0 ] * gradient_array[ : , 3 ], sel ) )
  w = aones( w )
  x = aget( gradient_array[ : , 0 ] * sin( aradians( gradient_array[ : , 1 ] ) ), sel )
  y = aget( gradient_array[ : , 0 ] * cos( aradians( gradient_array[ : , 1 ] ) ), sel )
  z = aget( gradient_array[ : , 2 ], sel )
  
  # create blank image
  try:
    fft_im_real = get_aips_file( uv.disk, 'FFT', 'ICL001', 0, 'MA' )
    fft_im_real_beam = get_facet_beam( fft_im_real )
    if fft_im_real_beam.exists():
      fft_im_real_beam.zap()
    if ( get_image_size( fft_im_real ) != [ fft_image_size, fft_image_size ] ):
      fft_im_real.zap()
      raise error( 'provoke remake of fft image' )
    fill_image( fft_im_real, do_all = True, value = 0. )
  except:
    fft_im_real = get_aips_file( uv.disk, 'FFT', 'ICL001', - 1, 'MA' )
    channel_count = get_channel_count( uv )
    call_aips_task( 'IMAGR', indata = uv, nchav = channel_count, 
        outdisk = fft_im_real.disk, outname = fft_im_real.name,
        outseq = fft_im_real.seq, outver = 0, allok = 0,
        cellsize = [ 1., 1. ], imsize = [ fft_image_size, fft_image_size ],
        do3dimag = 0, niter = 100000, flux = 100000., 
        fldsize = [ [ fft_image_size, fft_image_size ] ], nfield = 1,
        cmethod = '', minpatch = fft_image_size - 1, dotv = 0, overlap = 2,
        uvsize = [ fft_image_size, fft_image_size ], gain = 0.1 )
    fill_image( fft_im_real, do_all = True, value = 0. )
    fft_im_real.zap_table( 'CC', 0 )
    fft_im_real_beam = get_facet_beam( fft_im_real )
    fft_im_real_beam.zap()
  
  # create blank UV images
  fft_uv_real = get_aips_file( uv.disk, 'FFT', 'UVREAL', -1, 'MA' )
  fft_uv_imag = get_aips_file( uv.disk, 'FFT', 'UVIMAG', -1, 'MA' )
  if ( fft_uv_real.seq > fft_uv_imag.seq ):
    fft_uv_imag.seq = fft_uv_real.seq
  elif ( fft_uv_real.seq < fft_uv_imag.seq ):
    fft_uv_real.seq = fft_uv_imag.seq
  call_aips_task( 'FFT', indata = fft_im_real, opcode = 'MARE',
      outdisk = fft_uv_real.disk, outname = fft_uv_real.name, outseq = fft_uv_real.seq )
  
  # fill UV images with phase info
  if ( not fft_cell_size is None ):
    delta = 1. / fft_cell_size
  else:
    delta = abs( aget( gradient_array[ : , 0 ], sel ) ).max() / float( 
        ( fft_image_size / 2 ) - 10 )
  u = array( around( y / delta ), dtype = int32 ) + ( fft_image_size / 2 )
  v = array( around( x / delta ), dtype = int32 ) + ( fft_image_size / 2 )
  real_pixels = zeros( ( fft_image_size, fft_image_size ) , dtype = float32 )
  imag_pixels = zeros( ( fft_image_size, fft_image_size ) , dtype = float32 )
  wght_pixels = zeros( ( fft_image_size, fft_image_size ) , dtype = float32 )
  for i in range( len( sel ) ):
    if ( ( u[ i ] >= 0. ) and ( u[ i ] < fft_image_size ) and
         ( v[ i ] >= 0. ) and ( v[ i ] < fft_image_size ) ):
      c = r_phi_to_complex( [ 1., z[ i ] ] )
      w1 = wght_pixels[ u[ i ], v[ i ] ]
      w2 = w[ i ]
      wght_pixels[ u[ i ], v[ i ] ] = w1 + w2
      real_pixels[ u[ i ], v[ i ] ] = ( real_pixels[ u[ i ], v[ i ] ] * w1 + 
          c.real * w2 ) / ( w1 + w2 )
      imag_pixels[ u[ i ], v[ i ] ] = ( real_pixels[ u[ i ], v[ i ] ] * w1 + 
          c.imag * w2 ) / ( w1 + w2 )
  wiz_fft_uv_real = wizardry( fft_uv_real )
  plane = wiz_fft_uv_real[ 0 ]
  plane.pixels = transpose( real_pixels )
  plane.update()
  wiz_fft_uv_imag = wizardry( fft_uv_imag )
  plane = wiz_fft_uv_imag[ 0 ]
  plane.pixels = transpose( imag_pixels )
  plane.update()
  
  # FFT UV images to images
  fft_ma_real = get_aips_file( uv.disk, 'FFT', 'MAREAL', -1, 'MA' )
  fft_ma_imag = get_aips_file( uv.disk, 'FFT', 'MAIMAG', -1, 'MA' )
  if ( fft_ma_real.seq > fft_ma_imag.seq ):
    fft_ma_imag.seq = fft_ma_real.seq
  elif ( fft_ma_real.seq < fft_ma_imag.seq ):
    fft_ma_real.seq = fft_ma_imag.seq
  call_aips_task( 'FFT', indata = fft_uv_real, in2data = fft_uv_imag, opcode = 'UVCX',
      outdisk = fft_ma_real.disk, outname = fft_ma_real.name, outseq = fft_ma_real.seq )
  
  # determine gradient
  [ max_val, max_pos ] = get_image_maximum( fft_ma_real )
  pixels = get_image_pixels( fft_ma_real )
  sel = awhere( ( pixels != get_aips_magic_value() ) & 
      ( pixels > peak_ratio * max_val ) )
  rxy = sel + 1. - array( max_pos, dtype = sel.dtype )
  sel2 = awhere( add.reduce( rxy**2, axis = 1 ) < float( moment_radius**2 ) )
  sel = aget( sel, sel2 )
  mom = dot( aget( pixels, sel ), sel ) / aget( pixels, sel ).sum() + 1.
  max_pos = mom - array( get_pixel_reference( fft_ma_real ), dtype = mom.dtype )
  gradient = - max_pos / ( 2. * sqrt( 2. ) * delta )
  
  # check boundaries
  limits = [ warning * float( fft_image_size ), 
      ( 1 - warning ) * float( fft_image_size ) ]
  if ( ( min( mom ) < limits[ 0 ] ) or ( max( mom ) > limits[ 1 ] ) ):
    print 'WARNING: maximum detected near FFT grid boundaries'
    print limits, mom
  
  # cleanup
  fft_uv_real.zap()
  fft_uv_imag.zap()
  fft_ma_real.zap()
  fft_ma_imag.zap()
  
  return gradient

###############################################################################

def fit_ionospheric_mkl_model( uv, facets, facet_list = [], beta = 5. / 3., 
    r_0 = 1., order = 15, layer_heights = [ 200.e3 ], layer_weights = [ 1. ], 
    iterations = 4, equal_weights = True, normalize_weights = True, 
    e_steps = 1, double_precision = True, print_info = True,
    remove_gradient = False, extract_gradient = False,
    estimate_offsets = False, time_steps = [], include_airmass = True, 
    solution_version = 0, source_rms_limit = 1.e6,
    antenna_offset_limit = 1.e6, antenna_rms_limit = 1.e6,
    prop_time_limit = 1., prop_rms_limit = 0. ):
# prop_time_limit in minutes
  
  if double_precision:
    dtype = float64
  else:
    dtype = float32
  
  # select appropriate model functions
  # initialize model parameters
  model = 'mkl'
  layer_count = len( layer_heights )
  P = zeros( shape = ( order ), dtype = dtype )
  if ( len( layer_weights ) != layer_count ):
    raise error( 'layer heights and weights are different lengths' )
  
  if print_info:
    print 'starting model fits using model = ' + model
  
  # get all relevant data
  calibration_data = get_phase_calibration_data( uv, facets, facet_list = facet_list,
      solution_version = solution_version, print_info = print_info )
  time_table = calibration_data[ 0 ]
  center_table = calibration_data[ 1 ]
  peel_table = calibration_data[ 2 ]
  array_table = calibration_data[ 3 ]
  antenna_table = calibration_data[ 4 ]
  reference_array = calibration_data[ 5 ]
  phase_array = calibration_data[ 6 ]
  error_array = calibration_data[ 7 ]
  source_count = len( peel_table )
  antenna_count = len( antenna_table )
  
  fit_table = []
  pierce_table = []
  
  if estimate_offsets:
    solution_table = []

  prop_last_time = -1.e9
  prop_last_G = zeros( ( 2 ), dtype = dtype )
  
  # loop over time stamps
  if ( len( time_steps ) == 0 ):
    n_list = range( len( time_table ) )
  else:
    n_list = time_steps
  for n in n_list:

    if print_info:
      print '... time step n = ' + repr( n )
    
    repeat_fit = True
    rejected_sources = False
    rejected_antennas = False
    while repeat_fit:
      repeat_fit = False
      if ( not rejected_sources ):
        use_sources = [ True for k in range( source_count ) ]
      if ( not rejected_antennas ):
        use_antennas = [ True for i in range( antenna_count ) ]
      
      temp_pierce_table = []

      # containers for model input data
      ref_list = []
      ref_table = []
      phase_table = []
      error_table = []
      info_table = []
      Xpl_table = []
      pzal_table = []
      Bl_table = []

      if estimate_offsets:
        solution_row = [ [ time_table[ n ][ 0 ], 
            reference_array[ n ][ 0 ] + 1, 0., 0. ] ] + [ [ 0., 0., 0., 0. ] 
            for i in range( antenna_count ) ]

      # loop over height layers
      for l in range( layer_count ):
        height = layer_heights[ l ]
      
        # containers for model input data
        Xp_table = []
        pza_table = []
      
        # get pierce point coordinates
        pierce_table = calculate_pierce_coordinates( time_table[ n ], 
            center_table, peel_table, array_table, antenna_table, height = height, 
            iterations = iterations )
      
        # loop over pierce points
        j = 0
        for pierce_info in pierce_table:
          [ X, pza, [ k, i ] ] = pierce_info
          if ( ( error_array[ n ][ k ][ i ] > 0. ) and
              use_sources[ k ] and use_antennas[ i ] ):
          
            # store model fit input data
            Xp_table.append( X )
            pza_table.append( pza )
            if ( l == 0 ):
              phase_table.append( phase_array[ n ][ k ][ i ] )
              error_table.append( error_array[ n ][ k ][ i ] )
              if ( i == reference_array[ n ][ k ] ):
                ref_list.append( [ k, i, j ] )
              j = j + 1
              info_table.append( [ k, i, reference_array[ n ][ k ] ] )
          
            # temporarily store pierce point info
            temp_pierce_table.append( [ time_table[ n ][ 0 ], l + 1, k + 1, i + 1, X, pza ] )
      
        if ( l == 0 ):
          # loop over pierce points
          for pierce_info in pierce_table:
            [ X, pza, [ k, i ] ] = pierce_info
            if ( ( error_array[ n ][ k ][ i ] > 0. ) and
                use_sources[ k ] and use_antennas[ i ] ):
          
              # store index to reference antenna
              j = ( [ ( k == ref[ 0 ] ) for ref in ref_list ] ).index( True )
              ref_table.append( ref_list[ j ][ 2 ] )
      
        # skip times of poor data
        if ( len( Xp_table ) < order ):
          if print_info:
            print '...... skipping time step due to no / too few data points'
          fit_table.append( [ time_table[ n ][ 0 ], 0, azeros( P ).tolist(), 0. ] )
          for p in temp_pierce_table:
            pierce_table.append( p )
          if estimate_offsets:
            solution_table.append( solution_row )
          break
      
        # convert tables to arrays for fitting routines
        Xp_table = array( Xp_table, dtype = dtype )
        pza_table = array( pza_table, dtype = dtype )
        if ( not include_airmass ):
          pza_table = azeros( pza_table )
        if ( l == 0 ):
          ref_table = array( ref_table, dtype = int32 )
          phase_table = array( phase_table, dtype = dtype )
          error_table = array( error_table, dtype = dtype )
          if equal_weights:
            error_table = ones( shape( error_table ), dtype = dtype )
          info_table = array( info_table, dtype = int32 )
      
        # calculate structure matrix
        p_count = len( Xp_table )
        D_table = resize( Xp_table, ( p_count, p_count, 2 ) )
        D_table = transpose( D_table, ( 1, 0, 2 ) ) - D_table
        D_table = add.reduce( D_table**2, 2 )
        D_table = ( D_table / ( r_0**2 ) )**( beta / 2. )
      
        # calculate covariance matrix C
        # calculate partial product for interpolation B
        C_table = - D_table / 2.
        C_table = transpose( C_table - ( add.reduce( C_table, 0 ) / float( p_count ) ) )
        B_table = add.reduce( C_table, 0 ) / float( p_count )
        C_table = C_table - B_table
        C_table = ( C_table + transpose( C_table ) ) / 2.

        # incorporate airmass functions and layer weights
        # add layer C to total C
        A_table = resize( 1. / cos( aradians( pza_table ) ), ( p_count, p_count ) )
        A_table = A_table * transpose( A_table )
        if ( l == 0 ):
          Cl_table = C_table * A_table * ( layer_weights[ l ]**2 )
        else:
          Cl_table = Cl_table + C_table * A_table * ( layer_weights[ l ]**2 )

        # save tables per height layer
        Xpl_table.append( Xp_table )
        pzal_table.append( pza_table )
        Bl_table.append( B_table )

      # skip times of poor data
      if ( len( Xp_table ) < order ):
        continue

      # convert to arrays
      Xpl_table = array( Xpl_table, dtype = dtype )
      pzal_table = array( pzal_table, dtype = dtype )
      Bl_table = array( Bl_table, dtype = dtype )

      # eigenvalue decomposition
      # reforce symmetry
      # select subset of base vectors
      try:
        [ U_table, S, Ut_table ] = linalg.svd( Cl_table )
      except:
        if print_info:
          print '... SVD did not converge, skipping time interval'
        fit_table.append( [ time_table[ n ][ 0 ], 0, azeros( P ).tolist(), 0. ] )
        for p in temp_pierce_table:
          pierce_table.append( p )
        if estimate_offsets:
          solution_table.append( solution_row )
        continue
      U_table = ( U_table + transpose( Ut_table ) ) / 2.
      U_table = U_table[ : , 0 : order ]
      S = S[ 0 : order ]

      # see if we can use the gradient of the previous good fit as a starting point
      dtime = abs( time_table[ n ][ 0 ] - prop_last_time ) * 24. * 60.
      if ( dtime < prop_time_limit ):
        G = prop_last_G.copy()

      else:
        # estimate initial phase gradient on layer with heighest weight
        l = layer_weights.index( max( layer_weights ) )
        Xp_table = Xpl_table[ l ] / cos( aradians( pzal_table[ l ] ) ).reshape( 
            ( p_count, 1 ) )
        G = estimate_phase_gradient_through_fft( uv, Xp_table, ref_table, phase_table,
            error_table, fft_cell_size = 200., fft_image_size = 1024 )
#        if print_info:
#          print '... estimated gradient is %s' % repr( G.tolist() )
      
      # gradient fit to improve accuracy
      function_keywords = { 'phi_model' : phi_gradient_model, 
          'X_table' : Xpl_table[ l ], 'pza_table' : pzal_table[ l ], 
          'ref_table' : ref_table, 'phase_table' : phase_table,
          'error_table' : error_table, 'normalize_weights' : normalize_weights,
          'wrap_phases' : True }
      parameter_info = []
      for m in range( 2 ):
        par_info = { 'parname' : 'P_%d' % ( m ), 'value' : G[ m ], 
            'limits' : [ None, None ] }
        parameter_info.append( par_info )
      fit = mpfit( fit_phi_model, functkw = function_keywords, parinfo = parameter_info,
          quiet = True, autoderivative = True, debug = False, fastnorm = False, 
          nocovar = True, dblprec = double_precision )
      G = fit.params.copy()
      chi_G = sqrt( fit.fnorm )
      if print_info:
#        print '... fitted gradient is %s' % repr( G.tolist() )
        print '... gradient post-fit phase RMS is %s degrees' % repr( chi_G )
      
      # project gradient onto base vectors
      # distribute gradient over layers
      phi_table = azeros( phase_table )
      sel = ref_table.reshape( p_count, 1 )
      Xp_rel_table = Xp_table - aget( Xp_table, sel )
      Xp_rel_table = aput( Xp_rel_table, sel, aones( Xp_rel_table ) )
      for l in range( layer_count ):
        X_table = Xpl_table[ l ] / cos( aradians( pzal_table[ l ] ) ).reshape( 
            ( p_count, 1 ) )
        X_rel_table = X_table - aget( X_table, sel )
        ratio_table = sqrt( add.reduce( X_rel_table**2, 1 ) / 
            add.reduce( Xp_rel_table**2, 1 ) )
        ratio_table = aput( ratio_table, sel, azeros( ratio_table ) )
        sel2 = awhere( ratio_table != 0. )
        ratio = aget( ratio_table, sel2 ).mean()
        phi_table = phi_table + layer_weights[ l ] * dot( X_table, G / ratio )
      P0 = dot( transpose( U_table ), phi_table )

      if remove_gradient:
        PG = P0.copy()
        P0 = azeros( P0 )
        phi_G_table = dot( U_table, PG )
        phi_G_table = phi_G_table - aget( phi_G_table, sel )
        phase_table = amodulo( ( phase_table - phi_G_table ) + 180., 360. ) - 180.
      
      # define data passed to the model fit routines
      # define model parameter structure
      function_keywords = { 'U_table' : U_table, 'ref_table' : ref_table, 
          'phase_table' : phase_table, 'error_table' : error_table, 
          'normalize_weights' : normalize_weights }
    
      # determine lowest chi^2 by evolutionary steps
      P = P0.copy()
      chi2 = 1.e10
#      for e_step in range( 1, 1 + e_steps ):
      e_step = e_steps
      
      # each time start with same initial guess
      Pe = P0.copy()
      
      # determine evolution of fit
      P_fixed_table = ones( ( e_step, order ), dtype = bool )
      Sa = abs( S )
      Sa_max = max( Sa )
      Sa_min = min( Sa )
      Sa_max_sel = awhere( Sa == Sa_max )
      Sa_min_sel = awhere( Sa == Sa_min )
      P_fixed_table[ 0 ] = aput( P_fixed_table[ 0 ], Sa_max_sel, False )
      P_fixed_table[ e_step - 1 ] = aput( P_fixed_table[ e_step - 1 ], Sa_min_sel,
          False )
      e_delta = ( log( Sa_min ) - log( Sa_max ) ) / float( e_step )
      e_limits = exp( arange( e_step + 1, dtype = dtype ) * e_delta + log( Sa_max ) )
      for e in range( e_step ):
        for m in range( order ):
          if ( ( Sa[ m ] <=  e_limits[ e ] ) and ( Sa[ m ] > e_limits[ e + 1 ] ) ): 
            P_fixed_table[ e, m ] = False

      # evolutionary fit
      for e in range( e_step ):
        
        if alltrue( P_fixed_table[ e ] ):
          continue
        
        # specify fit parameters
        parameter_info = []
        for m in range( order ):
          par_info = { 'parname' : 'P_%d' % ( m + 1 ), 'value' : Pe[ m ], 
              'limits' : [ None, None ] }
          if P_fixed_table[ e, m ]:
            par_info[ 'limits' ] = [ Pe[ m ], Pe[ m ] ]
          parameter_info.append( par_info )
        
        # fit model to data points
        fit = mpfit( fit_phi_mkl_model, functkw = function_keywords, 
            parinfo = parameter_info, quiet = True, autoderivative = True, 
            debug = False, fastnorm = False, nocovar = True,
            dblprec = double_precision )
#        if print_info:
#          print '...... mpfit() returned status %d' % ( fit.status ) 
        if ( fit.errmsg != '' ):
          raise error( fit.errmsg )
        Pe = fit.params.copy()
      
        # store fit results
        if ( fit.fnorm < chi2 ):
#          best_e = e_step
          P = Pe.copy()
          chi = sqrt( fit.fnorm )
      
      if print_info:
#        print '... fitted parameters are %s' % repr( P.tolist() )
#        print '... used e-value is %s' % repr( best_e )
        print '... model post-fit phase RMS is %s degrees' % repr( chi )
      
      # check that model fit is better than gradient fit
      if ( chi > chi_G ):
        if print_info:
          print '... WARNING: model fit is worse than gradient fit'
          print '...... using gradient fit and allowing repeated fit'
        P = azeros( P )
        rejected_sources = False
        rejected_antennas = False
        prop_last_time = -1.e9

      # repeat bad fit if previous gradient was used
      if ( ( dtime < prop_time_limit ) and ( chi > prop_rms_limit ) ):
        if print_info:
          print '...... repeating fit using fresh gradient estimate'
        rejected_sources = False
        rejected_antennas = False
        repeat_fit = True
        prop_last_time = -1.e9
        continue
      
      # add removed gradient back
      if remove_gradient:
        P = P + PG
        phase_table = amodulo( ( phase_table + phi_G_table ) + 180., 360. ) - 180.

      # calculate model phases
      phi_table = dot( U_table, P )
      
      # extract gradient
      if extract_gradient:
        F_table = dot( U_table, P / S )
        phi_offset = phi_mkl_model( layer_weights,
            azeros( Xpl_table[ : , 0 : 1, : ] ), aones( pzal_table[ : , 0 : 1 ] ),
            Xpl_table, pzal_table, Bl_table, F_table, beta = beta, r_0 = r_0 )
        T1 = linalg.inv( dot( transpose( Xp_table ), Xp_table ) )
        T2 = dot( transpose( Xp_table ), dot( U_table, P ) - phi_offset )
        G = dot( T1, T2 )

      # cross check (phi_ip should match phi)
      # recalculate model phases at pierce points
#      phi_ip_table = phi_mkl_model( layer_weights, Xpl_table, pzal_table, Xpl_table, 
#          pzal_table, Bl_table, F_table, beta = beta, r_0 = r_0 )
    
      # split pierce points per source
      ref_list = []
      for ref in ref_table:
        if ( not ( ref in ref_list ) ):
          ref_list.append( ref )
    
      # estimate residual phase per antenna per source
      antenna_error_table = azeros( error_table )
      antenna_offset_table = azeros( error_table )
      reference_antennas = []
      for ref in ref_list:
        reference_antennas.append( info_table[ ref ][ 1 ] )
        sel = awhere( ref_table == ref )
        dphase_table = aget( phase_table, sel )
        dphase_table = resize( dphase_table, ( len( sel ), len( sel ) ) )
        dphase_table = transpose( dphase_table ) - dphase_table
        dphi_table = aget( phi_table, sel )
        dphi_table = resize( dphi_table, ( len( sel ), len( sel ) ) )
        dphi_table = transpose( dphi_table ) - dphi_table
        dweight_table = aget( error_table, sel )
        dweight_table = resize( dweight_table**2, ( len( sel ), len( sel ) ) )
        dweight_table = transpose( dweight_table ) + dweight_table
        dweight_table = 1. / sqrt( dweight_table )
        abs_err_table = add.reduce( dweight_table * 
            abs( amodulo( ( dphase_table - dphi_table ) + 180., 360. ) - 180. ), 1 )
        normalize_table = add.reduce( dweight_table, 1 ) - diagonal( dweight_table )
        abs_err_table = abs_err_table / normalize_table
        antenna_error_table = aput( antenna_error_table, sel, abs_err_table )
        offset_table = add.reduce( dweight_table * 
            ( amodulo( ( dphase_table - dphi_table ) + 180., 360. ) - 180. ), 1 )
        offset_table = offset_table / normalize_table
        antenna_offset_table = aput( antenna_offset_table, sel, offset_table )
    
      # estimate residual phase RMS per source
      source_rmss = zeros( ( source_count ), dtype = dtype )
      for k in range( source_count ):
        # select all phases related to single antenna
        sel = awhere( info_table[ : , 0 ] == k )
        # TODO: determine reasonable number of values for calculating the mean
        if ( len( sel ) < 3 ):
          continue
        # calculate mean phase offset per source
        phase_offsets = aget( antenna_offset_table, sel )
        rms_phase_offset = sqrt( ( ( amodulo( phase_offsets
            + 180., 360. ) - 180. )**2 ).mean() )
        source_rmss[ k ] = rms_phase_offset

      # check for source with excessive RMS
      if ( not rejected_sources ):
        sel = awhere( source_rmss > source_rms_limit )
        for s in sel:
          s = s[ 0 ]
          use_sources[ s ] = False
          rejected_sources = True
          if print_info:
            print '...... source %s has excessive RMS of %s degrees' % (
                repr( s + 1 ), repr( source_rmss[ s ] ) )
        if rejected_sources:
          if print_info:
            print '...... repeating fit'
          repeat_fit = True
          if ( chi > prop_rms_limit ):
            prop_last_time = -1.e9
          continue
      
      # estimate residual phase offset and stddev per antenna
      antenna_means = zeros( ( antenna_count ), dtype = dtype )
      antenna_stds = zeros( ( antenna_count ), dtype = dtype )
      for i in range( antenna_count ):
        # select all phases related to single antenna
        sel = awhere( info_table[ : , 1 ] == i )
        # TODO: determine reasonable number of values for calculating the antenna mean
        if ( len( sel ) < 3 ):
          continue
        # calculate mean phase offset per antenna
        phase_offsets = aget( antenna_offset_table, sel )
        mean_phase_offset = amean_phase( phase_offsets )
        std_phase_offset = sqrt( ( ( amodulo( ( phase_offsets - mean_phase_offset )
            + 180., 360. ) - 180. )**2 ).mean() )
        antenna_means[ i ] = mean_phase_offset
        antenna_stds[ i ] = std_phase_offset
        if estimate_offsets:
          # store phase error in solution table
          # note this automatically maps mean_phase in [ -180, 180 ] domain
          solution = r_phi_to_complex( [ 1., mean_phase_offset ] )
          solution_row[ i + 1 ] = [ solution.real, solution.imag, 0.,
              1. / std_phase_offset ]
      
      # check for antennas with excessive offset or stddev
      if ( not rejected_antennas ):
        sel1 = awhere( antenna_means > antenna_offset_limit )
        for s in sel1:
          s = s[ 0 ]
          if print_info:
            print '...... antenna %s has excessive offset of %s degrees' % (
                repr( s + 1 ), repr( antenna_means[ s ] ) )
          if ( s in reference_antennas ):
            if print_info:
              print '...... WARNING: reference antenna not dropped'
          else:
            use_antennas[ s ] = False
            rejected_antennas = True
        sel2 = awhere( antenna_stds > antenna_rms_limit )
        for s in sel2:
          s = s[ 0 ]
          if use_antennas[ s ]:
            if print_info:
              print '...... antenna %s has excessive stddev of %s degrees' % (
                  repr( s + 1 ), repr( antenna_stds[ s ] ) )
            if ( s in reference_antennas ):
              if print_info:
                print '...... WARNING: reference antenna not dropped'
            else:
              use_antennas[ s ] = False
              rejected_antennas = True
        if rejected_antennas:
          if print_info:
            print '...... repeating fit'
          repeat_fit = True
          if ( chi > prop_rms_limit ):
            prop_last_time = -1.e9
          continue

      # store fit results
      weight = 1. / chi
      fit_table.append( [ time_table[ n ][ 0 ], p_count, P.tolist(), weight ] )
      for p in temp_pierce_table:
        pierce_table.append( p )

      if estimate_offsets:
        # add solution row to solution table
        solution_table.append( solution_row )

      # see if we can use the gradient of the previous good fit as a starting point
      if ( chi < prop_rms_limit ):
        prop_last_G = G.copy()
        prop_last_time = time_table[ n ][ 0 ]

  # write fit results to model fit table
  write_ionospheric_mkl_model_fit_table( uv, layer_heights, layer_weights, fit_table,
      beta = beta, r_0 = r_0 )
  fit_version = uv.table_highver( 'NI' )
  write_ionospheric_mkl_pierce_table( uv, pierce_table, version = fit_version,
      iterations = iterations )
  
  # write solution table with antenna phase offsets to uv data
  if estimate_offsets:
    if print_info:
      print '... NOTE: writing solution table with antenna phase offsets to UV data'
    write_solution_table( uv, solution_table, out_version = 0 )
  
  return

###############################################################################

def fit_ionospheric_mkl_model_new( uv, facets, facet_list = [], beta = 5. / 3., 
    r_0 = 1., order = 15, layer_heights = [ 200.e3 ], layer_weights = [ 1. ], 
    iterations = 4, equal_weights = True, normalize_weights = True, 
    e_steps = 1, double_precision = True, print_info = True,
    extract_gradient = False, remove_polynomial = False,
    estimate_offsets = False, time_steps = [], include_airmass = True, 
    solution_version = 0, source_rms_limit = 1.e6,
    antenna_offset_limit = 1.e6, antenna_rms_limit = 1.e6,
    prop_time_limit = 1., prop_rms_limit = 0. ): # remove_gradient = False, 
# prop_time_limit in minutes
  
  if double_precision:
    dtype = float64
  else:
    dtype = float32
  
  # select appropriate model functions
  # initialize model parameters
  model = 'mkl'
  layer_count = len( layer_heights )
  P = zeros( shape = ( order ), dtype = dtype )
  if ( len( layer_weights ) != layer_count ):
    raise error( 'layer heights and weights are different lengths' )
  
  if print_info:
    print 'starting model fits using model = ' + model
  
  # get all relevant data
  calibration_data = get_phase_calibration_data( uv, facets, facet_list = facet_list,
      solution_version = solution_version, print_info = print_info )
  time_table = calibration_data[ 0 ]
  center_table = calibration_data[ 1 ]
  peel_table = calibration_data[ 2 ]
  array_table = calibration_data[ 3 ]
  antenna_table = calibration_data[ 4 ]
  reference_array = calibration_data[ 5 ]
  phase_array = calibration_data[ 6 ]
  error_array = calibration_data[ 7 ]
  source_count = len( peel_table )
  antenna_count = len( antenna_table )
  
  fit_table = []
  pierce_table = []
  
  if estimate_offsets:
    solution_table = []

  prop_last_time = -1.e9
  prop_last_G = zeros( ( 2 ), dtype = dtype )
  
  # loop over time stamps
  if ( len( time_steps ) == 0 ):
    n_list = range( len( time_table ) )
  else:
    n_list = time_steps
  for n in n_list:

    if print_info:
      print '... time step n = ' + repr( n )
    
    repeat_fit = True
    rejected_sources = False
    rejected_antennas = False
    while repeat_fit:
      repeat_fit = False
      if ( not rejected_sources ):
        use_sources = [ True for k in range( source_count ) ]
      if ( not rejected_antennas ):
        use_antennas = [ True for i in range( antenna_count ) ]
      
      temp_pierce_table = []

      # containers for model input data
      ref_list = []
      ref_table = []
      phase_table = []
      error_table = []
      info_table = []
      Xpl_table = []
      pzal_table = []
      Bl_table = []

      if estimate_offsets:
        solution_row = [ [ time_table[ n ][ 0 ], 
            reference_array[ n ][ 0 ] + 1, 0., 0. ] ] + [ [ 0., 0., 0., 0. ] 
            for i in range( antenna_count ) ]

      # loop over height layers
      l_max = layer_weights.index( max( layer_weights ) )
      for l in range( layer_count ):
        height = layer_heights[ l ]
      
        # containers for model input data
        Xp_table = []
        pza_table = []
      
        # get pierce point coordinates
        pierce_table = calculate_pierce_coordinates( time_table[ n ], 
            center_table, peel_table, array_table, antenna_table, height = height, 
            iterations = iterations )
      
        # loop over pierce points
        j = 0
        for pierce_info in pierce_table:
          [ X, pza, [ k, i ] ] = pierce_info
          if ( ( error_array[ n ][ k ][ i ] > 0. ) and
              use_sources[ k ] and use_antennas[ i ] ):
          
            # store model fit input data
            Xp_table.append( X )
            pza_table.append( pza )
            if ( l == 0 ):
              phase_table.append( phase_array[ n ][ k ][ i ] )
              error_table.append( error_array[ n ][ k ][ i ] )
              if ( i == reference_array[ n ][ k ] ):
                ref_list.append( [ k, i, j ] )
              j = j + 1
              info_table.append( [ k, i, reference_array[ n ][ k ] ] )
          
            # temporarily store pierce point info
            temp_pierce_table.append( [ time_table[ n ][ 0 ], l + 1, k + 1, i + 1, X, pza ] )
      
        if ( l == 0 ):
          # loop over pierce points
          for pierce_info in pierce_table:
            [ X, pza, [ k, i ] ] = pierce_info
            if ( ( error_array[ n ][ k ][ i ] > 0. ) and
                use_sources[ k ] and use_antennas[ i ] ):
          
              # store index to reference antenna
              j = ( [ ( k == ref[ 0 ] ) for ref in ref_list ] ).index( True )
              ref_table.append( ref_list[ j ][ 2 ] )
      
        # skip times of poor data
        skip = False
        if ( len( Xp_table ) <= order + 2 ):
          if print_info:
            print '...... skipping time step due to no / too few data points'
          fit_table.append( [ time_table[ n ][ 0 ], 0, azeros( P ).tolist(), 0. ] )
          for p in temp_pierce_table:
            pierce_table.append( p )
          if estimate_offsets:
            solution_table.append( solution_row )
          skip = True
          break
      
        # convert tables to arrays for fitting routines
        Xp_table = array( Xp_table, dtype = dtype )
        pza_table = array( pza_table, dtype = dtype )
        if ( not include_airmass ):
          pza_table = azeros( pza_table )
        if ( l == 0 ):
          ref_table = array( ref_table, dtype = int32 )
          phase_table = array( phase_table, dtype = dtype )
          error_table = array( error_table, dtype = dtype )
          if equal_weights:
            error_table = ones( shape( error_table ), dtype = dtype )
          info_table = array( info_table, dtype = int32 )
      
        # calculate structure matrix
        p_count = len( Xp_table )
        D_table = resize( Xp_table, ( p_count, p_count, 2 ) )
        D_table = transpose( D_table, ( 1, 0, 2 ) ) - D_table
        D_table = add.reduce( D_table**2, 2 )
        D_table = ( D_table / ( r_0**2 ) )**( beta / 2. )
      
        # calculate covariance matrix C
        # calculate partial product for interpolation B
        C_table = - D_table / 2.
        C_table = transpose( C_table - ( add.reduce( C_table, 0 ) / float( p_count ) ) )
        B_table = add.reduce( C_table, 0 ) / float( p_count )
        C_table = C_table - B_table
        C_table = ( C_table + transpose( C_table ) ) / 2.

        # incorporate airmass functions and layer weights
        # add layer C to total C
        A_table = resize( 1. / cos( aradians( pza_table ) ), ( p_count, p_count ) )
        A_table = A_table * transpose( A_table )
        if ( l == 0 ):
          Cl_table = C_table * A_table * ( layer_weights[ l ]**2 )
        else:
          Cl_table = Cl_table + C_table * A_table * ( layer_weights[ l ]**2 )
        if ( l == l_max ):
          R_table = dot( C_table, diag( 1. / cos( aradians( pza_table ) ) ) )

        # save tables per height layer
        Xpl_table.append( Xp_table )
        pzal_table.append( pza_table )
        Bl_table.append( B_table )

      # skip times of poor data
      if skip:
        continue

      # convert to arrays
      Xpl_table = array( Xpl_table, dtype = dtype )
      pzal_table = array( pzal_table, dtype = dtype )
      Bl_table = array( Bl_table, dtype = dtype )

      # eigenvalue decomposition
      # reforce symmetry
      # select subset of base vectors
      try:
        [ U_table, S, Ut_table ] = linalg.svd( Cl_table )
      except:
        if print_info:
          print '... SVD did not converge, skipping time interval'
        fit_table.append( [ time_table[ n ][ 0 ], 0, azeros( P ).tolist(), 0. ] )
        for p in temp_pierce_table:
          pierce_table.append( p )
        if estimate_offsets:
          solution_table.append( solution_row )
        continue
      U_table = ( U_table + transpose( Ut_table ) ) / 2.
      U_table = U_table[ : , 0 : order ]
      S = S[ 0 : order ]
      # calculate projection matrix for dominant layer
      R_table = dot( R_table, dot( U_table, diag( 1. / S ) ) )
      R_table = dot( linalg.inv( dot( transpose( R_table ), R_table ) ),
          transpose( R_table ) )
      
      # see if we can use the gradient of the previous good fit as a starting point
      dtime = abs( time_table[ n ][ 0 ] - prop_last_time ) * 24. * 60.
      Xp_table = Xpl_table[ l_max ] / cos( aradians( pzal_table[ l_max ] ) ).reshape( 
          ( p_count, 1 ) )
      if ( dtime < prop_time_limit ):
        G = prop_last_G.copy()
      else:
        # estimate initial phase gradient on layer with heighest weight
        G = estimate_phase_gradient_through_fft( uv, Xp_table, ref_table, phase_table,
            error_table, fft_cell_size = 200., fft_image_size = 1024 )
#        if print_info:
#          print '... estimated gradient is %s' % repr( G.tolist() )
      
#      if remove_gradient:
#        save_G = G.copy()
#        phi_G_table = dot( Xp_table, G )
#        sel = ref_table.reshape( p_count, 1 )
#        phi_G_table = phi_G_table - aget( phi_G_table, sel )
#        phase_table = amodulo( ( phase_table - phi_G_table ) + 180., 360. ) - 180.
#        G = azeros( G )
      
      # 2-step polynomial fit to improve accuracy
      function_keywords = { 'phi_model' : phi_poly_model, 
          'X_table' : Xpl_table[ l_max ], 'pza_table' : pzal_table[ l_max ], 
          'ref_table' : ref_table, 'phase_table' : phase_table,
          'error_table' : error_table, 'normalize_weights' : normalize_weights,
          'wrap_phases' : True }
      parameter_info = []
      for m in range( 5 ):
        if ( m < 2 ):
          par_info = { 'parname' : 'P_%d' % ( m ), 'value' : G[ m ], 
              'limits' : [ None, None ] }
        else:
          par_info = { 'parname' : 'P_%d' % ( m ), 'value' : 0., 
              'limits' : [ 0., 0. ] }
        parameter_info.append( par_info )
      fit = mpfit( fit_phi_model, functkw = function_keywords, parinfo = parameter_info,
          quiet = True, autoderivative = True, debug = False, fastnorm = False, 
          nocovar = True, dblprec = double_precision )
      pol = fit.params.copy()
      G = pol[ 0 : 2 ]
      chi_grad = sqrt( fit.fnorm )
      if print_info:
        print '... gradient post-fit phase RMS is %s degrees' % repr( chi_grad )
      
      function_keywords = { 'phi_model' : phi_poly_model, 
          'X_table' : Xpl_table[ l_max ], 'pza_table' : pzal_table[ l_max ], 
          'ref_table' : ref_table, 'phase_table' : phase_table,
          'error_table' : error_table, 'normalize_weights' : normalize_weights,
          'wrap_phases' : True }
      parameter_info = []
      for m in range( 5 ):
        if ( m < 2 ):
          par_info = { 'parname' : 'P_%d' % ( m ), 'value' : pol[ m ], 
              'limits' : [ pol[ m ], pol[ m ] ] }
        else:
          par_info = { 'parname' : 'P_%d' % ( m ), 'value' : pol[ m ], 
              'limits' : [ None, None ] }
        parameter_info.append( par_info )
      fit = mpfit( fit_phi_model, functkw = function_keywords, parinfo = parameter_info,
          quiet = True, autoderivative = True, debug = False, fastnorm = False, 
          nocovar = True, dblprec = double_precision )
      pol = fit.params.copy()
      chi_pol = sqrt( fit.fnorm )
      if print_info:
        print '... polynomial post-fit phase RMS is %s degrees' % repr( chi_pol )
      if ( chi_grad < chi_pol ):
        pol[ 2 : 5 ] = [ 0., 0., 0. ]
        chi_pol = chi_grad
      
      # project polynomial onto base vectors
#      if remove_gradient:
#        G = save_G.copy()
#        pol[ 0 : 2 ] = pol[ 0 : 2 ] + G
#        phase_table = amodulo( ( phase_table + phi_G_table ) + 180., 360. ) - 180.
#      phi_table = phi_poly_model( Xp_table, pol )
#      P0 = dot( transpose( U_table ), phi_table )
      phi_table = phi_poly_model( Xpl_table[ l_max ], pol ) / cos( 
          aradians( pzal_table[ l_max ] ) )
      P0 = dot( R_table, phi_table )
      sel = ref_table.reshape( p_count, 1 )
      phi_table = phi_table - aget( phi_table, sel )

#      # define data passed to the model fit routines
#      # define model parameter structure
#      function_keywords = { 'U_table' : U_table, 'ref_table' : ref_table, 
#          'phase_table' : phi_table, 'error_table' : error_table, 
#          'normalize_weights' : normalize_weights }
#      # specify fit parameters
#      parameter_info = []
#      for m in range( order ):
#        par_info = { 'parname' : 'P_%d' % ( m + 1 ), 'value' : P0[ m ], 
#            'limits' : [ None, None ] }
#        parameter_info.append( par_info )
#      # fit model to data points
#      fit = mpfit( fit_phi_mkl_model, functkw = function_keywords, 
#          parinfo = parameter_info, quiet = True, autoderivative = True, 
#          debug = False, fastnorm = False, nocovar = True,
#          dblprec = double_precision )
#      if ( fit.errmsg != '' ):
#        raise error( fit.errmsg )
#      P0 = fit.params.copy()
#      chi_pol_fit = sqrt( fit.fnorm )
#      if print_info:
#        print '... polynomial to model post-fit phase RMS is %s degrees' % repr( chi_pol_fit )
      
      if remove_polynomial:
        P_pol = P0.copy()
        P0 = azeros( P0 )
#        phi_pol_table = dot( U_table, P_pol )
#        sel = ref_table.reshape( p_count, 1 )
#        phi_pol_table = phi_pol_table - aget( phi_pol_table, sel )
        phi_pol_table = phi_table
        phase_table = amodulo( ( phase_table - phi_pol_table ) + 180., 360. ) - 180.
      
      # define data passed to the model fit routines
      # define model parameter structure
      function_keywords = { 'U_table' : U_table, 'ref_table' : ref_table, 
          'phase_table' : phase_table, 'error_table' : error_table, 
          'normalize_weights' : normalize_weights }
    
      # determine lowest chi^2 by evolutionary steps
      P = P0.copy()
      chi2 = 1.e10
#      for e_step in range( 1, 1 + e_steps ):
      e_step = e_steps
      
      # each time start with same initial guess
      Pe = P0.copy()
      
      # determine evolution of fit
      P_fixed_table = ones( ( e_step, order ), dtype = bool )
      Sa = abs( S )
      Sa_max = max( Sa )
      Sa_min = min( Sa )
      Sa_max_sel = awhere( Sa == Sa_max )
      Sa_min_sel = awhere( Sa == Sa_min )
      P_fixed_table[ 0 ] = aput( P_fixed_table[ 0 ], Sa_max_sel, False )
      P_fixed_table[ e_step - 1 ] = aput( P_fixed_table[ e_step - 1 ], Sa_min_sel,
          False )
      e_delta = ( log( Sa_min ) - log( Sa_max ) ) / float( e_step )
      e_limits = exp( arange( e_step + 1, dtype = dtype ) * e_delta + log( Sa_max ) )
      for e in range( e_step ):
        for m in range( order ):
          if ( ( Sa[ m ] <=  e_limits[ e ] ) and ( Sa[ m ] > e_limits[ e + 1 ] ) ): 
            P_fixed_table[ e, m ] = False

      # evolutionary fit
      for e in range( e_step ):
        
        if alltrue( P_fixed_table[ e ] ):
          continue
        
        # specify fit parameters
        parameter_info = []
        for m in range( order ):
          par_info = { 'parname' : 'P_%d' % ( m + 1 ), 'value' : Pe[ m ], 
              'limits' : [ None, None ] }
          if P_fixed_table[ e, m ]:
            par_info[ 'limits' ] = [ Pe[ m ], Pe[ m ] ]
          parameter_info.append( par_info )
        
        # fit model to data points
        fit = mpfit( fit_phi_mkl_model, functkw = function_keywords, 
            parinfo = parameter_info, quiet = True, autoderivative = True, 
            debug = False, fastnorm = False, nocovar = True,
            dblprec = double_precision )
#        if print_info:
#          print '...... mpfit() returned status %d' % ( fit.status ) 
        if ( fit.errmsg != '' ):
          raise error( fit.errmsg )
        Pe = fit.params.copy()
      
        # store fit results
        if ( fit.fnorm < chi2 ):
#          best_e = e_step
          P = Pe.copy()
          chi = sqrt( fit.fnorm )
      
      if print_info:
#        print '... fitted parameters are %s' % repr( P.tolist() )
#        print '... used e-value is %s' % repr( best_e )
        print '... model post-fit phase RMS is %s degrees' % repr( chi )
      
      # check that model fit is better than polynomial fit
      if ( chi > chi_pol ):
        if print_info:
          print '... WARNING: model fit is worse than polynomial fit'
#          print '...... using gradient fit and allowing repeated fit'
#        P = P0.copy()
#        rejected_sources = False
#        rejected_antennas = False
#        prop_last_time = -1.e9

      # repeat bad fit if previous gradient was used
      if ( ( dtime < prop_time_limit ) and ( chi > prop_rms_limit ) ):
        if print_info:
          print '...... repeating fit using fresh gradient estimate'
        rejected_sources = False
        rejected_antennas = False
        repeat_fit = True
        prop_last_time = -1.e9
        continue
      
      # add removed gradient back
      if remove_polynomial:
        P = P + P_pol
        phase_table = amodulo( ( phase_table + phi_pol_table ) + 180., 360. ) - 180.

      # calculate model phases
      phi_table = dot( U_table, P )
      
      # extract gradient
      if extract_gradient:
        F_table = dot( U_table, P / S )
        phi_offset = phi_mkl_model( layer_weights,
            azeros( Xpl_table[ : , 0 : 1, : ] ), aones( pzal_table[ : , 0 : 1 ] ),
            Xpl_table, pzal_table, Bl_table, F_table, beta = beta, r_0 = r_0 )
        T1 = linalg.inv( dot( transpose( Xp_table ), Xp_table ) )
        T2 = dot( transpose( Xp_table ), dot( U_table, P ) - phi_offset )
        G = dot( T1, T2 )

      # cross check (phi_ip should match phi)
      # recalculate model phases at pierce points
#      phi_ip_table = phi_mkl_model( layer_weights, Xpl_table, pzal_table, Xpl_table, 
#          pzal_table, Bl_table, F_table, beta = beta, r_0 = r_0 )
    
      # split pierce points per source
      ref_list = []
      for ref in ref_table:
        if ( not ( ref in ref_list ) ):
          ref_list.append( ref )
    
      # estimate residual phase per antenna per source
      antenna_error_table = azeros( error_table )
      antenna_offset_table = azeros( error_table )
      reference_antennas = []
      for ref in ref_list:
        reference_antennas.append( info_table[ ref ][ 1 ] )
        sel = awhere( ref_table == ref )
        dphase_table = aget( phase_table, sel )
        dphase_table = resize( dphase_table, ( len( sel ), len( sel ) ) )
        dphase_table = transpose( dphase_table ) - dphase_table
        dphi_table = aget( phi_table, sel )
        dphi_table = resize( dphi_table, ( len( sel ), len( sel ) ) )
        dphi_table = transpose( dphi_table ) - dphi_table
        dweight_table = aget( error_table, sel )
        dweight_table = resize( dweight_table**2, ( len( sel ), len( sel ) ) )
        dweight_table = transpose( dweight_table ) + dweight_table
        dweight_table = 1. / sqrt( dweight_table )
        abs_err_table = add.reduce( dweight_table * 
            abs( amodulo( ( dphase_table - dphi_table ) + 180., 360. ) - 180. ), 1 )
        normalize_table = add.reduce( dweight_table, 1 ) - diagonal( dweight_table )
        abs_err_table = abs_err_table / normalize_table
        antenna_error_table = aput( antenna_error_table, sel, abs_err_table )
        offset_table = add.reduce( dweight_table * 
            ( amodulo( ( dphase_table - dphi_table ) + 180., 360. ) - 180. ), 1 )
        offset_table = offset_table / normalize_table
        antenna_offset_table = aput( antenna_offset_table, sel, offset_table )
    
      # estimate residual phase RMS per source
      source_rmss = zeros( ( source_count ), dtype = dtype )
      for k in range( source_count ):
        # select all phases related to single antenna
        sel = awhere( info_table[ : , 0 ] == k )
        # TODO: determine reasonable number of values for calculating the mean
        if ( len( sel ) < 3 ):
          continue
        # calculate mean phase offset per source
        phase_offsets = aget( antenna_offset_table, sel )
        rms_phase_offset = sqrt( ( ( amodulo( phase_offsets
            + 180., 360. ) - 180. )**2 ).mean() )
        source_rmss[ k ] = rms_phase_offset

      # check for source with excessive RMS
      if ( not rejected_sources ):
        sel = awhere( source_rmss > source_rms_limit )
        for s in sel:
          s = s[ 0 ]
          use_sources[ s ] = False
          rejected_sources = True
          if print_info:
            print '...... source %s has excessive RMS of %s degrees' % (
                repr( s + 1 ), repr( source_rmss[ s ] ) )
        if rejected_sources:
          if print_info:
            print '...... repeating fit'
          repeat_fit = True
          if ( chi > prop_rms_limit ):
            prop_last_time = -1.e9
          continue
      
      # estimate residual phase offset and stddev per antenna
      antenna_means = zeros( ( antenna_count ), dtype = dtype )
      antenna_stds = zeros( ( antenna_count ), dtype = dtype )
      for i in range( antenna_count ):
        # select all phases related to single antenna
        sel = awhere( info_table[ : , 1 ] == i )
        # TODO: determine reasonable number of values for calculating the antenna mean
        if ( len( sel ) < 3 ):
          continue
        # calculate mean phase offset per antenna
        phase_offsets = aget( antenna_offset_table, sel )
        mean_phase_offset = amean_phase( phase_offsets )
        std_phase_offset = sqrt( ( ( amodulo( ( phase_offsets - mean_phase_offset )
            + 180., 360. ) - 180. )**2 ).mean() )
        antenna_means[ i ] = mean_phase_offset
        antenna_stds[ i ] = std_phase_offset
        if estimate_offsets:
          # store phase error in solution table
          # note this automatically maps mean_phase in [ -180, 180 ] domain
          solution = r_phi_to_complex( [ 1., mean_phase_offset ] )
          solution_row[ i + 1 ] = [ solution.real, solution.imag, 0.,
              1. / std_phase_offset ]
      
      # check for antennas with excessive offset or stddev
      if ( not rejected_antennas ):
        sel1 = awhere( antenna_means > antenna_offset_limit )
        for s in sel1:
          s = s[ 0 ]
          if print_info:
            print '...... antenna %s has excessive offset of %s degrees' % (
                repr( s + 1 ), repr( antenna_means[ s ] ) )
          if ( s in reference_antennas ):
            if print_info:
              print '...... WARNING: reference antenna not dropped'
          else:
            use_antennas[ s ] = False
            rejected_antennas = True
        sel2 = awhere( antenna_stds > antenna_rms_limit )
        for s in sel2:
          s = s[ 0 ]
          if use_antennas[ s ]:
            if print_info:
              print '...... antenna %s has excessive stddev of %s degrees' % (
                  repr( s + 1 ), repr( antenna_stds[ s ] ) )
            if ( s in reference_antennas ):
              if print_info:
                print '...... WARNING: reference antenna not dropped'
            else:
              use_antennas[ s ] = False
              rejected_antennas = True
        if rejected_antennas:
          if print_info:
            print '...... repeating fit'
          repeat_fit = True
          if ( chi > prop_rms_limit ):
            prop_last_time = -1.e9
          continue

      # store fit results
      weight = 1. / chi
      fit_table.append( [ time_table[ n ][ 0 ], p_count, P.tolist(), weight ] )
      for p in temp_pierce_table:
        pierce_table.append( p )

      if estimate_offsets:
        # add solution row to solution table
        solution_table.append( solution_row )

      # see if we can use the gradient of the previous good fit as a starting point
      if ( chi < prop_rms_limit ):
        prop_last_G = G.copy()
        prop_last_time = time_table[ n ][ 0 ]

  # write fit results to model fit table
  write_ionospheric_mkl_model_fit_table( uv, layer_heights, layer_weights, fit_table,
      beta = beta, r_0 = r_0 )
  fit_version = uv.table_highver( 'NI' )
  write_ionospheric_mkl_pierce_table( uv, pierce_table, version = fit_version,
      iterations = iterations )
  
  # write solution table with antenna phase offsets to uv data
  if estimate_offsets:
    if print_info:
      print '... NOTE: writing solution table with antenna phase offsets to UV data'
    write_solution_table( uv, solution_table, out_version = 0 )
  
  return

###############################################################################

def write_ionospheric_mkl_model_fit_table( uv, layer_heights, layer_weights,
    fit_table, version = 0, beta = 5. / 3., r_0 = 1., **keywords ):

  # create new NI table
  layer_count = len( layer_heights )
  order = len( fit_table[ 0 ][ 2 ] )
  new_ni_table = new_table( uv, 'NI', version, num_coef = order )
  row = new_table_row( new_ni_table )
  integration_time = restore_parameter( uv, 'integration_time' )
  time_interval = integration_time / ( 24. * 60. * 60. )

  # loop over 
  for fit_row in fit_table:
    [ time, p_count, coefs, weight ] = fit_row
    row.time = float32( time )
    row.time_interval = float32( time_interval )
    row.antenna_no = p_count
    row.source_id = 0
    row.subarray = 0
    row.weight = float32( weight )
    row.coef = [ float32( coef ) for coef in coefs ]
    new_ni_table.append( row )

  # add keywords
  new_ni_table.keywords[ 'MODEL' ] = 'mkl'
  new_ni_table.keywords[ 'REF_FREQ' ] = float32( get_central_frequency( uv ) )
  new_ni_table.keywords[ 'BETA' ] = float32( beta )
  new_ni_table.keywords[ 'R_0' ] = float32( r_0 )
  new_ni_table.keywords[ 'LAYERS' ] = layer_count
  for l in range( layer_count ):
    new_ni_table.keywords[ 'HEIGHT%d' % ( l + 1 ) ] = float32( layer_heights[ l ] )
    new_ni_table.keywords[ 'WEIGHT%d' % ( l + 1 ) ] = float32( layer_weights[ l ] )
  for key in keywords.keys():
    new_ni_table.keywords[ key ] = keywords[ key ]

  #close new table
  new_ni_table.close()

  return

###############################################################################

def write_ionospheric_mkl_pierce_table( uv, pierce_table, version = 0,
    iterations = 4, **keywords ):

  # create new OB table
  new_ob_table = new_table( uv, 'OB', version )
  row = new_table_row( new_ob_table )

  # loop over 
  for pierce_point in pierce_table:
    [ time, layer, source, antenna, X, za ] = pierce_point
    row.time = float32( time )
    row.orientation = float32( layer )
    row.subarray = source
    row.antenna_no = antenna
    row.orbxyz = [ float( x ) for x in X ] + [ float32( za ) ]
    new_ob_table.append( row )

  # add missing keywords and close new table
#  new_ob_table.keywords[ '...' ] = ...
  new_ob_table.keywords[ 'ITER' ] = iterations
  for key in keywords.keys():
    new_ob_table.keywords[ key ] = keywords[ key ]
  new_ob_table.close()

  return

###############################################################################

def generate_solution_tables_from_mkl_fit_table( uv, facets, reference_antenna = 0,
    fit_version = 0, solution_version = 0, include_delays = False, delay_sign = -1.,
    print_info = True, rejection_limit = None, rms_limit = None, facet_list = [],
    time_steps = [], rms_min = 0.1 ):
  
  wiz_uv = wizardry( uv )
  ni_table = wiz_uv.table( 'NI', fit_version )
  model = ni_table.keywords[ 'MODEL' ]
  if ( type( model ) == type( [] ) ):
    model = model[ 0 ]
  model = model.strip()
  if ( model != 'mkl' ):
    raise error( 'unknown model: %s' % ( model ) )
  order = int( ni_table.keywords[ 'NUM_COEF' ] )
  reference_frequency = float32( ni_table.keywords[ 'REF_FREQ' ] )
  beta = float32( ni_table.keywords[ 'BETA' ] )
  r_0 = float32( ni_table.keywords[ 'R_0' ] )
  layer_count = int( ni_table.keywords[ 'LAYERS' ] )
  layer_heights = []
  layer_weights = []
  for l in range( layer_count ):
    layer_heights.append( float32( ni_table.keywords[ 'HEIGHT%d' % ( l + 1 ) ] ) )
    layer_weights.append( float32( ni_table.keywords[ 'WEIGHT%d' % ( l + 1 ) ] ) )
  if ( len( facet_list ) == 0 ):
    used_facet_list = range( 1, 1 + restore_parameter( facets, 'facet_count' ) )
  else:
    used_facet_list = facet_list
  
  ob_table = wiz_uv.table( 'OB', ni_table.version )
  if ( reference_antenna == 0 ):
    ref_ant = 1
  else:
    ref_ant = reference_antenna
  r = ref_ant - 1
  if ( solution_version == 0 ):
    sn_version = facets.table_highver( 'SN' ) + 1
  else:
    sn_version = solution_version
  try:
    iterations = int( ob_table.keywords[ 'ITER' ] )
  except:
    iterations = 4
  
  # determine time grid for solutions
  time_list = get_time_list( uv )
  time_count = len( time_list )
  
  # read ionospheric fit table
  fit_time_list = []
  fit_count_list = []
  fit_coef_table = []
  fit_weight_list = []
  for row in ni_table:
    fit_time_list.append( float32( row.time ) )
    fit_count_list.append( row.antenna_no )
    fit_coef_table.append( [ float32( coef ) for coef in row.coef ] )
    fit_weight_list.append( float32( row.weight ) )
  fit_time_count = len( fit_time_list )
  fit_coef_table = array( fit_coef_table, dtype = float64 )
  
  # reject bad fits
  fit_weight_array = array( fit_weight_list, dtype = float64 )
  sel = awhere( fit_weight_array > 0. )
  fit_rms_array = 1. / aget( fit_weight_array, sel )
  sub_sel = awhere( fit_rms_array >= rms_min )
  sel = reshape( aget( sel.ravel(), sub_sel ), ( len( sub_sel ), 1 ) ) 
  if ( not rejection_limit is None ):
    fit_rms_array = 1. / aget( fit_weight_array, sel )
    len_rms = len( fit_rms_array ) + 1
    while( len( fit_rms_array ) != len_rms ):
      len_rms = len( fit_rms_array )
      rms_mean = fit_rms_array.mean()
      rms_std = sqrt( ( ( fit_rms_array - rms_mean )**2 ).mean() )
      sub_sel = awhere( fit_rms_array <= rms_mean + rejection_limit * rms_std )
      sel = reshape( aget( sel.ravel(), sub_sel ), ( len( sub_sel ), 1 ) ) 
      fit_rms_array = 1. / aget( fit_weight_array, sel )
    if print_info:
      print 'rejected %s of %s fits...' % ( repr( len( fit_weight_array ) - 
          len( sel ) ), repr( len( fit_weight_array ) ) )
      print 'phase rms mean = %s  stddev = %s' % ( repr( rms_mean ), 
          repr( rms_std ) )
  elif ( not rms_limit is None ):
    fit_rms_array = 1. / aget( fit_weight_array, sel )
    sub_sel = awhere( fit_rms_array <= rms_limit )
    sel = reshape( aget( sel.ravel(), sub_sel ), ( len( sub_sel ), 1 ) ) 
    rms_mean = aget( fit_rms_array, sub_sel ).mean()
    rms_std = sqrt( ( ( aget( fit_rms_array, sub_sel ) - rms_mean )**2 ).mean() )
    if print_info:
      print 'rejected %s of %s fits...' % ( repr( len( fit_weight_array ) - 
          len( sel ) ), repr( len( fit_weight_array ) ) )
      print 'phase rms mean = %s  stddev = %s' % ( repr( rms_mean ), 
          repr( rms_std ) )
  sel = sel.ravel().tolist()
  
  # read ionospheric pierce point table
  pierce_time_list = []
  pierce_X_table = []
  pierce_za_list = []
  pierce_index_table = []
  for row in ob_table:
    pierce_time_list.append( float32( row.time ) )
    pierce_X_table.append( [ float32( x ) for x in row.orbxyz[ 0 : 2 ] ] )
    pierce_za_list.append( float32( row.orbxyz[ 2 ] ) )
    pierce_index_table.append( [ int( round( row.orientation ) ) - 1, 
        row.subarray - 1, row.antenna_no - 1 ] )
  pierce_time_count = len( pierce_time_list )
  pierce_time_array = array( pierce_time_list, dtype = float64 )
  pierce_index_array = array( pierce_index_table, dtype = int64 )
  pierce_X_array = array( pierce_X_table, dtype = float64 )
  pierce_za_array = array( pierce_za_list, dtype = float64 )
  
  # get other relevant data
  calibration_data = get_phase_calibration_data( uv, facets, time_info = False, 
      source_info = True, antenna_info = True, calibration_info = False, 
      facet_list = used_facet_list, print_info = print_info )
  center_table = calibration_data[ 1 ]
  source_table = calibration_data[ 2 ]
  array_table = calibration_data[ 3 ]
  antenna_table = calibration_data[ 4 ]
  source_count = len( source_table )
  antenna_count = len( antenna_table )
  
  # generate time table
  fit_time_table = []
  fit_gst_list = get_gst_list( uv, fit_time_list )
  for n in range( fit_time_count ):
    fit_time_table.append( [ fit_time_list[ n ], fit_gst_list[ n ] ] )
  
  # generate reference antenna table
  reference_list = []
  for k in range( source_count ):
    reference_list.append( r )
  
  # generate solutions and write them to facets
  if print_info:
    print 'generating solutions ...'
  phase_table = - 360. * ones( ( time_count, source_count, antenna_count ), 
      dtype = float64 )
  weight_table = zeros( time_count, dtype = float64 )
  if include_delays:
    delay_table = zeros( ( time_count, source_count, antenna_count ), 
        dtype = float64 )
  
  for n in range( fit_time_count ):

    if ( len( time_steps ) > 0 ):
      if ( not n in time_steps ):
        continue

    if print_info:
      print '... time step n = ', n
    
    # only process non-rejected fits
    if ( not n in sel ):
      continue
    
    # save fit weight
    try:
      nn = time_list.index( fit_time_list[ n ] )
    except ValueError:
      raise error( 'no match for solution time in uv time' )
    weight_table[ nn ] = 1. / radians( 1. / fit_weight_list[ n ] )
    
    # get model fit parameters
    P = fit_coef_table[ n ]
    weight = fit_weight_list[ n ]
    
    if print_info:
      print '...... calculating base vectors'
    
    active_antennas = []
    Xpl_table = []
    pzal_table = []
    Bl_table = []
    for l in range( layer_count ):
      Xp_table = []
      pza_table = []
      
      # get pierce points 
#      for np in range( pierce_time_count ):
#        if ( ( pierce_time_list[ np ] == fit_time_list[ n ] ) and 
#             ( pierce_index_table[ np ][ 0 ] == l ) ):
#          Xp_table.append( pierce_X_table[ np ] )
#          pza_table.append( pierce_za_list[ np ] )
#          if ( not pierce_index_table[ np ][ 2 ] in active_antennas ):
#            active_antennas.append( pierce_index_table[ np ][ 2 ] )
      sel2 = awhere( ( pierce_time_array == fit_time_list[ n ] ) & 
          ( pierce_index_array[ : , 0 ] == l ) )
      Xp_table = Xp_table + aget( pierce_X_array, sel2 ).tolist()
      pza_table = pza_table + aget( pierce_za_array, sel2 ).tolist()
      for i in aget( pierce_index_array[ : , 2 ], sel2 ).tolist():
        if ( not i in active_antennas ):
           active_antennas.append( i )

      Xp_table = array( Xp_table, dtype = float64 )
      pza_table = array( pza_table, dtype = float64 )
      
      # calculate structure matrix
      p_count = len( Xp_table )
      if ( p_count != fit_count_list[ n ] ):
        raise error( 'pierce count does not match' )
      D_table = resize( Xp_table, ( p_count, p_count, 2 ) )
      D_table = transpose( D_table, ( 1, 0, 2 ) ) - D_table
      D_table = add.reduce( D_table**2, 2 )
      D_table = ( D_table / ( r_0**2 ) )**( beta / 2. )
      
      # calculate covariance matrix C
      # calculate partial product for interpolation B
      # reforce symmetry
      C_table = - D_table / 2.
      C_table = transpose( C_table - ( add.reduce( C_table, 0 ) / 
          float( p_count ) ) )
      B_table = add.reduce( C_table, 0 ) / float( p_count )
      C_table = C_table - B_table
      C_table = ( C_table + transpose( C_table ) ) / 2.
      
      # incorporate airmass functions and layer weights
      # add layer C to total C
      A_table = resize( 1. / cos( aradians( pza_table ) ), ( p_count, p_count ) )
      A_table = A_table * transpose( A_table )
      if ( l == 0 ):
        Cl_table = C_table * A_table * ( layer_weights[ l ]**2 )
      else:
        Cl_table = Cl_table + C_table * A_table * ( layer_weights[ l ]**2 )
      
      # save tables per height layer
      Xpl_table.append( Xp_table )
      pzal_table.append( pza_table )
      Bl_table.append( B_table )
    
    # convert to arrays
    Xpl_table = array( Xpl_table, dtype = float64 )
    pzal_table = array( pzal_table, dtype = float64 )
    Bl_table = array( Bl_table, dtype = float64 )
    
    # eigenvalue decomposition
    # reforce symmetry
    # select subset of base vectors
    [ U_table, S, Ut_table ] = linalg.svd( Cl_table )
    U_table = ( U_table + transpose( Ut_table ) ) / 2.
    U_table = U_table[ : , 0 : order ]
    S = S[ 0 : order ]
    
    # calculate interpolation matrix
    F_table = dot( U_table, P / S )
    
    if print_info:
      print '...... calculating pierce point coordinates'
    
    Xl_table = []
    zal_table = []
    ref_list = []
    ref_table = []
    for l in range( layer_count ):
      
      X_table = []
      za_table = []
      
      # get pierce point coordinates
      pierce_table = calculate_pierce_coordinates( fit_time_table[ n ], 
          center_table, source_table, array_table, antenna_table, 
          height = layer_heights[ l ], iterations = iterations )
      
      # put all new pierce points into one array
      j = 0
      for pierce_info in pierce_table:
        [ X, za, [ k, i ] ] = pierce_info
#        if ( i in active_antennas ):
        if True:
          X_table.append( X )
          za_table.append( za )
          if ( l == 0 ):
            if ( i == r ):
              ref_list.append( [ k, i, j ] )
            j = j + 1
        elif ( i == r ):
          raise error( 'no solution for reference antenna' )
      Xl_table.append( X_table )
      zal_table.append( za_table )
      
      # loop over pierce points
      # store index to reference antenna
      if ( l == 0 ):
        for pierce_info in pierce_table:
          [ X, pza, [ k, i ] ] = pierce_info
          j = ( [ ( k == ref[ 0 ] ) for ref in ref_list ] ).index( True )
          ref_table.append( ref_list[ j ][ 2 ] )
    
    Xl_table = array( Xl_table, dtype = float64 )
    zal_table = array( zal_table, dtype = float64 )
    
    if print_info:
      print '...... generating solutions'
    
    # calculate pierce point model solutions
    phi_table = phi_mkl_model( layer_weights, Xl_table, zal_table, Xpl_table, 
        pzal_table, Bl_table, F_table, beta = beta, r_0 = r_0 )
    
    # only store phase corrections of active antennas
    for j in range( len( pierce_table ) ):
      [ X, pza, [ k, i ] ] = pierce_table[ j ]
      if ( i in active_antennas ):
        [ ref_X, ref_pza, [ ref_k, ref_i ] ] = pierce_table[ ref_table[ j ] ]
        phi = phi_table[ j ]
        ref_phi = phi_table[ ref_table[ j ] ]
        phase_table[ nn, k, i ] = amodulo( ( phi - ref_phi ) + 180., 360. ) - 180.
        if include_delays:
          delay_table[ nn, k, i ] = ( delay_sign * ( phi / 360. ) / 
              reference_frequency )
  
  # write solutions to facets
  if print_info:
    print 'writing solution tables ...'
  
  delay = 0.
  for k in range( source_count ):
    solution_table = []
    for n in range( time_count ):
      
      # generate solution row per time step
      solution_row = [ [ time_list[ n ], ref_ant, 0., 0. ] ]
      weight = weight_table[ n ]
      for i in range( antenna_count ):
        if ( phase_table[ n, k, i ] < - 300. ):
          solution_row.append( [ 0., 0., 0., 0. ] ) # flag vis
        else:
          solution = r_phi_to_complex( [ 1., phase_table[ n, k, i ] ] )
          if include_delays:
            delay = delay_table[ n, k, i ]
          solution_row.append( [ solution.real, solution.imag, delay, weight ] )
      
      # add solution row to solution table
      solution_table.append( [ [ xx for xx in x ] for x in solution_row ] )
    
    # write solution table to facet
    if print_info:
      print '... writing solution table for source ', k + 1
    facet = get_facet( facets, used_facet_list[ k ] )
    write_solution_table( facet, solution_table, out_version = sn_version )
  
  return

###############################################################################

def fit_ionospheric_pmkl_model( uv, facets, facet_list = [], beta = 5./3., r_0 = 1., 
    layer_heights = [ 250.e3,350.e3 ], layer_weights = [ 0.5,0.5 ], order = 15, 
    iterations = 1, equal_weights = True, normalize_weights = True, 
    e_steps = 1, double_precision = True, print_info = True, 
    estimate_offsets = True, time_steps = [], include_airmass = True, 
    solution_version = 0, source_rms_limit = 45., first_order = False,
    antenna_offset_limit = 1.e6, antenna_rms_limit = 50., fix_gradient = True,
    prop_time_limit = 1., prop_rms_limit = 60., dof_factor = 3,
    exclude_antennas = [], density_scale = 1.e3, density_power = 1. ):
# prop_time_limit in minutes
  
  if double_precision:
    dtype = float64
  else:
    dtype = float32
  l_max = layer_weights.index( max( layer_weights ) )
  
  # select appropriate model functions
  # initialize model parameters
  model = 'pmkl'
  layer_count = len( layer_heights )
  P = zeros( shape = ( order ), dtype = dtype )
  if ( len( layer_weights ) != layer_count ):
    raise error( 'layer heights and weights are different lengths' )
  
  if print_info:
    print 'starting model fits using model = ' + model
  
  # get all relevant data
  calibration_data = get_phase_calibration_data( uv, facets, facet_list = facet_list,
      solution_version = solution_version, print_info = print_info )
  time_table = calibration_data[ 0 ]
  center_table = calibration_data[ 1 ]
  peel_table = calibration_data[ 2 ]
  array_table = calibration_data[ 3 ]
  antenna_table = calibration_data[ 4 ]
  reference_array = calibration_data[ 5 ]
  phase_array = calibration_data[ 6 ]
  error_array = calibration_data[ 7 ]
  source_count = len( peel_table )
  antenna_count = len( antenna_table )
  
  # select most common reference antenna
  reference_antenna = argmax( bincount( array( reference_array ).ravel() ) ) + 1
  
  # compensate for antenna density, downweight dense antennas
  weight_table = [ 1. for a in antenna_table ]
  if ( ( not density_scale is None ) and ( density_scale > 0. ) ):
    X_table = []
    for a in antenna_table:
      if ( a[ 0 ] != [ 0., 0., 0. ] ):
        X_table.append( a[ 0 ] )
    X_count = len( X_table )
    X_table = array( X_table, dtype = float64 )
    D_table = resize( X_table, ( X_count, X_count, 3 ) )
    D_table = transpose( D_table, ( 1, 0, 2 ) ) - D_table
    D_table = add.reduce( D_table**2, 2 )
#    W_table = exp( - D_table / density_scale**2 )
    W_table = D_table / density_scale**2
    W_table = aput( W_table, awhere( W_table > 100. ), 100. )
    W_table = exp( - W_table )
    W_table = 1. / add.reduce( W_table, 0 )
    if ( not density_power is None ):
      W_table = W_table**density_power
    j = 0
    for i in range( len( antenna_table ) ):
      if ( antenna_table[ i ][ 0 ] != [ 0., 0., 0. ] ):
        weight_table[ i ] = W_table[ j ]
        j = j + 1
  
  fit_table = []
  pierce_table = []
  
  if estimate_offsets:
    solution_table = []
  
  prop_last_time = -1.e9
#  prop_last_G = zeros( ( 2 ), dtype = dtype )
  prop_last_poly = zeros( ( 5 ), dtype = dtype )
  
  # loop over time stamps
  if ( len( time_steps ) == 0 ):
    n_list = range( len( time_table ) )
  else:
    n_list = time_steps
  for n in n_list:
    
    if print_info:
      print '... time step n = %d / %d' % ( n + 1, len( n_list ) )
    
    repeat_fit = True
    rejected_sources = False
    rejected_antennas = False
    while repeat_fit:
      skip = False
      repeat_fit = False
      if ( not rejected_sources ):
        use_sources = [ True for k in range( source_count ) ]
      if ( not rejected_antennas ):
        use_antennas = [ True for i in range( antenna_count ) ]
        for i in exclude_antennas:
          use_antennas[ i - 1 ] = False
      
      temp_pierce_table = []
      
      # containers for model input data
      ref_list = []
      ref_table = []
      phase_table = []
      error_table = []
      info_table = []
      Xpl_table = []
      pzal_table = []
      Bl_table = []
      
      if estimate_offsets:
        solution_row = [ [ time_table[ n ][ 0 ], reference_antenna, 0., 0. ] ] + [
            [ 0., 0., 0., 0. ] for i in range( antenna_count ) ]
      
      # loop over height layers
      for l in range( layer_count ):
        height = layer_heights[ l ]
        
        # containers for model input data
        Xp_table = []
        pza_table = []
        
        # get pierce point coordinates
        ipp_table = calculate_pierce_coordinates( time_table[ n ], 
            center_table, peel_table, array_table, antenna_table, height = height, 
            iterations = iterations )
        
        # loop over pierce points
        j = 0
        for ipp_info in ipp_table:
          [ X, pza, [ k, i ] ] = ipp_info
          if ( ( error_array[ n ][ k ][ i ] > 0. ) and
              use_sources[ k ] and use_antennas[ i ] ):
            
            # store model fit input data
            Xp_table.append( X )
            pza_table.append( pza )
            if ( l == 0 ):
              phase_table.append( phase_array[ n ][ k ][ i ] )
              if equal_weights:
                error_table.append( 1. / weight_table[ i ] )
              else:
                error_table.append( error_array[ n ][ k ][ i ] / weight_table[ i ] )
              if ( i == reference_array[ n ][ k ] ):
                ref_list.append( [ k, i, j ] )
              j = j + 1
              info_table.append( [ k, i, reference_array[ n ][ k ] ] )
            
            # temporarily store pierce point info
            temp_pierce_table.append( [ time_table[ n ][ 0 ], l + 1, k + 1, 
                i + 1, X, pza ] )
        
        if ( l == 0 ):
          # loop over pierce points
          for ipp_info in ipp_table:
            [ X, pza, [ k, i ] ] = ipp_info
            if ( ( error_array[ n ][ k ][ i ] > 0. ) and
                use_sources[ k ] and use_antennas[ i ] ):
              
              # store index to reference antenna
              j = ( [ ( k == ref[ 0 ] ) for ref in ref_list ] ).index( True )
              ref_table.append( ref_list[ j ][ 2 ] )
        
        # skip times of poor data
        if ( len( Xp_table ) <= dof_factor * order + 5 ):
          if print_info:
            print '...... skipping time step due to no / too few data points'
          fit_table.append( [ time_table[ n ][ 0 ], 0, 
              zeros( ( 5 ), dtype = dtype ), azeros( P ).tolist(), 0. ] )
          for p in temp_pierce_table:
            pierce_table.append( p )
          if estimate_offsets:
            solution_table.append( solution_row )
          skip = True
          break
        
        # convert tables to arrays for fitting routines
        Xp_table = array( Xp_table, dtype = dtype )
        pza_table = array( pza_table, dtype = dtype )
        if ( not include_airmass ):
          pza_table = azeros( pza_table )
        if ( l == 0 ):
          ref_table = array( ref_table, dtype = int32 )
          phase_table = array( phase_table, dtype = dtype )
          error_table = array( error_table, dtype = dtype )
#          if equal_weights:
#            error_table = ones( shape( error_table ), dtype = dtype )
          info_table = array( info_table, dtype = int32 )
        p_count = len( Xp_table )
        
        if ( order > 0 ):
          # calculate structure matrix
          D_table = resize( Xp_table, ( p_count, p_count, 2 ) )
          D_table = transpose( D_table, ( 1, 0, 2 ) ) - D_table
          D_table = add.reduce( D_table**2, 2 )
          D_table = ( D_table / ( r_0**2 ) )**( beta / 2. )
          
          # calculate covariance matrix C
          # calculate partial product for interpolation B
          C_table = - D_table / 2.
          C_table = transpose( C_table - ( add.reduce( C_table, 0 ) / float( p_count ) ) )
          B_table = add.reduce( C_table, 0 ) / float( p_count )
          C_table = C_table - B_table
          C_table = ( C_table + transpose( C_table ) ) / 2.
          
          # incorporate airmass functions and layer weights
          # add layer C to total C
          A_table = resize( 1. / cos( aradians( pza_table ) ), ( p_count, p_count ) )
          A_table = A_table * transpose( A_table )
          if ( l == 0 ):
            Cl_table = C_table * A_table * ( layer_weights[ l ]**2 )
          else:
            Cl_table = Cl_table + C_table * A_table * ( layer_weights[ l ]**2 )
          if ( l == l_max ):
            R_table = dot( C_table, diag( 1. / cos( aradians( pza_table ) ) ) )
        
        # save tables per height layer
        Xpl_table.append( Xp_table )
        pzal_table.append( pza_table )
        if ( order > 0 ):
          Bl_table.append( B_table )
      
      # skip times of poor data
      if skip:
        continue
      
      # convert to arrays
      Xpl_table = array( Xpl_table, dtype = dtype )
      pzal_table = array( pzal_table, dtype = dtype )
      if ( order > 0 ):
        Bl_table = array( Bl_table, dtype = dtype )
        
        # eigenvalue decomposition
        # reforce symmetry
        # select subset of base vectors
        try:
          [ U_table, S, Ut_table ] = linalg.svd( Cl_table )
        except:
          if print_info:
            print '... SVD did not converge, skipping time interval'
#          pdb.set_trace()
          fit_table.append( [ time_table[ n ][ 0 ], 0, 
              zeros( ( 5 ), dtype = dtype ), azeros( P ).tolist(), 0. ] )
          for p in temp_pierce_table:
            pierce_table.append( p )
          if estimate_offsets:
            solution_table.append( solution_row )
          skip = True
        else:
          U_table = ( U_table + transpose( Ut_table ) ) / 2.
          U_table = U_table[ : , 0 : order ]
          S = S[ 0 : order ]
          # calculate projection matrix for dominant layer
          R_table = dot( R_table, dot( U_table, diag( 1. / S ) ) )
          R_table = dot( linalg.inv( dot( transpose( R_table ), R_table ) ),
              transpose( R_table ) )
      
      # skip times of poor data
      if skip:
        continue
      
      # see if we can use the gradient of the previous good fit as a starting point
      dtime = abs( time_table[ n ][ 0 ] - prop_last_time ) * 24. * 60.
      Xp_table = Xpl_table[ l_max ] / cos( aradians( pzal_table[ l_max ] ) ).reshape( 
          ( p_count, 1 ) )
      if ( dtime < prop_time_limit ):
#        G = prop_last_G.copy()
        poly = prop_last_poly.copy()
        if first_order:
          poly[ 2 : 5 ] = 0.
      else:
        poly = azeros( prop_last_poly )
        # estimate initial phase gradient on layer with heighest weight
#        G = estimate_phase_gradient_through_fft( uv, Xp_table, ref_table, phase_table,
#            error_table, fft_cell_size = 200., fft_image_size = 1024 )
##        if print_info:
##          print '... estimated gradient is %s' % repr( G.tolist() )
        poly[ 0 : 2 ] = estimate_phase_gradient_through_fft( uv, Xp_table,
            ref_table, phase_table, error_table, fft_cell_size = 200.,
            fft_image_size = 1024 )
      
      # refine gradient fit
      function_keywords = { 'phi_model' : phi_poly_model, 
          'X_table' : Xpl_table[ l_max ], 'pza_table' : pzal_table[ l_max ], 
          'ref_table' : ref_table, 'phase_table' : phase_table,
          'error_table' : error_table, 'normalize_weights' : normalize_weights,
          'wrap_phases' : True }
      parameter_info = []
      for m in range( 5 ):
        if ( m < 2 ):
          par_info = { 'parname' : 'P_%d' % ( m ), 'value' : poly[ m ], 
              'limits' : [ None, None ] }
        else:
          par_info = { 'parname' : 'P_%d' % ( m ), 'value' : poly[ m ], 
              'limits' : [ poly[ m ], poly[ m ] ] }
        parameter_info.append( par_info )
      fit = mpfit( fit_phi_model, functkw = function_keywords, parinfo = parameter_info,
          quiet = True, autoderivative = True, debug = False, fastnorm = False, 
          nocovar = True, dblprec = double_precision )
      poly = fit.params.copy()
#      G = poly[ 0 : 2 ]
#      chi_grad = sqrt( fit.fnorm )
#      if print_info:
#        print '... gradient post-fit phase RMS is %s degrees' % repr( chi_grad )
      chi_poly = sqrt( fit.fnorm )
      
      if ( not first_order ):
        # refine higher order fit
        function_keywords = { 'phi_model' : phi_poly_model, 
            'X_table' : Xpl_table[ l_max ], 'pza_table' : pzal_table[ l_max ], 
            'ref_table' : ref_table, 'phase_table' : phase_table,
            'error_table' : error_table, 'normalize_weights' : normalize_weights,
            'wrap_phases' : True }
        parameter_info = []
        for m in range( 5 ):
          if ( fix_gradient and ( m < 2 ) ):
            par_info = { 'parname' : 'P_%d' % ( m ), 'value' : poly[ m ], 
                'limits' : [ poly[ m ], poly[ m ] ] }
          else:
            par_info = { 'parname' : 'P_%d' % ( m ), 'value' : poly[ m ], 
                'limits' : [ None, None ] }
          parameter_info.append( par_info )
        fit = mpfit( fit_phi_model, functkw = function_keywords, quiet = True, 
            autoderivative = True, debug = False, fastnorm = False, nocovar = True,
            dblprec = double_precision, parinfo = parameter_info )
        poly = fit.params.copy()
        chi_poly = sqrt( fit.fnorm )
      if print_info:
        print '... polynomial post-fit phase RMS is %s degrees' % repr( chi_poly )
      
#      if first_order: # ( first_order or ( chi_grad < chi_poly ) ):
#        poly[ 2 : 5 ] = [ 0., 0., 0. ]
#        chi_poly = chi_grad
      
      # remove polynomial model from data
      phi_poly_table = phi_poly_model( Xpl_table[ l_max ], poly ) / cos( 
          aradians( pzal_table[ l_max ] ) )
      sel = ref_table.reshape( p_count, 1 )
      phi_poly_table = phi_poly_table - aget( phi_poly_table, sel )
      if ( order > 0 ):
        phase_table = amodulo( ( phase_table - phi_poly_table ) + 180., 360. ) - 180.
        
        # define data passed to the model fit routines
        # define model parameter structure
        function_keywords = { 'U_table' : U_table, 'ref_table' : ref_table, 
            'phase_table' : phase_table, 'error_table' : error_table, 
            'normalize_weights' : normalize_weights }
        
        # specify fit parameters
        parameter_info = []
        for m in range( order ):
          par_info = { 'parname' : 'P_%d' % ( m + 1 ), 'value' : 0., 
              'limits' : [ None, None ] }
          parameter_info.append( par_info )
        
        # fit model to data points
        fit = mpfit( fit_phi_mkl_model, functkw = function_keywords, 
            parinfo = parameter_info, quiet = True, autoderivative = True, 
            debug = False, fastnorm = False, nocovar = True,
            dblprec = double_precision )
#        if print_info:
#          print '...... mpfit() returned status %d' % ( fit.status )
        if ( fit.errmsg != '' ):
#          raise error( fit.errmsg )
          if print_info:
            print 'WARNING: %s' % ( fit.errmsg )
          fit_table.append( [ time_table[ n ][ 0 ], 0, 
              zeros( ( 5 ), dtype = dtype ), azeros( P ).tolist(), 0. ] )
          for p in temp_pierce_table:
            pierce_table.append( p )
          if estimate_offsets:
            solution_table.append( solution_row )
          skip = True
          break
        
        # store fit results
        P = fit.params.copy()
        chi_pmkl = sqrt( fit.fnorm )
        
        if print_info:
#          print '... fitted parameters are %s' % repr( P.tolist() )
#          print '... used e-value is %s' % repr( best_e )
          print '... model post-fit phase RMS is %s degrees' % repr( chi_pmkl )
      
      # skip times of error in fit
      if skip:
        continue
      
      if ( order > 0 ):
        # check that model fit is better than polynomial fit
        if ( chi_pmkl < chi_poly ):
          chi = chi_pmkl
        else:
          if print_info:
            print '... WARNING: PMKL model fit is no improvement'
          P = azeros( P )
#          chi = chi_poly
          chi = chi_pmkl
      else:
        chi = chi_poly
      
      # repeat bad fit if previous gradient was used
      if ( ( dtime < prop_time_limit ) and ( chi > prop_rms_limit ) ):
        if print_info:
          print '...... repeating fit using fresh gradient estimate'
        rejected_sources = False
        rejected_antennas = False
        repeat_fit = True
        prop_last_time = -1.e9
        continue
      
      if ( order > 0 ):
        # add polynomial model back
        phase_table = amodulo( ( phase_table + phi_poly_table ) + 180., 360. ) - 180.
        
        # calculate model phases
        phi_mkl_table = dot( U_table, P )
        sel = ref_table.reshape( p_count, 1 )
        phi_mkl_table = phi_mkl_table - aget( phi_mkl_table, sel )
        phi_table = amodulo( ( phi_poly_table + phi_mkl_table ) + 180., 360. ) - 180.
      else:
        phi_table = amodulo( phi_poly_table + 180., 360. ) - 180.
      
      # split pierce points per source
      ref_list = []
      for ref in ref_table:
        if ( not ( ref in ref_list ) ):
          ref_list.append( ref )
      
      # estimate residual phase per antenna per source
      antenna_error_table = azeros( error_table )
      antenna_offset_table = azeros( error_table )
      reference_antennas = []
      for ref in ref_list:
        reference_antennas.append( info_table[ ref ][ 1 ] )
        sel = awhere( ref_table == ref )
        dphase_table = aget( phase_table, sel )
        dphase_table = resize( dphase_table, ( len( sel ), len( sel ) ) )
        dphase_table = transpose( dphase_table ) - dphase_table
        dphi_table = aget( phi_table, sel )
        dphi_table = resize( dphi_table, ( len( sel ), len( sel ) ) )
        dphi_table = transpose( dphi_table ) - dphi_table
        dweight_table = aget( error_table, sel )
        dweight_table = resize( dweight_table**2, ( len( sel ), len( sel ) ) )
        dweight_table = transpose( dweight_table ) + dweight_table
        dweight_table = 1. / sqrt( dweight_table )
        abs_err_table = add.reduce( dweight_table * 
            abs( amodulo( ( dphase_table - dphi_table ) + 180., 360. ) - 180. ), 1 )
        normalize_table = add.reduce( dweight_table, 1 ) - diagonal( dweight_table )
        if ( len( sel ) == 1 ):
          abs_err_table = azeros( abs_err_table )
        else:
          abs_err_table = abs_err_table / normalize_table
        antenna_error_table = aput( antenna_error_table, sel, abs_err_table )
        offset_table = add.reduce( dweight_table * 
            ( amodulo( ( dphase_table - dphi_table ) + 180., 360. ) - 180. ), 1 )
        if ( len( sel ) == 1 ):
          offset_table = azeros( offset_table )
        else:
          offset_table = offset_table / normalize_table
        antenna_offset_table = aput( antenna_offset_table, sel, offset_table )
      
      # estimate residual phase offset and stddev per antenna
      antenna_means = zeros( ( antenna_count ), dtype = dtype )
      antenna_stds = zeros( ( antenna_count ), dtype = dtype )
      for i in range( antenna_count ):
        # select all phases related to single antenna
        sel = awhere( info_table[ : , 1 ] == i )
        # TODO: determine reasonable number of values for calculating the antenna mean
        if ( len( sel ) < 3 ):
          continue
        # calculate mean phase offset per antenna
        phase_offsets = aget( antenna_offset_table, sel )
        mean_phase_offset = amean_phase( phase_offsets )
        std_phase_offset = sqrt( ( ( amodulo( ( phase_offsets - mean_phase_offset )
            + 180., 360. ) - 180. )**2 ).mean() )
        antenna_means[ i ] = mean_phase_offset
        antenna_stds[ i ] = std_phase_offset
        if estimate_offsets:
          # store phase error in solution table
          # note this automatically maps mean_phase in [ -180, 180 ] domain
          solution = r_phi_to_complex( [ 1., mean_phase_offset ] )
          solution_row[ i + 1 ] = [ solution.real, solution.imag, 0.,
              1. / std_phase_offset ]
      
      # check for antennas with excessive offset or stddev
      if ( not rejected_antennas ):
        sel1 = awhere( antenna_means > antenna_offset_limit )
        for s in sel1:
          s = s[ 0 ]
          if print_info:
            print '...... antenna %s has excessive offset of %s degrees' % (
                repr( s + 1 ), repr( antenna_means[ s ] ) )
          if ( s in reference_antennas ):
            if print_info:
              print '...... WARNING: reference antenna not dropped'
          else:
            use_antennas[ s ] = False
            rejected_antennas = True
        sel2 = awhere( antenna_stds > antenna_rms_limit )
        for s in sel2:
          s = s[ 0 ]
          if use_antennas[ s ]:
            if print_info:
              print '...... antenna %s has excessive stddev of %s degrees' % (
                  repr( s + 1 ), repr( antenna_stds[ s ] ) )
            if ( s in reference_antennas ):
              if print_info:
                print '...... WARNING: reference antenna not dropped'
            else:
              use_antennas[ s ] = False
              rejected_antennas = True
        if rejected_antennas:
          if print_info:
            print '...... repeating fit'
          repeat_fit = True
          if ( chi > prop_rms_limit ):
            prop_last_time = -1.e9
          continue
      
      # estimate residual phase RMS per source
      source_rmss = zeros( ( source_count ), dtype = dtype )
      for k in range( source_count ):
        # select all phases related to single antenna
        sel = awhere( info_table[ : , 0 ] == k )
        # TODO: determine reasonable number of values for calculating the mean
        if ( len( sel ) < 3 ):
          continue
        # calculate mean phase offset per source
        phase_offsets = aget( antenna_offset_table, sel )
        rms_phase_offset = sqrt( ( ( amodulo( phase_offsets
            + 180., 360. ) - 180. )**2 ).mean() )
        source_rmss[ k ] = rms_phase_offset
      
      # check for source with excessive RMS
      if ( not rejected_sources ):
        sel = awhere( source_rmss > source_rms_limit )
        for s in sel:
          s = s[ 0 ]
          use_sources[ s ] = False
          rejected_sources = True
          if print_info:
            print '...... source %s has excessive RMS of %s degrees' % (
                repr( facet_list[ s ] ), repr( source_rmss[ s ] ) )
        if rejected_sources:
          if print_info:
            print '...... repeating fit'
          repeat_fit = True
          if ( chi > prop_rms_limit ):
            prop_last_time = -1.e9
          continue
      
      # store fit results
      weight = 1. / chi
      fit_table.append( [ time_table[ n ][ 0 ], p_count, poly.tolist(),
          P.tolist(), weight ] )
      for p in temp_pierce_table:
        pierce_table.append( p )
      
      if estimate_offsets:
        # add solution row to solution table
        solution_table.append( solution_row )
      
#      # see if we can use the gradient of the previous good fit as a starting point
#      if ( chi < prop_rms_limit ):
#        prop_last_G = G.copy()
#        prop_last_time = time_table[ n ][ 0 ]
      
      # see if we can use the polynomial of the previous good fit as a starting point
      if ( chi < prop_rms_limit ):
        prop_last_poly = poly.copy()
        prop_last_time = time_table[ n ][ 0 ]
  
  # write fit results to model fit table
  write_ionospheric_pmkl_model_fit_table( uv, layer_heights, layer_weights,
       fit_table, beta = beta, r_0 = r_0 )
  fit_version = uv.table_highver( 'NI' )
  write_ionospheric_mkl_pierce_table( uv, pierce_table, version = fit_version,
      iterations = iterations )
  
  # write solution table with antenna phase offsets to uv data
  if estimate_offsets:
    if print_info:
      print '... NOTE: writing solution table with antenna phase offsets to UV data'
    write_solution_table( uv, solution_table, out_version = 0 )
  
  return

###############################################################################

def write_ionospheric_pmkl_model_fit_table( uv, layer_heights, layer_weights,
    fit_table, version = 0, beta = 5. / 3., r_0 = 1., **keywords ):

  # create new NI table
  layer_count = len( layer_heights )
  order = len( fit_table[ 0 ][ 2 ] ) + len( fit_table[ 0 ][ 3 ] )
  new_ni_table = new_table( uv, 'NI', version, num_coef = order )
  row = new_table_row( new_ni_table )
  integration_time = restore_parameter( uv, 'integration_time' )
  time_interval = integration_time / ( 24. * 60. * 60. )

  # loop over 
  for fit_row in fit_table:
    [ time, p_count, poly, coefs, weight ] = fit_row
    row.time = float32( time )
    row.time_interval = float32( time_interval )
    row.antenna_no = p_count
    row.source_id = 0
    row.subarray = 0
    row.weight = float32( weight )
    row.coef = [ float32( pol ) for pol in poly ] + [ float32( coef ) for coef in coefs ]
    new_ni_table.append( row )

  # add keywords
  new_ni_table.keywords[ 'MODEL' ] = 'pmkl'
  new_ni_table.keywords[ 'REF_FREQ' ] = float32( get_central_frequency( uv ) )
  new_ni_table.keywords[ 'BETA' ] = float32( beta )
  new_ni_table.keywords[ 'R_0' ] = float32( r_0 )
  new_ni_table.keywords[ 'LAYERS' ] = layer_count
  for l in range( layer_count ):
    new_ni_table.keywords[ 'HEIGHT%d' % ( l + 1 ) ] = float32( layer_heights[ l ] )
    new_ni_table.keywords[ 'WEIGHT%d' % ( l + 1 ) ] = float32( layer_weights[ l ] )
  for key in keywords.keys():
    new_ni_table.keywords[ key ] = keywords[ key ]

  #close new table
  new_ni_table.close()
  
  return

###############################################################################

def generate_solution_tables_from_pmkl_fit_table( uv, facets, fit_version = 0,
    reference_antenna = 0, solution_version = 0, include_delays = True,
    print_info = True, rejection_limit = None, rms_limit = 35., facet_list = [],
    time_steps = [], rms_min = 0.1, include_antennas = [], include_phases = True,
    include_delay_phases = True ):
  
  wiz_uv = wizardry( uv )
  ni_table = wiz_uv.table( 'NI', fit_version )
  model = ni_table.keywords[ 'MODEL' ]
  if ( type( model ) == type( [] ) ):
    model = model[ 0 ]
  model = model.strip()
  if ( model != 'pmkl' ):
    raise error( 'unknown model: %s' % ( model ) )
  order = int( ni_table.keywords[ 'NUM_COEF' ] ) - 5
  if include_delays:
    reference_frequency = get_frequency( uv )
    bandwidth = get_bandwidth( uv )
    solution_frequency = float32( ni_table.keywords[ 'REF_FREQ' ] )
    fractional_bandwidth = bandwidth / solution_frequency
  beta = float32( ni_table.keywords[ 'BETA' ] )
  r_0 = float32( ni_table.keywords[ 'R_0' ] )
  layer_count = int( ni_table.keywords[ 'LAYERS' ] )
  layer_heights = []
  layer_weights = []
  for l in range( layer_count ):
    layer_heights.append( float32( ni_table.keywords[ 'HEIGHT%d' % ( l + 1 ) ] ) )
    layer_weights.append( float32( ni_table.keywords[ 'WEIGHT%d' % ( l + 1 ) ] ) )
  l_max = layer_weights.index( max( layer_weights ) )
  if ( len( facet_list ) == 0 ):
    used_facet_list = range( 1, 1 + restore_parameter( facets, 'facet_count' ) )
  else:
    used_facet_list = facet_list
  
  ob_table = wiz_uv.table( 'OB', ni_table.version )
  if ( reference_antenna == 0 ):
    ref_ant = 1
  else:
    ref_ant = reference_antenna
  r = ref_ant - 1
  try:
    iterations = int( ob_table.keywords[ 'ITER' ] )
  except:
    iterations = 4
  
  # determine time grid for solutions
  time_list = get_time_list( uv )
  time_count = len( time_list )
  
  # read ionospheric fit table
  fit_time_list = []
  fit_count_list = []
  fit_poly_table = []
  fit_coef_table = []
  fit_weight_list = []
  for row in ni_table:
    fit_time_list.append( float32( row.time ) )
    fit_count_list.append( row.antenna_no )
    fit_poly_table.append( [ float32( pol ) for pol in row.coef[ 0 : 5 ] ] )
    fit_coef_table.append( [ float32( coef ) for coef in row.coef[ 5 : ] ] )
    fit_weight_list.append( float32( row.weight ) )
  fit_time_count = len( fit_time_list )
  fit_poly_table = array( fit_poly_table, dtype = float64 )
  fit_coef_table = array( fit_coef_table, dtype = float64 )
  
  # reject bad fits
  fit_weight_array = array( fit_weight_list, dtype = float64 )
  sel = awhere( fit_weight_array > 0. )
  fit_rms_array = 1. / aget( fit_weight_array, sel )
  sub_sel = awhere( fit_rms_array >= rms_min )
  sel = reshape( aget( sel.ravel(), sub_sel ), ( len( sub_sel ), 1 ) ) 
  if ( not rejection_limit is None ):
    fit_rms_array = 1. / aget( fit_weight_array, sel )
    len_rms = len( fit_rms_array ) + 1
    while( len( fit_rms_array ) != len_rms ):
      len_rms = len( fit_rms_array )
      rms_mean = fit_rms_array.mean()
      rms_std = sqrt( ( ( fit_rms_array - rms_mean )**2 ).mean() )
      sub_sel = awhere( fit_rms_array <= rms_mean + rejection_limit * rms_std )
      sel = reshape( aget( sel.ravel(), sub_sel ), ( len( sub_sel ), 1 ) ) 
      fit_rms_array = 1. / aget( fit_weight_array, sel )
    if print_info:
      print 'rejected %s of %s fits...' % ( repr( len( fit_weight_array ) - 
          len( sel ) ), repr( len( fit_weight_array ) ) )
      print 'phase rms mean = %s  stddev = %s' % ( repr( rms_mean ), 
          repr( rms_std ) )
  elif ( not rms_limit is None ):
    fit_rms_array = 1. / aget( fit_weight_array, sel )
    sub_sel = awhere( fit_rms_array <= rms_limit )
    sel = reshape( aget( sel.ravel(), sub_sel ), ( len( sub_sel ), 1 ) ) 
    rms_mean = aget( fit_rms_array, sub_sel ).mean()
    rms_std = sqrt( ( ( aget( fit_rms_array, sub_sel ) - rms_mean )**2 ).mean() )
    if print_info:
      print 'rejected %s of %s fits...' % ( repr( len( fit_weight_array ) - 
          len( sel ) ), repr( len( fit_weight_array ) ) )
      print 'phase rms mean = %s  stddev = %s' % ( repr( rms_mean ), 
          repr( rms_std ) )
  sel = sel.ravel().tolist()
  
  # read ionospheric pierce point table
  pierce_time_list = []
  pierce_X_table = []
  pierce_za_list = []
  pierce_index_table = []
  for row in ob_table:
    pierce_time_list.append( float32( row.time ) )
    pierce_X_table.append( [ float32( x ) for x in row.orbxyz[ 0 : 2 ] ] )
    pierce_za_list.append( float32( row.orbxyz[ 2 ] ) )
    pierce_index_table.append( [ int( round( row.orientation ) ) - 1, 
        row.subarray - 1, row.antenna_no - 1 ] )
  pierce_time_count = len( pierce_time_list )
  pierce_time_array = array( pierce_time_list, dtype = float64 )
  pierce_index_array = array( pierce_index_table, dtype = int64 )
  pierce_X_array = array( pierce_X_table, dtype = float64 )
  pierce_za_array = array( pierce_za_list, dtype = float64 )
  
  # get other relevant data
  calibration_data = get_phase_calibration_data( uv, facets, time_info = False, 
      source_info = True, antenna_info = True, calibration_info = False, 
      facet_list = used_facet_list, print_info = print_info )
  center_table = calibration_data[ 1 ]
  source_table = calibration_data[ 2 ]
  array_table = calibration_data[ 3 ]
  antenna_table = calibration_data[ 4 ]
  source_count = len( source_table )
  antenna_count = len( antenna_table )
  
  # generate time table
  fit_time_table = []
  fit_gst_list = get_gst_list( uv, fit_time_list )
  for n in range( fit_time_count ):
    fit_time_table.append( [ fit_time_list[ n ], fit_gst_list[ n ] ] )
  
  # generate reference antenna table
  reference_list = []
  for k in range( source_count ):
    reference_list.append( r )
  
  # generate solutions and write them to facets
  if print_info:
    print 'generating solutions ...'
  phase_table = - 360. * ones( ( time_count, source_count, antenna_count ), 
      dtype = float64 )
  weight_table = zeros( time_count, dtype = float64 )
  if include_delays:
    delay_table = zeros( ( time_count, source_count, antenna_count ), 
        dtype = float64 )
  
  for n in range( fit_time_count ):
    
    if ( len( time_steps ) > 0 ):
      if ( not n in time_steps ):
        continue
    
    if print_info:
      print '... time step n = %d / %d' % ( n + 1, fit_time_count )
    
    # only process non-rejected fits
    if ( not n in sel ):
      continue
    
    # save fit weight
    try:
      nn = time_list.index( fit_time_list[ n ] )
    except ValueError:
#      raise error( 'no match for solution time in uv time' )
      if print_info:
        print '...... no match for solution time in uv time'
      continue
    weight_table[ nn ] = 1. / radians( 1. / fit_weight_list[ n ] )
    
    # get model fit parameters
    poly = fit_poly_table[ n ]
    P = fit_coef_table[ n ]
    weight = fit_weight_list[ n ]
    
    if print_info:
      print '...... calculating base vectors'
    
    active_antennas = []
    Xpl_table = []
    pzal_table = []
    Bl_table = []
    for l in range( layer_count ):
      Xp_table = []
      pza_table = []
      
      # get pierce points 
#      for np in range( pierce_time_count ):
#        if ( ( pierce_time_list[ np ] == fit_time_list[ n ] ) and 
#             ( pierce_index_table[ np ][ 0 ] == l ) ):
#          Xp_table.append( pierce_X_table[ np ] )
#          pza_table.append( pierce_za_list[ np ] )
#          if ( not pierce_index_table[ np ][ 2 ] in active_antennas ):
#            active_antennas.append( pierce_index_table[ np ][ 2 ] )
      sel2 = awhere( ( pierce_time_array == fit_time_list[ n ] ) & 
          ( pierce_index_array[ : , 0 ] == l ) )
      Xp_table = Xp_table + aget( pierce_X_array, sel2 ).tolist()
      pza_table = pza_table + aget( pierce_za_array, sel2 ).tolist()
      for i in aget( pierce_index_array[ : , 2 ], sel2 ).tolist():
        if ( not i in active_antennas ):
           active_antennas.append( i )
      
      Xp_table = array( Xp_table, dtype = float64 )
      pza_table = array( pza_table, dtype = float64 )
      
      # calculate structure matrix
      p_count = len( Xp_table )
      if ( p_count != fit_count_list[ n ] ):
        raise error( 'pierce count does not match' )
      D_table = resize( Xp_table, ( p_count, p_count, 2 ) )
      D_table = transpose( D_table, ( 1, 0, 2 ) ) - D_table
      D_table = add.reduce( D_table**2, 2 )
      D_table = ( D_table / ( r_0**2 ) )**( beta / 2. )
      
      # calculate covariance matrix C
      # calculate partial product for interpolation B
      # reforce symmetry
      C_table = - D_table / 2.
      C_table = transpose( C_table - ( add.reduce( C_table, 0 ) / 
          float( p_count ) ) )
      B_table = add.reduce( C_table, 0 ) / float( p_count )
      C_table = C_table - B_table
      C_table = ( C_table + transpose( C_table ) ) / 2.
      
      # incorporate airmass functions and layer weights
      # add layer C to total C
      A_table = resize( 1. / cos( aradians( pza_table ) ), ( p_count, p_count ) )
      A_table = A_table * transpose( A_table )
      if ( l == 0 ):
        Cl_table = C_table * A_table * ( layer_weights[ l ]**2 )
      else:
        Cl_table = Cl_table + C_table * A_table * ( layer_weights[ l ]**2 )
      
      # save tables per height layer
      Xpl_table.append( Xp_table )
      pzal_table.append( pza_table )
      Bl_table.append( B_table )
    
    # convert to arrays
    Xpl_table = array( Xpl_table, dtype = float64 )
    pzal_table = array( pzal_table, dtype = float64 )
    Bl_table = array( Bl_table, dtype = float64 )
    
    # eigenvalue decomposition
    # reforce symmetry
    # select subset of base vectors
    try:
      [ U_table, S, Ut_table ] = linalg.svd( Cl_table )
    except:
      if print_info:
        print '...... SVD did not converge, skipping time interval'
      continue
    U_table = ( U_table + transpose( Ut_table ) ) / 2.
    U_table = U_table[ : , 0 : order ]
    S = S[ 0 : order ]
    
    # calculate interpolation matrix
    F_table = dot( U_table, P / S )
    
    if print_info:
      print '...... calculating pierce point coordinates'
    
    Xl_table = []
    zal_table = []
    ref_list = []
    ref_table = []
    for l in range( layer_count ):
      
      X_table = []
      za_table = []
      
      # get pierce point coordinates
      pierce_table = calculate_pierce_coordinates( fit_time_table[ n ], 
          center_table, source_table, array_table, antenna_table, 
          height = layer_heights[ l ], iterations = iterations )
      
      # put all new pierce points into one array
      j = 0
      for pierce_info in pierce_table:
        [ X, za, [ k, i ] ] = pierce_info
#        if ( i in active_antennas ):
        if True:
          X_table.append( X )
          za_table.append( za )
          if ( l == 0 ):
            if ( i == r ):
              ref_list.append( [ k, i, j ] )
            j = j + 1
        elif ( i == r ):
          raise error( 'no solution for reference antenna' )
      Xl_table.append( X_table )
      zal_table.append( za_table )
      
      # loop over pierce points
      # store index to reference antenna
      if ( l == 0 ):
        for pierce_info in pierce_table:
          [ X, pza, [ k, i ] ] = pierce_info
          j = ( [ ( k == ref[ 0 ] ) for ref in ref_list ] ).index( True )
          ref_table.append( ref_list[ j ][ 2 ] )
    
    Xl_table = array( Xl_table, dtype = float64 )
    zal_table = array( zal_table, dtype = float64 )
    
    if print_info:
      print '...... generating solutions'
    
    # calculate pierce point model solutions
    phi_poly_table = phi_poly_model( Xl_table[ l_max ], poly ) / cos( 
        aradians( zal_table[ l_max ] ) )
    phi_mkl_table = phi_mkl_model( layer_weights, Xl_table, zal_table, Xpl_table, 
        pzal_table, Bl_table, F_table, beta = beta, r_0 = r_0 )
    phi_table = phi_poly_table + phi_mkl_table
    
    # only store phase corrections of active antennas
    for j in range( len( pierce_table ) ):
      [ X, pza, [ k, i ] ] = pierce_table[ j ]
      if ( ( i in active_antennas ) or ( ( i + 1 ) in include_antennas ) ):
#        [ ref_X, ref_pza, [ ref_k, ref_i ] ] = pierce_table[ ref_table[ j ] ]
        if include_phases:
          phase_table[ nn, k, i ] = phi_table[ j ] - phi_table[ ref_table[ j ] ]
        else:
          phase_table[ nn, k, i ] = 0.
        if include_delays:
          dphase = phi_table[ j ] - phi_table[ ref_table[ j ] ]
          delay = - ( dphase / 360. ) / ( solution_frequency * 2. * pi )
          delay_table[ nn, k, i ] = delay
          if include_delay_phases:
            dphase = 360. * delay * ( solution_frequency - reference_frequency )
            phase_table[ nn, k, i ] = phase_table[ nn, k, i ] - dphase
        phase_table[ nn, k, i ] = amodulo( phase_table[ nn, k, i ] + 180., 360. ) - 180.
  
  # write solutions to facets
  if print_info:
    print 'writing solution tables ...'
  
  delay = 0.
  for k in range( source_count ):
    solution_table = []
    for n in range( time_count ):
      
      # generate solution row per time step
      solution_row = [ [ time_list[ n ], ref_ant, 0., 0. ] ]
      weight = weight_table[ n ]
      ref_solution = r_phi_to_complex( [ 1., phase_table[ n, k, r ] ] )
      if include_delays:
        ref_delay = delay_table[ n, k, r ]
      for i in range( antenna_count ):
        if ( ( phase_table[ n, k, i ] < -300. ) or ( phase_table[ n, k, r ] < -300. ) ):
          solution_row.append( [ 0., 0., 0., 0. ] ) # flag vis
        else:
          solution = r_phi_to_complex( [ 1., phase_table[ n, k, i ] ] ) / ref_solution
          if include_delays:
            delay = delay_table[ n, k, i ] - ref_delay
          solution_row.append( [ solution.real, solution.imag, delay, weight ] )
      
      # add solution row to solution table
      solution_table.append( [ [ xx for xx in x ] for x in solution_row ] )
    
    # write solution table to facet
    if print_info:
      print '... writing solution table for source ', used_facet_list[ k ]
    facet = get_facet( facets, used_facet_list[ k ] )
    if ( solution_version == 0 ):
      sn_version = facet.table_highver( 'SN' ) + 1
    else:
      sn_version = solution_version
    write_solution_table( facet, solution_table, out_version = sn_version )
  
  return

###############################################################################

def read_phase_calibration_data( uv, time_info = True, solution_version = 0,
    print_info = True, ignore_abundant_times = True ):

# TODO: possibly include delays into model fitting

  # get time info
  time_list = get_time_list( uv )
  time_count = len( time_list )
  time_table = None
  if print_info:
    print 'getting time info ...'
  gst_list = get_gst_list( uv, time_list = time_list )
  time_table = [ [ time_list[ n ], gst_list[ n ] ] for n in range( time_count ) ]

  # get source info
  center_table = None
  if print_info:
    print 'getting source info ...'
  obs_epoch = get_observing_epoch( uv )
  center_table = convert_radec_from_j2000( get_radec( uv ), obs_epoch )

  # get antenna info
  array_table = None
  antenna_table = None
  antenna_list = get_antenna_positions( uv )
  antenna_no_list = [ antenna[ 0 ] for antenna in antenna_list ]
  antenna_count = antenna_no_list[ - 1 ]
  if print_info:
    print 'getting antenna info ...'
  array_xyz = get_mean_antenna_position( uv ) # get_array_position( uv )
  array_table = [ array_xyz, xyz_to_geo_llh( array_xyz ) ]
  antenna_table = []
  for antenna_no in range( 1, 1 + antenna_count ):
    if ( antenna_no in antenna_no_list ):
      antenna_index = antenna_no_list.index( antenna_no )
      antenna_xyz = antenna_list[ antenna_index ][ 1 ]
      antenna_table.append( [ antenna_xyz, xyz_to_geo_llh( antenna_xyz ) ] )
    else:
      antenna_table.append( [ [ 0., 0., 0. ], [ 0., 0., 0. ] ] )
  
  # read calibration info from uv
  reference_table = None
  phase_table = None
  error_table = None
  if print_info:
    print 'getting calibration info ...'
  phase_array = 360. * ones( shape = ( time_count, antenna_count ), dtype = float64 )
  error_array = - 360. * ones( shape = ( time_count, antenna_count ), dtype = float64 )
  reference_array = zeros( shape = ( time_count ), dtype = int32 )
  if ( solution_version <= 0 ):
    sol_version = uv.table_highver( 'SN' ) + solution_version
  else:
    sol_version = solution_version
  solution_table = read_solution_table( uv, in_version = sol_version )
  n = 0
  w = 0
  for solution in solution_table:
    [ time, reference_antenna ] = solution[ 0 ][ 0 : 2 ]
    try:
      n = time_list.index( time )
    except:
      if ignore_abundant_times:
        w = w + 1
      else:
        raise error( 'time in solution table does not correspond to time ' + 
            'in observations' )
    reference_array[ n ] = reference_antenna - 1
    [ ref_gain_real, ref_gain_imag, ref_delay, ref_weight ] = solution[ 
        reference_antenna ][ 0 : 4 ]
    ref_gain = complex( ref_gain_real, ref_gain_imag )
    if ( ( ref_weight == 0. ) or ( ref_gain == complex( 0., 0. ) ) ):
#      raise error( 'reference antenna has no solution for time instance' )
      # skip time instance altogether
      continue
    [ ref_amp, ref_phase ] = complex_to_r_phi( ref_gain )
    for i in range( antenna_count ):
      [ gain_real, gain_imag, delay, weight ] = solution[ i + 1 ][ 0 : 4 ]
      gain = complex( gain_real, gain_imag )
      if ( ( weight == 0. ) or ( ( i != reference_array[ n ] ) and 
          ( gain == complex( 1., 0. ) ) ) or ( gain == complex( 0., 0. ) ) ):
        phase = 360.
        phase_error = - 360.
      else:
        [ amp, phase ] = complex_to_r_phi( gain )
        phase_error = degrees( 1. / weight )
      phase_array[ n ][ i ] = phase
      error_array[ n ][ i ] = phase_error
  if ( print_info and ( w > 0 ) ):
   print ( 'WARNING: %d time(s) in solution table do(es) not ' % ( w ) + 
       'correspond to time(s) in observations' )
  reference_table = reference_array.tolist()
  phase_table = phase_array.tolist()
  error_table = error_array.tolist()

  return [ time_table, center_table, array_table, antenna_table, 
      reference_table, phase_table, error_table ]

###############################################################################

def correlate_baseline_phases( uv, in_version = 0, max_period = 30.,
    min_period = 20., max_gap = 1., max_shift = 5., remove_average = True,
    time_step = 1, height = 300.e3, iterations = 1,
    print_info = True, tolerance = 0.1, antennas = [], min_points = 10 ):
# max_gap, max_period & min_period in minutes
  
  normalize_weights = True
  debug = False
  double_precision = True
  if double_precision:
    dtype = float64
  else:
    dtype = float32
  
  calibration_data = read_phase_calibration_data( uv, solution_version = in_version,
      print_info = print_info )
  time_array = array( calibration_data[ 0 ], dtype = dtype )
  center_table = calibration_data[ 1 ]
  array_table = calibration_data[ 2 ]
  antenna_table = calibration_data[ 3 ]
  reference_array = calibration_data[ 4 ]
  ph_array = array( calibration_data[ 5 ], dtype = dtype )
  er_array = array( calibration_data[ 6 ], dtype = dtype )
  time_count = len( time_array )
  antenna_count = len( antenna_table )
  if ( len( antennas ) > 0 ):
    ants = [ i for i in antennas ]
    ants.sort( cmp = lambda a, b: cmp( a, b ) )
  else:
    ants = range( 1, 1 + antenna_count )
  
  # calculate some time variables
  dtime = median( time_array[ 1 : , 0 ] - time_array[ : -1, 0 ] )
  time_tol = tolerance * dtime
  dtime = dtime * float( time_step )
  
  # identify scans or large gaps
  scan_list = []
  last_n = - 1
  first_n = 0
  for n in range( time_count ):
    sel = awhere( er_array[ n, : ] > 0. )
    if ( len( sel ) > 0 ):
      if ( last_n == - 1 ):
        scan_entry = [ n ]
        first_n = n
      elif ( time_array[ n, 0 ] - time_array[ last_n, 0 ] > 
          max_gap / ( 24. * 60. ) ):
        scan_entry.append( last_n )
        scan_list.append( scan_entry )
        scan_entry = [ n ]
        first_n = n
      elif ( time_array[ n, 0 ] - time_array[ first_n, 0 ] > 
          max_period / ( 24. * 60. ) ):
        scan_entry.append( n - 1 )
        scan_list.append( scan_entry )
        scan_entry = [ n ]
        first_n = n
      last_n = n
  scan_entry.append( last_n )
  scan_list.append( scan_entry )
  
  # merge short scans with neighbouring ones
  for s in range( len( scan_list ) ):
    scan = scan_list[ s ]
    if ( time_array[ scan[ 1 ], 0 ] - time_array[ scan[ 0 ], 0 ] < 
        min_period / ( 24. * 60. ) ):
      delta_low = 1.e7
      delta_high = 1.e7
      if ( s > 0 ):
        scan_low = scan_list[ s - 1 ]
        if ( scan_low == [ - 1.e7, - 1.e7 ] ):
          delta_low = 1.e7
        else:
          delta_low = time_array[ scan[ 0 ], 0 ] - time_array[ scan_low[ 1 ], 0 ]
      if ( s < len( scan_list ) - 1 ):
        scan_high = scan_list[ s + 1 ]
        if ( scan_high == [ - 1.e7, - 1.e7 ] ):
          delta_high = 1.e7
        else:
          delta_high = time_array[ scan_high[ 0 ], 0 ] - time_array[ scan[ 1 ], 0 ]
      if ( ( delta_low < 1.e5 ) or ( delta_high < 1.e5 ) ):
        if ( delta_low < delta_high ):
          if ( delta_low < max_gap  / ( 24. * 60. ) ):
            scan_list[ s - 1 ][ 1 ] = scan[ 1 ]
          scan_list[ s ] = [ - 1.e7, - 1.e7 ]
        else:
          if ( delta_high < max_gap  / ( 24. * 60. ) ):
            scan_list[ s + 1 ][ 0 ] = scan[ 0 ]
          scan_list[ s ] = [ - 1.e7, - 1.e7 ]
  while ( [ - 1.e7, - 1.e7 ] in scan_list ):
    scan_list.remove( [ - 1.e7, - 1.e7 ] )
  
  # assign time stamps to scans
  scan_count = len( scan_list )
  scan_sel = []
  for t in range( time_count ):
    scan_id = - 1
    for scan in scan_list:
      if ( t in range( scan[ 0 ], scan[ 1 ] + 1 ) ):
        scan_id = scan_list.index( scan )
        break
    scan_sel.append( scan_id )
  scan_sel = array( scan_sel, dtype = int32 ).reshape( len( scan_sel ), 1 )
  
  # remove entries within short, isolated scans
  time_list = []
  for scan in scan_list:
    time_list = time_list + range( scan[ 0 ], scan[ 1 ] + 1 )
  
  # unwrap and average phases
  uc_ph_array = azeros( ph_array )
  for i in range( antenna_count ):
    for s in range( scan_count ):
      scan = scan_list[ s ]
      sel = awhere( er_array[ scan[ 0 ] : scan[ 1 ] + 1, i ] > 0. )
      if ( len( sel ) > 0 ):
        p1 = aget( ph_array[ scan[ 0 ] : scan[ 1 ] + 1, i ], sel )
        if ( len( p1 ) < 20 ):
          p2 = adegrees( unwrap( aradians( p1 ) ) )
        else:
          p2 = aunwrap_phase( p1, alpha = 0.01 )
          if sometrue( isnan( p2 ) ):
            p2 = aunwrap_phase( p1, alpha = 0.001 )
          if ( sometrue( isnan( p2 ) ) or ( max( fabs( p2 ) ) > 1.e4 ) ):
            p2 = adegrees( unwrap( aradians( p1 ) ) )
        p2 = p2 - p2.mean()
        uc_ph_array[ scan[ 0 ] : scan[ 1 ] + 1, i ] = aput(
            uc_ph_array[ scan[ 0 ] : scan[ 1 ] + 1, i ], sel, p2 )
  
  # solve per scan
  results = []
  if print_info:
    print 'found %s scans' % ( repr( scan_count ) )
  for s in range( scan_count ):
    X_table = []
    za_table = []
    scan = scan_list[ s ]
    slen = scan[ 1 ] - scan[ 0 ] + 1
    if print_info:
      print '... processing scan %s' % ( repr( s + 1 ) )
    
    # get pierce point coordinates at center of scan
    stime_array = time_array[ scan[ 0 ] : scan[ 1 ] + 1, 0 ]
    center_time = ( stime_array[ 0 ] + stime_array[ -1 ] ) / 2.
    dtime_array = abs( stime_array - center_time )
    center_index = awhere( dtime_array == dtime_array.min() )[ 0, 0 ] + scan[ 0 ]
    pierce_table = calculate_pierce_coordinates( time_array[ center_index ],
        center_table, [ center_table ], array_table, antenna_table,
        height = height, iterations = iterations )
    
    # loop over pierce points
    for pierce_info in pierce_table:
      [ X, za, [ k, i ] ] = pierce_info
      # store model fit input data
      X_table.append( X )
      za_table.append( za )
    results.append( [ time_array[ center_index ], X_table, za_table ] )
    
    # loop over all possible lags
    time_offset = -dtime
    while ( time_offset < max_shift / ( 60. * 24. ) ):
      time_offset = time_offset + dtime
      stime_array_p = stime_array + ( time_offset / 2. )
      stime_array_m = stime_array - ( time_offset / 2. )
      dstime_array = resize( stime_array_m, ( slen, slen ) ) - transpose( 
          resize( stime_array_p, ( slen, slen ) ) )
      sel_pm = awhere( abs( dstime_array ) < time_tol )
      if ( len( sel_pm ) < min_points ):
        continue
      sel_p = sel_pm[ : , 0 : 1 ]
      sel_m = sel_pm[ : , 1 : 2 ]
      time_offset = time_offset + median( aget( dstime_array, sel_pm ) )
      
      # loop over all possible antenna combinations
      for a in ants:
        for b in ants:
          if ( b <= a ):
            continue
          for c in ants:
            for d in ants:
#              if ( ( d <= c ) or ( [ a, b ] == [ c, d ] ) ):
              if ( d <= c ):
                continue
              
              # check available phases
              ph_array_p = aget( uc_ph_array[ scan[ 0 ] : scan[ 1 ] + 1 ], sel_p )
              ph_array_m = aget( uc_ph_array[ scan[ 0 ] : scan[ 1 ] + 1 ], sel_m )
              er_array_p = aget( er_array[ scan[ 0 ] : scan[ 1 ] + 1 ], sel_p )
              er_array_m = aget( er_array[ scan[ 0 ] : scan[ 1 ] + 1 ], sel_m )
              sel_pm = awhere( ( er_array_p[ : , a - 1 ] > 0. ) & 
                  ( er_array_p[ : , b - 1 ] > 0. ) &
                  ( er_array_m[ : , c - 1 ] > 0. ) & 
                  ( er_array_m[ : , d - 1 ] > 0. ) )
              sel_mp = awhere( ( er_array_m[ : , a - 1 ] > 0. ) & 
                  ( er_array_m[ : , b - 1 ] > 0. ) &
                  ( er_array_p[ : , c - 1 ] > 0. ) & 
                  ( er_array_p[ : , d - 1 ] > 0. ) )
              
              # first process positive lag for ab and negative lag for cd
              if ( len( sel_pm ) > min_points ):
                ph_array_pab = ( aget( ph_array_p[ : , a - 1 ], sel_pm ) - 
                    aget( ph_array_p[ : , b - 1 ], sel_pm ) )
                ph_array_mcd = ( aget( ph_array_m[ : , c - 1 ], sel_pm ) - 
                    aget( ph_array_m[ : , d - 1 ], sel_pm ) )
                if remove_mean:
                  ph_array_pab = ph_array_pab - ph_array_pab.mean()
                  ph_array_mcd = ph_array_mcd - ph_array_mcd.mean()
                corr_pm = ( ph_array_pab * ph_array_mcd ).mean()
                results.append( [ a, b, c, d, time_offset, corr_pm ] )
              
              # then process positive lag for cd and negative lag for ab
              if ( ( len( sel_mp ) > min_points ) and ( time_offset > time_tol ) ):
                ph_array_mab = ( aget( ph_array_m[ : , a - 1 ], sel_mp ) - 
                    aget( ph_array_m[ : , b - 1 ], sel_mp ) )
                ph_array_pcd = ( aget( ph_array_p[ : , c - 1 ], sel_mp ) - 
                    aget( ph_array_p[ : , d - 1 ], sel_mp ) )
                if remove_mean:
                  ph_array_mab = ph_array_mab - ph_array_mab.mean()
                  ph_array_pcd = ph_array_pcd - ph_array_pcd.mean()
                corr_mp = ( ph_array_mab * ph_array_pcd ).mean()
                results.append( [ a, b, c, d, -time_offset, corr_mp ] )
  
  return results

###############################################################################

def correlate_antenna_phases( uv, in_version = 0, max_period = 30.,
    min_period = 20., max_gap = 1., max_shift = 5., remove_mean = True,
    time_step = 1, height = 300.e3, iterations = 1, min_points = 10,
    print_info = True, tolerance = 0.1, antennas = [] ):
# max_gap, max_period & min_period in minutes
  
  normalize_weights = True
  debug = False
  double_precision = True
  if double_precision:
    dtype = float64
  else:
    dtype = float32
  
  calibration_data = read_phase_calibration_data( uv, solution_version = in_version,
      print_info = print_info )
  time_array = array( calibration_data[ 0 ], dtype = dtype )
  center_table = calibration_data[ 1 ]
  array_table = calibration_data[ 2 ]
  antenna_table = calibration_data[ 3 ]
  reference_array = calibration_data[ 4 ]
  ph_array = array( calibration_data[ 5 ], dtype = dtype )
  er_array = array( calibration_data[ 6 ], dtype = dtype )
  time_count = len( time_array )
  antenna_count = len( antenna_table )
  if ( len( antennas ) > 0 ):
    ants = [ i for i in antennas ]
    ants.sort( cmp = lambda a, b: cmp( a, b ) )
  else:
    ants = range( 1, 1 + antenna_count )
  
  # calculate some time variables
  dtime = median( time_array[ 1 : , 0 ] - time_array[ : -1, 0 ] )
  time_tol = tolerance * dtime
  dtime = dtime * float( time_step )
  
  # identify scans or large gaps
  scan_list = []
  last_n = - 1
  first_n = 0
  for n in range( time_count ):
    sel = awhere( er_array[ n, : ] > 0. )
    if ( len( sel ) > 0 ):
      if ( last_n == - 1 ):
        scan_entry = [ n ]
        first_n = n
      elif ( time_array[ n, 0 ] - time_array[ last_n, 0 ] > 
          max_gap / ( 24. * 60. ) ):
        scan_entry.append( last_n )
        scan_list.append( scan_entry )
        scan_entry = [ n ]
        first_n = n
      elif ( time_array[ n, 0 ] - time_array[ first_n, 0 ] > 
          max_period / ( 24. * 60. ) ):
        scan_entry.append( n - 1 )
        scan_list.append( scan_entry )
        scan_entry = [ n ]
        first_n = n
      last_n = n
  scan_entry.append( last_n )
  scan_list.append( scan_entry )
  
  # merge short scans with neighbouring ones
  for s in range( len( scan_list ) ):
    scan = scan_list[ s ]
    if ( time_array[ scan[ 1 ], 0 ] - time_array[ scan[ 0 ], 0 ] < 
        min_period / ( 24. * 60. ) ):
      delta_low = 1.e7
      delta_high = 1.e7
      if ( s > 0 ):
        scan_low = scan_list[ s - 1 ]
        if ( scan_low == [ - 1.e7, - 1.e7 ] ):
          delta_low = 1.e7
        else:
          delta_low = time_array[ scan[ 0 ], 0 ] - time_array[ scan_low[ 1 ], 0 ]
      if ( s < len( scan_list ) - 1 ):
        scan_high = scan_list[ s + 1 ]
        if ( scan_high == [ - 1.e7, - 1.e7 ] ):
          delta_high = 1.e7
        else:
          delta_high = time_array[ scan_high[ 0 ], 0 ] - time_array[ scan[ 1 ], 0 ]
      if ( ( delta_low < 1.e5 ) or ( delta_high < 1.e5 ) ):
        if ( delta_low < delta_high ):
          if ( delta_low < max_gap  / ( 24. * 60. ) ):
            scan_list[ s - 1 ][ 1 ] = scan[ 1 ]
          scan_list[ s ] = [ - 1.e7, - 1.e7 ]
        else:
          if ( delta_high < max_gap  / ( 24. * 60. ) ):
            scan_list[ s + 1 ][ 0 ] = scan[ 0 ]
          scan_list[ s ] = [ - 1.e7, - 1.e7 ]
  while ( [ - 1.e7, - 1.e7 ] in scan_list ):
    scan_list.remove( [ - 1.e7, - 1.e7 ] )
  
  # assign time stamps to scans
  scan_count = len( scan_list )
  scan_sel = []
  for t in range( time_count ):
    scan_id = - 1
    for scan in scan_list:
      if ( t in range( scan[ 0 ], scan[ 1 ] + 1 ) ):
        scan_id = scan_list.index( scan )
        break
    scan_sel.append( scan_id )
  scan_sel = array( scan_sel, dtype = int32 ).reshape( len( scan_sel ), 1 )
  
  # remove entries within short, isolated scans
  time_list = []
  for scan in scan_list:
    time_list = time_list + range( scan[ 0 ], scan[ 1 ] + 1 )
  
  # unwrap and average phases
  uc_ph_array = azeros( ph_array )
  antenna_active = zeros( ( scan_count, antenna_count ), dtype = bool )
  for i in range( antenna_count ):
    for s in range( scan_count ):
      scan = scan_list[ s ]
      sel = awhere( er_array[ scan[ 0 ] : scan[ 1 ] + 1, i ] > 0. )
      if ( len( sel ) > 0 ):
        antenna_active[ s, i ] = True
        p1 = aget( ph_array[ scan[ 0 ] : scan[ 1 ] + 1, i ], sel )
        if ( len( p1 ) < 20 ):
          p2 = adegrees( unwrap( aradians( p1 ) ) )
        else:
          p2 = aunwrap_phase( p1, alpha = 0.01 )
          if sometrue( isnan( p2 ) ):
            p2 = aunwrap_phase( p1, alpha = 0.001 )
          if ( sometrue( isnan( p2 ) ) or ( max( fabs( p2 ) ) > 1.e4 ) ):
            p2 = adegrees( unwrap( aradians( p1 ) ) )
        p2 = p2 - p2.mean()
        uc_ph_array[ scan[ 0 ] : scan[ 1 ] + 1, i ] = aput(
            uc_ph_array[ scan[ 0 ] : scan[ 1 ] + 1, i ], sel, p2 )
  
  # solve per scan
  results = []
  if print_info:
    print 'found %s scans' % ( repr( scan_count ) )
  for s in range( scan_count ):
#    X_table = []
#    za_table = []
    X_table = zeros( ( antenna_count, 2 ), dtype = float64 )
    za_table = zeros( ( antenna_count ), dtype = float64 )
    scan = scan_list[ s ]
    slen = scan[ 1 ] - scan[ 0 ] + 1
    if print_info:
      print '... processing scan %s' % ( repr( s + 1 ) )
    
    # get pierce point coordinates at center of scan
    stime_array = time_array[ scan[ 0 ] : scan[ 1 ] + 1, 0 ]
    center_time = ( stime_array[ 0 ] + stime_array[ -1 ] ) / 2.
    dtime_array = abs( stime_array - center_time )
    center_index = awhere( dtime_array == dtime_array.min() )[ 0, 0 ] + scan[ 0 ]
    pierce_table = calculate_pierce_coordinates( time_array[ center_index ],
        center_table, [ center_table ], array_table, antenna_table,
        height = height, iterations = iterations )
    
    # loop over pierce points
    for pierce_info in pierce_table:
      [ X, za, [ k, i ] ] = pierce_info
      # store model fit input data
      X_table[ i ] = array( X, dtype = float64 )
      za_table[ i ] = za
    results.append( [ time_array[ center_index ], X_table, za_table ] )
    
    # loop over all possible lags
    time_offset = -dtime
    while ( time_offset < max_shift / ( 60. * 24. ) ):
      time_offset = time_offset + dtime
      stime_array_p = stime_array + ( time_offset / 2. )
      stime_array_m = stime_array - ( time_offset / 2. )
      dstime_array = resize( stime_array_m, ( slen, slen ) ) - transpose( 
          resize( stime_array_p, ( slen, slen ) ) )
      sel_pm = awhere( abs( dstime_array ) < time_tol )
      if ( len( sel_pm ) < min_points ):
        continue
      sel_p = sel_pm[ : , 0 : 1 ]
      sel_m = sel_pm[ : , 1 : 2 ]
      time_offset = time_offset + median( aget( dstime_array, sel_pm ) )
      
      # loop over all possible antenna combinations
      for a in ants:
        if ( not antenna_active[ s, a - 1 ] ):
          continue
        for b in ants:
          if ( ( not antenna_active[ s, b - 1 ] ) or ( b <= a ) ):
            continue
          if True:
              c = a
              d = b
              
              # check available phases
              ph_array_p = aget( uc_ph_array[ scan[ 0 ] : scan[ 1 ] + 1 ], sel_p )
              ph_array_m = aget( uc_ph_array[ scan[ 0 ] : scan[ 1 ] + 1 ], sel_m )
              er_array_p = aget( er_array[ scan[ 0 ] : scan[ 1 ] + 1 ], sel_p )
              er_array_m = aget( er_array[ scan[ 0 ] : scan[ 1 ] + 1 ], sel_m )
              sel_pm = awhere( ( er_array_p[ : , a - 1 ] > 0. ) & 
                  ( er_array_p[ : , b - 1 ] > 0. ) &
                  ( er_array_m[ : , c - 1 ] > 0. ) & 
                  ( er_array_m[ : , d - 1 ] > 0. ) )
              sel_mp = awhere( ( er_array_m[ : , a - 1 ] > 0. ) & 
                  ( er_array_m[ : , b - 1 ] > 0. ) &
                  ( er_array_p[ : , c - 1 ] > 0. ) & 
                  ( er_array_p[ : , d - 1 ] > 0. ) )
              
              # first process positive lag for ab and negative lag for cd
              if ( len( sel_pm ) > min_points ):
                ph_array_pab = ( aget( ph_array_p[ : , a - 1 ], sel_pm ) - 
                    aget( ph_array_p[ : , b - 1 ], sel_pm ) )
                ph_array_mcd = ( aget( ph_array_m[ : , c - 1 ], sel_pm ) - 
                    aget( ph_array_m[ : , d - 1 ], sel_pm ) )
                if remove_mean:
                  ph_array_pab = ph_array_pab - ph_array_pab.mean()
                  ph_array_mcd = ph_array_mcd - ph_array_mcd.mean()
                corr_pm = ( ( ph_array_pab - ph_array_mcd )**2 ).mean()
                results.append( [ a, b, c, d, time_offset, corr_pm ] )
              
              # then process positive lag for cd and negative lag for ab
              if ( ( len( sel_mp ) > min_points ) and ( time_offset > time_tol ) ):
                if ( [ c, d ] == [ a, b ] ):
                  corr_mp = corr_pm
                else:
                  ph_array_mab = ( aget( ph_array_m[ : , a - 1 ], sel_mp ) - 
                      aget( ph_array_m[ : , b - 1 ], sel_mp ) )
                  ph_array_pcd = ( aget( ph_array_p[ : , c - 1 ], sel_mp ) - 
                      aget( ph_array_p[ : , d - 1 ], sel_mp ) )
                  if remove_mean:
                    ph_array_mab = ph_array_mab - ph_array_mab.mean()
                    ph_array_pcd = ph_array_pcd - ph_array_pcd.mean()
                  corr_mp = ( ( ph_array_mab - ph_array_pcd )**2 ).mean()
                results.append( [ a, b, c, d, -time_offset, corr_mp ] )
  
  return results

###############################################################################

def correlate_antenna_phases2( uv, in_version = 0, max_period = 30.,
    min_period = 20., max_gap = 1., max_shift = 5., remove_mean = True,
    time_step = 1, height = 300.e3, iterations = 1, min_points = 10,
    print_info = True, tolerance = 0.1, antennas = [], baselines = [] ):
# max_gap, max_period & min_period in minutes
  
  normalize_weights = True
  debug = False
  double_precision = True
  if double_precision:
    dtype = float64
  else:
    dtype = float32
  
  calibration_data = read_phase_calibration_data( uv, solution_version = in_version,
      print_info = print_info )
  time_array = array( calibration_data[ 0 ], dtype = dtype )
  center_table = calibration_data[ 1 ]
  array_table = calibration_data[ 2 ]
  antenna_table = calibration_data[ 3 ]
  reference_array = calibration_data[ 4 ]
  ph_array = array( calibration_data[ 5 ], dtype = dtype )
  er_array = array( calibration_data[ 6 ], dtype = dtype )
  time_count = len( time_array )
  antenna_count = len( antenna_table )
  if ( len( antennas ) > 0 ):
    ants = [ i for i in antennas ]
    ants.sort( cmp = lambda a, b: cmp( a, b ) )
  else:
    ants = range( 1, 1 + antenna_count )
  if ( len( baselines ) > 0 ):
    bass = [ i for i in baselines ]
    bass.sort( cmp = lambda a, b: cmp( a, b ) )
  else:
    bass = range( 1, 1 + antenna_count )
  
  # calculate some time variables
  dtime = median( time_array[ 1 : , 0 ] - time_array[ : -1, 0 ] )
  time_tol = tolerance * dtime
  dtime = dtime * float( time_step )
  
  # identify scans or large gaps
  scan_list = []
  last_n = - 1
  first_n = 0
  for n in range( time_count ):
    sel = awhere( er_array[ n, : ] > 0. )
    if ( len( sel ) > 0 ):
      if ( last_n == - 1 ):
        scan_entry = [ n ]
        first_n = n
      elif ( time_array[ n, 0 ] - time_array[ last_n, 0 ] > 
          max_gap / ( 24. * 60. ) ):
        scan_entry.append( last_n )
        scan_list.append( scan_entry )
        scan_entry = [ n ]
        first_n = n
      elif ( time_array[ n, 0 ] - time_array[ first_n, 0 ] > 
          max_period / ( 24. * 60. ) ):
        scan_entry.append( n - 1 )
        scan_list.append( scan_entry )
        scan_entry = [ n ]
        first_n = n
      last_n = n
  scan_entry.append( last_n )
  scan_list.append( scan_entry )
  
  # merge short scans with neighbouring ones
  for s in range( len( scan_list ) ):
    scan = scan_list[ s ]
    if ( time_array[ scan[ 1 ], 0 ] - time_array[ scan[ 0 ], 0 ] < 
        min_period / ( 24. * 60. ) ):
      delta_low = 1.e7
      delta_high = 1.e7
      if ( s > 0 ):
        scan_low = scan_list[ s - 1 ]
        if ( scan_low == [ - 1.e7, - 1.e7 ] ):
          delta_low = 1.e7
        else:
          delta_low = time_array[ scan[ 0 ], 0 ] - time_array[ scan_low[ 1 ], 0 ]
      if ( s < len( scan_list ) - 1 ):
        scan_high = scan_list[ s + 1 ]
        if ( scan_high == [ - 1.e7, - 1.e7 ] ):
          delta_high = 1.e7
        else:
          delta_high = time_array[ scan_high[ 0 ], 0 ] - time_array[ scan[ 1 ], 0 ]
      if ( ( delta_low < 1.e5 ) or ( delta_high < 1.e5 ) ):
        if ( delta_low < delta_high ):
          if ( delta_low < max_gap  / ( 24. * 60. ) ):
            scan_list[ s - 1 ][ 1 ] = scan[ 1 ]
          scan_list[ s ] = [ - 1.e7, - 1.e7 ]
        else:
          if ( delta_high < max_gap  / ( 24. * 60. ) ):
            scan_list[ s + 1 ][ 0 ] = scan[ 0 ]
          scan_list[ s ] = [ - 1.e7, - 1.e7 ]
  while ( [ - 1.e7, - 1.e7 ] in scan_list ):
    scan_list.remove( [ - 1.e7, - 1.e7 ] )
  
  # assign time stamps to scans
  scan_count = len( scan_list )
  scan_sel = []
  for t in range( time_count ):
    scan_id = - 1
    for scan in scan_list:
      if ( t in range( scan[ 0 ], scan[ 1 ] + 1 ) ):
        scan_id = scan_list.index( scan )
        break
    scan_sel.append( scan_id )
  scan_sel = array( scan_sel, dtype = int32 ).reshape( len( scan_sel ), 1 )
  
  # remove entries within short, isolated scans
  time_list = []
  for scan in scan_list:
    time_list = time_list + range( scan[ 0 ], scan[ 1 ] + 1 )
  
  # unwrap and average phases
  uc_ph_array = azeros( ph_array )
  antenna_active = zeros( ( scan_count, antenna_count ), dtype = bool )
  for i in range( antenna_count ):
    for s in range( scan_count ):
      scan = scan_list[ s ]
      sel = awhere( er_array[ scan[ 0 ] : scan[ 1 ] + 1, i ] > 0. )
      if ( len( sel ) > 0 ):
        antenna_active[ s, i ] = True
        p1 = aget( ph_array[ scan[ 0 ] : scan[ 1 ] + 1, i ], sel )
        if ( len( p1 ) < 20 ):
          p2 = adegrees( unwrap( aradians( p1 ) ) )
        else:
          p2 = aunwrap_phase( p1, alpha = 0.01 )
          if sometrue( isnan( p2 ) ):
            p2 = aunwrap_phase( p1, alpha = 0.001 )
          if ( sometrue( isnan( p2 ) ) or ( max( fabs( p2 ) ) > 1.e4 ) ):
            p2 = adegrees( unwrap( aradians( p1 ) ) )
        p2 = p2 - p2.mean()
        uc_ph_array[ scan[ 0 ] : scan[ 1 ] + 1, i ] = aput(
            uc_ph_array[ scan[ 0 ] : scan[ 1 ] + 1, i ], sel, p2 )
  
  # solve per scan
  results = []
  if print_info:
    print 'found %s scans' % ( repr( scan_count ) )
  for s in range( scan_count ):
    X_table = []
    za_table = []
    scan = scan_list[ s ]
    slen = scan[ 1 ] - scan[ 0 ] + 1
    if print_info:
      print '... processing scan %s' % ( repr( s + 1 ) )
    
    # get pierce point coordinates at center of scan
    stime_array = time_array[ scan[ 0 ] : scan[ 1 ] + 1, 0 ]
    center_time = ( stime_array[ 0 ] + stime_array[ -1 ] ) / 2.
    dtime_array = abs( stime_array - center_time )
    center_index = awhere( dtime_array == dtime_array.min() )[ 0, 0 ] + scan[ 0 ]
    pierce_table = calculate_pierce_coordinates( time_array[ center_index ],
        center_table, [ center_table ], array_table, antenna_table,
        height = height, iterations = iterations )
    
    # loop over pierce points
    for pierce_info in pierce_table:
      [ X, za, [ k, i ] ] = pierce_info
      # store model fit input data
      X_table.append( X )
      za_table.append( za )
    results.append( [ time_array[ center_index ], X_table, za_table ] )
    
    # loop over all possible lags
    time_offset = -dtime
    while ( time_offset < max_shift / ( 60. * 24. ) ):
      time_offset = time_offset + dtime
      stime_array_p = stime_array + ( time_offset / 2. )
      stime_array_m = stime_array - ( time_offset / 2. )
      dstime_array = resize( stime_array_m, ( slen, slen ) ) - transpose( 
          resize( stime_array_p, ( slen, slen ) ) )
      sel_pm = awhere( abs( dstime_array ) < time_tol )
      if ( len( sel_pm ) < min_points ):
        continue
      sel_p = sel_pm[ : , 0 : 1 ]
      sel_m = sel_pm[ : , 1 : 2 ]
      time_offset = time_offset + median( aget( dstime_array, sel_pm ) )
      
      # loop over all possible antenna combinations
      for a in ants:
        if ( not antenna_active[ s, a - 1 ] ):
          continue
        for b in ants:
          if ( ( not antenna_active[ s, b - 1 ] ) or ( b <= a ) ):
            continue
          for c in bass:
            if ( ( not antenna_active[ s, c - 1 ] ) or ( c in [ a, b ] ) ):
              continue
            d = c
            if True:
              
              # check available phases
              ph_array_p = aget( uc_ph_array[ scan[ 0 ] : scan[ 1 ] + 1 ], sel_p )
              ph_array_m = aget( uc_ph_array[ scan[ 0 ] : scan[ 1 ] + 1 ], sel_m )
              er_array_p = aget( er_array[ scan[ 0 ] : scan[ 1 ] + 1 ], sel_p )
              er_array_m = aget( er_array[ scan[ 0 ] : scan[ 1 ] + 1 ], sel_m )
              sel_pm = awhere( ( er_array_p[ : , a - 1 ] > 0. ) & 
                  ( er_array_p[ : , b - 1 ] > 0. ) &
                  ( er_array_m[ : , c - 1 ] > 0. ) & 
                  ( er_array_m[ : , d - 1 ] > 0. ) )
              sel_mp = awhere( ( er_array_m[ : , a - 1 ] > 0. ) & 
                  ( er_array_m[ : , b - 1 ] > 0. ) &
                  ( er_array_p[ : , c - 1 ] > 0. ) & 
                  ( er_array_p[ : , d - 1 ] > 0. ) )
              
              # first process positive lag for ab and negative lag for cd
              if ( len( sel_pm ) > min_points ):

                ph_array_pab = ( aget( ph_array_p[ : , a - 1 ], sel_pm ) - 
                    aget( ph_array_p[ : , b - 1 ], sel_pm ) )
                ph_array_mcd = ( aget( ph_array_m[ : , c - 1 ], sel_pm ) - 
                    aget( ph_array_m[ : , d - 1 ], sel_pm ) )
                if remove_mean:
                  ph_array_pab = ph_array_pab - ph_array_pab.mean()
                  ph_array_mcd = ph_array_mcd - ph_array_mcd.mean()
                corr_pm = ( ( ph_array_pab - ph_array_mcd )**2 ).mean()
                results.append( [ a, b, c, d, time_offset, corr_pm ] )
              
              # then process positive lag for cd and negative lag for ab
              if ( ( len( sel_mp ) > min_points ) and ( time_offset > time_tol ) ):
                ph_array_mab = ( aget( ph_array_m[ : , a - 1 ], sel_mp ) - 
                    aget( ph_array_m[ : , b - 1 ], sel_mp ) )
                ph_array_pcd = ( aget( ph_array_p[ : , c - 1 ], sel_mp ) - 
                    aget( ph_array_p[ : , d - 1 ], sel_mp ) )
                if remove_mean:
                  ph_array_mab = ph_array_mab - ph_array_mab.mean()
                  ph_array_pcd = ph_array_pcd - ph_array_pcd.mean()
                corr_mp = ( ( ph_array_mab - ph_array_pcd )**2 ).mean()
                results.append( [ a, b, c, d, -time_offset, corr_mp ] )
  
  return results

###############################################################################

#def bandpass_filter_phases( phases, integration_time = 10., 
#    time_window = [ 1., 10. ], order = 5 ):
#
#  wn = [ integration_time / ( 30. * time_window[ 1 ] ),
#         integration_time / ( 30. * time_window[ 0 ] ) ]
#  [ b, a ] = ss.butter( order, wn, btype = 'bandpass', output = 'ba' )
#  samples = len( phases )
#  ph = array( phases )
#  x = zeros( ( samples + 2 * order ), dtype = ph.dtype )
#  x[ order : -order ] = array( ph )
#  y = ss.lfilter( b, a, x )

###############################################################################

def correlate_baseline_phases_new( uv, in_version = 0, max_period = 30.,
    min_period = 20., max_gap = 1., max_shift = 5., remove_mean = False,
    time_step = 1, height = 300.e3, iterations = 1, min_points = 10,
    print_info = True, tolerance = 0.25, antennas = [], span = 2,
    blength_range = [ 1.e3, 1.e6 ], bangle_range = [ 0., 60. ],
    clength_range = [], cangle_range = [],
    rratio_range = [ 0., 0.75 ], rangle_range = [ 0., 45. ] ):
# max_gap, max_period & min_period in minutes
# blength_range for individual baselines
# bangle_range for relative angle between baselines
# clength_range for distance between baselines
# cangle_range for angle of line between baselines
# rratio_range for relative length of baseline to line between baselines
# rangle_range for relative angle between baseline and line between baselines
  
  normalize_weights = True
  debug = False
  double_precision = True
  if double_precision:
    dtype = float64
  else:
    dtype = float32
  
  calibration_data = read_phase_calibration_data( uv, solution_version = in_version,
      print_info = print_info )
  time_array = array( calibration_data[ 0 ], dtype = dtype )
  center_table = calibration_data[ 1 ]
  array_table = calibration_data[ 2 ]
  antenna_table = calibration_data[ 3 ]
  reference_array = calibration_data[ 4 ]
  ph_array = array( calibration_data[ 5 ], dtype = dtype )
  er_array = array( calibration_data[ 6 ], dtype = dtype )
  time_count = len( time_array )
  antenna_count = len( antenna_table )
  if ( len( antennas ) > 0 ):
    ants = [ i for i in antennas ]
    ants.sort( cmp = lambda a, b: cmp( a, b ) )
  else:
    ants = range( 1, 1 + antenna_count )
  
  # calculate some time variables
  dtime = median( time_array[ 1 : , 0 ] - time_array[ : -1, 0 ] )
#  dtime = restore_parameter( uv, 'integration_time' ) / 86400.
  time_tol = tolerance * dtime
  dtime = dtime * float( time_step )
  
  # identify scans or large gaps
  scan_list = []
  last_n = - 1
  first_n = 0
  for n in range( time_count ):
    sel = awhere( er_array[ n, : ] > 0. )
    if ( len( sel ) > 0 ):
      if ( last_n == - 1 ):
        scan_entry = [ n ]
        first_n = n
      elif ( time_array[ n, 0 ] - time_array[ last_n, 0 ] > 
          max_gap / ( 24. * 60. ) ):
        scan_entry.append( last_n )
        scan_list.append( scan_entry )
        scan_entry = [ n ]
        first_n = n
      elif ( time_array[ n, 0 ] - time_array[ first_n, 0 ] > 
          max_period / ( 24. * 60. ) ):
        scan_entry.append( n - 1 )
        scan_list.append( scan_entry )
        scan_entry = [ n ]
        first_n = n
      last_n = n
  scan_entry.append( last_n )
  scan_list.append( scan_entry )
  
  # merge short scans with neighbouring ones
  for s in range( len( scan_list ) ):
    scan = scan_list[ s ]
    if ( time_array[ scan[ 1 ], 0 ] - time_array[ scan[ 0 ], 0 ] < 
        min_period / ( 24. * 60. ) ):
      delta_low = 1.e7
      delta_high = 1.e7
      if ( s > 0 ):
        scan_low = scan_list[ s - 1 ]
        if ( scan_low == [ - 1.e7, - 1.e7 ] ):
          delta_low = 1.e7
        else:
          delta_low = time_array[ scan[ 0 ], 0 ] - time_array[ scan_low[ 1 ], 0 ]
      if ( s < len( scan_list ) - 1 ):
        scan_high = scan_list[ s + 1 ]
        if ( scan_high == [ - 1.e7, - 1.e7 ] ):
          delta_high = 1.e7
        else:
          delta_high = time_array[ scan_high[ 0 ], 0 ] - time_array[ scan[ 1 ], 0 ]
      if ( ( delta_low < 1.e5 ) or ( delta_high < 1.e5 ) ):
        if ( delta_low < delta_high ):
          if ( delta_low < max_gap  / ( 24. * 60. ) ):
            scan_list[ s - 1 ][ 1 ] = scan[ 1 ]
          scan_list[ s ] = [ - 1.e7, - 1.e7 ]
        else:
          if ( delta_high < max_gap  / ( 24. * 60. ) ):
            scan_list[ s + 1 ][ 0 ] = scan[ 0 ]
          scan_list[ s ] = [ - 1.e7, - 1.e7 ]
  while ( [ - 1.e7, - 1.e7 ] in scan_list ):
    scan_list.remove( [ - 1.e7, - 1.e7 ] )
  
  # assign time stamps to scans
  scan_count = len( scan_list )
  scan_sel = []
  for t in range( time_count ):
    scan_id = - 1
    for scan in scan_list:
      if ( t in range( scan[ 0 ], scan[ 1 ] + 1 ) ):
        scan_id = scan_list.index( scan )
        break
    scan_sel.append( scan_id )
  scan_sel = array( scan_sel, dtype = int32 ).reshape( len( scan_sel ), 1 )
  
  # remove entries within short, isolated scans
  time_list = []
  for scan in scan_list:
    time_list = time_list + range( scan[ 0 ], scan[ 1 ] + 1 )
  if print_info:
    print 'found %s scans' % ( repr( scan_count ) )
  
  # unwrap and average phases
  uc_ph_array = azeros( ph_array )
  antenna_active = zeros( ( scan_count, antenna_count ), dtype = bool )
  for i in range( antenna_count ):
    for s in range( scan_count ):
      scan = scan_list[ s ]
      sel = awhere( er_array[ scan[ 0 ] : scan[ 1 ] + 1, i ] > 0. )
      if ( len( sel ) > 0 ):
        antenna_active[ s, i ] = True
        p1 = aget( ph_array[ scan[ 0 ] : scan[ 1 ] + 1, i ], sel )
        if ( len( p1 ) < 20 ):
          p2 = adegrees( unwrap( aradians( p1 ) ) )
        else:
          p2 = aunwrap_phase( p1, alpha = 0.01 )
          if sometrue( isnan( p2 ) ):
            p2 = aunwrap_phase( p1, alpha = 0.001 )
          if ( sometrue( isnan( p2 ) ) or ( max( fabs( p2 ) ) > 1.e4 ) ):
            p2 = adegrees( unwrap( aradians( p1 ) ) )
        p2 = p2 - p2.mean()
        uc_ph_array[ scan[ 0 ] : scan[ 1 ] + 1, i ] = aput(
            uc_ph_array[ scan[ 0 ] : scan[ 1 ] + 1, i ], sel, p2 )
  
  # interpolate small gaps
  if ( span > 0 ):
    for i in range( antenna_count ):
      for s in range( scan_count ):
        scan = scan_list[ s ]
        if antenna_active[ s, i ]:
          for n in range( scan[ 0 ], scan[ 1 ] + 1 ):
            if ( er_array[ n, i ] <= 0. ):
              sel = awhere( abs( time_array[ scan[ 0 ] : scan[ 1 ] + 1, 0 ] - 
                  time_array[ n, 0 ] ) < float( span ) * dtime + time_tol )
              if ( len( sel ) > 0 ):
                sel2 = awhere( 
                    aget( er_array[ scan[ 0 ] : scan[ 1 ] + 1, i ], sel ) > 0. )
                if ( len( sel2 ) > 0 ):
                  uc_ph_array[ n, i ] = aget( aget( uc_ph_array[ scan[ 0 ] : 
                      scan[ 1 ] + 1, i ], sel ), sel2 ).mean()
                  er_array[ n, i ] = aget( aget( er_array[ scan[ 0 ] : 
                      scan[ 1 ] + 1, i ], sel ), sel2 ).mean()
  
  # check baselines against input criteria
  results = []
  baselines = []
  if print_info:
    print 'selecting baselines'
  for a in ants:
    if ( not antenna_active[ s, a - 1 ] ):
      continue
    x_a = array( antenna_table[ a - 1 ][ 0 ] )
    for b in ants:
      if ( ( not antenna_active[ s, b - 1 ] ) or ( b <= a ) ):
        continue
      x_b = array( antenna_table[ b - 1 ][ 0 ] )
      l_ab = sqrt( add.reduce( ( x_a - x_b )**2 ) )
      if ( len( blength_range ) == 2 ):
        if ( ( l_ab < blength_range[ 0 ] ) or ( l_ab > blength_range[ 1 ] ) ):
          continue
      x_ab = ( x_a + x_b ) / 2.
      llh_ab = xyz_to_geo_llh( x_ab )
      for c in ants:
        if ( ( not antenna_active[ s, c - 1 ] ) or ( c <= a ) or ( c == b ) ):
          continue
        x_c = array( antenna_table[ c - 1 ][ 0 ] )
        for d in ants:
          if ( ( not antenna_active[ s, d - 1 ] ) or ( d <= c ) or ( d == b ) ):
            continue
          x_d = array( antenna_table[ d - 1 ][ 0 ] )
          l_cd = sqrt( add.reduce( ( x_c - x_d )**2 ) )
          if ( len( blength_range ) == 2 ):
            if ( ( l_cd < blength_range[ 0 ] ) or ( l_cd > blength_range[ 1 ] ) ):
              continue
          x_cd = ( x_c + x_d ) / 2.
          llh_cd = xyz_to_geo_llh( x_cd )
          cos_ra_abcd = dot( ( x_a - x_b ), ( x_c - x_d ) ) / ( l_ab * l_cd )
          ra_abcd = degrees( acos( max( min( cos_ra_abcd, 1. ), -1. ) ) )
          if ( ra_abcd > 90. ):
            ara_abcd = 180. - ra_abcd
          else:
            ara_abcd = ra_abcd
          if ( len( bangle_range ) == 2 ):
            if ( ( ara_abcd < bangle_range[ 0 ] ) or ( ara_abcd > bangle_range[ 1 ] ) ):
              continue
          l_abcd = sqrt( add.reduce( ( x_ab - x_cd )**2 ) )
          if ( len( clength_range ) == 2 ):
            if ( ( l_abcd < clength_range[ 0 ] ) or ( l_abcd > clength_range[ 1 ] ) ):
              continue
          if ( len( rratio_range ) == 2 ):
            if ( ( l_ab / l_abcd < rratio_range[ 0 ] ) or 
                ( l_ab / l_abcd > rratio_range[ 1 ] ) or
                ( l_cd / l_abcd < rratio_range[ 0 ] ) or 
                ( l_cd / l_abcd > rratio_range[ 1 ] ) ):
              continue
          [ r_abcd, a_abcd ] = calculate_angular_separation( llh_cd[ 0 : 2 ],
              llh_ab[ 0 : 2 ] )
          if ( len( cangle_range ) == 2 ):
            da_abcd0 = amodulo( ( a_abcd - cangle_range[ 0 ] ) + 180., 360. ) - 180.
            da_abcd1 = amodulo( ( a_abcd - cangle_range[ 1 ] ) + 180., 360. ) - 180.
            if ( ( da_abcd0 < 0. ) or ( da_abcd1 > 0. ) ):
              continue
          cos_ra_ababcd = dot( ( x_a - x_b ), ( x_ab - x_cd ) ) / ( l_ab * l_abcd )
          ra_ababcd = degrees( acos( max( min( cos_ra_ababcd, 1. ), -1. ) ) )
          if ( ra_ababcd > 90. ):
            ara_ababcd = 180. - ra_ababcd
          else:
            ara_ababcd = ra_ababcd
          if ( len( rangle_range ) == 2 ):
            if ( ( ara_ababcd < rangle_range[ 0 ] ) or ( ara_ababcd > rangle_range[ 1 ] ) ):
              continue
          cos_ra_cdabcd = dot( ( x_c - x_d ), ( x_ab - x_cd ) ) / ( l_cd * l_abcd )
          ra_cdabcd = degrees( acos( max( min( cos_ra_cdabcd, 1. ), -1. ) ) )
          if ( ra_cdabcd > 90. ):
            ara_cdabcd = 180. - ara_cdabcd
          else:
            ara_cdabcd = ra_cdabcd
          if ( len( rangle_range ) == 2 ):
            if ( ( ara_cdabcd < rangle_range[ 0 ] ) or ( ara_cdabcd > rangle_range[ 1 ] ) ):
              continue
          baselines.append( [ a, b, c, d ] )
          results.append( [ a, b, c, d, l_ab, l_cd, ra_abcd, l_abcd, a_abcd,
              ra_ababcd, ra_cdabcd ] )
  
  # solve per scan
  for s in range( scan_count ):
    X_table = zeros( ( antenna_count, 2 ), dtype = float64 )
    za_table = zeros( ( antenna_count ), dtype = float64 )
    scan = scan_list[ s ]
    slen = scan[ 1 ] - scan[ 0 ] + 1
    if print_info:
      print '... processing scan %s' % ( repr( s + 1 ) )
    
    # get pierce point coordinates at center of scan
    stime_array = time_array[ scan[ 0 ] : scan[ 1 ] + 1, 0 ]
    center_time = ( stime_array[ 0 ] + stime_array[ -1 ] ) / 2.
    dtime_array = abs( stime_array - center_time )
    center_index = awhere( dtime_array == dtime_array.min() )[ 0, 0 ] + scan[ 0 ]
    pierce_table = calculate_pierce_coordinates( time_array[ center_index ],
        center_table, [ center_table ], array_table, antenna_table,
        height = height, iterations = iterations )
    
    # loop over pierce points
    for pierce_info in pierce_table:
      [ X, za, [ k, i ] ] = pierce_info
      # store model fit input data
      X_table[ i ] = array( X, dtype = float64 )
      za_table[ i ] = za
    results.append( [ time_array[ center_index ], X_table, za_table ] )
    
    # loop over all possible lags
    time_offset = -dtime
    while ( time_offset < max_shift / ( 60. * 24. ) ):
      time_offset = time_offset + dtime
      if ( print_info ):
        print '...... time lag = %s' % ( repr( time_to_dhms( time_offset ) ) )
#      pdb.set_trace()
      stime_array_p = stime_array + ( time_offset / 2. )
      stime_array_m = stime_array - ( time_offset / 2. )
      dstime_array = resize( stime_array_m, ( slen, slen ) ) - transpose( 
          resize( stime_array_p, ( slen, slen ) ) )
      sel_pm = awhere( abs( dstime_array ) < time_tol )
      if ( len( sel_pm ) < min_points ):
        if ( print_info ):
          print '......... skipping lag, too little data'
        continue
      sel_p = sel_pm[ : , 0 : 1 ]
      sel_m = sel_pm[ : , 1 : 2 ]
#      time_offset = time_offset + median( aget( dstime_array, sel_pm ) )
#      sel_dstime_array = aget( dstime_array, sel_pm )
#      sel_min = awhere( abs( sel_dstime_array ) == abs( sel_dstime_array ).min() )[ 0 : 1 ]
#      pdb.set_trace()
#      time_offset = time_offset + aget( sel_dstime_array, sel_min )
      
      # loop over all possible antenna combinations
      for bb in baselines:
        [ a, b, c, d ] = bb
        
        # check available phases
        ph_array_p = aget( uc_ph_array[ scan[ 0 ] : scan[ 1 ] + 1 ], sel_p )
        ph_array_m = aget( uc_ph_array[ scan[ 0 ] : scan[ 1 ] + 1 ], sel_m )
        er_array_p = aget( er_array[ scan[ 0 ] : scan[ 1 ] + 1 ], sel_p )
        er_array_m = aget( er_array[ scan[ 0 ] : scan[ 1 ] + 1 ], sel_m )
        sel_pm = awhere( ( er_array_p[ : , a - 1 ] > 0. ) & 
            ( er_array_p[ : , b - 1 ] > 0. ) &
            ( er_array_m[ : , c - 1 ] > 0. ) & 
            ( er_array_m[ : , d - 1 ] > 0. ) )
        sel_mp = awhere( ( er_array_m[ : , a - 1 ] > 0. ) & 
            ( er_array_m[ : , b - 1 ] > 0. ) &
            ( er_array_p[ : , c - 1 ] > 0. ) & 
            ( er_array_p[ : , d - 1 ] > 0. ) )
        
        # first process positive lag for ab and negative lag for cd
        if ( len( sel_pm ) > min_points ):
          ph_array_pab = ( aget( ph_array_p[ : , a - 1 ], sel_pm ) - 
              aget( ph_array_p[ : , b - 1 ], sel_pm ) )
          ph_array_mcd = ( aget( ph_array_m[ : , c - 1 ], sel_pm ) - 
              aget( ph_array_m[ : , d - 1 ], sel_pm ) )
          if remove_mean:
            ph_array_pab = ph_array_pab - ph_array_pab.mean()
            ph_array_mcd = ph_array_mcd - ph_array_mcd.mean()
          corr_pm = ( ph_array_pab * ph_array_mcd ).mean()
          results.append( [ a, b, c, d, time_offset, corr_pm ] )
#          pdb.set_trace()
          if ( time_offset < time_tol ):
            plot( aget( aget( stime_array, sel_p ), sel_pm ), ph_array_pab, 'b' )
            plot( aget( aget( stime_array, sel_p ), sel_pm ), ph_array_mcd, 'r' )
            show()
        
        # then process positive lag for cd and negative lag for ab
        if ( ( len( sel_mp ) > min_points ) and ( time_offset > time_tol ) ):
          ph_array_mab = ( aget( ph_array_m[ : , a - 1 ], sel_mp ) - 
              aget( ph_array_m[ : , b - 1 ], sel_mp ) )
          ph_array_pcd = ( aget( ph_array_p[ : , c - 1 ], sel_mp ) - 
              aget( ph_array_p[ : , d - 1 ], sel_mp ) )
          if remove_mean:
            ph_array_mab = ph_array_mab - ph_array_mab.mean()
            ph_array_pcd = ph_array_pcd - ph_array_pcd.mean()
          corr_mp = ( ph_array_mab * ph_array_pcd ).mean()
          results.append( [ a, b, c, d, -time_offset, corr_mp ] )
  
  return results

###############################################################################

def correlate_baseline_phases_newest( uv, in_version = 0, window = 5.,
    max_lag = 5., max_gap = 1., remove_mean = False,
    time_step = 1, height = 300.e3, iterations = 1, min_points = 10,
    print_info = True, tolerance = 0.25, antennas = [], span = 3,
    blength_range = [ 1.e3, 1.e6 ], bangle_range = [ 0., 60. ],
    clength_range = [], cangle_range = [],
    rratio_range = [ 0., 0.75 ], rangle_range = [ 0., 45. ] ):
# window, max_lag & max_gap in minutes
# blength_range for individual baselines
# bangle_range for relative angle between baselines
# clength_range for distance between baselines
# cangle_range for angle of line between baselines
# rratio_range for relative length of baseline to line between baselines
# rangle_range for relative angle between baseline and line between baselines
  
  normalize_weights = True
  debug = False
  double_precision = True
  if double_precision:
    dtype = float64
  else:
    dtype = float32
  
  calibration_data = read_phase_calibration_data( uv, solution_version = in_version,
      print_info = print_info )
  time_array = array( calibration_data[ 0 ], dtype = dtype )
  center_table = calibration_data[ 1 ]
  array_table = calibration_data[ 2 ]
  antenna_table = calibration_data[ 3 ]
  reference_array = calibration_data[ 4 ]
  ph_array = array( calibration_data[ 5 ], dtype = dtype )
  er_array = array( calibration_data[ 6 ], dtype = dtype )
  time_count = len( time_array )
  antenna_count = len( antenna_table )
  if ( len( antennas ) > 0 ):
    ants = [ i for i in antennas ]
    ants.sort( cmp = lambda a, b: cmp( a, b ) )
  else:
    ants = range( 1, 1 + antenna_count )
  
  # calculate some time variables
  dtime = median( time_array[ 1 : , 0 ] - time_array[ : -1, 0 ] )
#  dtime = restore_parameter( uv, 'integration_time' ) / 86400.
  time_tol = tolerance * dtime
  dtime = dtime * float( time_step )
  max_scan = window + max_lag
  min_scan = window + 0.5 * max_lag
  
  # identify scans or large gaps
  scan_list = []
  last_n = - 1
  first_n = 0
  for n in range( time_count ):
    sel = awhere( er_array[ n, : ] > 0. )
    if ( len( sel ) > 0 ):
      if ( last_n == - 1 ):
        scan_entry = [ n ]
        first_n = n
      elif ( time_array[ n, 0 ] - time_array[ last_n, 0 ] > 
          max_gap / ( 24. * 60. ) ):
        scan_entry.append( last_n )
        scan_list.append( scan_entry )
        scan_entry = [ n ]
        first_n = n
      elif ( time_array[ n, 0 ] - time_array[ first_n, 0 ] > 
          max_scan / ( 24. * 60. ) ):
        scan_entry.append( n - 1 )
        scan_list.append( scan_entry )
        scan_entry = [ n ]
        first_n = n
      last_n = n
  scan_entry.append( last_n )
  scan_list.append( scan_entry )
  
  # merge short scans with neighbouring ones
  for s in range( len( scan_list ) ):
    scan = scan_list[ s ]
    if ( time_array[ scan[ 1 ], 0 ] - time_array[ scan[ 0 ], 0 ] < 
        min_scan / ( 24. * 60. ) ):
      delta_low = 1.e7
      delta_high = 1.e7
      if ( s > 0 ):
        scan_low = scan_list[ s - 1 ]
        if ( scan_low == [ - 1.e7, - 1.e7 ] ):
          delta_low = 1.e7
        else:
          delta_low = time_array[ scan[ 0 ], 0 ] - time_array[ scan_low[ 1 ], 0 ]
      if ( s < len( scan_list ) - 1 ):
        scan_high = scan_list[ s + 1 ]
        if ( scan_high == [ - 1.e7, - 1.e7 ] ):
          delta_high = 1.e7
        else:
          delta_high = time_array[ scan_high[ 0 ], 0 ] - time_array[ scan[ 1 ], 0 ]
      if ( ( delta_low < 1.e5 ) or ( delta_high < 1.e5 ) ):
        if ( delta_low < delta_high ):
          if ( delta_low < max_gap  / ( 24. * 60. ) ):
            scan_list[ s - 1 ][ 1 ] = scan[ 1 ]
          scan_list[ s ] = [ - 1.e7, - 1.e7 ]
        else:
          if ( delta_high < max_gap  / ( 24. * 60. ) ):
            scan_list[ s + 1 ][ 0 ] = scan[ 0 ]
          scan_list[ s ] = [ - 1.e7, - 1.e7 ]
  while ( [ - 1.e7, - 1.e7 ] in scan_list ):
    scan_list.remove( [ - 1.e7, - 1.e7 ] )
  
  # assign time stamps to scans
  scan_count = len( scan_list )
  scan_sel = []
  for t in range( time_count ):
    scan_id = - 1
    for scan in scan_list:
      if ( t in range( scan[ 0 ], scan[ 1 ] + 1 ) ):
        scan_id = scan_list.index( scan )
        break
    scan_sel.append( scan_id )
  scan_sel = array( scan_sel, dtype = int32 ).reshape( len( scan_sel ), 1 )
  
  # remove entries within short, isolated scans
  time_list = []
  for scan in scan_list:
    time_list = time_list + range( scan[ 0 ], scan[ 1 ] + 1 )
  if print_info:
    print 'found %s scans' % ( repr( scan_count ) )
  
  # unwrap and average phases per scan
  uc_ph_array = azeros( ph_array )
  antenna_active = zeros( ( scan_count, antenna_count ), dtype = bool )
  for i in range( antenna_count ):
    for s in range( scan_count ):
      scan = scan_list[ s ]
      sel = awhere( er_array[ scan[ 0 ] : scan[ 1 ] + 1, i ] > 0. )
      if ( len( sel ) > 0 ):
        antenna_active[ s, i ] = True
        p1 = aget( ph_array[ scan[ 0 ] : scan[ 1 ] + 1, i ], sel )
        if ( len( p1 ) < 20 ):
          p2 = adegrees( unwrap( aradians( p1 ) ) )
        else:
          p2 = aunwrap_phase( p1, alpha = 0.01 )
          if sometrue( isnan( p2 ) ):
            p2 = aunwrap_phase( p1, alpha = 0.001 )
          if ( sometrue( isnan( p2 ) ) or ( max( fabs( p2 ) ) > 1.e4 ) ):
            p2 = adegrees( unwrap( aradians( p1 ) ) )
        p2 = p2 - p2.mean()
        uc_ph_array[ scan[ 0 ] : scan[ 1 ] + 1, i ] = aput(
            uc_ph_array[ scan[ 0 ] : scan[ 1 ] + 1, i ], sel, p2 )
  
  # interpolate small gaps
  new_uc_ph_array = uc_ph_array.copy()
  new_er_array = er_array.copy()
  sel = awhere( new_er_array <= 0. )
  new_er_array = aput( new_er_array, sel, 0. )
#  lllist = []
  if ( span > 0 ):
    for i in range( antenna_count ):
      for s in range( scan_count ):
        scan = scan_list[ s ]
        if antenna_active[ s, i ]:
          for n in range( scan[ 0 ], scan[ 1 ] + 1 ):
            if ( er_array[ n, i ] <= 0. ):
              sel = awhere( abs( time_array[ scan[ 0 ] : scan[ 1 ] + 1, 0 ] - 
                  time_array[ n, 0 ] ) < float( span ) * dtime + time_tol )
              if ( len( sel ) > 0 ):
                sel2 = awhere( 
                    aget( er_array[ scan[ 0 ] : scan[ 1 ] + 1, i ], sel ) > 0. )
                if ( len( sel2 ) > 0 ):
                  new_uc_ph_array[ n, i ] = aget( aget( uc_ph_array[ scan[ 0 ] : 
                      scan[ 1 ] + 1, i ], sel ), sel2 ).mean()
                  new_er_array[ n, i ] = aget( aget( er_array[ scan[ 0 ] : 
                      scan[ 1 ] + 1, i ], sel ), sel2 ).mean()
#                  if ( not [ s, i ] in lllist ):
#                    lllist.append( [ s, i ] )
#  for ll in lllist:
#    [ s, i ] = ll
#    clf()
#    plot( time_array[ scan_list[ s ][ 0 ] : scan_list[ s ][ 1 ] + 1, 0 ],
#        uc_ph_array[ scan_list[ s ][ 0 ] : scan_list[ s ][ 1 ] + 1, i ], 'b' )
#    plot( time_array[ scan_list[ s ][ 0 ] : scan_list[ s ][ 1 ] + 1, 0 ],
#        new_uc_ph_array[ scan_list[ s ][ 0 ] : scan_list[ s ][ 1 ] + 1, i ], 'r' )
#    title( 'scan = %s, antenna = %s' % ( repr( s ), repr( i+1 ) ) )
#    show()
  uc_ph_array = new_uc_ph_array
  er_array = new_er_array
  
  # check baselines against input criteria
  results = []
  baselines = []
  if print_info:
    print 'selecting baselines'
  for a in ants:
    if ( not antenna_active[ s, a - 1 ] ):
      continue
    x_a = array( antenna_table[ a - 1 ][ 0 ] )
    for b in ants:
      if ( ( not antenna_active[ s, b - 1 ] ) or ( b <= a ) ):
        continue
      x_b = array( antenna_table[ b - 1 ][ 0 ] )
      l_ab = sqrt( add.reduce( ( x_a - x_b )**2 ) )
      if ( len( blength_range ) == 2 ):
        if ( ( l_ab < blength_range[ 0 ] ) or ( l_ab > blength_range[ 1 ] ) ):
          continue
      x_ab = ( x_a + x_b ) / 2.
      llh_ab = xyz_to_geo_llh( x_ab )
      for c in ants:
#        if ( ( not antenna_active[ s, c - 1 ] ) or ( c <= a ) or ( c == b ) ):
        if ( ( not antenna_active[ s, c - 1 ] ) or ( c < a ) ):
          continue
        x_c = array( antenna_table[ c - 1 ][ 0 ] )
        for d in ants:
#          if ( ( not antenna_active[ s, d - 1 ] ) or ( d <= c ) or ( d == b ) ):
          if ( ( not antenna_active[ s, d - 1 ] ) or ( d <= c ) or 
              ( [ c, d ] == [ a, b ] ) ):
            continue
          x_d = array( antenna_table[ d - 1 ][ 0 ] )
          l_cd = sqrt( add.reduce( ( x_c - x_d )**2 ) )
          if ( len( blength_range ) == 2 ):
            if ( ( l_cd < blength_range[ 0 ] ) or ( l_cd > blength_range[ 1 ] ) ):
              continue
          x_cd = ( x_c + x_d ) / 2.
          llh_cd = xyz_to_geo_llh( x_cd )
          cos_ra_abcd = dot( ( x_a - x_b ), ( x_c - x_d ) ) / ( l_ab * l_cd )
          ra_abcd = degrees( acos( max( min( cos_ra_abcd, 1. ), -1. ) ) )
          if ( ra_abcd > 90. ):
            ara_abcd = 180. - ra_abcd
          else:
            ara_abcd = ra_abcd
          if ( len( bangle_range ) == 2 ):
            if ( ( ara_abcd < bangle_range[ 0 ] ) or ( ara_abcd > bangle_range[ 1 ] ) ):
              continue
          l_abcd = sqrt( add.reduce( ( x_ab - x_cd )**2 ) )
          if ( len( clength_range ) == 2 ):
            if ( ( l_abcd < clength_range[ 0 ] ) or ( l_abcd > clength_range[ 1 ] ) ):
              continue
          if ( len( rratio_range ) == 2 ):
            if ( ( l_ab / l_abcd < rratio_range[ 0 ] ) or 
                ( l_ab / l_abcd > rratio_range[ 1 ] ) or
                ( l_cd / l_abcd < rratio_range[ 0 ] ) or 
                ( l_cd / l_abcd > rratio_range[ 1 ] ) ):
              continue
          [ r_abcd, a_abcd ] = calculate_angular_separation( llh_cd[ 0 : 2 ],
              llh_ab[ 0 : 2 ] )
          if ( len( cangle_range ) == 2 ):
            da_abcd0 = amodulo( ( a_abcd - cangle_range[ 0 ] ) + 180., 360. ) - 180.
            da_abcd1 = amodulo( ( a_abcd - cangle_range[ 1 ] ) + 180., 360. ) - 180.
            if ( ( da_abcd0 < 0. ) or ( da_abcd1 > 0. ) ):
              continue
          cos_ra_ababcd = dot( ( x_a - x_b ), ( x_ab - x_cd ) ) / ( l_ab * l_abcd )
          ra_ababcd = degrees( acos( max( min( cos_ra_ababcd, 1. ), -1. ) ) )
          if ( ra_ababcd > 90. ):
            ara_ababcd = 180. - ra_ababcd
          else:
            ara_ababcd = ra_ababcd
          if ( len( rangle_range ) == 2 ):
            if ( ( ara_ababcd < rangle_range[ 0 ] ) or ( ara_ababcd > rangle_range[ 1 ] ) ):
              continue
          cos_ra_cdabcd = dot( ( x_c - x_d ), ( x_ab - x_cd ) ) / ( l_cd * l_abcd )
          ra_cdabcd = degrees( acos( max( min( cos_ra_cdabcd, 1. ), -1. ) ) )
          if ( ra_cdabcd > 90. ):
            ara_cdabcd = 180. - ra_cdabcd
          else:
            ara_cdabcd = ra_cdabcd
          if ( len( rangle_range ) == 2 ):
            if ( ( ara_cdabcd < rangle_range[ 0 ] ) or ( ara_cdabcd > rangle_range[ 1 ] ) ):
              continue
          baselines.append( [ a, b, c, d ] )
          results.append( [ a, b, c, d, l_ab, l_cd, ra_abcd, l_abcd, a_abcd,
              ra_ababcd, ra_cdabcd ] )
  
  # correlate per scan
  for s in range( scan_count ):
    X_table = zeros( ( antenna_count, 2 ), dtype = float64 )
    za_table = zeros( ( antenna_count ), dtype = float64 )
    scan = scan_list[ s ]
    slen = scan[ 1 ] - scan[ 0 ] + 1
    if print_info:
      print '... processing scan %s' % ( repr( s + 1 ) )
    
    # get pierce point coordinates at center of scan
    stime_array = time_array[ scan[ 0 ] : scan[ 1 ] + 1, 0 ]
    center_time = ( stime_array[ 0 ] + stime_array[ -1 ] ) / 2.
    dtime_array = abs( stime_array - center_time )
    center_index = awhere( dtime_array == dtime_array.min() )[ 0, 0 ] + scan[ 0 ]
    pierce_table = calculate_pierce_coordinates( time_array[ center_index ],
        center_table, [ center_table ], array_table, antenna_table,
        height = height, iterations = iterations )
    
    # loop over pierce points
    for pierce_info in pierce_table:
      [ X, za, [ k, i ] ] = pierce_info
      # store model fit input data
      X_table[ i ] = array( X, dtype = float64 )
      za_table[ i ] = za
    results.append( [ time_array[ center_index ], X_table, za_table ] )
    
    # loop over all possible lags
    center_time = time_array[ center_index, 0 ]
    time_offset = -dtime
#    pdb.set_trace()
    while ( time_offset < max_lag / ( 60. * 24. ) ):
      time_offset = time_offset + dtime
      if ( print_info ):
        print '...... time lag = %s' % ( repr( time_to_dhms( time_offset ) ) )
      stime_array_p = stime_array + ( time_offset / 2. )
      stime_array_m = stime_array - ( time_offset / 2. )
      dstime_array = resize( stime_array_m, ( slen, slen ) ) - transpose( 
          resize( stime_array_p, ( slen, slen ) ) )
      sel_pm = awhere( abs( dstime_array ) < time_tol )
      if ( len( sel_pm ) < min_points ):
        if ( print_info ):
          print '......... skipping lag, too little data'
        continue
      sel_p = sel_pm[ : , 0 : 1 ]
      sel_m = sel_pm[ : , 1 : 2 ]
      
      # narrow selection down to time window
      adstime_array_p = abs( aget( stime_array_p, sel_p ) - center_time )
      adstime_array_m = abs( aget( stime_array_m, sel_m ) - center_time )
      sel_pm = awhere( ( adstime_array_p < window / ( 2. * 60. * 24. ) ) & 
          ( adstime_array_m < window / ( 2. * 60. * 24. ) ) )
      if ( len( sel_pm ) < min_points ):
        if ( print_info ):
          print '......... skipping lag, too little data in windows'
        continue
      sel_p = aget( sel_p, sel_pm )
      sel_m = aget( sel_m, sel_pm )
      
      # loop over all possible antenna combinations
      for bb in baselines:
        [ a, b, c, d ] = bb
        
        # check available phases
        ph_array_p = aget( uc_ph_array[ scan[ 0 ] : scan[ 1 ] + 1 ], sel_p )
        ph_array_m = aget( uc_ph_array[ scan[ 0 ] : scan[ 1 ] + 1 ], sel_m )
        er_array_p = aget( er_array[ scan[ 0 ] : scan[ 1 ] + 1 ], sel_p )
        er_array_m = aget( er_array[ scan[ 0 ] : scan[ 1 ] + 1 ], sel_m )
        sel_pm = awhere( ( er_array_p[ : , a - 1 ] > 0. ) & 
            ( er_array_p[ : , b - 1 ] > 0. ) &
            ( er_array_m[ : , c - 1 ] > 0. ) & 
            ( er_array_m[ : , d - 1 ] > 0. ) )
        sel_mp = awhere( ( er_array_m[ : , a - 1 ] > 0. ) & 
            ( er_array_m[ : , b - 1 ] > 0. ) &
            ( er_array_p[ : , c - 1 ] > 0. ) & 
            ( er_array_p[ : , d - 1 ] > 0. ) )
        
        # first process positive lag for ab and negative lag for cd
        if ( len( sel_pm ) > min_points ):
          ph_array_pab = ( aget( ph_array_p[ : , a - 1 ], sel_pm ) - 
              aget( ph_array_p[ : , b - 1 ], sel_pm ) )
          ph_array_mcd = ( aget( ph_array_m[ : , c - 1 ], sel_pm ) - 
              aget( ph_array_m[ : , d - 1 ], sel_pm ) )
          if remove_mean:
            ph_array_pab = ph_array_pab - ph_array_pab.mean()
            ph_array_mcd = ph_array_mcd - ph_array_mcd.mean()
          corr_pm = ( ph_array_pab * ph_array_mcd ).mean()
          results.append( [ a, b, c, d, time_offset, corr_pm ] )
#          pdb.set_trace()
#          if ( time_offset < time_tol ):
#            plot( aget( aget( stime_array, sel_p ), sel_pm ), ph_array_pab, 'b' )
#            plot( aget( aget( stime_array, sel_p ), sel_pm ), ph_array_mcd, 'r' )
#            title( 'antennas = %s, corr = %s' % ( repr( [ a, b, c, d ] ), repr( corr_pm ) ) )
#            show()
        
        # then process positive lag for cd and negative lag for ab
        if ( ( len( sel_mp ) > min_points ) and ( time_offset > time_tol ) ):
          ph_array_mab = ( aget( ph_array_m[ : , a - 1 ], sel_mp ) - 
              aget( ph_array_m[ : , b - 1 ], sel_mp ) )
          ph_array_pcd = ( aget( ph_array_p[ : , c - 1 ], sel_mp ) - 
              aget( ph_array_p[ : , d - 1 ], sel_mp ) )
          if remove_mean:
            ph_array_mab = ph_array_mab - ph_array_mab.mean()
            ph_array_pcd = ph_array_pcd - ph_array_pcd.mean()
          corr_mp = ( ph_array_mab * ph_array_pcd ).mean()
          results.append( [ a, b, c, d, -time_offset, corr_mp ] )
  
  return results

###############################################################################

def fit_cosine( P, dojac = None, phases = None, amplitudes = None, errors = None ):
  if ( ( not dojac is None ) or ( phases is None ) or ( amplitudes is None ) or
      ( errors is None ) ):
    return - 1, None, None
  model = cos( radians( phases + P[ 0 ] ) ) / P[ 1 ]
  damp = amplitudes - model
#  chi_array = sign( damp ) * ( abs( damp )**0.3 ) # * ( 1 + cos( radians( P[ 0 ] ) )**2 )
  chi_array = sign( damp ) * ( abs( damp ) ) / errors
  return 0, chi_array, None

def fit_velocity( phases, amplitudes, errors, print_info = True ):
  function_keywords = { 'phases' : phases, 'amplitudes' : amplitudes,
      'errors' : errors }
  parameter_info = [ 
      { 'parname' : 'direction', 'value' : 0., 'limits' : [ -360., 360. ] },
      { 'parname' : 'speed', 'value' : 100., 'limits' : [ 0., 1.e12 ] } ]
  fit = mpfit( fit_cosine, functkw = function_keywords,
      parinfo = parameter_info, quiet = True, autoderivative = True,
      debug = False, fastnorm = False, nocovar = False, dblprec = True )
  if ( not ( fit.status in [ 1, 2, 3 ] ) ):
    return [ None, None ]
  [ direction, speed ] = fit.params.copy()
  [ dir_error, speed_error ] = sqrt( abs( diagonal( fit.covar ) ) )
  if print_info:
    print 'post-fit error is %s' % repr( chi )
    print 'speed is %s' % repr( speed )
    print 'direction is %s' % repr( direction )
  return [ [ speed, speed_error ], [ direction, dir_error ] ]

###############################################################################

def fit_ionospheric_pmkl2_model( uv, facets, facet_list = [], beta = 5. / 3., 
    r_0 = 1., order = 15, layer_heights = [ 200.e3 ], layer_weights = [ 1. ], 
    iterations = 4, equal_weights = True, normalize_weights = True, 
    e_steps = 1, double_precision = True, print_info = True, 
    estimate_offsets = False, time_steps = [], include_airmass = True, 
    solution_version = 0, source_rms_limit = 1.e6, first_order = False,
    antenna_offset_limit = 1.e6, antenna_rms_limit = 1.e6, fix_gradient = True,
    prop_time_limit = 1., prop_rms_limit = 0., dof_factor = 3,
    exclude_antennas = [], density_scale = None, density_power = None ):
# prop_time_limit in minutes
  
  if double_precision:
    dtype = float64
  else:
    dtype = float32
  l_max = layer_weights.index( max( layer_weights ) )
  
  # select appropriate model functions
  # initialize model parameters
  model = 'pmkl2'
  layer_count = len( layer_heights )
  if first_order:
    new_order = order - 2
  else:
    new_order = order - 5
  if ( new_order < 0 ):
    raise error( 'order is too low to fit polynomial' )
  P = zeros( shape = ( new_order ), dtype = dtype )
  if ( len( layer_weights ) != layer_count ):
    raise error( 'layer heights and weights are different lengths' )
  
  if print_info:
    print 'starting model fits using model = ' + model
  
  # get all relevant data
  calibration_data = get_phase_calibration_data( uv, facets, facet_list = facet_list,
      solution_version = solution_version, print_info = print_info )
  time_table = calibration_data[ 0 ]
  center_table = calibration_data[ 1 ]
  peel_table = calibration_data[ 2 ]
  array_table = calibration_data[ 3 ]
  antenna_table = calibration_data[ 4 ]
  reference_array = calibration_data[ 5 ]
  phase_array = calibration_data[ 6 ]
  error_array = calibration_data[ 7 ]
  source_count = len( peel_table )
  antenna_count = len( antenna_table )
  
  # compensate for antenna density, downweight dense antennas
  weight_table = [ 1. for a in antenna_table ]
  if ( ( not density_scale is None ) and ( density_scale > 0. ) ):
    X_table = []
    for a in antenna_table:
      if ( a[ 0 ] != [ 0., 0., 0. ] ):
        X_table.append( a[ 0 ] )
    X_count = len( X_table )
    X_table = array( X_table, dtype = float64 )
    D_table = resize( X_table, ( X_count, X_count, 3 ) )
    D_table = transpose( D_table, ( 1, 0, 2 ) ) - D_table
    D_table = add.reduce( D_table**2, 2 )
#    W_table = exp( - D_table / density_scale**2 )
    W_table = D_table / density_scale**2
    W_table = aput( W_table, awhere( W_table > 100. ), 100. )
    W_table = exp( - W_table )
    W_table = 1. / add.reduce( W_table, 0 )
    if ( not density_power is None ):
      W_table = W_table**density_power
    j = 0
    for i in range( len( antenna_table ) ):
      if ( antenna_table[ i ][ 0 ] != [ 0., 0., 0. ] ):
        weight_table[ i ] = W_table[ j ]
        j = j + 1
  
  fit_table = []
  pierce_table = []
  
  if estimate_offsets:
    solution_table = []

  prop_last_time = -1.e9
  prop_last_poly = zeros( ( 5 ), dtype = dtype )
  
  # loop over time stamps
  if ( len( time_steps ) == 0 ):
    n_list = range( len( time_table ) )
  else:
    n_list = time_steps
  for n in n_list:

    if print_info:
      print '... time step n = %d / %d' % ( n + 1, len( n_list ) )
    
    repeat_fit = True
    rejected_sources = False
    rejected_antennas = False
    while repeat_fit:
      skip = False
      repeat_fit = False
      if ( not rejected_sources ):
        use_sources = [ True for k in range( source_count ) ]
      if ( not rejected_antennas ):
        use_antennas = [ True for i in range( antenna_count ) ]
        for i in exclude_antennas:
          use_antennas[ i - 1 ] = False
      
      temp_pierce_table = []

      # containers for model input data
      ref_list = []
      ref_table = []
      phase_table = []
      error_table = []
      info_table = []
      Xpl_table = []
      pzal_table = []
#      Bl_table = []

      if estimate_offsets:
        solution_row = [ [ time_table[ n ][ 0 ], 
            reference_array[ n ][ 0 ] + 1, 0., 0. ] ] + [ [ 0., 0., 0., 0. ] 
            for i in range( antenna_count ) ]

      # loop over height layers
      for l in range( layer_count ):
        height = layer_heights[ l ]
      
        # containers for model input data
        Xp_table = []
        pza_table = []
      
        # get pierce point coordinates
        ipp_table = calculate_pierce_coordinates( time_table[ n ], 
            center_table, peel_table, array_table, antenna_table, height = height, 
            iterations = iterations )
      
        # loop over pierce points
        j = 0
        for ipp_info in ipp_table:
          [ X, pza, [ k, i ] ] = ipp_info
          if ( ( error_array[ n ][ k ][ i ] > 0. ) and
              use_sources[ k ] and use_antennas[ i ] ):
          
            # store model fit input data
            Xp_table.append( X )
            pza_table.append( pza )
            if ( l == 0 ):
              phase_table.append( phase_array[ n ][ k ][ i ] )
              if equal_weights:
                error_table.append( 1. / weight_table[ i ] )
              else:
                error_table.append( error_array[ n ][ k ][ i ] / weight_table[ i ] )
              if ( i == reference_array[ n ][ k ] ):
                ref_list.append( [ k, i, j ] )
              j = j + 1
              info_table.append( [ k, i, reference_array[ n ][ k ] ] )
          
            # temporarily store pierce point info
            temp_pierce_table.append( [ time_table[ n ][ 0 ], l + 1, k + 1, 
                i + 1, X, pza ] )
        
        if ( l == 0 ):
          # loop over pierce points
          for ipp_info in ipp_table:
            [ X, pza, [ k, i ] ] = ipp_info
            if ( ( error_array[ n ][ k ][ i ] > 0. ) and
                use_sources[ k ] and use_antennas[ i ] ):
              
              # store index to reference antenna
              j = ( [ ( k == ref[ 0 ] ) for ref in ref_list ] ).index( True )
              ref_table.append( ref_list[ j ][ 2 ] )
        
        # skip times of poor data
        if ( len( Xp_table ) <= dof_factor * order ):
          if print_info:
            print '...... skipping time step due to no / too few data points'
          if first_order:
            fit_table.append( [ time_table[ n ][ 0 ], 0, 
                zeros( ( 2 ), dtype = dtype ), azeros( P ).tolist(), 0. ] )
          else:
            fit_table.append( [ time_table[ n ][ 0 ], 0, 
                zeros( ( 5 ), dtype = dtype ), azeros( P ).tolist(), 0. ] )
          for p in temp_pierce_table:
            pierce_table.append( p )
          if estimate_offsets:
            solution_table.append( solution_row )
          skip = True
          break
      
        # convert tables to arrays for fitting routines
        Xp_table = array( Xp_table, dtype = dtype )
        pza_table = array( pza_table, dtype = dtype )
        if ( not include_airmass ):
          pza_table = azeros( pza_table )
        if ( l == 0 ):
          ref_table = array( ref_table, dtype = int32 )
          phase_table = array( phase_table, dtype = dtype )
          error_table = array( error_table, dtype = dtype )
          info_table = array( info_table, dtype = int32 )
        p_count = len( Xp_table )
      
        if ( new_order > 0 ):
          # calculate structure matrix
          D_table = resize( Xp_table, ( p_count, p_count, 2 ) )
          D_table = transpose( D_table, ( 1, 0, 2 ) ) - D_table
          D_table = add.reduce( D_table**2, 2 )
          D_table = ( D_table / ( r_0**2 ) )**( beta / 2. )
        
          # calculate covariance matrix C
          C_table = - D_table / 2.
          C_table = transpose( C_table - ( add.reduce( C_table, 0 ) / float( p_count ) ) )
          # calculate partial product for interpolation B
#          B_table = add.reduce( C_table, 0 ) / float( p_count )
#          C_table = C_table - B_table
          C_table = C_table - add.reduce( C_table, 0 ) / float( p_count )
          C_table = ( C_table + transpose( C_table ) ) / 2.
        
          # incorporate airmass functions and layer weights
          # add layer C to total C
          A_table = resize( 1. / cos( aradians( pza_table ) ), ( p_count, p_count ) )
          A_table = A_table * transpose( A_table )
          if ( l == 0 ):
            Cl_table = C_table * A_table * ( layer_weights[ l ]**2 )
          else:
            Cl_table = Cl_table + C_table * A_table * ( layer_weights[ l ]**2 )
# test          if ( l == l_max ):
# test            R_table = dot( C_table, diag( 1. / cos( aradians( pza_table ) ) ) )
          
        # save tables per height layer
        Xpl_table.append( Xp_table )
        pzal_table.append( pza_table )
#        if ( new_order > 0 ):
#          Bl_table.append( B_table )
        
      # skip times of poor data
      if skip:
        continue
        
      # convert to arrays
      Xpl_table = array( Xpl_table, dtype = dtype )
      pzal_table = array( pzal_table, dtype = dtype )
      if ( new_order > 0 ):
#        Bl_table = array( Bl_table, dtype = dtype )
        
        # project out polynomial fit
        Xp_table = Xpl_table[ l_max ]
        pza_table = pzal_table[ l_max ]
        if first_order:
          Yt_table = array( [ Xp_table[ : , 0 ], Xp_table[ : , 1 ] ] )
        else:
          Yt_table = array( [ Xp_table[ : , 0 ], Xp_table[ : , 1 ], Xp_table[ : , 0 ]**2,
              Xp_table[ : , 0 ] * Xp_table[ : , 1 ], Xp_table[ : , 1 ]**2 ] )
        Yt_table = array( [ Xp_table[ : , 0 ], Xp_table[ : , 1 ] ] )
        Yt_table = Yt_table / cos( aradians( pza_table ) )
        Y_table = transpose( Yt_table )
        Yt_table = dot( linalg.inv( dot( Yt_table, Y_table ) ), Yt_table )
        Y_table = dot( Y_table, Yt_table )
        Y_table = diag( ones( ( p_count ), dtype = Y_table.dtype ) ) - Y_table
        # here, transpose( Y_table ) = Y_table
        Cl_table = dot( Y_table, dot( Cl_table, Y_table ) )
        
        # eigenvalue decomposition
        # reforce symmetry
        # select subset of base vectors
        try:
          [ U_table, S, Ut_table ] = linalg.svd( Cl_table )
        except:
          if print_info:
            print '... SVD did not converge, skipping time interval'
          if first_order:
            fit_table.append( [ time_table[ n ][ 0 ], 0, 
                zeros( ( 2 ), dtype = dtype ), azeros( P ).tolist(), 0. ] )
          else:
            fit_table.append( [ time_table[ n ][ 0 ], 0, 
                zeros( ( 5 ), dtype = dtype ), azeros( P ).tolist(), 0. ] )
          for p in temp_pierce_table:
            pierce_table.append( p )
          if estimate_offsets:
            solution_table.append( solution_row )
          continue
        U_table = ( U_table + transpose( Ut_table ) ) / 2.
        U_table = U_table[ : , 0 : new_order ]
        S = S[ 0 : new_order ]
# test       R_table = dot( R_table, dot( U_table, diag( 1. / S ) ) )
# test       R_table = dot( linalg.inv( dot( transpose( R_table ), R_table ) ),
# test           transpose( R_table ) )
      
      # see if we can use the gradient of the previous good fit as a starting point
      dtime = abs( time_table[ n ][ 0 ] - prop_last_time ) * 24. * 60.
      Xp_table = Xpl_table[ l_max ] / cos( aradians( pzal_table[ l_max ] ) ).reshape( 
          ( p_count, 1 ) )
      if ( dtime < prop_time_limit ):
        poly = prop_last_poly.copy()
        if first_order:
          poly[ 2 : 5 ] = 0.
      else:
        poly = azeros( prop_last_poly )
        # estimate initial phase gradient on layer with heighest weight
        poly[ 0 : 2 ] = estimate_phase_gradient_through_fft( uv, Xp_table,
            ref_table, phase_table, error_table, fft_cell_size = 200.,
            fft_image_size = 1024 )
      
      # refine gradient fit
      function_keywords = { 'phi_model' : phi_poly_model, 
          'X_table' : Xpl_table[ l_max ], 'pza_table' : pzal_table[ l_max ], 
          'ref_table' : ref_table, 'phase_table' : phase_table,
          'error_table' : error_table, 'normalize_weights' : normalize_weights,
          'wrap_phases' : True }
      parameter_info = []
      for m in range( 5 ):
        if ( m < 2 ):
          par_info = { 'parname' : 'P_%d' % ( m ), 'value' : poly[ m ], 
              'limits' : [ None, None ] }
        else:
          par_info = { 'parname' : 'P_%d' % ( m ), 'value' : poly[ m ], 
              'limits' : [ poly[ m ], poly[ m ] ] }
        parameter_info.append( par_info )
      fit = mpfit( fit_phi_model, functkw = function_keywords, parinfo = parameter_info,
          quiet = True, autoderivative = True, debug = False, fastnorm = False, 
          nocovar = True, dblprec = double_precision )
      poly = fit.params.copy()
      chi_poly = sqrt( fit.fnorm )
      
      if ( not first_order ):
        # refine higher order fit
        function_keywords = { 'phi_model' : phi_poly_model, 
            'X_table' : Xpl_table[ l_max ], 'pza_table' : pzal_table[ l_max ], 
            'ref_table' : ref_table, 'phase_table' : phase_table,
            'error_table' : error_table, 'normalize_weights' : normalize_weights,
            'wrap_phases' : True }
        parameter_info = []
        for m in range( 5 ):
          if ( fix_gradient and ( m < 2 ) ):
            par_info = { 'parname' : 'P_%d' % ( m ), 'value' : poly[ m ], 
                'limits' : [ poly[ m ], poly[ m ] ] }
          else:
            par_info = { 'parname' : 'P_%d' % ( m ), 'value' : poly[ m ], 
                'limits' : [ None, None ] }
          parameter_info.append( par_info )
        fit = mpfit( fit_phi_model, functkw = function_keywords, quiet = True, 
            autoderivative = True, debug = False, fastnorm = False, nocovar = True,
            dblprec = double_precision, parinfo = parameter_info )
        poly = fit.params.copy()
        chi_poly = sqrt( fit.fnorm )
      if print_info:
        print '... polynomial post-fit phase RMS is %s degrees' % repr( chi_poly )
      
      # remove polynomial model from data
      phi_poly_table = phi_poly_model( Xpl_table[ l_max ], poly ) / cos( 
          aradians( pzal_table[ l_max ] ) )
      sel = ref_table.reshape( p_count, 1 )
      phi_poly_table = phi_poly_table - aget( phi_poly_table, sel )
      if ( new_order > 0 ):
        phase_table = amodulo( ( phase_table - phi_poly_table ) + 180., 360. ) - 180.
        
        # define data passed to the model fit routines
        # define model parameter structure
        function_keywords = { 'U_table' : U_table, 'ref_table' : ref_table, 
            'phase_table' : phase_table, 'error_table' : error_table, 
            'normalize_weights' : normalize_weights }
        
        # specify fit parameters
        parameter_info = []
        for m in range( new_order ):
          par_info = { 'parname' : 'P_%d' % ( m + 1 ), 'value' : 0., 
              'limits' : [ None, None ] }
          parameter_info.append( par_info )
        
        # fit model to data points
        fit = mpfit( fit_phi_mkl_model, functkw = function_keywords, 
            parinfo = parameter_info, quiet = True, autoderivative = True, 
            debug = False, fastnorm = False, nocovar = True,
            dblprec = double_precision )
        if ( fit.errmsg != '' ):
          if print_info:
            print 'WARNING: %s' % ( fit.errmsg )
          if first_order:
            fit_table.append( [ time_table[ n ][ 0 ], 0, 
                zeros( ( 2 ), dtype = dtype ), azeros( P ).tolist(), 0. ] )
          else:
            fit_table.append( [ time_table[ n ][ 0 ], 0, 
                zeros( ( 5 ), dtype = dtype ), azeros( P ).tolist(), 0. ] )
          for p in temp_pierce_table:
            pierce_table.append( p )
          if estimate_offsets:
            solution_table.append( solution_row )
          skip = True
          break

        # store fit results
        P = fit.params.copy()
        chi_pmkl = sqrt( fit.fnorm )
      
        if print_info:
          print '... model post-fit phase RMS is %s degrees' % repr( chi_pmkl )
      
      # skip times of error in fit
      if skip:
        continue
      
      if ( new_order > 0 ):
        # check that model fit is better than polynomial fit
        if ( chi_pmkl < chi_poly ):
          chi = chi_pmkl
        else:
          if print_info:
            print '... WARNING: %s model fit is no improvement' % ( model )
          P = azeros( P )
          chi = chi_pmkl
      else:
        chi = chi_poly

      # repeat bad fit if previous gradient was used
      if ( ( dtime < prop_time_limit ) and ( chi > prop_rms_limit ) ):
        if print_info:
          print '...... repeating fit using fresh gradient estimate'
        rejected_sources = False
        rejected_antennas = False
        repeat_fit = True
        prop_last_time = -1.e9
        continue
      
      if ( new_order > 0 ):
        # add polynomial model back
        phase_table = amodulo( ( phase_table + phi_poly_table ) + 180., 360. ) - 180.
        
        # calculate model phases
        phi_mkl_table = dot( U_table, P )
        sel = ref_table.reshape( p_count, 1 )
        phi_mkl_table = phi_mkl_table - aget( phi_mkl_table, sel )
        phi_table = amodulo( ( phi_poly_table + phi_mkl_table ) + 180., 360. ) - 180.
      else:
        phi_table = amodulo( phi_poly_table + 180., 360. ) - 180.
      
      # split pierce points per source
      ref_list = []
      for ref in ref_table:
        if ( not ( ref in ref_list ) ):
          ref_list.append( ref )
      
      # estimate residual phase per antenna per source
      antenna_error_table = azeros( error_table )
      antenna_offset_table = azeros( error_table )
      reference_antennas = []
      for ref in ref_list:
        reference_antennas.append( info_table[ ref ][ 1 ] )
        sel = awhere( ref_table == ref )
        dphase_table = aget( phase_table, sel )
        dphase_table = resize( dphase_table, ( len( sel ), len( sel ) ) )
        dphase_table = transpose( dphase_table ) - dphase_table
        dphi_table = aget( phi_table, sel )
        dphi_table = resize( dphi_table, ( len( sel ), len( sel ) ) )
        dphi_table = transpose( dphi_table ) - dphi_table
        dweight_table = aget( error_table, sel )
        dweight_table = resize( dweight_table**2, ( len( sel ), len( sel ) ) )
        dweight_table = transpose( dweight_table ) + dweight_table
        dweight_table = 1. / sqrt( dweight_table )
        abs_err_table = add.reduce( dweight_table * 
            abs( amodulo( ( dphase_table - dphi_table ) + 180., 360. ) - 180. ), 1 )
        normalize_table = add.reduce( dweight_table, 1 ) - diagonal( dweight_table )
        abs_err_table = abs_err_table / normalize_table
        antenna_error_table = aput( antenna_error_table, sel, abs_err_table )
        offset_table = add.reduce( dweight_table * 
            ( amodulo( ( dphase_table - dphi_table ) + 180., 360. ) - 180. ), 1 )
        offset_table = offset_table / normalize_table
        antenna_offset_table = aput( antenna_offset_table, sel, offset_table )
      
      # estimate residual phase offset and stddev per antenna
      antenna_means = zeros( ( antenna_count ), dtype = dtype )
      antenna_stds = zeros( ( antenna_count ), dtype = dtype )
      for i in range( antenna_count ):
        # select all phases related to single antenna
        sel = awhere( info_table[ : , 1 ] == i )
        # TODO: determine reasonable number of values for calculating the antenna mean
        if ( len( sel ) < 3 ):
          continue
        # calculate mean phase offset per antenna
        phase_offsets = aget( antenna_offset_table, sel )
        mean_phase_offset = amean_phase( phase_offsets )
        std_phase_offset = sqrt( ( ( amodulo( ( phase_offsets - mean_phase_offset )
            + 180., 360. ) - 180. )**2 ).mean() )
        antenna_means[ i ] = mean_phase_offset
        antenna_stds[ i ] = std_phase_offset
        if estimate_offsets:
          # store phase error in solution table
          # note this automatically maps mean_phase in [ -180, 180 ] domain
          solution = r_phi_to_complex( [ 1., mean_phase_offset ] )
          solution_row[ i + 1 ] = [ solution.real, solution.imag, 0.,
              1. / std_phase_offset ]
      
      # check for antennas with excessive offset or stddev
      if ( not rejected_antennas ):
        sel1 = awhere( antenna_means > antenna_offset_limit )
        for s in sel1:
          s = s[ 0 ]
          if print_info:
            print '...... antenna %s has excessive offset of %s degrees' % (
                repr( s + 1 ), repr( antenna_means[ s ] ) )
          if ( s in reference_antennas ):
            if print_info:
              print '...... WARNING: reference antenna not dropped'
          else:
            use_antennas[ s ] = False
            rejected_antennas = True
        sel2 = awhere( antenna_stds > antenna_rms_limit )
        for s in sel2:
          s = s[ 0 ]
          if use_antennas[ s ]:
            if print_info:
              print '...... antenna %s has excessive stddev of %s degrees' % (
                  repr( s + 1 ), repr( antenna_stds[ s ] ) )
            if ( s in reference_antennas ):
              if print_info:
                print '...... WARNING: reference antenna not dropped'
            else:
              use_antennas[ s ] = False
              rejected_antennas = True
        if rejected_antennas:
          if print_info:
            print '...... repeating fit'
          repeat_fit = True
          if ( chi > prop_rms_limit ):
            prop_last_time = -1.e9
          continue
      
      # estimate residual phase RMS per source
      source_rmss = zeros( ( source_count ), dtype = dtype )
      for k in range( source_count ):
        # select all phases related to single antenna
        sel = awhere( info_table[ : , 0 ] == k )
        # TODO: determine reasonable number of values for calculating the mean
        if ( len( sel ) < 3 ):
          continue
        # calculate mean phase offset per source
        phase_offsets = aget( antenna_offset_table, sel )
        rms_phase_offset = sqrt( ( ( amodulo( phase_offsets
            + 180., 360. ) - 180. )**2 ).mean() )
        source_rmss[ k ] = rms_phase_offset
      
      # check for source with excessive RMS
      if ( not rejected_sources ):
        sel = awhere( source_rmss > source_rms_limit )
        for s in sel:
          s = s[ 0 ]
          use_sources[ s ] = False
          rejected_sources = True
          if print_info:
            print '...... source %s has excessive RMS of %s degrees' % (
                repr( facet_list[ s ] ), repr( source_rmss[ s ] ) )
        if rejected_sources:
          if print_info:
            print '...... repeating fit'
          repeat_fit = True
          if ( chi > prop_rms_limit ):
            prop_last_time = -1.e9
          continue
      
      # store fit results
      weight = 1. / chi
      if first_order:
        fit_table.append( [ time_table[ n ][ 0 ], p_count, poly[ 0 : 2 ].tolist(),
            P.tolist(), weight ] )
      else:
        fit_table.append( [ time_table[ n ][ 0 ], p_count, poly.tolist(),
            P.tolist(), weight ] )
      for p in temp_pierce_table:
        pierce_table.append( p )
      
      if estimate_offsets:
        # add solution row to solution table
        solution_table.append( solution_row )
      
      # see if we can use the polynomial of the previous good fit as a starting point
      if ( chi < prop_rms_limit ):
        prop_last_poly = poly.copy()
        prop_last_time = time_table[ n ][ 0 ]
  
  # write fit results to model fit table
  write_ionospheric_pmkl2_model_fit_table( uv, layer_heights, layer_weights,
       fit_table, beta = beta, r_0 = r_0 )
  fit_version = uv.table_highver( 'NI' )
  write_ionospheric_mkl_pierce_table( uv, pierce_table, version = fit_version,
      iterations = iterations )
  
  # write solution table with antenna phase offsets to uv data
  if estimate_offsets:
    if print_info:
      print '... NOTE: writing solution table with antenna phase offsets to UV data'
    write_solution_table( uv, solution_table, out_version = 0 )
  
  return

###############################################################################

def write_ionospheric_pmkl2_model_fit_table( uv, layer_heights, layer_weights,
    fit_table, version = 0, beta = 5. / 3., r_0 = 1., **keywords ):

  # create new NI table
  layer_count = len( layer_heights )
  poly_order = len( fit_table[ 0 ][ 2 ] )
  mkl_order = len( fit_table[ 0 ][ 3 ] )
  order = poly_order + mkl_order
  new_ni_table = new_table( uv, 'NI', version, num_coef = order )
  row = new_table_row( new_ni_table )
  integration_time = restore_parameter( uv, 'integration_time' )
  time_interval = integration_time / ( 24. * 60. * 60. )

  # loop over 
  for fit_row in fit_table:
    [ time, p_count, poly, coefs, weight ] = fit_row
    row.time = float32( time )
    row.time_interval = float32( time_interval )
    row.antenna_no = p_count
    row.source_id = 0
    row.subarray = 0
    row.weight = float32( weight )
    row.coef = [ float32( pol ) for pol in poly ] + [ float32( coef ) for coef in coefs ]
    new_ni_table.append( row )

  # add keywords
  new_ni_table.keywords[ 'MODEL' ] = 'pmkl2'
  new_ni_table.keywords[ 'REF_FREQ' ] = float32( get_central_frequency( uv ) )
  new_ni_table.keywords[ 'BETA' ] = float32( beta )
  new_ni_table.keywords[ 'R_0' ] = float32( r_0 )
  new_ni_table.keywords[ 'POLY' ] = poly_order
  new_ni_table.keywords[ 'LAYERS' ] = layer_count
  for l in range( layer_count ):
    new_ni_table.keywords[ 'HEIGHT%d' % ( l + 1 ) ] = float32( layer_heights[ l ] )
    new_ni_table.keywords[ 'WEIGHT%d' % ( l + 1 ) ] = float32( layer_weights[ l ] )
  for key in keywords.keys():
    new_ni_table.keywords[ key ] = keywords[ key ]

  #close new table
  new_ni_table.close()
  
  return

###############################################################################

def generate_solution_tables_from_pmkl2_fit_table( uv, facets, fit_version = 0,
    reference_antenna = 0, solution_version = 0, include_delays = False,
    print_info = True, rejection_limit = None, rms_limit = None, facet_list = [],
    time_steps = [], rms_min = 0.1, include_antennas = [], include_phases = True ):
  
  wiz_uv = wizardry( uv )
  ni_table = wiz_uv.table( 'NI', fit_version )
  model = ni_table.keywords[ 'MODEL' ]
  if ( type( model ) == type( [] ) ):
    model = model[ 0 ]
  model = model.strip()
  if ( model != 'pmkl2' ):
    raise error( 'unexpected model: %s' % ( model ) )
  poly_order = int( ni_table.keywords[ 'POLY' ] )
  order = int( ni_table.keywords[ 'NUM_COEF' ] ) - poly_order
  if include_delays:
    reference_frequency = get_frequency( uv )
    bandwidth = get_bandwidth( uv )
    solution_frequency = float32( ni_table.keywords[ 'REF_FREQ' ] )
    fractional_bandwidth = bandwidth / solution_frequency
  beta = float32( ni_table.keywords[ 'BETA' ] )
  r_0 = float32( ni_table.keywords[ 'R_0' ] )
  layer_count = int( ni_table.keywords[ 'LAYERS' ] )
  layer_heights = []
  layer_weights = []
  for l in range( layer_count ):
    layer_heights.append( float32( ni_table.keywords[ 'HEIGHT%d' % ( l + 1 ) ] ) )
    layer_weights.append( float32( ni_table.keywords[ 'WEIGHT%d' % ( l + 1 ) ] ) )
  l_max = layer_weights.index( max( layer_weights ) )
  if ( len( facet_list ) == 0 ):
    used_facet_list = range( 1, 1 + restore_parameter( facets, 'facet_count' ) )
  else:
    used_facet_list = facet_list
  
  ob_table = wiz_uv.table( 'OB', ni_table.version )
  if ( reference_antenna == 0 ):
    ref_ant = 1
  else:
    ref_ant = reference_antenna
  r = ref_ant - 1
  try:
    iterations = int( ob_table.keywords[ 'ITER' ] )
  except:
    iterations = 4
  
  # determine time grid for solutions
  time_list = get_time_list( uv )
  time_count = len( time_list )
  
  # read ionospheric fit table
  fit_time_list = []
  fit_count_list = []
  fit_poly_table = []
  fit_coef_table = []
  fit_weight_list = []
  for row in ni_table:
    fit_time_list.append( float32( row.time ) )
    fit_count_list.append( row.antenna_no )
    fit_poly_table.append( [ float32( pol ) for pol in row.coef[ 0 : poly_order ] ] )
    fit_coef_table.append( [ float32( coef ) for coef in row.coef[ poly_order : ] ] )
    fit_weight_list.append( float32( row.weight ) )
  fit_time_count = len( fit_time_list )
  fit_poly_table = array( fit_poly_table, dtype = float64 )
  fit_coef_table = array( fit_coef_table, dtype = float64 )
  
  # reject bad fits
  fit_weight_array = array( fit_weight_list, dtype = float64 )
  sel = awhere( fit_weight_array > 0. )
  fit_rms_array = 1. / aget( fit_weight_array, sel )
  sub_sel = awhere( fit_rms_array >= rms_min )
  sel = reshape( aget( sel.ravel(), sub_sel ), ( len( sub_sel ), 1 ) ) 
  if ( not rejection_limit is None ):
    fit_rms_array = 1. / aget( fit_weight_array, sel )
    len_rms = len( fit_rms_array ) + 1
    while( len( fit_rms_array ) != len_rms ):
      len_rms = len( fit_rms_array )
      rms_mean = fit_rms_array.mean()
      rms_std = sqrt( ( ( fit_rms_array - rms_mean )**2 ).mean() )
      sub_sel = awhere( fit_rms_array <= rms_mean + rejection_limit * rms_std )
      sel = reshape( aget( sel.ravel(), sub_sel ), ( len( sub_sel ), 1 ) ) 
      fit_rms_array = 1. / aget( fit_weight_array, sel )
    if print_info:
      print 'rejected %s of %s fits...' % ( repr( len( fit_weight_array ) - 
          len( sel ) ), repr( len( fit_weight_array ) ) )
      print 'phase rms mean = %s  stddev = %s' % ( repr( rms_mean ), 
          repr( rms_std ) )
  elif ( not rms_limit is None ):
    fit_rms_array = 1. / aget( fit_weight_array, sel )
    sub_sel = awhere( fit_rms_array <= rms_limit )
    sel = reshape( aget( sel.ravel(), sub_sel ), ( len( sub_sel ), 1 ) ) 
    rms_mean = aget( fit_rms_array, sub_sel ).mean()
    rms_std = sqrt( ( ( aget( fit_rms_array, sub_sel ) - rms_mean )**2 ).mean() )
    if print_info:
      print 'rejected %s of %s fits...' % ( repr( len( fit_weight_array ) - 
          len( sel ) ), repr( len( fit_weight_array ) ) )
      print 'phase rms mean = %s  stddev = %s' % ( repr( rms_mean ), 
          repr( rms_std ) )
  sel = sel.ravel().tolist()
  
  # read ionospheric pierce point table
  pierce_time_list = []
  pierce_X_table = []
  pierce_za_list = []
  pierce_index_table = []
  for row in ob_table:
    pierce_time_list.append( float32( row.time ) )
    pierce_X_table.append( [ float32( x ) for x in row.orbxyz[ 0 : 2 ] ] )
    pierce_za_list.append( float32( row.orbxyz[ 2 ] ) )
    pierce_index_table.append( [ int( round( row.orientation ) ) - 1, 
        row.subarray - 1, row.antenna_no - 1 ] )
  pierce_time_count = len( pierce_time_list )
  pierce_time_array = array( pierce_time_list, dtype = float64 )
  pierce_index_array = array( pierce_index_table, dtype = int64 )
  pierce_X_array = array( pierce_X_table, dtype = float64 )
  pierce_za_array = array( pierce_za_list, dtype = float64 )
  
  # get other relevant data
  calibration_data = get_phase_calibration_data( uv, facets, time_info = False, 
      source_info = True, antenna_info = True, calibration_info = False, 
      facet_list = used_facet_list, print_info = print_info )
  center_table = calibration_data[ 1 ]
  source_table = calibration_data[ 2 ]
  array_table = calibration_data[ 3 ]
  antenna_table = calibration_data[ 4 ]
  source_count = len( source_table )
  antenna_count = len( antenna_table )
  
  # generate time table
  fit_time_table = []
  fit_gst_list = get_gst_list( uv, fit_time_list )
  for n in range( fit_time_count ):
    fit_time_table.append( [ fit_time_list[ n ], fit_gst_list[ n ] ] )
  
  # generate reference antenna table
  reference_list = []
  for k in range( source_count ):
    reference_list.append( r )
  
  # generate solutions and write them to facets
  if print_info:
    print 'generating solutions ...'
  phase_table = - 360. * ones( ( time_count, source_count, antenna_count ), 
      dtype = float64 )
  weight_table = zeros( time_count, dtype = float64 )
  if include_delays:
    delay_table = zeros( ( time_count, source_count, antenna_count ), 
        dtype = float64 )
  
  for n in range( fit_time_count ):
    
    if ( len( time_steps ) > 0 ):
      if ( not n in time_steps ):
        continue
    
    if print_info:
      print '... time step n = %d / %d' % ( n + 1, fit_time_count )
    
    # only process non-rejected fits
    if ( not n in sel ):
      continue
    
    # save fit weight
    try:
      nn = time_list.index( fit_time_list[ n ] )
    except ValueError:
#      raise error( 'no match for solution time in uv time' )
      if print_info:
        print '...... no match for solution time in uv time'
      continue
    weight_table[ nn ] = 1. / radians( 1. / fit_weight_list[ n ] )
    
    # get model fit parameters
    poly = fit_poly_table[ n ]
    if ( poly_order == 2 ):
      poly = array( poly.tolist() + [ 0.,0.,0. ], dtype = poly.dtype )
    P = fit_coef_table[ n ]
    weight = fit_weight_list[ n ]
    
    if print_info:
      print '...... calculating base vectors'
    
    active_antennas = []
    Xpl_table = []
    pzal_table = []
    Bl_table = []
    for l in range( layer_count ):
      Xp_table = []
      pza_table = []
      
      # get pierce points 
#      for np in range( pierce_time_count ):
#        if ( ( pierce_time_list[ np ] == fit_time_list[ n ] ) and 
#             ( pierce_index_table[ np ][ 0 ] == l ) ):
#          Xp_table.append( pierce_X_table[ np ] )
#          pza_table.append( pierce_za_list[ np ] )
#          if ( not pierce_index_table[ np ][ 2 ] in active_antennas ):
#            active_antennas.append( pierce_index_table[ np ][ 2 ] )
      sel2 = awhere( ( pierce_time_array == fit_time_list[ n ] ) & 
          ( pierce_index_array[ : , 0 ] == l ) )
      Xp_table = Xp_table + aget( pierce_X_array, sel2 ).tolist()
      pza_table = pza_table + aget( pierce_za_array, sel2 ).tolist()
      for i in aget( pierce_index_array[ : , 2 ], sel2 ).tolist():
        if ( not i in active_antennas ):
           active_antennas.append( i )
      
      Xp_table = array( Xp_table, dtype = float64 )
      pza_table = array( pza_table, dtype = float64 )
      
      # calculate structure matrix
      p_count = len( Xp_table )
      if ( p_count != fit_count_list[ n ] ):
        raise error( 'pierce count does not match' )
      D_table = resize( Xp_table, ( p_count, p_count, 2 ) )
      D_table = transpose( D_table, ( 1, 0, 2 ) ) - D_table
      D_table = add.reduce( D_table**2, 2 )
      D_table = ( D_table / ( r_0**2 ) )**( beta / 2. )
      
      # calculate covariance matrix C
      # calculate partial product for interpolation B
      # reforce symmetry
      C_table = - D_table / 2.
      C_table = transpose( C_table - ( add.reduce( C_table, 0 ) / 
          float( p_count ) ) )
      B_table = add.reduce( C_table, 0 ) / float( p_count )
      C_table = C_table - B_table
      C_table = ( C_table + transpose( C_table ) ) / 2.
      
      # incorporate airmass functions and layer weights
      # add layer C to total C
      A_table = resize( 1. / cos( aradians( pza_table ) ), ( p_count, p_count ) )
      A_table = A_table * transpose( A_table )
      if ( l == 0 ):
        Cl_table = C_table * A_table * ( layer_weights[ l ]**2 )
      else:
        Cl_table = Cl_table + C_table * A_table * ( layer_weights[ l ]**2 )
      
      # save tables per height layer
      Xpl_table.append( Xp_table )
      pzal_table.append( pza_table )
      Bl_table.append( B_table )
    
    # convert to arrays
    Xpl_table = array( Xpl_table, dtype = float64 )
    pzal_table = array( pzal_table, dtype = float64 )
    Bl_table = array( Bl_table, dtype = float64 )
    
    # project out polynomial fit
    Xp_table = Xpl_table[ l_max ]
    pza_table = pzal_table[ l_max ]
    if ( poly_order == 2 ):
      Yt_table = array( [ Xp_table[ : , 0 ], Xp_table[ : , 1 ] ] )
    else:
      Yt_table = array( [ Xp_table[ : , 0 ], Xp_table[ : , 1 ], Xp_table[ : , 0 ]**2,
          Xp_table[ : , 0 ] * Xp_table[ : , 1 ], Xp_table[ : , 1 ]**2 ] )
    Yt_table = array( [ Xp_table[ : , 0 ], Xp_table[ : , 1 ] ] )
    Yt_table = Yt_table / cos( aradians( pza_table ) )
    Y_table = transpose( Yt_table )
    Yt_table = dot( linalg.inv( dot( Yt_table, Y_table ) ), Yt_table )
    Y_table = dot( Y_table, Yt_table )
    Y_table = diag( ones( ( p_count ), dtype = Y_table.dtype ) ) - Y_table
    # here, transpose( Y_table ) = Y_table
    Cl_table = dot( Y_table, dot( Cl_table, Y_table ) )
    
    # eigenvalue decomposition
    # reforce symmetry
    # select subset of base vectors
    [ U_table, S, Ut_table ] = linalg.svd( Cl_table )
    U_table = ( U_table + transpose( Ut_table ) ) / 2.
    U_table = U_table[ : , 0 : order ]
    S = S[ 0 : order ]
    
    # calculate interpolation matrix
    F_table = dot( U_table, P / S )
    
    if print_info:
      print '...... calculating pierce point coordinates'
    
    Xl_table = []
    zal_table = []
    ref_list = []
    ref_table = []
    for l in range( layer_count ):
      
      X_table = []
      za_table = []
      
      # get pierce point coordinates
      pierce_table = calculate_pierce_coordinates( fit_time_table[ n ], 
          center_table, source_table, array_table, antenna_table, 
          height = layer_heights[ l ], iterations = iterations )
      
      # put all new pierce points into one array
      j = 0
      for pierce_info in pierce_table:
        [ X, za, [ k, i ] ] = pierce_info
#        if ( i in active_antennas ):
        if True:
          X_table.append( X )
          za_table.append( za )
          if ( l == 0 ):
            if ( i == r ):
              ref_list.append( [ k, i, j ] )
            j = j + 1
        elif ( i == r ):
          raise error( 'no solution for reference antenna' )
      Xl_table.append( X_table )
      zal_table.append( za_table )
      
      # loop over pierce points
      # store index to reference antenna
      if ( l == 0 ):
        for pierce_info in pierce_table:
          [ X, pza, [ k, i ] ] = pierce_info
          j = ( [ ( k == ref[ 0 ] ) for ref in ref_list ] ).index( True )
          ref_table.append( ref_list[ j ][ 2 ] )
    
    Xl_table = array( Xl_table, dtype = float64 )
    zal_table = array( zal_table, dtype = float64 )
    
    if print_info:
      print '...... generating solutions'
    
    # calculate pierce point model solutions
    phi_poly_table = phi_poly_model( Xl_table[ l_max ], poly ) / cos( 
        aradians( zal_table[ l_max ] ) )
    phi_mkl_table = phi_mkl_model( layer_weights, Xl_table, zal_table, Xpl_table, 
        pzal_table, Bl_table, F_table, beta = beta, r_0 = r_0 )
    phi_table = phi_poly_table + phi_mkl_table
    
    # only store phase corrections of active antennas
    for j in range( len( pierce_table ) ):
      [ X, pza, [ k, i ] ] = pierce_table[ j ]
      if ( ( i in active_antennas ) or ( ( i + 1 ) in include_antennas ) ):
#        [ ref_X, ref_pza, [ ref_k, ref_i ] ] = pierce_table[ ref_table[ j ] ]
        if include_phases:
          phase_table[ nn, k, i ] = phi_table[ j ] - phi_table[ ref_table[ j ] ]
        else:
          phase_table[ nn, k, i ] = 0.
        if include_delays:
          dphase = phi_table[ j ] - phi_table[ ref_table[ j ] ]
          delay = - ( dphase / 360. ) / ( solution_frequency * 2. * pi )
          delay_table[ nn, k, i ] = delay
          dphase = 360. * delay * ( solution_frequency - reference_frequency )
          phase_table[ nn, k, i ] = phase_table[ nn, k, i ] - dphase
        phase_table[ nn, k, i ] = amodulo( phase_table[ nn, k, i ] + 180., 360. ) - 180.
  
  # write solutions to facets
  if print_info:
    print 'writing solution tables ...'
  
  delay = 0.
  for k in range( source_count ):
    solution_table = []
    for n in range( time_count ):
      
      # generate solution row per time step
      solution_row = [ [ time_list[ n ], ref_ant, 0., 0. ] ]
      weight = weight_table[ n ]
      ref_solution = r_phi_to_complex( [ 1., phase_table[ n, k, r ] ] )
      if include_delays:
        ref_delay = delay_table[ n, k, r ]
      for i in range( antenna_count ):
        if ( ( phase_table[ n, k, i ] < -300. ) or ( phase_table[ n, k, r ] < -300. ) ):
          solution_row.append( [ 0., 0., 0., 0. ] ) # flag vis
        else:
          solution = r_phi_to_complex( [ 1., phase_table[ n, k, i ] ] ) / ref_solution
          if include_delays:
            delay = delay_table[ n, k, i ] - ref_delay
          solution_row.append( [ solution.real, solution.imag, delay, weight ] )
      
      # add solution row to solution table
      solution_table.append( [ [ xx for xx in x ] for x in solution_row ] )
    
    # write solution table to facet
    if print_info:
      print '... writing solution table for source ', used_facet_list[ k ]
    facet = get_facet( facets, used_facet_list[ k ] )
    if ( solution_version == 0 ):
      sn_version = facet.table_highver( 'SN' ) + 1
    else:
      sn_version = solution_version
    write_solution_table( facet, solution_table, out_version = sn_version )
  
  return

###############################################################################

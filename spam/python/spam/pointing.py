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
try:
  import _pointing
except:
  __pointing = False
else:
  __pointing = True

###############################################################################

def fit_beam_offset_model( P, dojac = None, beam_offset_model = None, beam_parameters = None,
    rp_table = None, gain_zero_table = None, amplitude_table = None, error_table = None,
    normalize_weights = None ):
# X_table contains the source positions, with the field center rotated to zero RADEC

  # handle input parameters
  if ( ( not dojac is None ) or ( beam_offset_model is None ) or ( beam_parameters is None ) or 
      ( rp_table is None ) or ( gain_zero_table is None ) or ( amplitude_table is None ) or 
      ( error_table is None ) or ( normalize_weights is None ) ):
    return - 1, None, None

  # calculate gain ratios at source positions
  drdp = array( [ sqrt( ( P * P ).sum() ), degrees( atan2( P[ 0 ], P[ 1 ] ) ) ], dtype = P.dtype )
  gain_table = beam_offset_model( rp_table, drdp, beam_parameters )
  gain_table = gain_zero_table / gain_table

  # calculate chi2 terms
  chi_list = []
  if normalize_weights:
    normalize_list = []
  gain_count = len( gain_table )
  dlog_amplitude_table = log( amplitude_table )
  dlog_amplitude_table = resize( dlog_amplitude_table, ( gain_count, gain_count ) )
  dlog_amplitude_table = transpose( dlog_amplitude_table ) - dlog_amplitude_table
  dlog_gain_table = log( gain_table )
  dlog_gain_table = resize( dlog_gain_table, ( gain_count, gain_count ) )
  dlog_gain_table = transpose( dlog_gain_table ) - dlog_gain_table
  dweight_table = resize( error_table**2, ( gain_count, gain_count ) )
  dweight_table = transpose( dweight_table ) + dweight_table
  dweight_table = 1. / sqrt( dweight_table )
  chi_table = dweight_table * ( dlog_amplitude_table - dlog_gain_table )
  dummy = [ [ chi_list.append( chi_table[ b, a ] )
      for a in range( b + 1, gain_count ) ] for b in range( gain_count ) ]
  if normalize_weights:
    dummy = [ [ normalize_list.append( dweight_table[ b, a ] )
        for a in range( b + 1, gain_count ) ] for b in range( gain_count ) ]

  # make (normalized) chi2 array
  chi_array = array( chi_list, dtype = chi_table.dtype )
  if normalize_weights:
    normalize_array = array( normalize_list, dtype = dweight_table.dtype )
    factor = 1. / sqrt( ( normalize_array**2 ).sum() )
    chi_array = factor * chi_array

  return 0, chi_array, None

###############################################################################

# AIPS PBPARM beam model

def aips_beam_model( rp_table, drdp, beam_parameters ):
  freq = beam_parameters[ 0 ]
  pbparms = beam_parameters[ 1 : ].tolist()
  gain_table = []
  for rp in rp_table:
    radec = calculate_offset_position( [ 0., 0. ], rp[ 0 ], rp[ 1 ] )
    radec = calculate_offset_position( radec, drdp[ 0 ], drdp[ 1 ] )
    [ r, p ] = calculate_angular_separation( [ 0., 0. ], radec )
    gain_table.append( calculate_pbparm_attenuation( freq, r, pbparms ) )
  return array( gain_table, dtype = rp_table.dtype )

###############################################################################

def get_amplitude_calibration_data( uv, facets, facet_list = [], time_info = True,
    source_info = True, calibration_info = True, print_info = False, solution_version = - 1 ):

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
    center_table = get_radec( uv )
    source_table = [ get_radec( get_facet( facets, source ) ) for source in source_list ]

  # get antenna info
  antenna_list = get_antenna_positions( uv )
  antenna_list = [ antenna[ 0 ] for antenna in antenna_list ]
  antenna_count = antenna_list[ - 1 ]

  # read calibration info from peel facets
  amplitude_table = None
  error_table = None
  if calibration_info:
    if print_info:
      print 'getting calibration info ...'
    amplitude_array = - 1. * ones( shape = ( time_count, source_count, antenna_count ), dtype = float64 )
    error_array = - 1. * ones( shape = ( time_count, source_count, antenna_count ), dtype = float64 )
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
      for solution in solution_table:
        [ time, reference_antenna, dummy, dummy2 ] = solution[ 0 ]
        try:
          n = time_list.index( time )
        except:
          raise error( 'time in solution table does not correspond to time in observations' )
        for i in range( antenna_count ):
          [ gain_real, gain_imag, delay, weight ] = solution[ i + 1 ]
          gain = complex( gain_real, gain_imag )
          if ( ( weight == 0. ) or 
              ( ( i != reference_antenna - 1 ) and ( gain == complex( 1., 0. ) ) ) or 
              ( gain == complex( 0., 0. ) ) ):
            amp = - 1.
            amp_error = - 1.
          else:
            [ amp, phase ] = complex_to_r_phi( gain )
            amp_error = 1. / weight # TODO: need something better than this
          amplitude_array[ n ][ k ][ i ] = amp
          error_array[ n ][ k ][ i ] = amp_error
      k = k + 1
    amplitude_table = amplitude_array.tolist()
    error_table = error_array.tolist()

  return [ time_table, center_table, source_table, amplitude_table, error_table ]

###############################################################################

def fit_pointing_model( uv, facets, facet_list = [], time_steps = [], equal_weights = False,
    normalize_weights = True, movie_file_name = '', format = 'mpg', xy_range = [ - 6., 6. ],
    double_precision = True, beam_cutoff = 0.3, solution_version = 0, fit_version = 0,
    print_info = False ):

  if double_precision:
    dtype = float64
  else:
    dtype = float32

  make_movie = False
  if ( movie_file_name != '' ):
    make_movie = True
    image_file_name_list = []

  # select appropriate model functions
  # initialize model parameters
  model = 'aips_beam'
  beam_offset_model = aips_beam_model
  P = zeros( shape = ( 2 ), dtype = dtype )
  pbparm3 = restore_parameter( uv, 'pbparm3' )
  pbparm4 = restore_parameter( uv, 'pbparm4' )
  pbparm5 = restore_parameter( uv, 'pbparm5' )
  pbparm6 = restore_parameter( uv, 'pbparm6' )
  pbparm7 = restore_parameter( uv, 'pbparm7' )
  beam_parameters = array( [ get_central_frequency( uv ), beam_cutoff, 1., 
      pbparm3, pbparm4, pbparm5, pbparm6, pbparm7 ], dtype = dtype )

  if print_info:
    print 'starting model fits using model = ' + model

  # get all relevant data
  calibration_data = get_amplitude_calibration_data( uv, facets, facet_list = facet_list,
      solution_version = solution_version, print_info = print_info )
  time_table = calibration_data[ 0 ]
  center_table = calibration_data[ 1 ]
  source_table = calibration_data[ 2 ]
  amplitude_array = calibration_data[ 3 ]
  error_array = calibration_data[ 4 ]
  source_count = len( amplitude_array[ 0 ] )
  antenna_count = len( amplitude_array[ 0 ][ 0 ] )

  fit_table = []

  # rotate source positions to center RADEC = [ 0., 0. ]
  # calculate zero offset gains
  # only use sources that are within the beam cutoff
  rp_all_table = []
  for k in range( source_count ):
    rp = calculate_angular_separation( center_table, source_table[ k ] )
    rp_all_table.append( rp )
  rp_all_array = array( rp_all_table, dtype = dtype )
  gain_zero_all_array = beam_offset_model( rp_all_array, P, beam_parameters )
  source_sel = awhere( gain_zero_all_array > 0. ).ravel().tolist()
  source_count = len( source_sel )
  if print_info:
    print 'using %s sources' % ( repr( source_count ) )

  # loop over time stamps
  if ( len( time_steps ) == 0 ):
    n_list = range( len( time_table ) )
  else:
    n_list = time_steps
  for n in n_list:

    if print_info:
      print '... time step n = ' + repr( n + 1 )

    for i in range( antenna_count ):

      if print_info:
        print '...... antenna i = ' + repr( i + 1 )

      # containers for model input data
      rp_table = []
      gain_zero_table = []
      amplitude_table = []
      error_table = []

      # containers for plot data
      if make_movie:
        x = []
        y = []
        z = []

      for k in source_sel:
        if ( error_array[ n ][ k ][ i ] > 0. ):
          rp_table.append( rp_all_array[ k ] )
          gain_zero_table.append( gain_zero_all_array[ k ] )
          amplitude_table.append( amplitude_array[ n ][ k ][ i ] )
          error_table.append( error_array[ n ][ k ][ i ] )

        # store plot data
        if make_movie:
          rp = rp_all_table[ k ]
          x.append( rp[ 0 ] * sin( radians( rp[ 1 ] ) ) )
          y.append( rp[ 0 ] * cos( radians( rp[ 1 ] ) ) )
          z.append( amplitude_array[ n ][ k ][ i ] )

      # skip times of poor data
      if ( len( rp_table ) < 3 ):
        if print_info:
          print '......... skipping due to no / too few data points'
        fit_table.append( [ time_table[ n ][ 0 ], i, azeros( P ).tolist(), 0. ] )
        continue

      # convert tables to arrays for fitting routines
      rp_table = array( rp_table, dtype = dtype )
      gain_zero_table = array( gain_zero_table, dtype = dtype )
      amplitude_table = array( amplitude_table, dtype = dtype )
      error_table = array( error_table, dtype = dtype )
      if equal_weights:
        error_table = ones( shape( error_table ), dtype = dtype )

      # set initial guess to zero
      P = azeros( P )

      # define data passed to the model fit routines
      # define model parameter structure
      function_keywords = { 'beam_offset_model' : beam_offset_model, 'beam_parameters' : beam_parameters,
          'rp_table' : rp_table, 'gain_zero_table' : gain_zero_table, 'amplitude_table' : amplitude_table,
          'error_table' : error_table, 'normalize_weights' : normalize_weights }

      # specify fit parameters
      parameter_info = []
      for m in range( len( P ) ):
        par_info = { 'parname' : 'P_%d' % ( m ), 'value' : P[ m ], 'limits' : [ None, None ] }
        parameter_info.append( par_info )

      # fit model to data points
      fit = mpfit( fit_beam_offset_model, functkw = function_keywords, parinfo = parameter_info,
          quiet = True, autoderivative = True, debug = False, fastnorm = False, nocovar = True,
          dblprec = double_precision )
#      if print_info:
#        print '...... mpfit() returned status %d' % ( fit.status ) 
      if ( fit.errmsg != '' ):
        raise error( fit.errmsg )
      P = fit.params.copy()
      chi2 = fit.fnorm

      # store fit results
      weight = 1. / sqrt( chi2 )
      fit_table.append( [ time_table[ n ][ 0 ], i, P.tolist(), weight ] )
      if print_info:
        print '......... fitted parameters are %s' % repr( P.tolist() )
        print '......... chi-2 value is %s' % repr( chi2 )

  # write fit results to model fit table
  write_beam_offset_fit_table( uv, fit_table, model = model, out_version = fit_version )

  # make movie from image files
  if make_movie:
    make_movie_from_images( image_file_name_list, movie_file_name, format = format )

  return










if False:

    # calculate mean absolute difference between measured phase and model phase per antenna per source
    antenna_error_table = azeros( error_table )
    for ref in ref_list:
      index_list = []
      for j in range( len( ref_table ) ):
        if ( ref_table[ j ] == ref ):
          index_list.append( j )
      ref_count = len( index_list )
      dphase_table = array( [ phase_table[ index_list[ a ] ] for a in range( ref_count ) ],
          dtype = phase_table.dtype )
      dphase_table = resize( dphase_table, ( ref_count, ref_count ) )
      dphase_table = transpose( dphase_table ) - dphase_table
      dphi_table = array( [ phi_table[ index_list[ a ] ] for a in range( ref_count ) ],
          dtype = phi_table.dtype )
      dphi_table = resize( dphi_table, ( ref_count, ref_count ) )
      dphi_table = transpose( dphi_table ) - dphi_table
      dweight_table = array( [ error_table[ index_list[ a ] ] for a in range( ref_count ) ],
          dtype = error_table.dtype )
      dweight_table = resize( dweight_table**2, ( ref_count, ref_count ) )
      dweight_table = transpose( dweight_table ) + dweight_table
      dweight_table = 1. / sqrt( dweight_table )
      chi_table = add.reduce( dweight_table * 
          abs( amodulo( ( dphase_table - dphi_table ) + 180., 360. ) - 180. ), 1 )
      normalize_table = add.reduce( dweight_table, 1 ) - diagonal( dweight_table )
      chi_table = chi_table / normalize_table
      for a in range( ref_count ):
        antenna_error_table[ index_list[ a ] ] = chi_table[ a ]

    # estimate mean error per antenna
    if remove_average:
      # determine average error per antenna
      for i in range( len( antenna_table ) ):
        # select all phases related to single antenna
        sel = awhere( info_table[ : , 1 ] == i )
        # TODO: determine reasonable number of values for calculating the antenna mean
        if ( len( sel ) == 0 ):
          continue

        # determine error average
        # better not use rms, as we're looking for systematic large errors
        mean_error = aget( antenna_error_table, sel ).mean()


# TODO: include movie making in code above


    if make_movie:

      if ( not plot_gradient ):
        Q = array( [ 0., 0. ], dtype = dtype )
        R = P
      else:
        # extract gradient
        F_table = dot( U_table, Si * P )
        offset = phi_klmap_model( array( [ 0., 0. ], dtype = dtype ), Xp_table, B_table, F_table,
            beta = beta, r_0 = r_0 )
        Q = dot( U_table, P ) - offset
        Q = dot( transpose( Xp_table ), Q )
        Q2 = dot( transpose( Xp_table ), Xp_table )
        Q = dot( linalg.inv( Q2 ), Q )
        # remove gradient from base parameters
        R = dot( Xp_table, Q )
#        R = dot( Ut_table, R )[ 0 : order ]
        R = dot( transpose( U_table ), R )
        R = P - R
        # calculate new interpolation matrix
        F_table = dot( U_table, Si * R )

        if print_info:
          print '... extracted gradient is %s' % repr( Q.tolist() )

      # calculate Earth referenced coordinates of pierce points of array center towards pointing center
      [ center_pxyz, center_pza ] = calculate_pierce_point( array_table[ 0 ], center_table, 
          time_table[ n ][ 1 ], height = height )
      if ( not include_airmass ):
        center_pza = 0.
      center_p_geo_llh = xyz_to_geo_llh( center_pxyz )

      # determine grid of ionospheric phases from model
      delta_xy = ( xy_range[ 1 ] - xy_range[ 0 ] ) / 100.
      xx = arange( xy_range[ 0 ], xy_range[ 1 ] + delta_xy, delta_xy, dtype = Xp_table.dtype )
      yy = arange( xy_range[ 0 ], xy_range[ 1 ] + delta_xy, delta_xy, dtype = Xp_table.dtype )
      X_table = array( [ [ [ xx[ a ], yy[ b ] ]
        for a in range( len( xx ) ) ] for b in range( len( yy ) ) ], dtype = Xp_table.dtype )
      X_table = reshape( X_table, ( len( yy ) * len( xx ), 2 ) )
      PHI = phi_klmap_model( X_table, Xp_table, B_table, F_table, beta = beta, r_0 = r_0 )
      PHI = reshape( PHI, ( len( yy ), len( xx ) ) )
      if plot_gradient:
        PHI[ 0 : 11, 0 : 11 ] = array( [ [ Q[ 0 ] * xx[ 45 + a ] + Q[ 1 ] * yy[ 45 + b ]
            for a in range( 11 ) ] for b in range( 11 ) ], dtype = phase_table.dtype )
      MODEL_PHASE = amodulo( ( PHI / cos( radians( center_pza ) ) ) + 180., 360. ) - 180.

      # plot phase grid below measured phases
      clf()
      imshow( MODEL_PHASE, extent = ( xx.min(), xx.max(), yy.min(), yy.max() ),
          interpolation = 'nearest', vmin = - 180., vmax = 180., origin = 'lower' )
      axis( xy_range + xy_range, 'scaled' )
      hsv()
      hold( True )

      # loop over pierce points
      for j in range( len( Xp_table ) ):

        # get measured antenna phases
        phase = phase_table[ j ]
        error = error_table[ j ]
        X = Xp_table[ j ]
        pza = pza_table[ j ]
        ref_phase = phase_table[ ref_table[ j ] ]
        ref_X = Xp_table[ ref_table[ j ] ]
        ref_pza = pza_table[ ref_table[ j ] ]

        # remove reference phase and phase gradient from measured phases
        phi = ( X * Q ).sum() / cos( radians( pza ) )
        ref_phi = ( ref_X * Q ).sum() / cos( radians( ref_pza ) )
        phase = amodulo( ( ( phase - ref_phase ) - ( phi - ref_phi ) ) + 180., 360. ) - 180.

        # remove offset to background
        offset = phi_klmap_model( ref_X, Xp_table, B_table, F_table, beta = beta, r_0 = r_0 )
        offset = offset / cos( radians( center_pza ) )
        z[ j ] = amodulo( ( phase + offset ) + 180., 360. ) - 180.

      # scale point sizes with antenna errors
      error_range = [ 0., 60. ] # degrees
      scale_range = [ 5., 40 ] # pixels
      factor = ( scale_range[ 1 ] - scale_range[ 0 ] ) / ( error_range[ 1 ] - error_range[ 0 ] )
      e = antenna_error_table
      e = aput( e, awhere( ( e > 0. ) & ( e < error_range[ 0 ] ) ), error_range[ 0 ] )
      e = aput( e, awhere( e > error_range[ 1 ] ), error_range[ 1 ] )
      e = factor * ( e - error_range[ 0 ] ) + scale_range[ 0 ]
      e = aput( e, awhere( e < scale_range[ 0 ] ), 0. )

      # plot pierce points with phases
      scatter( x, y, e, z, vmin = - 180., vmax = 180. )
      axis( [ xy_range[ 0 ], xy_range[ 1 ], xy_range[ 0 ], xy_range[ 1 ] ], 'scaled' )
      hsv()
      xlabel( r'$\Delta$longitude [deg]' )
      ylabel( r'$\Delta$latitude [deg]' )
      title( r'n = %5d,  $\sigma_{\textrm{phase}}$ = %7.3f deg' % ( n, 1. / weight ) )
      cb = colorbar()
      cb.ax.set_ylabel( r'phase [deg]' )
#      show()

      # save plot to image file
      frame_file_name = movie_file_name + '_%04d.png' % ( n )
      if file_exists( frame_file_name ):
        remove_file( frame_file_name )
      savefig( frame_file_name, dpi = 100 )

      # store image filename
      image_file_name_list.append( frame_file_name )



###############################################################################

def write_beam_offset_fit_table( uv, fit_table, model = 'polynomial', 
    out_version = 0, **keywords ):

  # create new NI table
  num_coef = len( fit_table[ 0 ][ 2 ] )
  new_ni_table = new_table( uv, 'NI', out_version, num_coef = num_coef )
  row = new_table_row( new_ni_table )
  integration_time = restore_parameter( uv, 'integration_time' )
  time_interval = integration_time / ( 24. * 60. * 60. )

  # loop over 
  for fit_row in fit_table:
    [ time, i, coefs, weight ] = fit_row
    row.time = time
    row.time_interval = time_interval
    row.antenna_no = i + 1
    row.source_id = 0
    row.subarray = 0
    row.weight = weight
    row.coef = coefs
    new_ni_table.append( row )

  # add additional keywords and close new table
  new_ni_table.keywords[ 'MODEL' ] = model
  for key in keywords.keys():
    new_ni_table.keywords[ key ] = keywords[ key ]
  new_ni_table.close()

  return

###############################################################################

def generate_solution_tables_from_beam_offset_fit_table( uv, facets, fit_version = 0,
    solution_version = 0, print_info = False, rejection_limit = None, beam_cutoff = 0.1,
    reference_antenna = 0 ):

  wiz_uv = wizardry( uv )
  ni_table = wiz_uv.table( 'NI', fit_version )
  num_coef = ni_table.keywords[ 'NUM_COEF' ]
  model = ni_table.keywords[ 'MODEL' ]
  model = model.strip()
  if ( model != 'aips_beam'[ 0 : 8 ] ):
    raise error( 'unknown model: %s' % ( model ) )
  beam_offset_model = aips_beam_model
  if ( solution_version == 0 ):
    sn_version = facets.table_highver( 'SN' ) + 1
  else:
    sn_version = solution_version
  delay = 0.
  if ( reference_antenna == 0 ):
    ref_ant = 1
  else:
    ref_ant = reference_antenna

  pbparm3 = restore_parameter( uv, 'pbparm3' )
  pbparm4 = restore_parameter( uv, 'pbparm4' )
  pbparm5 = restore_parameter( uv, 'pbparm5' )
  pbparm6 = restore_parameter( uv, 'pbparm6' )
  pbparm7 = restore_parameter( uv, 'pbparm7' )
  beam_parameters = array( [ get_central_frequency( uv ), beam_cutoff, 1., 
      pbparm3, pbparm4, pbparm5, pbparm6, pbparm7 ], dtype = float64 )

  # determine time grid for solutions
  time_list = get_time_list( uv )
  time_count = len( time_list )

  # read beam offset fit table
  fit_time_list = []
  fit_antenna_list = []
  fit_coef_table = []
  fit_weight_list = []
  for row in ni_table:
    fit_time_list.append( row.time )
    fit_antenna_list.append( row.antenna_no - 1 )
    fit_coef_table.append( row.coef )
    fit_weight_list.append( row.weight )
  fit_coef_array = array( fit_coef_table, dtype = float64 )
  fit_weight_array = array( fit_weight_list )
#  fit_time_count = len( fit_time_list )
  antenna_count = max( fit_antenna_list ) + 1

  # reject bad fits
  sel = awhere( fit_weight_array > 0. )
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
      print 'rejected %s of %s fits...' % ( repr( len( fit_weight_array ) - len( sel ) ),
          repr( len( fit_weight_array ) ) )
      print 'phase rms mean = %s  stddev = %s' % ( repr( rms_mean ), repr( rms_std ) )
  sel = sel.ravel().tolist()

  # get other relevant data
  calibration_data = get_amplitude_calibration_data( uv, facets, time_info = False, source_info = True,
    calibration_info = False, print_info = print_info )
  center_table = calibration_data[ 1 ]
  source_table = calibration_data[ 2 ]
  source_count = len( source_table )

  # generate solutions and write them to facets
  if print_info:
    print 'generating solutions ...'
  amplitude_table = zeros( ( time_count, source_count, antenna_count ), dtype = float64 )
  weight_table = zeros( ( time_count, antenna_count ), dtype = float64 )

  # rotate source positions to center RADEC = [ 0., 0. ]
  # calculate zero offset gains
  rp_table = []
  for k in range( source_count ):
    rp = calculate_angular_separation( center_table, source_table[ k ] )
    rp_table.append( rp )
  rp_array = array( rp_table, dtype = float64 )
  gain_zero_array = beam_offset_model( rp_array, zeros( ( 2 ), dtype = float64 ), beam_parameters )
  # note that gains below the cutoff have gain_zero = 0
  # these will get gain = 1. in the end

  last_time = - 1000000.
  for f in range( len( fit_time_list ) ):

    # only process non-rejected fits
    if ( f in sel ):

      fit_time = fit_time_list[ f ]
      if ( fit_time != last_time ):
        last_time = fit_time
        try:
          nn = time_list.index( fit_time )
        except ValueError:
          raise error( 'no match for solution time in uv time' )
        if print_info:
          print '... time step n = ', nn

      if print_info:
        print '...... generating solutions'

      # get model fit parameters
      i = fit_antenna_list[ f ]
      P = fit_coef_array[ f ]
      w = fit_weight_array[ f ]

      # calculate gain ratios at source positions
      drdp = array( [ sqrt( ( P * P ).sum() ), degrees( atan2( P[ 0 ], P[ 1 ] ) ) ], dtype = P.dtype )
      gain_array = beam_offset_model( rp_array, drdp, beam_parameters )

      # replace gains below cutoff with 1.
      gain_sel = awhere( ( gain_zero_array == 0. ) | ( gain_array == 0. ) )
      gain_array = gain_zero_array / gain_array
      gain_array = aput( gain_array, gain_sel, 1. )

      amplitude_table[ nn, : , i ] = gain_array
      weight_table[ nn, i ] = 1. / radians( 1. / w )

  # write solutions to facets
  if print_info:
    print 'writing solution tables ...'

  for k in range( source_count ):
    solution_table = []
    for n in range( time_count ):

      # generate solution row per time step
      solution_row = [ [ time_list[ n ], ref_ant, 0., 0. ] ]
      for i in range( antenna_count ):
        if ( amplitude_table[ n, k, i ] == 0. ):
          solution_row.append( [ 1., 0., 0., 0. ] ) # gain 1, but with zero weight
        else:
          solution_row.append( [ amplitude_table[ n, k, i ], 0., 0., weight_table[ n, i ] ] )

      # add solution row to solution table
      solution_table.append( solution_row )

    # write solution table to facet
    if print_info:
      print '... writing solution table for source k = ', k
    facet = get_facet( facets, k + 1 )
    write_solution_table( facet, solution_table, out_version = sn_version )

  return

###############################################################################


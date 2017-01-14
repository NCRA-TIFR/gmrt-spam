###############################################################################

# import Python modules
from sys import *
from os import *
from math import *
import pdb

# import 3rd party modules
from numpy import *

# import user modules
from files import *
from aips import *
from solutions import *
from sphere import *
from parameter import *
from skymodel import *
from error import *
from flag import *
from calibrate import *

###############################################################################

def prepare_uv_data( uv, stokes = 'I', solution_version = None, flag_version = None, 
    diameter = None, field_factor = 1.6 ):
  
  # TODO: put in primary beam parameters at different frequencies for GMRT, VLA and WSRT (and others ??)
  prep_uv = get_aips_file( uv.disk, uv.name, 'PREP', - 1, 'UV' )
  
  # do some data checks
  uv = get_aips_file( uv.disk, uv.name, uv.klass, uv.seq, 'UV' )
  if ( not uv.exists() ):
    raise error( 'AIPS UV data %s not found' % ( str( uv ) ) )
  if ( uv.header.object == 'MULTI' ):
    raise error( 'only single source UV data is currently supported' )
  try:
    if_index = uv.header.ctype.index( 'IF' )
  except:
    pass
  else:
    if ( uv.header.naxis[ if_index ] != 1 ):
      raise error( 'only single IF UV data is currently supported' )
#  if uv.header.naxis[ uv.header.ctype.index( 'STOKES' ) ] != 1:
#    raise error( 'only single stokes UV data is currently supported' )
  
  # convert epoch B1950. to J2000.
  epoch = get_epoch( uv )
  if ( epoch == 1950. ):
    call_aips_task( 'UVFIX', indata = uv, outdata = prep_uv )
  elif ( epoch == 2000. ):
    call_aips_task( 'MOVE', indata = uv, outdata = prep_uv, userid = get_aips_userid() )
  else:
    raise error( 'AIPS UV data has unknown epoch %s' % ( repr( epoch ) ) )
  
  # apply flags
  if ( not flag_version is None ):
    if table_exists( prep_uv, 'FG', flag_version ):
      flag_uv = apply_flag_table( prep_uv, version = flag_version )
      prep_uv.zap()
      flag_uv.rename( name = prep_uv.name, klass = prep_uv.klass, seq = prep_uv.seq )
  
  # apply solutions
  if ( not solution_version is None ):
    if table_exists( prep_uv, 'SN', solution_version ):
      cal_uv = apply_solution_table( prep_uv, version = solution_version )
      prep_uv.zap()
      cal_uv.rename( name = prep_uv.name, klass = prep_uv.klass, seq = prep_uv.seq )
  
  # convert to stokes
  if ( stokes != '' ):
    stokes_uv = get_aips_file( uv.disk, uv.name, stokes, - 1, 'UV' )
    call_aips_task( 'SPLIT', indata = prep_uv, stokes = stokes, douvcomp = 0,
        outdisk = stokes_uv.disk, outclass = stokes_uv.klass, outseq = stokes_uv.seq )
    prep_uv.zap()
    stokes_uv.rename( name = prep_uv.name, klass = prep_uv.klass, seq = prep_uv.seq )
  
  # sort data to time-baseline order
  sort_uv = get_aips_file( uv.disk, uv.name, 'SORT', - 1, 'UV' )
  call_aips_task( 'UVSRT', indata = prep_uv, outdata = sort_uv, sort = 'TB' )
  prep_uv.zap()
  sort_uv.rename( name = prep_uv.name, klass = prep_uv.klass, seq = prep_uv.seq )
  
  # determine several useful parameters
  scale_table = array( get_frequency_list( prep_uv ) ) / get_frequency( prep_uv )
  uvw_table = []
  time_table = []
  weight_table = []
  wiz_uv = wizardry( prep_uv )
  for group in wiz_uv:
    uvw_table.append( [ uvw for uvw in group.uvw ] )
    if time_table == []:
      time_table.append( group.time )
    elif group.time != time_table[ -1 ]:
      time_table.append( group.time )
    vis = array( group.visibility, dtype = float32 )
    weights = transpose( vis )[ 2 ].ravel().sum()
    weight_table.append( float( weights ) )
  uvw_table = array( uvw_table, dtype = float32 )
  time_table = array( time_table, dtype = float32 )
  weight_table = array( weight_table, dtype = float32 )
  delta_time_table = time_table[ 1 : ] - time_table[ : -1 ]
  integration_time = 24. * 60. * 60. * median( delta_time_table ) # seconds
  uv2_table = array( [ ( uv**2 ).sum() for uv in uvw_table[ : , 0 : 2 ] ], dtype = float32 )
  uv_max = scale_table.max() * sqrt( uv2_table.max() ) # lambda
  uv_min = scale_table.min() * sqrt( uv2_table.min() ) # lambda
  w_table = array( uvw_table[ : , 2 ], dtype = float32 )
  w_max = scale_table.max() * ( abs( w_table ) ).max()
  time_count = len( time_table )
  time_min = time_table[ 0 ]
  time_max = time_table[ -1 ]
  visibility_count = int( len( uvw_table ) * get_channel_count( prep_uv ) )
  sum_weights = weight_table.sum()
#  beam_size = float( around( 1.22 * 3600. * degrees( 1. / uv_max ), 1 ) ) # arcsec
  beam_size = float( around( 1. * 3600. * degrees( 1. / uv_max ), 1 ) ) # arcsec
  cell_size = around( beam_size / 4., 1 ) # arcsec / pixel
  
  # get array specific parameters
  if ( not diameter is None ):
    dish_diameter = diameter
    pbparm = [ 0., 0., 0., 0., 0., 0., 0. ]
  elif ( prep_uv.header.telescop.find( 'GMRT' ) != -1  ):
    dish_diameter = 45.
    frequency = get_central_frequency( prep_uv )
    if ( frequency > 140.e6 ) and ( frequency < 170.e6 ):
      pbparm = [ 0.3, 1., -4.04, 76.2, -68.8, 22.03, 0. ]
    elif ( frequency > 200.e6 ) and ( frequency < 280.e6 ):
      pbparm = [ 0.3, 1., -3.366, 46.159, -29.963, 7.529, 0. ]
    elif ( frequency > 300.e6 ) and ( frequency < 350.e6 ):
      pbparm = [ 0.3, 1., -3.397, 47.192, -30.931, 7.803, 0. ] 
    elif ( frequency > 30.e6 ) and ( frequency < 80.e6 ):
      pbparm = [ 0.3, 1., -4., 60., -50., 15., 0. ]
    elif ( frequency > 580.e6 ) and ( frequency < 650.e6 ):
      pbparm = [ 0.3, 1., -3.486, 47.749, -35.203, 10.399, 0. ] 
    elif ( frequency > 1.e9 ) and ( frequency < 1.8e9 ):
      pbparm = [ 0.3, 1., -2.27961, 21.4611, -9.7929, 1.80153, 0. ] 
    else:
      raise error( 'beam shape of telescope at %d MHz is unknown' % (
          int( frequency / 1.e6 ) ) )
  elif ( prep_uv.header.telescop.find( 'VLA' ) != -1  ):
    dish_diameter = 25.
    frequency = get_central_frequency( prep_uv )
    if ( frequency > 65.e6 ) and ( frequency < 85.e6 ):
      pbparm = [ 0.3, 1., -0.897, 2.71, -0.242, 0., 0. ]
    elif ( frequency > 200.e6 ) and ( frequency < 500.e6 ):
      pbparm = [ 0.3, 1., -0.935, 3.23, -0.378, 0., 0. ]
    else:
      raise error( 'beam shape of telescope at %d MHz is unknown' % (
          int( frequency / 1.e6 ) ) )
  elif ( prep_uv.header.telescop.find( 'WSRT' ) != -1  ):
    dish_diameter = 25.
    frequency = get_central_frequency( prep_uv )
    # primary beam shape = cos(C*nu*r)^6, nu in GHz and r in degrees
    if ( frequency > 1.30e9 ) and ( frequency < 1.60e9 ):
      # C = 68.
      pbparm = [ 0.3, 1., -1.070370370, 5.091967688, -1.423131628, 0.2608246907,
          -0.03324696029 ]
    else:
      raise error( 'beam shape of telescope at %d MHz is unknown' % (
          int( frequency / 1.e6 ) ) )
  elif ( ( prep_uv.header.telescop.find( 'LOFAR36' ) != -1  ) or 
         ( prep_uv.header.telescop.find( 'LOFAR49' ) != -1  ) ):
    dish_diameter = 25.
    pbparm = [ 0., 0., 0., 0., 0., 0., 0. ]
  elif ( prep_uv.header.telescop.find( 'LOFAR' ) != -1  ):
#    dish_diameter = 50.
#    pbparm = [ 0., 0., 0., 0., 0., 0., 0. ]
    frequency = get_central_frequency( prep_uv )
    if ( frequency > 10.e6 ) and ( frequency < 90.e6 ):
      dish_diameter = 80.
      pbparm = [ 0.3, 1., -4.04, 76.2, -68.8, 22.03, 0. ]
    elif ( frequency > 110.e6 ) and ( frequency < 300.e6 ):
      dish_diameter = 30.
      pbparm = [ 0.3, 1., -0.935, 3.23, -0.378, 0., 0. ]
    else:
      raise error( 'beam shape of telescope at %d MHz is unknown' % (
          int( frequency / 1.e6 ) ) )
  else:
    raise error( 'telescope array %s not supported' % ( prep_uv.header.telescop.strip() ) )
  field_size = float( around( field_factor * 
      degrees( get_central_wavelength( prep_uv ) / dish_diameter ), 2 ) )

  # store parameters with UV data
  store_parameter( prep_uv, 'visibility_count', visibility_count )
  store_parameter( prep_uv, 'sum_weights', sum_weights )
  store_parameter( prep_uv, 'uv_min', uv_min )
  store_parameter( prep_uv, 'uv_max', uv_max )
  store_parameter( prep_uv, 'w_max', w_max )
  store_parameter( prep_uv, 'time_count', time_count )
  store_parameter( prep_uv, 'integration_time', integration_time )
  store_parameter( prep_uv, 'time_min', time_min )
  store_parameter( prep_uv, 'time_max', time_max )
  store_parameter( prep_uv, 'beam_size', beam_size )
  store_parameter( prep_uv, 'cell_size', cell_size )
  store_parameter( prep_uv, 'dish_diameter', dish_diameter )
  store_parameter( prep_uv, 'field_size', field_size )
  store_parameter( prep_uv, 'pbparm3', pbparm[ 2 ] )
  store_parameter( prep_uv, 'pbparm4', pbparm[ 3 ] )
  store_parameter( prep_uv, 'pbparm5', pbparm[ 4 ] )
  store_parameter( prep_uv, 'pbparm6', pbparm[ 5 ] )
  store_parameter( prep_uv, 'pbparm7', pbparm[ 6 ] )
  
  return prep_uv

###############################################################################

def calibrate_pb( uv, facets, sigma = 0., signal_to_noise = 10.,
    reference_antenna = 0, conversion_method = 'DFT', do_amplitude = False,
    amplitude_interval = 5., phase_interval_min = 1. / 60.,
    apply_solutions = False, calib_params = {}, interpolation_method = '',
    snr_limit = 2., print_info = True, sidelobe_rejection = 0.,
    frequency_correction = False, normalize_gains = True ):
  
  print 'WARNING: calibrate_pb() is obsolete, please use calibrate_model() instead'
  
  # determine total flux in model
  facet_count = restore_parameter( facets, 'facet_count' )
  total_model_flux = 0.
  for i in range( 1, facet_count + 1 ):
    facet = get_facet( facets, i )
    total_model_flux = total_model_flux + get_model_flux( facet )
  
  # calculate noise per integration interval
  time_count = len( get_time_list( uv ) )
  try:
    cpb_noise = restore_parameter( uv, 'cpb_noise' )
  except:
    if ( sigma > 0. ):
      raise error( 'cpb_noise is not yet defined' )
    cpb_noise = 0.
  noise_per_interval = cpb_noise * sqrt( float( time_count ) )
  
  # calculate solution interval
  integration_time = restore_parameter( uv, 'integration_time' ) / 60.
  sn_per_interval = total_model_flux / noise_per_interval
  phase_interval = ( signal_to_noise / sn_per_interval ) * integration_time
  interval_count = ceil( signal_to_noise / sn_per_interval )
  if ( phase_interval < phase_interval_min ):
    solution_interval = phase_interval_min
  else:
    solution_interval = integration_time * interval_count
  if ( not do_amplitude ):
    store_parameter( uv, 'solution_interval', solution_interval )
    if print_info:
#      print '... interval_count = %s' % ( repr( interval_count ) )
      print '... calibration solution interval = %s sec' % ( 
          repr( 60. * solution_interval ) )
  
  # calibrate UV data
  calibrate_model( uv, facets, reference_antenna = reference_antenna,
      phase_interval = solution_interval, do_amplitude = do_amplitude,
      amplitude_interval = amplitude_interval, keep_solutions = True, 
      interpolation_method = interpolation_method, snr_limit = snr_limit,
      apply_solutions = apply_solutions, conversion_method = conversion_method,
      calib_params = calib_params, sigma = sigma, sidelobe_rejection = sidelobe_rejection,
      frequency_correction = frequency_correction, normalize_gains = normalize_gains )
  
  return

###############################################################################

def fit_gaussian_to_peak( facet, pos = [], offset_ratio_max = 0.5, double_ratio_min = 0.25,
    return_double_fit = False ):

  fit_results = []

  if ( pos == [] ):
    [ max_flux, [ max_x, max_y ] ] = get_image_extremum( facet, force_positive = True )
  else:
    [ max_x, max_y ] = [ int( around( pos[ 0 ] ) ), int( around( pos[ 1 ] ) ) ]
    max_flux = get_pixel_value( facet, [ max_x, max_y ] )
  facet_size = get_image_size( facet )
  [ beam_bmaj_pix, beam_bmin_pix, beam_bpa_pix ] = convert_beam_size( facet, to_pixel = True )
  x_min = max( [ 1, max_x - 32 ] )
  x_max = min( [ facet_size[ 0 ], max_x + 31 ] )
  y_min = max( [ 1, max_y - 32 ] )
  y_max = min( [ facet_size[ 1 ], max_y + 31 ] )
  good_fit = True

  # check if fit is done on beam
  if ( [ beam_bmaj_pix, beam_bmin_pix, beam_bpa_pix ] == [ 0.,0.,0. ] ):
    try:
      [ fmax, fpos, fwidth, domax, dopos, dowidth ] = call_aips_task( 'JMFIT', 
          indata = facet, ngauss = 1, niter = 200, gmax = [ max_flux, 0,0,0 ],
          blc = [ x_min, y_min, 0,0,0,0,0 ], trc = [ x_max, y_max, 0,0,0,0,0 ],
          gpos = [ [ max_x, max_y ], [ 0,0 ], [ 0,0 ], [ 0,0 ] ],
          gwidth = [ [ 4.,4.,0. ], [ 0,0,0 ], [ 0,0,0 ], [ 0,0,0 ] ],
          domax = [ 1,0,0,0 ], dopos = [ [ 1,1 ], [ 0,0 ], [ 0,0 ], [ 0,0 ] ], 
          dowidth = [ [ 1,1,1 ], [ 0,0,0 ], [ 0,0,0 ], [ 0,0,0 ] ], 
          dooutput = -1, domodel = -1, outvers = -1, offset = 0.,
          outputs = [ 'fmax', 'fpos', 'fwidth', 'domax', 'dopos', 'dowidth' ] )
    except:
      return fit_results
    
    peak_flux = fmax[ 0 ]
    peak_flux_err = domax[ 0 ]
    [ peak_x, peak_y ] = fpos[ 0 ]
    [ peak_x_err, peak_y_err ] = dopos[ 0 ]
    [ bmaj_pix, bmin_pix, bpa_pix ] = fwidth[ 0 ]
    [ bmaj_err_pix, bmin_err_pix, bpa_err_pix ] = dowidth[ 0 ]
    # sanity check on fit
    if ( ( peak_x < 0.5 ) or ( peak_x > facet_size[ 0 ] + 0.5 ) or 
        ( peak_y < 0.5 ) or ( peak_y > facet_size[ 1 ] + 0.5 ) ):
      return fit_results
    # convert bmaj, bmin and bpa to arcsec
    # TODO: use errors 
    [ int_bmaj, int_bmin, int_bpa ] = convert_beam_size( facet, 
        beam = [ bmaj_pix, bmin_pix, bpa_pix ] )
    fit_results.append( [ [ peak_x, peak_y ], peak_flux, [ int_bmaj, int_bmin, int_bpa ] ] )
    return fit_results

  # fit gaussian to peak
  try:
    [ fmax, fpos, fwidth, domax, dopos, dowidth ] = call_aips_task( 'JMFIT', 
        indata = facet, ngauss = 1, niter = 200, gmax = [ max_flux, 0, 0, 0 ], 
        blc = [ x_min, y_min, 0,0,0,0,0 ], trc = [ x_max, y_max, 0,0,0,0,0 ],
        gpos = [ [ max_x, max_y ], [ 0,0 ], [ 0,0 ], [ 0,0 ] ],
        gwidth = [ [ beam_bmaj_pix, beam_bmin_pix, beam_bpa_pix ], [ 0,0,0 ], [ 0,0,0 ], [ 0,0,0 ] ],
        domax = [ 1,0,0,0 ], dopos = [ [ 1,1 ], [ 0,0 ], [ 0,0 ], [ 0,0 ] ], 
        dowidth = [ [ 1,1,1 ], [ 0,0,0 ], [ 0,0,0 ], [ 0,0,0 ] ], 
        dooutput = -1, domodel = -1, outvers = -1, offset = 0.,
        outputs = [ 'fmax', 'fpos', 'fwidth', 'domax', 'dopos', 'dowidth' ] )
  except:
    return fit_results
  
  peak_flux = fmax[ 0 ]
  peak_flux_err = domax[ 0 ]
  [ peak_x, peak_y ] = fpos[ 0 ]
  [ peak_x_err, peak_y_err ] = dopos[ 0 ]
  [ bmaj_pix, bmin_pix, bpa_pix ] = fwidth[ 0 ]
  [ bmaj_err_pix, bmin_err_pix, bpa_err_pix ] = dowidth[ 0 ]
  # sanity check on fit
  if ( ( peak_x < 0.5 ) or ( peak_x > facet_size[ 0 ] + 0.5 ) or 
      ( peak_y < 0.5 ) or ( peak_y > facet_size[ 1 ] + 0.5 ) ):
    return fit_results
  r2 = ( peak_x - max_x )**2 + ( peak_y - max_y )**2
  if  ( r2 < ( offset_ratio_max * beam_bmaj_pix )**2 ):
    # convert bmaj, bmin and bpa to arcsec
    # TODO: use errors 
    [ int_bmaj, int_bmin, int_bpa ] = convert_beam_size( facet, 
        beam = [ bmaj_pix, bmin_pix, bpa_pix ] )
    fit_results.append( [ [ peak_x, peak_y ], peak_flux, [ int_bmaj, int_bmin, int_bpa ] ] )
    return fit_results
  
  # fitted peak is offset a long way from the peak flux pixel
  # maybe we're dealing with a close double here
  if ( peak_flux < ( ( 1. + double_ratio_min ) / 2. ) * max_flux ):
    return fit_results
  
  # try fitting a double gaussian
  try:
    [ fmax, fpos, fwidth, domax, dopos, dowidth ] = call_aips_task( 'JMFIT', 
        indata = facet, ngauss = 2, gmax = [ max_flux, peak_flux, 0,0 ], niter = 200,
        blc = [ x_min, y_min, 0,0,0,0,0 ], trc = [ x_max, y_max, 0,0,0,0,0 ], 
        gpos = [ [ max_x, max_y ], [ peak_x, peak_y ], [ 0,0 ], [ 0,0 ] ],
        gwidth = [ [ beam_bmaj_pix, beam_bmin_pix, beam_bpa_pix ], 
        [ beam_bmaj_pix, beam_bmin_pix, beam_bpa_pix ], [ 0,0,0 ], [ 0,0,0 ] ],
        domax = [ 1,1,0,0 ], dopos = [ [ 1,1 ], [ 1,1 ], [ 0,0 ], [ 0,0 ] ], 
        dowidth = [ [ 1,1,1 ], [ 1,1,1 ], [ 0,0,0 ], [ 0,0,0 ] ], 
        dooutput = -1, domodel = -1, outvers = -1, offset = 0.,
        outputs = [ 'fmax', 'fpos', 'fwidth', 'domax', 'dopos', 'dowidth' ] )
  except:
    return fit_results
  
  # use the component closest to the peak flux pixel
  [ x1, y1 ] = fpos[ 0 ]
  [ x2, y2 ] = fpos[ 1 ]
  r1_2 = ( x1 - max_x )**2 + ( y1 - max_y )**2
  r2_2 = ( x2 - max_x )**2 + ( y2 - max_y )**2
  # sanity check on fit
  if ( ( x1 < 0.5 ) or ( x1 > facet_size[ 0 ] + 0.5 ) or 
      ( y1 < 0.5 ) or ( y1 > facet_size[ 1 ] + 0.5 ) ):
    r1_2 = 1.e6
  if ( ( x2 < 0.5 ) or ( x2 > facet_size[ 0 ] + 0.5 ) or 
      ( y2 < 0.5 ) or ( y2 > facet_size[ 1 ] + 0.5 ) ):
    r2_2 = 1.e6
  if ( r1_2 <= r2_2 ):
    rc_2 = r1_2
    ro_2 = r2_2
    closest_index = 0
    other_index = 1
  else:
    rc_2 = r2_2
    ro_2 = r1_2
    closest_index = 1
    other_index = 0
  [ peak_x, peak_y ] = fpos[ closest_index ]
  peak_flux = fmax[ closest_index ]
  peak_flux_err = domax[ closest_index ]
  [ peak_x_err, peak_y_err ] = dopos[ closest_index ]
  [ bmaj_pix, bmin_pix, bpa_pix ] = fwidth[ closest_index ]
  [ bmaj_err_pix, bmin_err_pix, bpa_err_pix ] = dowidth[ closest_index ]
  # sanity check on fit
  if  ( rc_2 > ( offset_ratio_max * beam_bmaj_pix )**2 ):
    return fit_results
  # convert bmaj, bmin and bpa to arcsec
  # TODO: use errors 
  [ int_bmaj, int_bmin, int_bpa ] = convert_beam_size( facet, 
      beam = [ bmaj_pix, bmin_pix, bpa_pix ] )
  fit_results.append( [ [ peak_x, peak_y ], peak_flux, [ int_bmaj, int_bmin, int_bpa ] ] )
  
  if return_double_fit:
    peak_flux = fmax[ other_index ]
    peak_flux_err = domax[ other_index ]
    [ peak_x, peak_y ] = fpos[ other_index ]
    [ peak_x_err, peak_y_err ] = dopos[ other_index ]
    [ bmaj_pix, bmin_pix, bpa_pix ] = fwidth[ other_index ]
    [ bmaj_err_pix, bmin_err_pix, bpa_err_pix ] = dowidth[ other_index ]
    # convert bmaj, bmin and bpa to arcsec
    # TODO: use errors 
    [ int_bmaj, int_bmin, int_bpa ] = convert_beam_size( facet, 
        beam = [ bmaj_pix, bmin_pix, bpa_pix ] )
    # sanity check on fit
    if  ( ro_2 < ( offset_ratio_max * beam_bmaj_pix )**2 ):
      fit_results.append( [ [ peak_x, peak_y ], peak_flux, [ int_bmaj, int_bmin, int_bpa ] ] )
      fit_results.sort( cmp = lambda a, b: cmp( b[ 1 ], a[ 1 ] ) )
  
  return fit_results

###############################################################################

def add_clean_boxes_old( uv, facets, sigma_min, facet_file_name = '', 
    peak_flux_ratio_max = None, clean_box_radius = 10, keep_boxes = True,
    new_boxes_per_facet_max = 100, facet_based_noise = False, 
    print_info = True, allow_negatives = True ):
  
  if ( facet_file_name == '' ):
    used_facet_file_name = restore_parameter( facets, 'facet_file_name' )
  else:
    used_facet_file_name = facet_file_name
  
  # if requested, remove clean box definitions from current boxfile
  facet_count = restore_parameter( facets, 'facet_count' )
  facet_list = range( 1, facet_count + 1 )
  temp_facet_file_name = used_facet_file_name + '.TEMP'
  if file_exists( temp_facet_file_name ):
    remove_file( temp_facet_file_name )
  if keep_boxes:
    copy_file( used_facet_file_name, temp_facet_file_name )
  else:
    extract_facet_definitions( used_facet_file_name, facet_list,
        temp_facet_file_name, include_clean_boxes = False )
  
  # blank non-search area and determine peak in local noise
  # determine minimum flux for sources
  peak_flux_max = 0.
  peak_rms = restore_parameter( uv, 'cpb_noise' )
  res_facets = get_aips_file( facets.disk, 'RES', facets.klass, - 1, 'MA' )
  for i in range( 1, facet_count + 1 ):
    facet_i = get_facet( facets, i )
    res_facet_i = get_facet( res_facets, i )
    call_aips_task( 'MOVE', indata = facet_i, outdata = res_facet_i,
        userid = get_aips_userid(), opcode = '' )
    fill_facet( res_facet_i, facet_file_name = temp_facet_file_name,
        do_edge_circle = True )
    if allow_negatives:
      [ peak_flux, pos ] = get_image_extremum( res_facet_i, force_positive = False )
      peak_flux = abs( peak_flux )
    else:
      [ peak_flux, pos ] = get_image_extremum( res_facet_i, force_positive = True )
    if ( peak_flux > peak_flux_max ):
      peak_flux_max = peak_flux
    if facet_based_noise:
      rms = get_image_rms( res_facet_i )
      [ avg, noise ] = call_aips_task( 'IMEAN', indata = res_facet_i,
          pixavg = 0., pixstd = rms, pixrange = [ -5. * rms, 5. * rms  ], 
          outputs = [ 'pixavg', 'pixstd' ] )
      if ( ( noise <= 0. ) or ( noise > 2. * rms ) ):
        if print_info:
          print '... WARNING: histogram noise fit failed, using image RMS instead'
        noise = rms
      if ( noise > peak_rms ):
        peak_rms = noise
  flux_min = sigma_min * peak_rms
  if ( peak_flux_ratio_max is None ):
    peak_flux_min = flux_min
  else:
    peak_flux_min = max( [ flux_min, peak_flux_max / peak_flux_ratio_max ] )
  if print_info:
    print 'clean box flux min = ' + repr( peak_flux_min )
  store_parameter( facets, 'clean_box_flux_min', peak_flux_min )
  
  # search for positive peaks and add clean boxes to them
  facet_has_box = [ False for i in range( facet_count ) ]
  temp2_facet_file_name = used_facet_file_name + '.TEMP2'
  if ( peak_flux_max >= peak_flux_min ):
    for i in range( 1, 1 + facet_count ):
      res_facet_i = get_facet( res_facets, i )
#      [ peak_flux, pos ] = get_image_maximum( res_facet_i )
      if allow_negatives:
        [ peak_flux, pos ] = get_image_extremum( res_facet_i, force_positive = False )
        peak_flux = abs( peak_flux )
      else:
        [ peak_flux, pos ] = get_image_extremum( res_facet_i, force_positive = True )
      new_boxes_count = 0
      while ( ( peak_flux >= peak_flux_min ) and 
          ( new_boxes_count < new_boxes_per_facet_max ) ):
        # prepare filling of peak region
        extract_facet_definitions( temp_facet_file_name, facet_list,
            temp2_facet_file_name, include_clean_boxes = False )
        add_circular_clean_box( temp2_facet_file_name, i, pos, clean_box_radius )
        # check if source is closest to center of current facet
        peak_radec = calculate_source_radec( res_facet_i, pos )
        overlapping_facet_list = [ i ] + get_facet_overlap( res_facet_i )
        [ [ main_i, main_pos ] ] = find_source_facets( facets, peak_radec,
            primary_facet_only = True, facet_list = overlapping_facet_list )
        while ( main_i < i ):
          res_facet_j = get_facet( res_facets, main_i )
          pixel_value = get_pixel_value( res_facet_j, main_pos, to_float = False )
          if ( pixel_value == get_aips_magic_value() ):
            break
          else:
            overlapping_facet_list.remove( main_i )
            [ [ main_i, main_pos ] ] = find_source_facets( facets, peak_radec, 
                primary_facet_only = True, facet_list = overlapping_facet_list )
        if ( i == main_i ):
          if print_info:
            print ( '... adding clean box to facet ' + repr( i ) + 
                ' at position ' + repr( pos ) )
          facet_has_box[ i - 1 ] = True
          new_boxes_count = new_boxes_count + 1
          add_circular_clean_box( temp_facet_file_name, i, pos, clean_box_radius )
          fill_facet( res_facet_i, facet_file_name = temp2_facet_file_name,
              do_edge_circle = False )
        else:
          fill_facet( res_facet_i, facet_file_name = temp2_facet_file_name,
              do_edge_circle = False, value = 0. )
#        [ peak_flux, pos ] = get_image_maximum( res_facet_i )
        if allow_negatives:
          [ peak_flux, pos ] = get_image_extremum( res_facet_i, force_positive = False )
          peak_flux = abs( peak_flux )
        else:
          [ peak_flux, pos ] = get_image_extremum( res_facet_i, force_positive = True )
  
  # delete residual facets
  remove_facets( res_facets )
  remove_file( temp2_facet_file_name )
  
  # write dummy clean boxes to facets without clean box
  if keep_boxes:
    for i in range( 1, facet_count + 1 ):
      if facet_has_box[ i - 1 ]:
        remove_circular_clean_box( temp_facet_file_name, i, [ 5, 5 ], 1 )
        remove_rectangular_clean_box( temp_facet_file_name, i, [ 0, 0 ], [ 0, 0 ] )
  else: # ( not keep_boxes ):
    for i in range( 1, facet_count + 1 ):
      if ( not facet_has_box[ i - 1 ] ):
        add_rectangular_clean_box( temp_facet_file_name, i, [ 0, 0 ], [ 0, 0 ] )
  
  # replace old facet file with new one
  move_file( temp_facet_file_name, used_facet_file_name )
  
  return

###############################################################################

def add_clean_boxes( uv, facets, facet_list = [], facet_file_name = '', 
    clean_box_radius = 5, keep_boxes = False, new_boxes_per_facet_max = 100,
    facet_based_boxes = True, sidelobe_rejection = 1.05, print_info = True,
    box_sigma = 5., convolve_factor = 1. ):
  
  # process input parameters
  if ( facet_file_name == '' ):
    used_facet_file_name = restore_parameter( facets, 'facet_file_name' )
  else:
    used_facet_file_name = facet_file_name
  used_facet_file_name = path.expandvars( used_facet_file_name )
  if ( len( facet_list ) == 0 ):
    used_facet_list = range( 1, 1 + restore_parameter( facets, 'facet_count' ) )
  else:
    used_facet_list = facet_list
  
  # make work copies of facets and facet_file
  res_facet_file_name = used_facet_file_name + '.NEW'
  extract_facet_definitions( used_facet_file_name, used_facet_list,
      res_facet_file_name, include_clean_boxes = keep_boxes,
      renumber_facets = False )
  res_facets = get_aips_file( facets.disk, 'RES', facets.klass, - 1, 'MA' )
  for i in used_facet_list:
    facet = get_facet( facets, i )
    res_facet = get_facet( res_facets, i )
    call_aips_task( 'MOVE', indata = facet, outdata = res_facet,
        userid = get_aips_userid() )
  fill_facets( res_facets, used_facet_list, facet_file_name = res_facet_file_name,
      do_boxes = True, do_edge_circle = True )
  
  # convolve facets to capture large-scale emission
  if ( convolve_factor > 1. ):
    temp_facets = get_aips_file( facets.disk, 'TEMP', facets.klass, -1, 'MA' )
    for i in used_facet_list:
      facet = get_facet( res_facets, i )
      temp_facet = get_facet( temp_facets, i )
      [ bmaj, bmin, bpa ] = get_beam_size( facet )
      bmaj = convolve_factor * bmaj
      bmin = convolve_factor * bmin
      call_aips_task( 'CONVL', indata = facet, outdata = temp_facet,
          opcode = 'GAUS', bmaj = bmaj, bmin = bmin, bpa = bpa, doblank = 1 )
      facet.zap()
      temp_facet.rename( name = facet.name, klass = facet.klass, seq = facet.seq )
  
  # determine clean box thresholds
  global_noise = restore_parameter( uv, 'cpb_noise' )
  facet_threshold_list = []
  for i in used_facet_list:
    facet = get_facet( res_facets, i )
    if facet_based_boxes:
      rms = get_image_rms( facet )
      [ avg, noise ] = call_aips_task( 'IMEAN', indata = facet,
          pixavg = 0., pixstd = rms, pixrange = [ -5. * rms, 5. * rms  ], 
          outputs = [ 'pixavg', 'pixstd' ] )
      if ( ( noise <= 0. ) or ( noise > 2. * rms ) ):
        if print_info:
          print '... WARNING: histogram noise fit failed, using image RMS instead'
        noise = rms
    else:
      noise = global_noise
    threshold = abs( get_image_minimum( facet )[ 0 ] ) * sidelobe_rejection
    threshold = max( threshold, box_sigma * noise )
    facet_threshold_list.append( threshold )
  if ( not facet_based_boxes ):
    threshold = max( facet_threshold_list )
    for i in range( len( used_facet_list ) ):
      facet_threshold_list[ i ] = threshold
  if print_info:
    print 'clean box flux min = ' + repr( min( facet_threshold_list ) )
  store_parameter( facets, 'clean_box_flux_min', min( facet_threshold_list ) )
  
  # search for positive peaks and add clean boxes to them
  facet_has_box = [ False for i in used_facet_list ]
  new_box_list = []
  for j in range( len( used_facet_list ) ):
    i = used_facet_list[ j ]
    facet = get_facet( res_facets, i )
    [ peak_flux, peak_pos ] = get_image_maximum( facet )
    new_boxes_count = 0
    while ( ( peak_flux >= facet_threshold_list[ j ] ) and 
        ( new_boxes_count < new_boxes_per_facet_max ) ):
      
      # add clean box & fill peak region
      add_circular_clean_box( res_facet_file_name, i, peak_pos, clean_box_radius )
      fill_facets( res_facets, [ i ], do_edge_circle = False, do_boxes = True,
          facet_file_name = res_facet_file_name )
      
      # check if peak is closest to center of current facet
      # if not, remove clean box
      peak_radec = calculate_source_radec( facet, peak_pos )
      try:
        overlapping_facet_list = [ i ] + get_facet_overlap( facet )
        [ [ prim_i, prim_pos ] ] = find_source_facets( facets, peak_radec,
              primary_facet_only = True, facet_list = overlapping_facet_list )
      except:
        prim_i = i
      if ( i == prim_i ):
        if print_info:
          print ( '... adding clean box to facet ' + repr( i ) + 
              ' at position ' + repr( peak_pos ) )
        facet_has_box[ j ] = True
        new_boxes_count = new_boxes_count + 1
      else:
        remove_circular_clean_box( res_facet_file_name, i, peak_pos, clean_box_radius )
      # look for next peak
      [ peak_flux, peak_pos ] = get_image_maximum( facet )
  
  # write dummy clean boxes to facets without clean box
  for j in range( len( used_facet_list ) ):
    i = used_facet_list[ j ]
    if ( keep_boxes and facet_has_box[ j ] ):
      remove_circular_clean_box( res_facet_file_name, i, [ 5, 5 ], 1 )
      remove_rectangular_clean_box( res_facet_file_name, i, [ 0, 0 ], [ 0, 0 ] )
    elif ( not ( keep_boxes or facet_has_box[ j ] ) ):
      add_rectangular_clean_box( res_facet_file_name, i, [ 0, 0 ], [ 0, 0 ] )
  
  # replace old facet file with new one
  old_facet_file_name = used_facet_file_name + '.OLD'
  move_file( used_facet_file_name, old_facet_file_name )
  if ( len( facet_list ) == 0 ):
    move_file( res_facet_file_name, used_facet_file_name )
  else:
    temp_facet_file_name = used_facet_file_name + '.TEMP'
    extract_facet_definitions( old_facet_file_name, used_facet_list,
        temp_facet_file_name, include_clean_boxes = True,
        renumber_facets = False, invert_selection = True )
    merge_facet_definitions( temp_facet_file_name, res_facet_file_name,
        used_facet_file_name, renumber_facets = False )
    remove_file( res_facet_file_name )
    remove_file( temp_facet_file_name )
  
  # cleanup
  remove_facets( res_facets )
  
  return

###############################################################################

def add_centered_pb_facets( uv, pb_facets, sigma_min, blank_radius = 10,
    restore_components = False, apply_solutions = False, print_info = True,
    new_facet_count_max = None, center_facets = None, gain = 0.1, factor = 0.,
    imagr_params = {}, facet_list = [], clean_radius = 10, min_separation = 0.,
    keep_model_components = False ):
# IMPORTANT: supply exactly the same imaging parameters as used for production
# of the pb_facets
# min_separation in arcsec
  
  # restore parameters
  pb_facet_count = restore_parameter( pb_facets, 'facet_count' )
  pb_facet_file_name = restore_parameter( pb_facets, 'facet_file_name' )
  pb_facet_file_name = path.expandvars( pb_facet_file_name )
  pb_facet_size = get_image_size( pb_facets )
  try:
    old_added_facet_count = restore_parameter( pb_facets, 'added_facet_count' )
  except:
    old_added_facet_count = 0
  old_pb_facet_count = pb_facet_count - old_added_facet_count
  temp_facet_file_name = pb_facet_file_name + '.TEMP'
  
  # save clean box information
  clean_box_list = get_clean_boxes( pb_facet_file_name )
  clean_box_radec_list = []
  for clean_box in clean_box_list:
    [ i, a, b, c, d ] = clean_box
    facet = get_facet( pb_facets, i )
    if ( a == -1 ): # circular area
      radec = calculate_source_radec( facet, [ c, d ] )
    else: # rectangular area
      radec = calculate_source_radec( facet,
          [ float( a + c ) / 2., float( b + d ) / 2. ] )
    # link clean box to facets
    if ( i in range( 1, 1 + old_pb_facet_count ) ):
      clean_box_radec_list.append( [ i, radec ] )
    else:
      clean_box_radec_list.append( [ -1, radec ] )
  
  # create list of facets to add
  new_facet_list = []
  new_center_facet_list = []
  if ( not center_facets is None ):
    center_facet_count = restore_parameter( center_facets, 'facet_count' )
    center_facet_file_name = restore_parameter( center_facets, 'facet_file_name' )
    center_facet_file_name = path.expandvars( center_facet_file_name )
    if ( len( facet_list ) > 0 ):
      center_facet_list = facet_list
    else:
      center_facet_list = range( 1, 1 + center_facet_count )
    for i in center_facet_list:
      facet = get_facet( center_facets, i )
      if ( min_separation > 0. ):
        results = find_source_facets( pb_facets, get_radec( facet ),
            primary_facet_only = False )
        close_facet_found = False
        for [ j, pos ] in results:
          if ( j > old_pb_facet_count ):
            continue
          pb_facet = get_facet( pb_facets, j )
          [ r, p ] = calculate_angular_separation( get_radec( facet ),
              get_radec( pb_facet ) )
          if ( r < min_separation / 3600. ):
            close_facet_found = True
            break
        if close_facet_found:
          continue
      new_center_facet_list.append( i )
      pos = get_pixel_reference( facet )
      fit_peak = get_pixel_value( facet, pos )
      peak_radec = get_radec( facet )
      new_facet_list.append( [ fit_peak, peak_radec ] )
    # add clean box information
    center_box_list = get_clean_boxes( center_facet_file_name,
        facet_list = new_center_facet_list )
    for clean_box in center_box_list:
      clean_box_list.append( clean_box )
      [ i, a, b, c, d ] = clean_box
      facet = get_facet( center_facets, i )
      if ( a == -1 ): # circular area
        radec = calculate_source_radec( facet, [ c, d ] )
      else: # rectangular area
        radec = calculate_source_radec( facet, 
            [ float( a + c ) / 2., float( b + d ) / 2. ] )
      clean_box_radec_list.append( [ -1, radec ] )
  else:
    # restore model components
    if restore_components:
      facets = restore_model_components( pb_facets, imagr_params = imagr_params )
    else:
      facets = pb_facets
    # blank non-search area and determine peak in local noise
    res_facets = get_aips_file( pb_facets.disk, 'RES', pb_facets.klass, - 1, 'MA' )
    peak_rms = restore_parameter( uv, 'cpb_noise' )
    for i in range( 1, old_pb_facet_count + 1 ):
      if ( len( facet_list ) > 0 ):
        if ( not i in facet_list ):
          continue
      facet_i = get_facet( facets, i )
      res_facet_i = get_facet( res_facets, i )
      call_aips_task( 'MOVE', indata = facet_i, outdata = res_facet_i,
          userid = get_aips_userid(), opcode = '' )
      fill_image( res_facet_i, do_edge_circle = True )
      if restore_components:
        facet_i.zap()
    flux_min = sigma_min * peak_rms
    if print_info:
      print '... minimum source flux = %s Jy' % ( repr( flux_min ) )
    # search for positive peaks
    extract_facet_definitions( pb_facet_file_name, range( 1, 1 + old_pb_facet_count ),
      temp_facet_file_name, include_clean_boxes = False )
    for i in range( 1, old_pb_facet_count + 1 ):
      res_facet_i = get_facet( res_facets, i )
      [ peak_flux, pos ] = get_image_extremum( res_facet_i, force_positive = True )
      while ( peak_flux >= flux_min ):
        fit_results = fit_gaussian_to_peak( res_facet_i, pos = pos,
            return_double_fit = True )
        if ( len( fit_results ) == 0 ):
          fit_pos = [ float( p ) for p in pos ]
          fit_peak = peak_flux
        else:
          [ fit_pos, fit_peak, fit_beam ] = fit_results[ 0 ]
        # blank source area in residual facet
        add_circular_clean_box( temp_facet_file_name, i, pos, blank_radius )
        fill_facets( res_facets, [ i ], facet_file_name = temp_facet_file_name,
            do_boxes = True )
        # check if source is closest to center of current facet
        peak_radec = calculate_source_radec( res_facet_i, fit_pos )
        overlapping_facet_list = [ i ] + get_facet_overlap( res_facet_i )
        for j in range( old_pb_facet_count + 1, pb_facet_count + 1 ):
          if ( j in overlapping_facet_list ):
            overlapping_facet_list.remove( j )
        [ [ main_i, main_pos ] ] = find_source_facets( pb_facets, peak_radec,
            primary_facet_only = True, facet_list = overlapping_facet_list )
        if ( i == main_i ):
          # save new facet radec
          new_facet_list.append( [ fit_peak, peak_radec ] )
        else:
          # double check that peak was already associated with previous main facet
          while ( main_i < i ):
            res_facet_j = get_facet( res_facets, main_i )
            pixel_value = get_pixel_value( res_facet_j, main_pos, to_float = False )
            if ( pixel_value == 0. ):
              break
            else:
              # if not, find next main facet
              overlapping_facet_list.remove( main_i )
              [ [ main_i, main_pos ] ] = find_source_facets( pb_facets, peak_radec, 
                  primary_facet_only = True, facet_list = overlapping_facet_list )
              if ( i == main_i ):
                new_facet_list.append( [ fit_peak, peak_radec ] )
          if ( i != main_i ):
            remove_circular_clean_box( temp_facet_file_name, i, pos, blank_radius )
        # mark the actual clean area with zero
        fill_facets( res_facets, [ i ], facet_file_name = temp_facet_file_name,
            do_boxes = True, value = 0. )
        # find next peak in residual facet
        [ peak_flux, pos ] = get_image_extremum( res_facet_i, force_positive = True )
    if print_info:
      print '... found %s sources above minimum flux' % ( repr( len( new_facet_list ) ) )
    # delete residual facets
    remove_facets( res_facets )
    remove_file( temp_facet_file_name )
    # remove restored facets
    if restore_components:
      remove_facets( facets )
  
  # remove previously added facets
  new_facet_file_name = pb_facet_file_name + '.NEW'
  extract_facet_definitions( pb_facet_file_name, range( 1, 1 + old_pb_facet_count ),
      new_facet_file_name, include_clean_boxes = False )
  for i in range( 1 + old_pb_facet_count, 1 + pb_facet_count ):
    pb_facet_i = get_facet( pb_facets, i )
    if pb_facet_i.exists():
      pb_facet_i.zap()
    pb_beam_i = get_facet_beam( pb_facet_i )
    if pb_beam_i.exists():
      pb_beam_i.zap()
  if print_info:
    print '... removed %s previously added facets' % ( repr( old_added_facet_count ) )
  
  # add new facets
  new_facet_list.sort( cmp = lambda a, b: cmp( b[ 0 ], a[ 0 ] ) )
  new_added_facet_count = len( new_facet_list )
  if ( not new_facet_count_max is None ):
    new_added_facet_count = min( [ new_added_facet_count, new_facet_count_max ] )
  for i in range( new_added_facet_count ):
    peak_radec = new_facet_list[ i ][ 1 ]
    new_i = 1 + old_pb_facet_count + i
    add_facet( new_facet_file_name, peak_radec, pb_facet_size, facet_id = new_i )
  new_pb_facet_count = old_pb_facet_count + new_added_facet_count
  # generate new facets images and beams
  if ( new_added_facet_count > 0 ):
    # set calibration parameters
    if apply_solutions:
      solution_switch = 100
      solution_version = uv.table_highver( 'SN' )
    else:
      solution_switch = -1
      solution_version = -1
    # extract new facet definitions
    temp_facet_list = range( 1 + old_pb_facet_count, 1 + new_pb_facet_count )
    extract_facet_definitions( new_facet_file_name, temp_facet_list,
        temp_facet_file_name )
    # make beams and dirty images of new facets
    temp_facets = get_aips_file( pb_facets.disk, 'TEMP', 'ICL001', - 1, 'MA' )
    cell_size = restore_parameter( uv, 'cell_size' )
    uv_size = restore_parameter( uv, 'pb_image_size' )
    channel_count = get_channel_count( uv )
    model_version = pb_facets.table_highver( 'CC' )
    call_aips_task( 'IMAGR', indata = uv, nchav = channel_count, flagver = -1,
        docalib = solution_switch, gainuse = solution_version, gain = gain,
        outdisk = temp_facets.disk, outname = temp_facets.name, maxpixel = 0,
        outseq = temp_facets.seq, outver = model_version, in2disk = uv.disk,
        imsize = pb_facet_size, do3dimag = 1, cellsize = [ cell_size, cell_size ],
        niter = 1000000, flux = 100000., dotv = 0, boxfile = temp_facet_file_name,
        overlap = 2, nfield = new_added_facet_count, allokay = 0, factor = factor,
        imagrprm = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 ],
        cmethod = '', minpatch = pb_facet_size[ 0 ] - 1,
        uvsize = [ uv_size, uv_size ], **imagr_params )
    # copy model components
    if ( keep_model_components and ( not center_facets is None ) ):
      for i in range( 1, 1 + new_added_facet_count ):
        center_facet = get_facet( center_facets, new_center_facet_list[ i - 1 ] )
        new_facet = get_facet( temp_facets, i )
        if table_exists( new_facet, 'CC', model_version ):
          new_facet.zap_table( 'CC', model_version )
        call_aips_task( 'TACOP', indata = center_facet, outdata = new_facet,
            inext = 'CC', invers = 0, outvers = model_version, ncount = 1 )
    # rename facets and beams to place them behind the already existing facets
    for i in range( 1, 1 + new_added_facet_count ):
      temp_facet_i = get_facet( temp_facets, i )
      temp_beam_i = get_facet_beam( temp_facet_i )
      new_i = new_pb_facet_count - new_added_facet_count + i
      facet_i = get_facet( pb_facets, new_i )
      beam_i = get_facet_beam( facet_i )
      if facet_i.exists():
        facet_i.clrstat()
        facet_i.zap()
      if beam_i.exists():
        beam_i.clrstat()
        beam_i.zap()
      temp_facet_i.rename( name = facet_i.name, klass = facet_i.klass, seq = facet_i.seq )
      temp_beam_i.rename( name = beam_i.name, klass = beam_i.klass, seq = beam_i.seq )
    # update parameters
    store_parameter( pb_facets, 'facet_count', new_pb_facet_count )
    store_parameter( pb_facets, 'added_facet_count', new_added_facet_count )
    store_parameter( uv, 'pb_facet_count', new_pb_facet_count )
    store_parameter( uv, 'added_pb_facet_count', new_added_facet_count )
  # update facet overlap info
  determine_facet_overlap( pb_facets )
  if print_info:
    print '... added %s centered facets' % ( repr( new_added_facet_count ) )
  
  # store clean box information in new facet configuration
  facet_has_clean_box = [ False for i in range( new_pb_facet_count ) ]
  for k in range( len( clean_box_list ) ):
    [ i, a, b, c, d ] = clean_box_list[ k ]
    [ j, radec ] = clean_box_radec_list[ k ]
    # find target facet
    if ( j > 0 ):
      facet = get_facet( pb_facets, j )
      overlap_list = [ j ] + get_facet_overlap( facet )
      source_facet_list = find_source_facets( pb_facets, radec, primary_facet_only = True,
          facet_list = overlap_list )
    else:
      source_facet_list = find_source_facets( pb_facets, radec, primary_facet_only = True )
    if ( len( source_facet_list ) == 0 ):
      # clean box lies outside all facets, so discard
      if print_info:
        print "WARNING: discarding clean box at radec %s" % ( 
            repr( degdeg_to_hmsdms( radec ) ) )
      continue
    [ main_i, main_pos ] = source_facet_list[ 0 ]
    if ( main_i != j ):
      # recalculate clean box coordinates in new facet
      # TODO: should we adjust clean box size to compensate for fractional pixel shift ?
      if ( a == -1 ): # circular area
        [ c, d ] = [ int( around( main_pos[ 0 ] ) ), int( around( main_pos[ 1 ] ) ) ]
      else:
        [ dx, dy ] = [ abs( float( c - a ) / 2. ), abs( float( d - b ) / 2. ) ]
        [ a, b ] = [ int( around( main_pos[ 0 ] - dx ) ), 
            int( around( main_pos[ 1 ] - dy ) ) ]
        [ c, d ] = [ int( around( main_pos[ 0 ] + dx ) ),
            int( around( main_pos[ 1 ] + dy ) ) ]
    if ( a == -1 ): # circular area
      add_circular_clean_box( new_facet_file_name, main_i, [ c, d ], b )
    else: # rectangular area
      add_rectangular_clean_box( new_facet_file_name, main_i, [ a, b ], [ c, d ] )
    if ( not facet_has_clean_box[ main_i - 1 ] ):
      facet_has_clean_box[ main_i - 1 ] = True
  
  # add dummy clean boxes
  for i in range( 1, 1 + new_pb_facet_count ):
    if ( not facet_has_clean_box[ i - 1 ] ):
      add_rectangular_clean_box( new_facet_file_name, i, [ 0, 0 ], [ 0, 0 ] )
  
  # replace facet file
  move_file( pb_facet_file_name, pb_facet_file_name + '.OLD' )
  remove_file( temp_facet_file_name )
  move_file( new_facet_file_name, pb_facet_file_name )
  
  return new_added_facet_count

###############################################################################

def re_center_pb_facets( uv, facets, facet_list = [], facet_file_name = '',
    search_radius = 10., max_shift = 8., gain = 0.1, factor = 0., imagr_params = {} ):
# shifts facets by a fraction of a pixel to center sources
# clean boxes stay where they are
# TODO: should we update the facet overlap info (probably not)?

  pb_facet_count = restore_parameter( uv, 'pb_facet_count' )
  try:
    added_pb_facet_count = restore_parameter( uv, 'added_pb_facet_count' )
  except:
    return
  first_pb_facet_count = pb_facet_count - added_pb_facet_count
  first_facet_list = range( 1, 1 + first_pb_facet_count )
  added_facet_list = range( 1 + first_pb_facet_count, 1 + pb_facet_count )
  if ( len( facet_list ) > 0 ):
    for i in range( 1 + first_pb_facet_count, 1 + pb_facet_count ):
      if ( not i in facet_list ):
        added_facet_list.remove( i )
  
  # create temporary copy of facet file
  if ( facet_file_name == '' ):
    used_facet_file_name = restore_parameter( facets, 'facet_file_name' )
  else:
    used_facet_file_name = facet_file_name
  used_facet_file_name = path.expandvars( used_facet_file_name )
  temp_facet_file_name = used_facet_file_name + '.TEMP'
  extract_facet_definitions( used_facet_file_name, added_facet_list,
      temp_facet_file_name, include_clean_boxes = False )

  # blank non-search area
  res_facets = get_aips_file( facets.disk, 'RES', facets.klass, - 1, 'MA' )
  i = 0
  facets_changed = False
  for j in added_facet_list:
    i = i + 1
    facet_j = get_facet( facets, j )
    res_facet_i = get_facet( res_facets, i )
    call_aips_task( 'MOVE', indata = facet_j, outdata = res_facet_i, 
        userid = get_aips_userid(), opcode = '' )
    pos = get_pixel_reference( res_facet_i )
    add_circular_clean_box( temp_facet_file_name, i, pos, search_radius )
    fill_facets( res_facets, [ i ], facet_file_name = temp_facet_file_name,
        do_boxes = True, invert = True )
    fit_results = fit_gaussian_to_peak( res_facet_i, pos = pos, return_double_fit = True )
    if ( len( fit_results ) == 0 ):
      [ peak, peak_pos ] = get_image_extremum( res_facet_i, force_positive = True )
      if ( peak_pos == pos ):
        res_facet_i.zap()
        continue
      fit_pos = [ float( p ) for p in peak_pos ]
    else:
      [ fit_pos, fit_peak, fit_beam ] = fit_results[ 0 ]
    shift2 = ( fit_pos[ 0 ] - pos[ 0 ] )**2 + ( fit_pos[ 1 ] - pos[ 1 ] )**2
    if ( shift2 <= max_shift**2 ):
      peak_radec = calculate_source_radec( res_facet_i, fit_pos )
      facet_size = get_image_size( res_facet_i )
      replace_facet( used_facet_file_name, j, facet_size, peak_radec, keep_boxes = True )
      facets_changed = True
    res_facet_i.zap()

  # generate new facets images and beams
  if facets_changed:

    # extract new facet definitions
    extract_facet_definitions( used_facet_file_name, added_facet_list,
        temp_facet_file_name, include_clean_boxes = False )

    # make dirty images of new facets
    temp_facets = get_aips_file( facets.disk, 'TEMP', 'ICL001', - 1, 'MA' )
    cell_size = restore_parameter( uv, 'cell_size' )
    uv_size = restore_parameter( uv, 'pb_image_size' )
    channel_count = get_channel_count( uv )
    facet_size = get_image_size( facets )
    call_aips_task( 'IMAGR', indata = uv, nchav = channel_count, docalib = -1,
        gainuse = -1, outdisk = temp_facets.disk, outname = temp_facets.name,
        outseq = temp_facets.seq, outver = 0, in2disk = uv.disk, flagver = -1,
        cellsize = [ cell_size, cell_size ], imsize = facet_size, do3dimag = 1,
        niter = 1000000, flux = 100000., boxfile = temp_facet_file_name, dotv = 0,
        overlap = 2, nfield = len( added_facet_list ), maxpixel = 0, factor = factor,
        imagrprm = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 ],
        cmethod = '', minpatch = facet_size[ 0 ] - 1, gain = gain, allokay = 0,
        uvsize = [ uv_size, uv_size ], **imagr_params )

    i = 0
    for j in added_facet_list:
      i = i + 1
      temp_facet_i = get_facet( temp_facets, i )
      temp_beam_i = get_facet_beam( temp_facet_i )
      facet_j = get_facet( facets, j )
      beam_j = get_facet_beam( facet_j )

      # if present, copy solution tables (keep version number)
      if table_exists( facet_j, 'SN', 0 ):
        solution_version_high = facet_j.table_highver( 'SN' )
        for solution_version in range( 1, 1 + solution_version_high ):
          call_aips_task( 'TACOP', indata = facet_j, inext = 'SN',
              invers = solution_version, ncount = 1,
              outdata = temp_facet_i, outvers = solution_version )

      # if present, copy parameter tables (keep version number)
      if table_exists( facet_j, 'PS', 0 ):
        parameter_version_high = facet_j.table_highver( 'PS' )
        for parameter_version in range( 1, 1 + parameter_version_high ):
          call_aips_task( 'TACOP', indata = facet_j, inext = 'PS',
              invers = parameter_version, ncount = 1,
              outdata = temp_facet_i, outvers = parameter_version )

      # remove old facets and beams
      facet_j.zap()
      beam_j.zap()

      # rename facets and beams to place them behind the already existing facets
      temp_facet_i.rename( name = facet_j.name, klass = facet_j.klass, seq = facet_j.seq )
      temp_beam_i.rename( name = beam_j.name, klass = beam_j.klass, seq = beam_j.seq )

  remove_file( temp_facet_file_name )

  return

###############################################################################

def restore_model_components( facets, facet_list = [], model_version = 0,
    cross_restore = True, convolve = False, imagr_params = {}, subtract = False ):
  
  # process input params
  r_params = {}
  restoring_beam = False
  if ( ( 'bmaj' in imagr_params ) or ( 'bmin' in imagr_params ) or 
      ( 'bpa' in imagr_params ) ):
    if ( ( 'bmaj' in imagr_params ) and ( 'bmin' in imagr_params ) and 
        ( 'bpa' in imagr_params ) ):
      r_params[ 'bmaj' ] = imagr_params[ 'bmaj' ]
      r_params[ 'bmin' ] = imagr_params[ 'bmin' ]
      r_params[ 'bpa' ] = imagr_params[ 'bpa' ]
      restoring_beam = True
    else:
      raise error( 'bmaj/bmin/bpa need to be specified together' )
  rst_facets = get_aips_file( facets.disk, 'RST', facets.klass, - 1, 'MA' )
  if ( len( facet_list ) > 0 ):
    rst_facet_list = facet_list
  else:
    rst_facet_list = range( 1, 1 + restore_parameter( facets, 'facet_count' ) )
  if subtract:
    mode = 'SUB'
  else:
    mode = 'ADD'
  
  # make working copies of facets and remove tables
  for i in rst_facet_list:
    facet = get_facet( facets, i )
    rst_facet = get_facet( rst_facets, i )
    call_aips_task( 'MOVE', indata = facet, outdata = rst_facet,
        userid = get_aips_userid() )
    while table_exists( rst_facet, 'CC', 0 ):
      rst_facet.zap_table( 'CC', 0 )
    while table_exists( rst_facet, 'SN', 0 ):
      rst_facet.zap_table( 'SN', 0 )
  
  # convolve background to restoring beam
  temp_facets = get_aips_file( facets.disk, 'TEMP', facets.klass, -1, 'MA' )
  if ( convolve and restoring_beam ):
    for i in rst_facet_list:
      
      # measure beam size and store in header
      beam = get_facet_beam( get_facet( facets, i ) )
      [ fit_pos, fit_peak, fit_beam ] = ( fit_gaussian_to_peak( beam ) )[ 0 ]
      set_beam_size( rst_facet, fit_beam )
      
      # convolve background
      temp_facet = get_facet( temp_facets, i )
      call_aips_task( 'CONVL', indata = rst_facet, outdata = temp_facet, 
          opcode = 'GAUS', bmaj = r_params[ 'bmaj' ], bmin = r_params[ 'bmin' ],
          bpa = r_params[ 'bpa' ], factor = 0. )
      rst_facet.zap()
      temp_facet.rename( name = rst_facet.name, klass = rst_facet.klass,
          seq = rst_facet.seq )
  
  # restore model components to facets
  model_facet_list = []
  for i in rst_facet_list:
    facet = get_facet( facets, i )
    if ( not model_table_empty( facet, model_version = model_version ) ):
      model_facet_list.append( i )
  for i in rst_facet_list:
    if ( not i in model_facet_list ):
      continue
    rst_facet = get_facet( rst_facets, i )
    temp_facet = get_facet( temp_facets, i )
    facet = get_facet( facets, i )
    call_aips_task( 'CCRES', indata = rst_facet, in2data = facet,
        invers = model_version, outdata = temp_facet, optype = mode, **r_params )
    temp_facet.zap_table( 'CC', 1 )
    ### CCRES bug fix start
    table_file_name = get_aips_file_name( temp_facet, table = 'CC', version = 1 )
    if file_exists( table_file_name ):
      remove_file( table_file_name )
    ### CCRES bug fix end
    rst_facet.zap()
    temp_facet.rename( name = rst_facet.name, klass = rst_facet.klass,
        seq = rst_facet.seq )
  
  # cross-restoration of model components to facets (not into CC table !)
  if cross_restore:
    for i in rst_facet_list:
      facet = get_facet( facets, i )
      cross_facet_list = get_facet_overlap( facet )
      for j in cross_facet_list:
        if ( not j in model_facet_list ):
          continue
        rst_facet = get_facet( rst_facets, i )
        temp_facet = get_facet( temp_facets, i )
        facet = get_facet( facets, j )
        call_aips_task( 'CCRES', indata = rst_facet, in2data = facet,
            invers = model_version, outdata = temp_facet, optype = mode, **r_params )
        temp_facet.zap_table( 'CC', 1 )
        ### CCRES bug fix start
        table_file_name = get_aips_file_name( temp_facet, table = 'CC', version = 1 )
        if file_exists( table_file_name ):
          remove_file( table_file_name )
        ### CCRES bug fix end
        rst_facet.zap()
        temp_facet.rename( name = rst_facet.name, klass = rst_facet.klass,
            seq = rst_facet.seq )
  
  # make sure we put the right model table back at the right place
  # also restore solution tables if present
  for i in rst_facet_list:
    facet = get_facet( facets, i )
    rst_facet = get_facet( rst_facets, i )
    for j in range( 1, 1 + facet.table_highver( 'CC' ) ):
      if table_exists( facet, 'CC', j ):
        call_aips_task( 'TACOP', indata = facet, inext = 'CC', invers = j, 
            ncount = 1, outdata = rst_facet, outvers = j )
    for j in range( 1, 1 + facet.table_highver( 'SN' ) ):
      if table_exists( facet, 'SN', j ):
        call_aips_task( 'TACOP', indata = facet, inext = 'SN', invers = j, 
            ncount = 1, outdata = rst_facet, outvers = j )
  
  return get_facet( rst_facets, rst_facet_list[ 0 ] )

###############################################################################

def re_image_clean_facet( uv, facets, facet_id, clean_flux, subtract = True,
    facet_file_name = '', frequency_correction = False, conversion_method = 'DFT', 
    model_version = 0, do_sdi_clean = False, remake_beam = True, imagr_params = {},
    solution_version = 0 ):

  # process inputs
  if ( facet_file_name == '' ):
    sel_facet_file_name = restore_parameter( facets, 'facet_file_name' )
  else:
    sel_facet_file_name = facet_file_name
  sel_facet_file_name = path.expandvars( sel_facet_file_name )
  if frequency_correction:
    dish_diameter = restore_parameter( uv, 'dish_diameter' )
  else:
    dish_diameter = 0.
  facet = get_facet( facets, facet_id )
  beam = get_facet_beam( facet )
  if table_exists( facet, 'CC', model_version ):
    raise error( 'merging of model tables not allowed here' )

  # get more imaging parameters
  cell_size = get_pixel_size( facet, make_absolute = True )
  image_size = get_image_size( facet )
  uv_size = restore_parameter( uv, 'pb_image_size' )
  channel_count = get_channel_count( uv )
  if do_sdi_clean:
    sdi_param = 0.1
  else:
    sdi_param = 0.

  # copy facet and beam
  work_facet_file_name = sel_facet_file_name + '.%03d' % ( facet_id )
  extract_facet_definitions( sel_facet_file_name, [ facet_id ], work_facet_file_name )
  work_facet = get_aips_file( facet.disk, 'WORK', 'ICL001', -1, 'MA' )
  work_beam = get_facet_beam( work_facet )

  # re-image facet, while applying SN to UV
  call_aips_task( 'TACOP', indata = facet, inext = 'SN', invers = solution_version,
      ncount = 1, outdata = uv, outvers = 0 )

  if remake_beam:
    call_aips_task( 'IMAGR', indata = uv, nchav = channel_count, docalib = 100,
        gainuse = 0, outdisk = work_facet.disk, outname = work_facet.name,
        outseq = work_facet.seq, outver = model_version, in2disk = uv.disk,
        cellsize = cell_size, imsize = image_size, do3dimag = 1, niter = 1000000,
        flux = 0.99 * clean_flux, boxfile = work_facet_file_name, flagver = -1,
        cmethod = conversion_method, minpatch = image_size[ 0 ] - 1, maxpixel = 0, 
        overlap = 2, nfield = 1, allokay = 0, uvsize = [ uv_size, uv_size ], dotv = 0,
        imagrprm = [ dish_diameter, 0,0, sdi_param, 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0 ],
        bcomp = [ 0 for j in range( 64 ) ], **imagr_params )
    for table in facet.tables:
      table_version = table[ 0 ]
      table_type = ( table[ 1 ].split() )[ 1 ]
      if ( not table_type in [ 'HI', 'CG' ] ):
        call_aips_task( 'TACOP', indata = facet, outdata = work_facet, ncount = 1,
            inext = table_type, invers = table_version, outvers = table_version )

  else: # not remake_beam
    call_aips_task( 'MOVE', indata = facet, outdata = work_facet,
        userid = get_aips_userid() )
    call_aips_task( 'MOVE', indata = beam, outdata = work_beam,
        userid = get_aips_userid() )
    call_aips_task( 'IMAGR', indata = uv, nchav = channel_count, docalib = 100,
        gainuse = 0, outdisk = work_facet.disk, outname = work_facet.name,
        outseq = work_facet.seq, outver = model_version, in2disk = uv.disk,
        cellsize = cell_size, imsize = image_size, do3dimag = 1, niter = 1000000,
        flux = 0.95 * clean_flux, boxfile = work_facet_file_name, maxpixel = 0,
        cmethod = conversion_method, minpatch = image_size[ 0 ] - 1, dotv = 0,
        overlap = 2, nfield = 1, allokay = 1, uvsize = [ uv_size, uv_size ],
        imagrprm = [ dish_diameter, 0,0, sdi_param, 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0 ],
        bcomp = [ 0 for j in range( 64 ) ], flagver = -1, **imagr_params )
  
  # clean up
  uv.zap_table( 'SN', 0 )
  remove_file( work_facet_file_name )
  
  # remove cleaned flux
  model_component_count = get_model_component_count( work_facet, model_version )
  if ( subtract and ( model_component_count > 0 ) ):
    
    # subtract model from image
    scale_model_flux( work_facet, -1., model_version = model_version )
    rst_facet = restore_model_components( work_facet, facet_list = [ 1 ],
        cross_restore = False, imagr_params = imagr_params )
    work_facet.zap()
    rst_facet.rename( name = work_facet.name, klass = work_facet.klass,
        seq = work_facet.seq )
    work_facet.zap_table( 'CC', 0 )
    
    # subtract model from visibilities
    sub_uv = subtract_model( uv, work_facet, facet_list = [ 1 ],
        apply_solutions = True, solution_version = solution_version,
        keep_solutions = True, model_version = model_version,
        conversion_method = conversion_method, flag_solutions = False,
        frequency_correction = frequency_correction )
    uv.zap()
    sub_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  
  # replace original facet and beam
  facet.zap()
  if beam.exists():
    beam.zap()
  work_facet.rename( name = facet.name, klass = facet.klass, seq = facet.seq )
  work_beam.rename( name = beam.name, klass = beam.klass, seq = beam.seq )

  return

###############################################################################

def re_image_clean_pb_facets( uv, pb_facets, facet_list = [], sigma = 2.,
    conversion_method = 'DFT', apply_solutions = True, do_sdi_clean = False,
    restore_components = True, center_facets = False, print_info = True,
    frequency_correction = False, solution_version = 0, remake_beams = True,
    imagr_params = {}, model_version = -1, db_contrast = 0.25,
    shrink_factor = 0.975 ):
  
  # TODO: define proper flux limit for switching between facets
  
  # process parameters
  if frequency_correction:
    dish_diameter = restore_parameter( uv, 'dish_diameter' )
  else:
    dish_diameter = 0.
  cell_size = restore_parameter( uv, 'cell_size' )
  cell_size = [ cell_size, cell_size ]
  facet_size = restore_parameter( uv, 'pb_facet_size' )
  facet_size = [ facet_size, facet_size ]
  channel_count = get_channel_count( uv )
  uv_size = restore_parameter( uv, 'pb_image_size' )
  cpb_noise = restore_parameter( uv, 'cpb_noise' )
  flux = sigma * cpb_noise
  if do_sdi_clean:
    sdi_param = 0.1
  else:
    sdi_param = 0.
  
  # make copies of selected facets
  # rename selected facets to form consecutive series
  pb_facet_file_name = restore_parameter( pb_facets, 'facet_file_name' )
  pb_facet_file_name_e = path.expandvars( pb_facet_file_name )
  pb_facet_count = restore_parameter( pb_facets, 'facet_count' )
  sel_facets = get_aips_file( pb_facets.disk, 'SEL', 'ICL001', -1, 'MA' )
  if ( len( facet_list ) > 0 ):
    used_facet_list = facet_list
    sel_facet_count = len( facet_list )
  else:
    sel_facet_count = pb_facet_count
    used_facet_list = range( 1, 1 + sel_facet_count )
  sel_facet_list = range( 1, 1 + sel_facet_count )

  # check state of execution
  missing_facets = False
  missing_solutions = False
  missing_models = False
  restart_imaging = ( model_version >= 0 )
  facet_based_solutions = False
  for i in used_facet_list:
    pb_facet = get_facet( pb_facets, i )
    if ( not pb_facet.exists() ):
      missing_facets = True
      continue
    if table_exists( pb_facet, 'SN', solution_version ):
      facet_based_solutions = True
    else:
      missing_solutions = True
    if ( not table_exists( pb_facet, 'CC', model_version ) ):
      missing_models = True
  if ( apply_solutions and facet_based_solutions and 
      ( missing_facets or missing_solutions ) ):
    raise error( 'missing facet solution tables' )
  if ( restart_imaging and ( missing_facets or missing_models ) ):
    raise error( 'missing facet model tables' )
  
  # re-center added facets
  if center_facets:
    re_center_pb_facets( uv, pb_facets, facet_file_name = facet_file_name,
        facet_list = facet_list, imagr_params = imagr_params )

  # make work copies of facets
  sel_facet_file_name = pb_facet_file_name_e + '.SEL'
  extract_facet_definitions( pb_facet_file_name_e, used_facet_list, sel_facet_file_name )
  sel_facets = get_aips_file( pb_facets.disk, 'SEL', 'ICL001', -1, 'MA' )
  for i in sel_facet_list:
    j = used_facet_list[ i - 1 ]
    pb_facet = get_facet( pb_facets, j )
    pb_beam = get_facet_beam( pb_facet )
    facet = get_facet( sel_facets, i )
    beam = get_facet_beam( facet )
    call_aips_task( 'MOVE', indata = pb_facet, outdata = facet,
        userid = get_aips_userid() )
    call_aips_task( 'MOVE', indata = pb_beam, outdata = beam,
        userid = get_aips_userid() )
    for table in facet.tables:
      table_version = table[ 0 ]
      table_type = ( table[ 1 ].split() )[ 1 ]
      if ( not table_type in [ 'HI','CG','PS' ] ):
        facet.zap_table( table_type, table_version )
    if restart_imaging:
      call_aips_task( 'TACOP', indata = pb_facet, outdata = facet, inext = 'CC',
          invers = model_version, outvers = 1, ncount = 1 )
      new_model_version = 2
    else:
      new_model_version = 1
    if facet_based_solutions:
      call_aips_task( 'TACOP', indata = pb_facet, outdata = facet, inext = 'SN',
          invers = solution_version, outvers = 1, ncount = 1 )
  store_parameter( sel_facets, 'facet_file_name', sel_facet_file_name )
  if ( sel_facet_count < pb_facet_count ):
    store_parameter( sel_facets, 'facet_count', sel_facet_count )
    determine_facet_overlap( sel_facets )
  
  # solutions available per facet
  if facet_based_solutions:
    
    # make work copy of visibility data
    # flag bad solutions
    for i in sel_facet_list:
      if ( i == 1 ):
        flag_version = -1
      else:
        flag_version = 0
      facet = get_facet( sel_facets, i )
      flag_bad_solutions( uv, uvim = facet, apply_flags = False,
          flag_version = flag_version, solution_version = 1 )
    work_uv = apply_flag_table( uv )
    uv.zap_table( 'FG', 0 )
    # remove larger tables to speed up things
    while table_exists( work_uv, 'NI', 0 ):
      work_uv.zap_table( 'NI', 0 )
    while table_exists( work_uv, 'OB', 0 ):
      work_uv.zap_table( 'OB', 0 )
    while table_exists( work_uv, 'SN', 0 ):
      work_uv.zap_table( 'SN', 0 )
    while table_exists( work_uv, 'FG', 0 ):
      work_uv.zap_table( 'FG', 0 )
    
    # check for restart
    if restart_imaging:
      temp_uv = subtract_model( work_uv, sel_facets, apply_solutions = True,
          conversion_method = conversion_method, keep_solutions = True,
          frequency_correction = frequency_correction, flag_solutions = False,
          solution_version = 1, model_version = 1 )
      work_uv.zap()
      temp_uv.rename( name = work_uv.name, klass = work_uv.klass, seq = work_uv.seq )
    
    # dirty image all facets
    for i in sel_facet_list:
      re_image_clean_facet( work_uv, sel_facets, i, 100000., subtract = False,
          facet_file_name = sel_facet_file_name, remake_beam = remake_beams,
          do_sdi_clean = do_sdi_clean, frequency_correction = frequency_correction,
          model_version = new_model_version, conversion_method = conversion_method,
          solution_version = 1, imagr_params = imagr_params )
      if restart_imaging:
        facet = get_facet( sel_facets, i )
        facet.zap_table( 'CC', 2 )
    
    # determine list of extrema
    extremum_list = [ 0. ]
    temp_facets = get_aips_file( pb_facets.disk, 'TEMP', 'ICL001', -1, 'MA' )
    for i in sel_facet_list:
      facet = get_facet( sel_facets, i )
      temp_facet = get_facet( temp_facets, i )
      call_aips_task( 'MOVE', indata = facet, outdata = temp_facet,
          userid = get_aips_userid() )
    fill_facets( temp_facets, sel_facet_list, facet_file_name = sel_facet_file_name,
        do_boxes = True, invert = True, do_edge_circle = True )
    for i in sel_facet_list:
      temp_facet = get_facet( temp_facets, i )
      facet_extremum = get_image_extremum( temp_facet )
      if ( facet_extremum is None ):
        extremum_list.append( -100000. )
      else:
        extremum_list.append( abs( facet_extremum[ 0 ] ) )
      temp_facet.zap()
    
    # clean and re-image facets
    total_components = 0
    total_flux = 0.
    old_i = -1
    old_clean_flux = 1.e9
    while ( max( extremum_list ) > flux ):
      
      # get facet with max_abs_residual
      # determine facet min clean flux level
      i = extremum_list.index( max( extremum_list ) )
      clean_flux = float( max( [ flux, db_contrast * max( extremum_list ) ] ) )
      if print_info:
        print '... processing facet %s with absolute peak flux %s Jy' % (
            repr( i ), repr( max( extremum_list ) ) )
        print '... cleaning down to %s Jy' % ( repr( clean_flux ) )
      
      # check for endless loops
      if ( ( i == old_i ) and ( clean_flux == old_clean_flux ) ):
        if print_info:
          print '...... WARNING: detecting endless loop, skipping facet '
        new_i = extremum_list.index( max( extremum_list[ 0 : i ] + 
            extremum_list[ i + 1 : ] ) )
        extremum_list[ i ] = shrink_factor * extremum_list[ new_i ]
        continue
      old_i = i
      old_clean_flux = clean_flux
      
      # re-image and clean facet
      re_image_clean_facet( work_uv, sel_facets, i, clean_flux, subtract = True,
          facet_file_name = sel_facet_file_name, remake_beam = False,
          do_sdi_clean = do_sdi_clean, frequency_correction = frequency_correction,
          model_version = 2, conversion_method = conversion_method,
          solution_version = 1, imagr_params = imagr_params )
      
      # some bookkeeping
      facet = get_facet( sel_facets, i )
      dcomponents = get_model_component_count( facet, model_version = 2 )
      dflux = get_model_flux( facet, model_version = 2 )
      total_components = total_components + dcomponents
      total_flux = total_flux + dflux
      if print_info:
        print '... total # components = %s' % ( repr( total_components ) )
        print '... total cleaned flux = %s' % ( repr( total_flux ) )
      
      # merge clean components
      combine_model_tables( facet, model_versions = [ 1,2 ] )
      facet.zap_table( 'CC', 1 )
      facet.zap_table( 'CC', 2 )
      call_aips_task( 'TACOP', indata = facet, outdata = facet, inext = 'CC',
          invers = 3, outvers = 1, ncount = 1 )
      facet.zap_table( 'CC', 3 )
      
      # determine new extremum
      temp_facet = get_facet( temp_facets, i )
      call_aips_task( 'MOVE', indata = facet, outdata = temp_facet,
          userid = get_aips_userid() )
      fill_facets( temp_facets, [ i ], facet_file_name = sel_facet_file_name,
          do_boxes = True, invert = True, do_edge_circle = True )
      [ facet_extremum, pos ] = get_image_extremum( temp_facet )
      extremum_list[ i ] = abs( facet_extremum )
      temp_facet.zap()
      if print_info:
        sel = awhere( array( extremum_list ) > flux )
        print '... %s facets remain to be cleaned' % ( len( sel ) )
    
    # when cleaning is done, re-image all residual facets with solutions applied
    if ( total_components > 0 ):
      for i in sel_facet_list:
        re_image_clean_facet( work_uv, sel_facets, i, 100000., subtract = False,
            facet_file_name = sel_facet_file_name, remake_beam = remake_beams,
            do_sdi_clean = do_sdi_clean, frequency_correction = frequency_correction,
            model_version = 2, conversion_method = conversion_method,
            solution_version = 1, imagr_params = imagr_params )
        facet = get_facet( sel_facets, i )
        facet.zap_table( 'CC', 2 )
    
    # restore clean components to facets
    if restore_components:
      rst_facets = restore_model_components( sel_facets, facet_list = sel_facet_list, 
          imagr_params = imagr_params, model_version = 1 )
      for i in sel_facet_list:
        facet = get_facet( sel_facets, i )
        rst_facet = get_facet( rst_facets, i )
        facet.zap()
        rst_facet.rename( name = facet.name, klass = facet.klass, seq = facet.seq )
  
  else: # not facet_based_solutions
    
    # set calibration parameters
    # make work copy of visibility data
    if apply_solutions:
      docalib = 100
      work_uv = flag_bad_solutions( uv, solution_version = solution_version,
          flag_version = -1, apply_flags = True )
    else:
      docalib = -1
      work_uv = get_aips_file( uv.disk, uv.name, 'WORK', -1, 'UV' )
      call_aips_task( 'MOVE', indata = uv, outdata = work_uv, userid = get_aips_userid() )
    
    # check for restart
    if restart_imaging:
      temp_uv = subtract_model( work_uv, sel_facets, keep_solutions = True,
          conversion_method = conversion_method, apply_solutions = apply_calibration,
          frequency_correction = frequency_correction, flag_solutions = False,
          solution_version = solution_version, model_version = 1 )
      work_uv.zap()
      temp_uv.rename( name = work_uv.name, klass = work_uv.klass, seq = work_uv.seq )
      
    # re-image all facets in one run
    if ( remake_beams or missing_facets ):
      temp_facets = get_aips_file( pb_facets.disk, 'TEMP', 'ICL001', -1, 'MA' )
      call_aips_task( 'IMAGR', indata = work_uv, nchav = channel_count,
          docalib = docalib, gainuse = solution_version, outdisk = temp_facets.disk,
          outname = temp_facets.name, outseq = temp_facets.seq, do3dimag = 1,
          outver = new_model_version, imsize = facet_size, in2disk = uv.disk, 
          cmethod = conversion_method, minpatch = facet_size[ 0 ] - 1, allokay = 0, 
          flux = 0.99 * flux, boxfile = sel_facet_file_name, dotv = 0, overlap = 2,
          nfield = sel_facet_count, cellsize = cell_size, flagver = -1, maxpixel = 0, 
          imagrprm = [ dish_diameter, 0,0, sdi_param, 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0 ], 
          bcomp = [ 0 for j in range( 64 ) ], uvsize = [ uv_size, uv_size ],
          niter = 1000000, **imagr_params )
      
      # replace facets
      for i in sel_facet_list:
        facet = get_facet( sel_facets, i )
        beam = get_facet_beam( facet )
        temp_facet = get_facet( temp_facets, i )
        temp_beam = get_facet_beam( temp_facet )
        call_aips_task( 'TACOP', indata = facet, outdata = temp_facet,
            inext = 'PS', ncount = 1 )
        if restart_imaging:
          call_aips_task( 'TACOP', indata = facet, outdata = temp_facet,
              inext = 'CC', invers = 1, outvers = 1, ncount = 1 )
        if facet.exists():
          facet.zap()
        if beam.exists():
          beam.zap()
        temp_facet.rename( name = facet.name, klass = facet.klass, seq = facet.seq )
        temp_beam.rename( name = beam.name, klass = beam.klass, seq = beam.seq )
    
    else: # not remake_beams
      call_aips_task( 'IMAGR', indata = work_uv, nchav = channel_count,
          docalib = docalib, gainuse = solution_version, outdisk = sel_facets.disk,
          outname = sel_facets.name, outseq = sel_facets.seq, outver = new_model_version,
          imsize = facet_size, in2disk = uv.disk, do3dimag = 1, niter = 1000000,
          flux = 0.95 * flux, boxfile = sel_facet_file_name, cellsize = cell_size,
          cmethod = conversion_method, minpatch = facet_size[ 0 ] - 1, flagver = -1,
          allokay = 1, dotv = 0, overlap = 2, nfield = sel_facet_count, maxpixel = 0,
          imagrprm = [ dish_diameter, 0,0, sdi_param, 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0 ],
          bcomp = [ 0 for j in range( 64 ) ], uvsize = [ uv_size, uv_size ],
          **imagr_params )
    
    # restore model components
    if ( not restore_components ):
      for i in sel_facet_list:
        facet = get_facet( sel_facets, i )
        scale_model_flux( facet, -1. )
      rst_facets = restore_model_components( sel_facets, facet_list = sel_facet_list, 
          imagr_params = imagr_params )
    elif restart_imaging:
      rst_facets = restore_model_components( sel_facets, facet_list = sel_facet_list, 
          imagr_params = imagr_params, model_version = 1 )
    
    # merge models, replace facets
    if ( ( not restore_components ) or restart_imaging ):
      for i in sel_facet_list:
        facet = get_facet( sel_facets, i )
        rst_facet = get_facet( rst_facets, i )
        if ( not restore_components ):
          facet.zap_table( 'CC', 3 )
        while table_exists( rst_facet, 'CC', 0 ):
          rst_facet.zap_table( 'CC', 0 )
        if restart_imaging:
          combine_model_tables( facet, model_versions = [ 1,2 ] )
          call_aips_task( 'TACOP', indata = facet, outdata = rst_facet, inext = 'CC', 
              invers = 3, outvers = 1, ncount = 1 )
        facet.zap()
        rst_facet.rename( name = facet.name, klass = facet.klass, seq = facet.seq )
  
  # clean up
  work_uv.zap()

  # merge facets with originals
  for i in sel_facet_list:
    j = used_facet_list[ i - 1 ]
    facet = get_facet( sel_facets, i )
    beam = get_facet_beam( facet )
    pb_facet = get_facet( pb_facets, j )
    pb_beam = get_facet_beam( pb_facet )
    if table_exists( facet, 'PS', 0 ):
      facet.zap_table( 'PS', 0 )
    for table in pb_facet.tables:
      table_version = table[ 0 ]
      table_type = ( table[ 1 ].split() )[ 1 ]
      if ( not table_type in [ 'HI','CG','CC','SN' ] ):
        call_aips_task( 'TACOP', indata = pb_facet, outdata = facet,
            inext = table_type, invers = table_version, outvers = 0, ncount = 1 )
    pb_facet.zap()
    pb_beam.zap()
    facet.rename( name = pb_facet.name, klass = pb_facet.klass, seq = pb_facet.seq )
    beam.rename( name = pb_beam.name, klass = pb_beam.klass, seq = pb_beam.seq )
  remove_file( sel_facet_file_name )
  if missing_facets:
    determine_facet_overlap( pb_facets )
  
  # determine noise level
  store_parameter( pb_facets, 'clean_flux_min', flux )
  cpb_noise = measure_cpb_noise( uv, pb_facets )
  
  return

###############################################################################

def selfcal_pb_facets( uv, pb_facets, selfcal_cycle_min = 2, phase_interval = 0.,
    reference_antenna = 0, improvement_limit = 0.05, try_final_amplitude = False,
    amplitude_interval = 5., interpolation_method = '', conversion_method = 'DFT',
    calib_params = {}, frequency_correction = False, normalize_gains = True,
    snr_limit = 2., add_boxes = True, box_sigma = 5., facet_based_boxes = True,
    keep_boxes = False, sidelobe_rejection = 1.05, clean_box_radius = 5, 
    do_sdi_clean = False, imagr_params = {}, remake_beams = True, sigma = 2.,
    restore_components = True, center_facets = False, print_info = True ):
  
  cpb_noise = restore_parameter( uv, 'cpb_noise' )
  noise_improvement = 1.
  selfcal_i = 0
  while ( selfcal_i < selfcal_cycle_min ) or ( noise_improvement > improvement_limit ):
    selfcal_i = selfcal_i + 1
    last_cpb_noise = cpb_noise

    # selfcal and image pb
    if ( add_boxes and ( selfcal_i > 1 ) ):
      add_clean_boxes( uv, pb_facets, sidelobe_rejection = sidelobe_rejection,
          facet_based_boxes = facet_based_boxes, clean_box_radius = clean_box_radius, 
           print_info = print_info, box_sigma = box_sigma, keep_boxes = keep_boxes )
    calibrate_model( uv, pb_facets, print_info = print_info, snr_limit = snr_limit,
        reference_antenna = reference_antenna, conversion_method = conversion_method,
        frequency_correction = frequency_correction, phase_interval = phase_interval,
        calib_params = calib_params, apply_solutions = False,
        interpolation_method = interpolation_method )
    re_image_clean_pb_facets( uv, pb_facets, sigma = sigma, imagr_params = imagr_params,
        do_sdi_clean = do_sdi_clean, frequency_correction = frequency_correction, 
        conversion_method = conversion_method, restore_components = restore_components,
        apply_solutions = True, center_facets = ( center_facets and ( selfcal_i > 1 ) ),
        print_info = print_info, remake_beams = ( remake_beams and ( selfcal_i == 1 ) ) )

    # determine noise reduction
    cpb_noise = restore_parameter( uv, 'cpb_noise' )
    noise_improvement = 1. - ( cpb_noise / last_cpb_noise )

    # only keep solutions when noise improvement was detected
    solution_version = uv.table_highver( 'SN' )
    if ( selfcal_i == 1 ):
      last_solution_version = solution_version
    elif ( noise_improvement >= 0. ):
      last_solution_version = solution_version
    else:
      while ( uv.table_highver( 'SN' ) > last_solution_version ):
        uv.zap_table( 'SN', 0 )
      solution_version = last_solution_version
  
  if try_final_amplitude:
    selfcal_i = selfcal_i + 1
    last_cpb_noise = cpb_noise
    
    # selfcal and image pb
    ref_ant = get_reference_antenna( uv  )
    calibrate_model( uv, pb_facets, print_info = print_info, snr_limit = snr_limit,
        reference_antenna = ref_ant, conversion_method = conversion_method,
        frequency_correction = frequency_correction, phase_interval = phase_interval,
        calib_params = calib_params, apply_solutions = True, do_amplitude = True,
        amplitude_interval = amplitude_interval, keep_flags = False,
        interpolation_method = interpolation_method, normalize_gains = normalize_gains )
    combine_solutions( uv, in_version_1 = solution_version )
    re_image_clean_pb_facets( uv, pb_facets, sigma = sigma, imagr_params = imagr_params,
        do_sdi_clean = do_sdi_clean, frequency_correction = frequency_correction, 
        conversion_method = conversion_method, restore_components = restore_components,
        apply_solutions = True, center_facets = center_facets, remake_beams = False,
        print_info = print_info )
    
    # determine noise reduction
    cpb_noise = restore_parameter( temp_uv, 'cpb_noise' )
    noise_improvement = 1. - ( cpb_noise / last_cpb_noise )
    
    # only keep solutions when noise improvement was detected
    if ( noise_improvement < 0. ):
      while ( uv.table_highver( 'SN' ) > last_solution_version ):
        uv.zap_table( 'SN', 0 )
  
  # only re-image in case of bad degradation
  if ( noise_improvement < -improvement_limit ):
    re_image_clean_pb_facets( uv, pb_facets, sigma = sigma, imagr_params = imagr_params,
        do_sdi_clean = do_sdi_clean, frequency_correction = frequency_correction, 
        conversion_method = conversion_method, restore_components = restore_components,
        apply_solutions = True, center_facets = center_facets, remake_beams = False,
        print_info = print_info )
  
  return

###############################################################################

def merge_uv( up_uv, down_uv, keep_solutions = True ):

  new_uv = get_aips_file( up_uv.disk, up_uv.name, 'MERGE', - 1, 'UV' )
  call_aips_task( 'DBCON', indata = up_uv, in2data = down_uv, outdata = new_uv, doarray = 1 )

  # restore solution table
  while table_exists( new_uv, 'SN', 0 ):
    new_uv.zap_table( 'SN', 0 )
  if ( keep_solutions and ( not up_uv is None ) ):
    if table_exists( up_uv, 'SN', 0 ):
      call_aips_task( 'TACOP', indata = up_uv, inext = 'SN', invers = 0, ncount = 1,
          outdata = new_uv, outvers = 0 )

  return new_uv

###############################################################################

def split_uv_on_time_range( uv, time_rts, keep_solutions = True, time_list = [] ):
  
  # retrieve relevant time ranges
  time_rise = time_rts[ 0 ]
  time_transit = time_rts[ 1 ]
  time_set = time_rts[ 2 ]
  if ( len( time_list ) > 0 ):
    time_array = array( time_list )
    time_min = time_array.min()
    time_max = time_array.max()
  else:
    time_min = restore_parameter( uv, 'time_min' )
    time_max = restore_parameter( uv, 'time_max' )
  
  # determine what time ranges source is visible during observing run
  ud_list = []
  if ( time_rise == time_set ):
    if ( time_rise == time_transit ): # source never visible
      ud_list.append( [ 'DOWN', time_min, time_max ] )
    else: # source always visible
      ud_list.append( [ 'UP', time_min, time_max ] )
  else:
    time_up_before_min = time_min - ( ( time_min - time_rise ) % 1. )
    time_down_before_min = time_min - ( ( time_min - time_set ) % 1. )
    time_up_after_max = time_max + ( ( time_rise - time_max ) % 1. )
    time_down_after_max = time_max + ( ( time_set - time_max ) % 1. )
    # source is visible at start of observing run
    if ( time_up_before_min > time_down_before_min ): 
      time_up = time_up_before_min
      time_down = time_down_before_min + 1.
      while ( time_up < time_max ):
        ud_list.append( [ 'UP', max( [ time_up, time_min ] ), 
            min( [ time_down, time_max ] ) ] )
        time_up = time_up + 1.
        if ( time_down < time_max ):
          ud_list.append( [ 'DOWN', time_down, min( [ time_up, time_max ] ) ] )
        time_down = time_down + 1.
    # source is not visible at start of observing run
    else: 
      time_down = time_down_before_min
      time_up = time_up_before_min + 1.
      while ( time_down < time_max ):
        ud_list.append( [ 'DOWN', max( [ time_down, time_min ] ), 
            min( [ time_up, time_max ] ) ] )
        time_down = time_down + 1.
        if ( time_up < time_max ):
          ud_list.append( [ 'UP', time_up, min( [ time_down, time_max ] ) ] )
        time_up = time_up + 1.
  
  # split all up-times and down-times
  up_count = 0
  down_count = 0
  for row in ud_list:
    if ( len( time_list ) > 0 ):
      sel = awhere( ( time_array >= row[ 1 ] ) & ( time_array <= row[ 2 ] ) )
      if ( len( sel ) == 0 ):
        continue
    if ( row[ 0 ] == 'UP' ):
      up_uv = get_aips_file( uv.disk, uv.name, 'UP%02d' % ( up_count + 1 ), -1,
          'UV' )
      try:
        call_aips_task( 'UVCOP', indata = uv, outdata = up_uv, 
            uvcopprm = [ 1,0,0,0,0,3,0 ], timerang = time_to_dhms( row[ 1 ] ) +
            time_to_dhms( row[ 2 ] ) )
      except: # no data selected
        pass
      else:
        up_count = up_count + 1
    else: # row[ 0 ] == 'DOWN'
      down_uv = get_aips_file( uv.disk, uv.name, 'DOWN%02d' % ( down_count + 1 ),
          -1, 'UV' )
      try:
        call_aips_task( 'UVCOP', indata = uv, outdata = down_uv,
            uvcopprm = [ 1,0,0,0,0,3,0 ], timerang = time_to_dhms( row[ 1 ] ) +
            time_to_dhms( row[ 2 ] ) )
      except: # no data selected
        pass
      else:
        down_count = down_count + 1
  
  # combine up-times
  if ( up_count > 0 ):
    new_up_uv = get_aips_file( uv.disk, uv.name, 'UP', - 1, 'UV' )
    up_uv = get_aips_file( uv.disk, uv.name, 'UP01', 0, 'UV' )
    up_uv.rename( name = new_up_uv.name, klass = new_up_uv.klass,
        seq = new_up_uv.seq )
    if ( up_count > 1 ):
      for i in range( 2, up_count + 1 ):
        next_up_uv = get_aips_file( up_uv.disk, up_uv.name, 'UP%02d' % ( i ),
            up_uv.seq, 'UV' )
        merged_uv = merge_uv( new_up_uv, next_up_uv )
        new_up_uv.zap()
        next_up_uv.zap()
        merged_uv.rename( name = new_up_uv.name, klass = new_up_uv.klass,
            seq = new_up_uv.seq )
  else:
    new_up_uv = None
  
  # combine down-times
  if ( down_count > 0 ):
    new_down_uv = get_aips_file( uv.disk, uv.name, 'DOWN', - 1, 'UV' )
    down_uv = get_aips_file( uv.disk, uv.name, 'DOWN01', 0, 'UV' )
    down_uv.rename( name = new_down_uv.name, klass = new_down_uv.klass,
        seq = new_down_uv.seq )
    if ( down_count > 1 ):
      for i in range( 2, down_count + 1 ):
        next_down_uv = get_aips_file( new_down_uv.disk, down_uv.name, 
            'DOWN%02d' % ( i ), down_uv.seq, 'UV' )
        merged_uv = merge_uv( new_down_uv, next_down_uv )
        new_down_uv.zap()
        next_down_uv.zap()
        merged_uv.rename( name = new_down_uv.name, klass = new_down_uv.klass,
            seq = new_down_uv.seq )
  else:
    new_down_uv = None
  
  # remove solution tables, as they are messed up
  if ( not new_up_uv is None ):
    while table_exists( new_up_uv, 'SN', 0 ):
      new_up_uv.zap_table( 'SN', 0 )
  if ( not new_down_uv is None ):
    while table_exists( new_down_uv, 'SN', 0 ):
      new_down_uv.zap_table( 'SN', 0 )
  
  # restore original solution table
  if ( keep_solutions and ( table_exists( uv, 'SN', 0 ) ) ):
    if ( not new_up_uv is None ):
      call_aips_task( 'TACOP', indata = uv, inext = 'SN', invers = 0, ncount = 1,
          outdata = new_up_uv, outvers = 0 )
    if ( not new_down_uv is None ):
      call_aips_task( 'TACOP', indata = uv, inext = 'SN', invers = 0, ncount = 1,
          outdata = new_down_uv, outvers = 0 )
  
  return [ new_up_uv, new_down_uv ]

###############################################################################

def subtract_as_model( uv, facets, facet_list = [], sigma = 0., apply_solutions = False,
    keep_solutions = False, conversion_method = 'DFT', model_version = 0 ):

  if ( len( facet_list ) == 0 ):
    used_facet_count = restore_parameter( facets, 'facet_count' )
    used_facet_list = range( 1, 1 + used_facet_count )
  else:
    used_facet_count = len( facet_list )
    used_facet_list = [ i for i in facet_list ]
  time_list = get_time_list( uv )

  # subtract source models over correct time ranges
  sub_uv = get_aips_file( uv.disk, uv.name, 'SUB', - 1, 'UV' )
  call_aips_task( 'MOVE', indata = uv, outdata = sub_uv, userid = get_aips_userid() )
  for i in used_facet_list:
    facet_i = get_facet( facets, i )
    if ( abs( get_model_flux( facet_i ) ) > 0. ):
      rts_times = calculate_rise_transit_set_times( sub_uv, radec = get_radec( facet_i ) )
      [ up_uv, down_uv ] = split_uv_on_time_range( sub_uv, rts_times, 
          time_list = time_list )
      if ( not up_uv is None ):
        sub_up_uv = subtract_model( up_uv, facets, facet_list = [ i ], sigma = sigma,
            apply_solutions = apply_solutions, keep_solutions = keep_solutions,
            conversion_method = conversion_method, model_version = model_version,
            frequency_correction = False )
        sub_uv.zap()
        up_uv.zap()
        if ( not down_uv is None ):
          sub_up_down_uv = merge_uv( sub_up_uv, down_uv )
          down_uv.zap()
          sub_up_uv.zap()
          sub_up_down_uv.rename( name = sub_uv.name, klass = sub_uv.klass, seq = sub_uv.seq )
        else:
          sub_up_uv.rename( name = sub_uv.name, klass = sub_uv.klass, seq = sub_uv.seq )
      elif ( not down_uv is None ):
        down_uv.zap()
      
  return sub_uv

###############################################################################

def add_as_model( uv, facets, facet_list = [], sigma = 0., apply_solutions = False,
    keep_solutions = False, conversion_method = 'DFT', model_version = 0 ):

  if ( len( facet_list ) == 0 ):
    used_facet_count = restore_parameter( facets, 'facet_count' )
    used_facet_list = range( 1, 1 + used_facet_count )
  else:
    used_facet_count = len( facet_list )
    used_facet_list = [ i for i in facet_list ]
  time_list = get_time_list( uv )

  # subtract source models over correct time ranges
  add_uv = get_aips_file( uv.disk, uv.name, 'ADD', - 1, 'UV' )
  call_aips_task( 'MOVE', indata = uv, outdata = add_uv, userid = get_aips_userid() )
  for i in used_facet_list:
    facet_i = get_facet( facets, i )
    if ( abs( get_model_flux( facet_i ) ) > 0. ):
      rts_times = calculate_rise_transit_set_times( add_uv, radec = get_radec( facet_i ) )
      [ up_uv, down_uv ] = split_uv_on_time_range( add_uv, rts_times,
          time_list = time_list )
      if ( not up_uv is None ):
        add_up_uv = add_model( up_uv, facets, facet_list = [ i ], sigma = sigma,
            apply_solutions = apply_solutions, keep_solutions = keep_solutions,
            conversion_method = conversion_method, model_version = model_version,
            frequency_correction = False )
        add_uv.zap()
        up_uv.zap()
        if ( not down_uv is None ):
          add_up_down_uv = merge_uv( add_up_uv, down_uv )
          down_uv.zap()
          add_up_uv.zap()
          add_up_down_uv.rename( name = add_uv.name, klass = add_uv.klass, seq = add_uv.seq )
        else:
          add_up_uv.rename( name = add_uv.name, klass = add_uv.klass, seq = add_uv.seq )
      elif ( not down_uv is None ):
        down_uv.zap()

  return add_uv

###############################################################################

def image_clean_o_facets( uv, apply_solutions = False, sigma = 3., box_sigma = 5.,
    facet_based_boxes = True, sidelobe_rejection = 1.05, conversion_method = 'DFT', 
    do_sdi_clean = False, imagr_params = {}, print_info = True, clean_box_radius = 5 ):
  
  # make outlier facets
  # assume that outlier facets are above horizon during observations
  cell_size = restore_parameter( uv, 'cell_size' )
  uv_size = restore_parameter( uv, 'pb_image_size' )
  o_facet_count = restore_parameter( uv, 'o_facet_count' )
  o_facet_size = restore_parameter( uv, 'o_facet_size' )
  o_facet_file_name = restore_parameter( uv, 'o_facet_file_name' )
  o_facet_file_name_e = path.expandvars( o_facet_file_name )
  channel_count = get_channel_count( uv )
  o_facets = get_aips_file( uv.disk, 'O', 'ICL001', -1, 'MA' )
  cpb_noise = restore_parameter( uv, 'cpb_noise' )
  flux = sigma * cpb_noise
  if apply_solutions:
    solution_switch = 100
    solution_version = uv.table_highver( 'SN' )
  else:
    solution_switch = -1
    solution_version = -1
  if do_sdi_clean:
    sdi_param = 0.1
  else:
    sdi_param = 0.
  
  # image dirty facets
  call_aips_task( 'IMAGR', indata = uv, nchav = channel_count, dotv = 0, overlap = 2,
      docalib = solution_switch, gainuse = solution_version, nfield = o_facet_count,
      outdisk = o_facets.disk, outname = o_facets.name, outseq = o_facets.seq,
      do3dimag = 1, allokay = 0, outver = 0, maxpixel = 0, minpatch = o_facet_size - 1,
      cellsize = [ cell_size, cell_size ], imsize = [ o_facet_size, o_facet_size ],
      niter = 1000000, flux = 100000., boxfile = o_facet_file_name_e, flagver = -1,
      in2disk = uv.disk, cmethod = conversion_method, bcomp = [ 0 for j in range( 64 ) ], 
      imagrprm = [ 0,0,0, sdi_param, 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0 ], 
      uvsize = [ uv_size, uv_size ], **imagr_params )
  
  # store parameters
  store_parameter( o_facets, 'facet_count', o_facet_count )
  store_parameter( o_facets, 'facet_file_name', o_facet_file_name )
  
  # automatically generate clean boxes
  determine_facet_overlap( o_facets )
  add_clean_boxes( uv, o_facets, sidelobe_rejection = sidelobe_rejection,
      facet_based_boxes = facet_based_boxes, clean_box_radius = clean_box_radius, 
      print_info = print_info, box_sigma = box_sigma )
  
  # remove old model tables
  for i in range( 1, o_facet_count + 1 ):
    o_facet_i = get_facet( o_facets, i )
    o_facet_i.zap_table( 'CC', 0 )
  # re-image and clean facets
  call_aips_task( 'IMAGR', indata = uv, nchav = channel_count, docalib = solution_switch,
      gainuse = solution_version, outdisk = o_facets.disk, outname = o_facets.name,
      outseq = o_facets.seq, do3dimag = 1, nfield = o_facet_count, allokay = 1,
      outver = 0, cellsize = [ cell_size, cell_size ], boxfile = o_facet_file_name_e,
      imsize = [ o_facet_size, o_facet_size ], niter = 1000000, flux = 0.95 * flux,
      in2disk = uv.disk, cmethod = conversion_method, minpatch = o_facet_size - 1,
      dotv = 0, overlap = 2, bcomp = [ 0 for j in range( 64 ) ], maxpixel = 0,
      imagrprm = [ 0,0,0, sdi_param, 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0 ], flagver = -1,
      uvsize = [ uv_size, uv_size ], **imagr_params )
  
  return o_facets
  
###############################################################################

def image_clean_s_facets( uv, apply_solutions = False, sigma = 3., box_sigma = 5.,
    facet_based_boxes = True, sidelobe_rejection = 1.05, conversion_method = 'GRID',
    imagr_params = {}, do_sdi_clean = False, print_info = True,  clean_box_radius = 5 ):
  
  # make initial Sun facets
  cell_size = restore_parameter( uv, 'cell_size' )
  uv_size = restore_parameter( uv, 'pb_image_size' )
  s_facet_count = restore_parameter( uv, 's_facet_count' )
  if s_facet_count > 1:
    raise error( 'currently only one Sun facet supported' )
  s_facet_size = restore_parameter( uv, 's_facet_size' )
  s_facet_file_name = restore_parameter( uv, 's_facet_file_name' )
  s_facet_file_name_e = path.expandvars( s_facet_file_name ) 
  channel_count = get_channel_count( uv )
  s_facets = get_aips_file( uv.disk, 'S', 'ICL001', -1, 'MA' )
  cpb_noise = restore_parameter( uv, 'cpb_noise' )
  flux = sigma * cpb_noise
  if apply_solutions:
    solution_switch = 100
    solution_version = 0 # uv.table_highver( 'SN' )
  else:
    solution_switch = -1
    solution_version = -1
  if do_sdi_clean:
    sdi_param = 0.1
  else:
    sdi_param = 0.
  time_list = get_time_list( uv )
  
  # image dirty facets
  call_aips_task( 'IMAGR', indata = uv, nchav = channel_count, docalib = solution_switch,
      gainuse = solution_version, outdisk = s_facets.disk, outname = s_facets.name,
      outseq = s_facets.seq, do3dimag = 1, nfield = s_facet_count, allokay = 0,
      outver = 0, cellsize = [ cell_size, cell_size ], boxfile = s_facet_file_name_e,
      imsize = [ s_facet_size, s_facet_size ], niter = 1000000, flux = 1000000.,
      in2disk = uv.disk, cmethod = conversion_method, minpatch = s_facet_size - 1,
      dotv = 0, overlap = 2, bcomp = [ 0 for j in range( 64 ) ], maxpixel = 0,
      imagrprm = [ 0,0,0, sdi_param, 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0 ], flagver = -1,
      uvsize = [ uv_size, uv_size ], **imagr_params )
  
  # store parameters
  store_parameter( s_facets, 'facet_count', s_facet_count )
  store_parameter( s_facets, 'facet_file_name', s_facet_file_name )
  determine_facet_overlap( s_facets )
  
  # remove empty model tables
  for i in range( 1, s_facet_count + 1 ):
    s_facet_i = get_facet( s_facets, i )
    s_facet_i.zap_table( 'CC', 0 )
  
  # re-image and clean Sun facet using UV data only while above the horizon
  rts_times = calculate_rise_transit_set_times( uv, radec = get_radec( s_facets ) )
  [ uv_up, uv_down ] = split_uv_on_time_range( uv, rts_times,
      time_list = time_list )
  if ( not uv_up is None ):
    call_aips_task( 'IMAGR', indata = uv_up, nchav = channel_count, flagver = -1,
        docalib = solution_switch, gainuse = solution_version, do3dimag = 1,
        outdisk = s_facets.disk, outname = s_facets.name, outseq = s_facets.seq,
        nfield = s_facet_count, allokay = 1, outver = 0, niter = 1000000,
        cellsize = [ cell_size, cell_size ], imsize = [ s_facet_size, s_facet_size ],
        flux = 1000000., boxfile = s_facet_file_name_e, in2disk = uv_up.disk,
        cmethod = conversion_method, minpatch = s_facet_size - 1, dotv = 0,
        overlap = 2, bcomp = [ 0 for j in range( 64 ) ], maxpixel = 0,
        imagrprm = [ 0,0,0, sdi_param, 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0 ],
        uvsize = [ uv_size, uv_size ], **imagr_params )
    uv_up.zap()
  else:
    pass # leave the facet as is
  if ( not uv_down is None ):
    uv_down.zap()
  
  # automatically generate clean boxes
  add_clean_boxes( uv, s_facets, sidelobe_rejection = sidelobe_rejection,
      facet_based_boxes = facet_based_boxes, clean_box_radius = clean_box_radius, 
       print_info = print_info, box_sigma = box_sigma )
  
  # re-image and clean Sun facet using UV data only while above the horizon
  rts_times = calculate_rise_transit_set_times( uv, radec = get_radec( s_facets ) )
  [ uv_up, uv_down ] = split_uv_on_time_range( uv, rts_times,
      time_list = time_list )
  if ( not uv_up is None ):
    # remove old model tables
    for i in range( 1, s_facet_count + 1 ):
      s_facet_i = get_facet( s_facets, i )
      s_facet_i.zap_table( 'CC', 0 )
    call_aips_task( 'IMAGR', indata = uv_up, nchav = channel_count, flagver = -1,
        docalib = solution_switch, gainuse = solution_version, outver = 0,
        outdisk = s_facets.disk, outname = s_facets.name, outseq = s_facets.seq, 
        cellsize = [ cell_size, cell_size ], imsize = [ s_facet_size, s_facet_size ],
        do3dimag = 1, niter = 1000000, flux = 0.95 * flux, boxfile = s_facet_file_name_e,
        in2disk = uv_up.disk,  nfield = s_facet_count, cmethod = conversion_method,
        minpatch = s_facet_size - 1, dotv = 0, overlap = 2, allokay = 1, maxpixel = 0, 
        imagrprm = [ 0,0,0, sdi_param, 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0 ],
        uvsize = [ uv_size, uv_size ], **imagr_params )
    uv_up.zap()
  else:
    pass # leave the facet as is
  if ( not uv_down is None ):
    uv_down.zap()
  
  return s_facets

###############################################################################

def image_s_facets_time_series( uv, apply_solutions = False, min_fraction = 0.5,
    imagr_params = {}, print_info = True, time_step = 1. ):
  
  # get imaging parameters
  cell_size = restore_parameter( uv, 'cell_size' )
  uv_size = restore_parameter( uv, 'pb_image_size' )
  s_facet_count = restore_parameter( uv, 's_facet_count' )
  if s_facet_count > 1:
    raise error( 'currently only one Sun facet supported' )
  integration_time = restore_parameter( uv, 'integration_time' )
  s_facet_size = restore_parameter( uv, 's_facet_size' )
  s_facet_file_name = restore_parameter( uv, 's_facet_file_name' )
  s_facet_file_name_e = path.expandvars( s_facet_file_name ) 
  channel_count = get_channel_count( uv )
  s_facets = get_aips_file( uv.disk, 'S', 'ITM001', -1, 'MA' )
#  cpb_noise = restore_parameter( uv, 'cpb_noise' )
#  flux = sigma * cpb_noise
  if apply_solutions:
    solution_switch = 100
    solution_version = 0 # uv.table_highver( 'SN' )
  else:
    solution_switch = -1
    solution_version = -1
  conversion_method = 'GRID'
  sdi_param = 0.
  time_array = array( get_time_list( uv ), dtype = float64 )
  if print_info:
    print 'imaging Sun every %s minutes over time range %s - %s' % (
        repr( time_step ), time_to_string( time_array[ 0 ] ), 
        time_to_string( time_array[ -1 ] ) )
  
  # make time series of images
  t = 0
  i = 0
  dtime = time_step / 1440.
  ditime = integration_time / 86400.
  min_sel = int( ceil( min_fraction * dtime / ditime ) )
  temp_facets = get_aips_file( uv.disk, 'TEMP', 'ICL001', -1, 'MA' )
  while ( t < len( time_array ) ):
    time_min = time_array[ t ] - 0.5 * ditime
    time_max = time_array[ t ] + dtime - 0.5 * ditime
    sel = awhere( ( time_array > time_min ) & ( time_array < time_max ) )
    t = t + max( 1, len( sel ) )
    if ( len( sel ) < min_sel ):
      if print_info:
        print '... skipping time range %s - %s' % (
            time_to_string( time_min ), time_to_string( time_max ) )
      continue
    sel_time_array = aget( time_array, sel )
    time_min = sel_time_array[ 0 ] - 0.5 * ditime
    time_max = sel_time_array[ -1 ] + 0.5 * ditime
    if print_info:
      print '... imaging time range %s - %s' % (
          time_to_string( time_min ), time_to_string( time_max ) )
  
    # image dirty facets
    temp_facet = get_facet( temp_facets, 1 )
    temp_beam = get_facet_beam( temp_facet )
    try:
      call_aips_task( 'IMAGR', indata = uv, nchav = channel_count, docalib = solution_switch,
          gainuse = solution_version, outdisk = temp_facet.disk, outname = temp_facet.name,
          outseq = temp_facet.seq, do3dimag = 1, nfield = s_facet_count, allokay = 0,
          outver = 0, cellsize = [ cell_size, cell_size ], boxfile = s_facet_file_name_e,
          imsize = [ s_facet_size, s_facet_size ], niter = 1000000, flux = 1000000.,
          in2disk = uv.disk, cmethod = conversion_method, minpatch = s_facet_size - 1,
          dotv = 0, overlap = 2, bcomp = [ 0 for j in range( 64 ) ], maxpixel = 0,
          imagrprm = [ 0,0,0, sdi_param, 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0 ], flagver = -1,
          timerang = timetime_to_dhmsdhms( [ time_min, time_max ] ),
          uvsize = [ uv_size, uv_size ], **imagr_params )
    except RuntimeError:
      continue
    
    i = i + 1
    s_facet = get_facet( s_facets, i )
#    s_beam = get_facet_beam( s_facet )
    temp_facet.rename( name = s_facet.name, klass = s_facet.klass, seq = s_facet.seq )
#    temp_beam.rename( name = s_beam.name, klass = s_beam.klass, seq = s_beam.seq )
    temp_beam.zap()
    
    # store parameters
    store_parameter( s_facet, 'time_center', ( time_max + time_min ) / 2. )
    store_parameter( s_facet, 'time_step', ( time_max - time_min ) * 1440. )
  
  store_parameter( s_facets, 'facet_count', i )
  
  return s_facets

###############################################################################

def image_clean_a_facets( uv, apply_solutions = False, sigma = 3., box_sigma = 5.,
    do_sdi_clean = False, conversion_method = 'GRID', facet_based_boxes = True,
    imagr_params = {}, sidelobe_rejection = 1.05, clean_box_radius = 5,
    print_info = True ):
  
  # make initial A-team facets (no cleaning)
  cell_size = restore_parameter( uv, 'cell_size' )
  uv_size = restore_parameter( uv, 'pb_image_size' )
  a_facet_count = restore_parameter( uv, 'a_facet_count' )
  a_facet_size = restore_parameter( uv, 'a_facet_size' )
  a_facet_file_name = restore_parameter( uv, 'a_facet_file_name' )
  a_facet_file_name_e = path.expandvars( a_facet_file_name )
  channel_count = get_channel_count( uv )
  a_facets = get_aips_file( uv.disk, 'A', 'ICL001', - 1, 'MA' )
  cpb_noise = restore_parameter( uv, 'cpb_noise' )
  flux = sigma * cpb_noise
  if apply_solutions:
    solution_switch = 100
    solution_version = 0 # uv.table_highver( 'SN' )
  else:
    solution_switch = -1
    solution_version = -1
  if do_sdi_clean:
    sdi_param = 0.1
  else:
    sdi_param = 0.
  time_list = get_time_list( uv )
  
  # image dirty facets
  # this call also adds empty model tables to the facet
  call_aips_task( 'IMAGR', indata = uv, nchav = channel_count, docalib = solution_switch,
      gainuse = solution_version, outdisk = a_facets.disk, outname = a_facets.name,
      outseq = a_facets.seq, outver = 0, cellsize = [ cell_size, cell_size ],
      imsize = [ a_facet_size, a_facet_size ], do3dimag = 1, niter = 1000000,
      boxfile = a_facet_file_name_e, in2disk = uv.disk, cmethod = conversion_method,
      minpatch = a_facet_size - 1, dotv = 0, overlap = 2, nfield = a_facet_count,
      imagrprm = [ 0,0,0, sdi_param, 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0 ], allokay = 0, 
      bcomp = [ 0 for j in range( 64 ) ], maxpixel = 0, flagver = -1, flux = 100000.,
      uvsize = [ uv_size, uv_size ], **imagr_params )
  
  # store parameters
  store_parameter( a_facets, 'facet_count', a_facet_count )
  store_parameter( a_facets, 'facet_file_name', a_facet_file_name )
  determine_facet_overlap( a_facets )
  
  # re-image A-team facets using UV data only while source above the horizon
  a_list = []
  new_facets = get_aips_file( a_facets.disk, 'NEW_A', a_facets.klass, - 1, 'MA' )
  for i in range( 1, a_facet_count + 1 ):
    a_facet = get_facet( a_facets, i )
    a_beam = get_facet_beam( a_facet )
    rts_times = calculate_rise_transit_set_times( uv, radec = get_radec( a_facet ) )
    [ uv_up, uv_down ] = split_uv_on_time_range( uv, rts_times, time_list = time_list )
    if ( not uv_up is None ):
      new_facet = get_facet( new_facets, 1 )
      new_beam = get_facet_beam( new_facet )
      new_facet_file_name = a_facet_file_name_e + '.%03d' % ( i )
      extract_facet_definitions( a_facet_file_name_e, [ i ], new_facet_file_name )
      call_aips_task( 'IMAGR', indata = uv_up, nchav = channel_count, outver = 0,
          docalib = solution_switch, gainuse = solution_version, allokay = 0,
          outdisk = new_facet.disk, outname = new_facet.name, outseq = new_facet.seq,
          cellsize = [ cell_size, cell_size ], imsize = [ a_facet_size, a_facet_size ],
          do3dimag = 1, niter = 1000000, flux = 100000., boxfile = new_facet_file_name,
          in2disk = uv_up.disk, cmethod = conversion_method, nfield = 1, flagver = -1,
          minpatch = a_facet_size - 1, dotv = 0, overlap = 2, maxpixel = 0,
          imagrprm = [ 0,0,0, sdi_param, 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0 ],
          bcomp = [ 0 for j in range( 64 ) ], uvsize = [ uv_size, uv_size ],
          **imagr_params )
      remove_file( new_facet_file_name )
      call_aips_task( 'TACOP', indata = a_facet, outdata = new_facet, inext = 'PS',
          ncount = 1 )
      a_facet.zap()
      a_beam.zap()
      new_facet.rename( name = a_facet.name, klass = a_facet.klass, seq = a_facet.seq )
      new_beam.rename( name = a_beam.name, klass = a_beam.klass, seq = a_beam.seq )
      uv_up.zap()
      peak_flux = get_image_maximum( a_facet )[ 0 ]
      a_list.append( [ i, peak_flux ] )
    else:
      pass # leave the facet as is
    if ( not uv_down is None ):
      uv_down.zap()
  a_list.sort( cmp = lambda a, b: cmp( b[ 1 ], a[ 1 ] ) )
  
  # automatically generate clean boxes
  add_clean_boxes( uv, a_facets, sidelobe_rejection = sidelobe_rejection,
      facet_based_boxes = facet_based_boxes, clean_box_radius = clean_box_radius, 
       print_info = print_info, box_sigma = box_sigma )
  
  # re-image and clean A-team facets using UV data only while above the horizon
  for a in a_list:
    [ i, peak_flux ] = a
    a_facet = get_facet( a_facets, i )
    a_beam = get_facet_beam( a_facet )
    rts_times = calculate_rise_transit_set_times( uv, radec = get_radec( a_facet ) )
    [ uv_up, uv_down ] = split_uv_on_time_range( uv, rts_times,
        time_list = time_list )
    if ( not uv_up is None ):
      new_facet = get_facet( a_facets, i )
      new_beam = get_facet_beam( new_facet )
      new_facet.rename( name = new_facets.name, klass = new_facets.klass,
          seq = new_facets.seq )
      new_beam.rename( name = new_facets.name,
          klass = get_facet_beam( new_facet ).klass,
          seq = new_facets.seq )
      new_facet_file_name = a_facet_file_name_e + '.%03d' % ( i )
      extract_facet_definitions( a_facet_file_name_e, [ i ], new_facet_file_name )
      new_facet.zap_table( 'CC', 0 )
      call_aips_task( 'IMAGR', indata = uv_up, nchav = channel_count, outver = 0,
          docalib = solution_switch, gainuse = solution_version, do3dimag = 1,
          outdisk = new_facet.disk, outname = new_facet.name, outseq = new_facet.seq,
          cellsize = [ cell_size, cell_size ], imsize = [ a_facet_size, a_facet_size ],
          niter = 1000000, flux = 0.95 * flux, in2disk = uv_up.disk, flagver = -1,
          boxfile = new_facet_file_name, cmethod = conversion_method, nfield = 1,
          minpatch = a_facet_size - 1, dotv = 0, overlap = 2, allokay = 1, maxpixel = 0,
          imagrprm = [ 0,0,0, sdi_param, 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0 ],
          bcomp = [ 0 for j in range( 64 ) ], uvsize = [ uv_size, uv_size ],
          **imagr_params )
      new_facet.rename( name = a_facet.name, klass = a_facet.klass, seq = a_facet.seq )
      new_beam.rename( name = a_beam.name, klass = a_beam.klass, seq = a_beam.seq )
      uv_up.zap()
      remove_file( new_facet_file_name )
    else:
      pass # leave the facet as is
    if ( not uv_down is None ):
      uv_down.zap()
  
  return a_facets

###############################################################################

def image_clean_pb_facets( uv, sigma = 2., improvement_limit = 0.05,
    apply_solutions = True, do_sdi_clean = False, restore_components = True,
    conversion_method = 'DFT', add_boxes = True, box_sigma = 5., keep_boxes = False,
    frequency_correction = False, imagr_params = {}, facet_based_boxes = True,
    sidelobe_rejection = 1.05, clean_box_radius = 5, print_info = True ):

  # image facets (but no cleaning yet)
  if frequency_correction:
    dish_diameter = restore_parameter( uv, 'dish_diameter' )
  else:
    dish_diameter = 0.
  cell_size = restore_parameter( uv, 'cell_size' )
  uv_size = restore_parameter( uv, 'pb_image_size' )
  pb_facet_count = restore_parameter( uv, 'pb_facet_count' )
  pb_facet_size = restore_parameter( uv, 'pb_facet_size' )
  pb_facet_file_name = restore_parameter( uv, 'pb_facet_file_name' )
  pb_facet_file_name_e = path.expandvars( pb_facet_file_name )
  channel_count = get_channel_count( uv )
  pb_facets = get_aips_file( uv.disk, 'PB', 'ICL001', -1, 'MA' )
  work_uv = get_aips_file( uv.disk, uv.name, 'IMAGR', -1, 'UV' )
  if apply_solutions:
    solution_switch = 100
    solution_version = uv.table_highver( 'SN' )
  else:
    solution_switch = -1
    solution_version = -1
  if do_sdi_clean:
    sdi_param = 0.1
  else:
    sdi_param = 0.

  # image dirty facets
  call_aips_task( 'IMAGR', indata = uv, nchav = channel_count, flagver = -1,
      docalib = solution_switch, gainuse = solution_version, in2data = work_uv,
      outdisk = pb_facets.disk, outname = pb_facets.name, outseq = pb_facets.seq,
      outver = 1, cellsize = [ cell_size, cell_size ], flux = 100000., dotv = 0,
      do3dimag = 1, niter = 1000000, imsize = [ pb_facet_size, pb_facet_size ],
      boxfile = pb_facet_file_name_e, cmethod = conversion_method, allokay = 0,
      minpatch = pb_facet_size - 1, overlap = 2, nfield = pb_facet_count,
      imagrprm = [ dish_diameter,0,0,sdi_param,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0 ],
      bcomp = [ 100000 for i in range( 64 ) ], maxpixel = 0, 
      uvsize = [ uv_size, uv_size ], **imagr_params )

  # store parameters to facets
  # determine facet overlap
  # determine initial noise
  store_parameter( pb_facets, 'facet_count', pb_facet_count )
  store_parameter( pb_facets, 'facet_file_name', pb_facet_file_name )
  determine_facet_overlap( pb_facets )
  cpb_noise = measure_cpb_noise( uv, pb_facets )

  # clean facets in steps to converge to noise level
  model_version = pb_facets.table_highver( 'CC' )
  model_versions = [ model_version, model_version + 1 ]
  clean_i = 0
  clean_rounds_min = 3
  noise_reduction = 1.
  boxes = keep_boxes
  while ( ( clean_i < clean_rounds_min ) or ( noise_reduction > improvement_limit ) ):
    clean_i = clean_i + 1
    last_cpb_noise = cpb_noise
    
    # combine previous model tables
    if ( clean_i > 1 ):
      if print_info:
        print '... combining previous clean model tables'
      for i in range( 1, 1 + pb_facet_count ):
        facet = get_facet( pb_facets, i )
        combine_model_tables( facet, model_versions = model_versions )
        for j in model_versions:
          facet.zap_table( 'CC', j )
        call_aips_task( 'TACOP', indata = facet, outdata = facet, inext = 'CC',
            invers = 0, outvers = model_version, ncount = 1 )
        facet.zap_table( 'CC', 0 )
    
    # automatically generate clean boxes
    if add_boxes:
      add_clean_boxes( uv, pb_facets, clean_box_radius = clean_box_radius,
          keep_boxes = boxes, facet_based_boxes = facet_based_boxes,
          sidelobe_rejection = sidelobe_rejection, print_info = print_info,
          box_sigma = box_sigma )
      boxes = True
    
    # start/continue clean on facets
    if print_info:
      print '... cleaning down to %s Jy' % ( repr( sigma * cpb_noise ) )
    call_aips_task( 'IMAGR', indata = uv, nchav = channel_count, flagver = -1,
        docalib = solution_switch, gainuse = solution_version, in2data = work_uv,
        outdisk = pb_facets.disk, outname = pb_facets.name, outseq = pb_facets.seq,
        outver = model_version + 1, cellsize = [ cell_size, cell_size ],
        flux = 0.95 * sigma * cpb_noise, dotv = 0, do3dimag = 1, niter = 1000000,
        imsize = [ pb_facet_size, pb_facet_size ], boxfile = pb_facet_file_name_e,
        cmethod = conversion_method, minpatch = pb_facet_size - 1, overlap = 2,
        nfield = pb_facet_count, allokay = 2, bcomp = [ 100000 for i in range( 64 ) ],
        imagrprm = [ dish_diameter,0,0,sdi_param,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0 ],
        maxpixel = 0, uvsize = [ uv_size, uv_size ], **imagr_params )
    store_parameter( pb_facets, 'clean_flux_min', sigma * cpb_noise )
    
    # determine noise reduction
    cpb_noise = measure_cpb_noise( uv, pb_facets )
    noise_reduction = 1. - ( cpb_noise / last_cpb_noise )
    if print_info:
      print '... relative noise reduction = %s percent' % ( 
          repr( 100. * noise_reduction ) )
  
  # restore model components to image
  # flux from the last imaging run was restored by IMAGR itself
  if restore_components:
    if print_info:
      print '... restoring cleaned flux in image'
    restore_version = model_version
  else:
    if print_info:
      print '... removing cleaned flux from image'
    for i in range( 1, pb_facet_count + 1 ):
      facet = get_facet( pb_facets, i )
      scale_model_flux( facet, -1., model_version + 1 )
    restore_version = 0
  rst_facets = restore_model_components( pb_facets, imagr_params = imagr_params,
      model_version = restore_version )
  for i in range( 1, 1 + pb_facet_count ):
    facet = get_facet( rst_facets, i )
    while table_exists( facet, 'CC', 0 ):
      facet.zap_table( 'CC', 0 )
  
  # combine final model tables
  if print_info:
    print '... combining final clean model tables'
  for i in range( 1, 1 + pb_facet_count ):
    pb_facet = get_facet( pb_facets, i )
    combine_model_tables( pb_facet, model_versions = model_versions )
    rst_facet = get_facet( rst_facets, i )
    call_aips_task( 'TACOP', indata = pb_facet, outdata = rst_facet, inext = 'CC',
        invers = 0, outvers = model_version, ncount = 1 )
    pb_facet.zap()
    rst_facet.rename( name = pb_facet.name, klass = pb_facet.klass, seq = pb_facet.seq )
    
  # cleanup
  work_uv.zap()
  store_parameter( pb_facets, 'clean_flux_min', sigma * cpb_noise )

  return pb_facets

###############################################################################

def image_clean_pbo_facets( uv, sigma = 2., improvement_limit = 0.05,
    apply_solutions = True, do_sdi_clean = False, restore_components = True,
    conversion_method = 'DFT', add_boxes = True, box_sigma = 5., keep_boxes = False,
    frequency_correction = False, imagr_params = {}, facet_based_boxes = True,
    sidelobe_rejection = 1.05, clean_box_radius = 5, print_info = True ):
  
  # make copy of UV data
  temp_uv = get_aips_file( uv.disk, uv.name, 'TEMP', -1, 'UV' )
  call_aips_task( 'MOVE', indata = uv, outdata = temp_uv, userid = get_aips_userid() )
  
  # combine primary beam and nearby outlier facet definitions
  pb_facet_file_name = restore_parameter( uv, 'pb_facet_file_name' )
  pb_facet_file_name_e = path.expandvars( pb_facet_file_name )
  o_facet_file_name = restore_parameter( uv, 'o_facet_file_name' )
  o_facet_file_name_e = path.expandvars( o_facet_file_name )
  pbo_facet_file_name = pb_facet_file_name_e + '.PBO'
  merge_facet_definitions( pb_facet_file_name_e, o_facet_file_name_e,
      pbo_facet_file_name )
  pb_facet_count = restore_parameter( uv, 'pb_facet_count' )
  o_facet_count = restore_parameter( uv, 'o_facet_count' )
  pbo_facet_count = pb_facet_count + o_facet_count
  
  # image primary beam and nearby outlier facets
  store_parameter( temp_uv, 'pb_facet_count', pbo_facet_count )
  store_parameter( temp_uv, 'pb_facet_file_name', pbo_facet_file_name )
  pbo_facets = image_clean_pb_facets( temp_uv, sigma = sigma,
      improvement_limit = improvement_limit, apply_solutions = apply_solutions,
      do_sdi_clean = do_sdi_clean, restore_components = restore_components,
      keep_boxes = keep_boxes, conversion_method = conversion_method,
      add_boxes = add_boxes, box_sigma = box_sigma, print_info = print_info,
      frequency_correction = frequency_correction, imagr_params = imagr_params,
      facet_based_boxes = facet_based_boxes, clean_box_radius = clean_box_radius,
      sidelobe_rejection = sidelobe_rejection )
  
  # split primary beam and nearby outlier facet definitions
  # automatically move o_sources to pb_facets that have overlap with pb_facets
  pb_facet_list = range( 1, 1 + pb_facet_count )
  o_facet_list = range( 1 + pb_facet_count, 1 + pb_facet_count + o_facet_count )
  o_to_pb_facet_list = []
  for i in o_facet_list:
    o_facet = get_facet( pbo_facets, i )
    if ( not model_table_empty( o_facet ) ):
      facet_overlap = get_facet_overlap( o_facet )
      for j in facet_overlap:
        if ( j in pb_facet_list ):
          o_to_pb_facet_list.append( i )
          break
  for i in o_to_pb_facet_list:
    o_facet_list.remove( i )
    pb_facet_list.append( i )
  remove_file( pb_facet_file_name_e )
  extract_facet_definitions( pbo_facet_file_name, pb_facet_list, pb_facet_file_name_e )
  remove_file( o_facet_file_name_e )
  if ( len( o_facet_list ) > 0 ):
    extract_facet_definitions( pbo_facet_file_name, o_facet_list, o_facet_file_name_e )
  remove_file( pbo_facet_file_name )
  
  # split primary beam and nearby outlier facets
  pb_facets = pbo_facets
  o_facets = get_aips_file( pbo_facets.disk, 'O', pbo_facets.klass, -1, 'MA' )
  pb_facet_count = 0
  o_facet_count = 0
  for i in range( 1, 1 + pbo_facet_count ):
    pbo_facet = get_facet( pbo_facets, i )
    pbo_beam = get_facet_beam( pbo_facet )
    if ( i in pb_facet_list ):
      pb_facet_count = pb_facet_count + 1
      pb_facet = get_facet( pb_facets, pb_facet_count )
      pb_beam = get_facet_beam( pb_facet )
      if ( aips_file_name_to_string( pb_facet ) != 
          aips_file_name_to_string( pbo_facet ) ):
        pbo_beam.rename( name = pb_beam.name, klass = pb_beam.klass,
            seq = pb_beam.seq )
        pbo_facet.rename( name = pb_facet.name, klass = pb_facet.klass,
            seq = pb_facet.seq )
    else: # ( i in o_facet_list )
      o_facet_count = o_facet_count + 1
      o_facet = get_facet( o_facets, o_facet_count )
      o_beam = get_facet_beam( o_facet )
      pbo_beam.rename( name = o_beam.name, klass = o_beam.klass,
          seq = o_beam.seq )
      pbo_facet.rename( name = o_facet.name, klass = o_facet.klass,
          seq = o_facet.seq )
  store_parameter( uv, 'pb_facet_count', pb_facet_count )
  store_parameter( pb_facets, 'facet_count', pb_facet_count )
  store_parameter( pb_facets, 'facet_file_name', pb_facet_file_name )
  store_parameter( uv, 'o_facet_count', o_facet_count )
  if ( o_facet_count > 0 ):
    store_parameter( o_facets, 'facet_count', o_facet_count )
    store_parameter( o_facets, 'facet_file_name', o_facet_file_name )
  else:
    o_facets = None

  # recalculate overlaps
  determine_facet_overlap( pb_facets )
  determine_facet_overlap( o_facets )
  
  # save noise level and delete copy of UV
  cpb_noise = restore_parameter( temp_uv, 'cpb_noise' )
  store_parameter( uv, 'cpb_noise', cpb_noise )
  temp_uv.zap()
  
  return [ pb_facets, o_facets ]

###############################################################################

def apply_pb_correction( uv, image, cutoff = 0.3, invert = False ):
  pbparms = [ cutoff, 1 ] + get_pbparms( uv )
  radec = get_radec( uv )
  coordinates = degdeg_to_hmsdms( radec )
  if ( radec[ 1 ] < 0. ):
    coordinates[ 3 ] = -abs( coordinates[ 3 ] )
    coordinates[ 4 ] = -abs( coordinates[ 4 ] )
    coordinates[ 5 ] = -abs( coordinates[ 5 ] )
  if invert:
    do_inverse = 1
  else:
    do_inverse = -1
  pbcor_image = get_aips_file( image.disk, image.name, 'PBCOR', -1, 'MA' )
  call_aips_task( 'PBCOR', indata = image, outdata = pbcor_image,
      pbparm = pbparms, coordina = coordinates, doinvers = do_inverse )
  return pbcor_image

###############################################################################

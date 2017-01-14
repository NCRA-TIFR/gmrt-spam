###############################################################################

# import Python modules
from sys import *
from os import *
from datetime import *
from math import *
import pdb
import scipy.interpolate as interpolate

# import user modules
from files import *
from aips import *
from acalc import *
from parameter import *
from mpfit import *
from solutions import *
from plot import *
from error import *
from ionosphere import *
from unwrap import *
from image import *

###############################################################################

def calibrate_instrumental_flux_and_bandpass( uv, flux, channel = 0, 
    reference_antenna = 0, iterations = 5, image_size = 256, imagr_params = {},
    calib_params = {} ):

# uv contains single calibrator, all auto-polarizations
# note that cross-polarizations (LR,RL,XY,YX) are not calibrated
# this function attaches an SN and BP table to uv,
# to be applied to both calibrator and target before channel collapse (SPLIT)

  # split off single channel
  if ( channel == 0 ):
    used_channel = int( get_channel_count( uv ) / 2 )
  else:
    used_channel = channel
  uv_1ch = get_aips_file( uv.disk, uv.name, '1CHAN', - 1, 'UV' )
  call_aips_task( 'SPLIT', indata = uv, bchan = used_channel, echan = used_channel,
      outdisk = uv_1ch.disk, outclass = uv_1ch.klass, outseq = uv_1ch.seq, douvcomp = 0 )

  # analyze data
  prep_uv = prepare_uv_data( uv_1ch )
  call_aips_task( 'TACOP', indata = prep_uv, outdata = uv_1ch, inext = 'PS', invers = 0, ncount = 1 )
  prep_uv.zap()

  # cover primary beam area with number of facets
  define_pb_facets( uv_1ch, max_smearing = 0.2 ) # primary beam

  # build source model
  source_list = [ [ get_radec( uv_1ch ), flux ] ]
  pbm_facets = make_pb_model_from_source_list( uv_1ch, source_list )

  # amplitude calibrate against source model on longest time scale
  # apply to full data
  cal_uv_1ch = get_aips_file( uv.disk, uv.name, '1CHCAL', - 1, 'UV' )
  dtime = ( restore_parameter( uv_1ch, 'time_max' ) - 
      restore_parameter( uv_1ch, 'time_min' ) )
  call_aips_task( 'CALIB', indata = uv_1ch, in2data = pbm_facets, invers = 0, ncomp = [ 0 ],
      flux = 0., nmaps = 1, cmodel = 'COMP', outdata = cal_uv_1ch, cmethod = 'DFT', snver = 0,
      refant = reference_antenna, solint = 3. * 24. * 60. * dtime, soltype = 'L1', solmode = 'A&P',
      aparm = [ 4, 0, 0, 0, 0, 0, 1, 0, 0 ], cparm = [ 0, 0, 0, 0, 0, 0 ], **calib_params )
  call_aips_task( 'TACOP', indata = uv_1ch, outdata = uv, inext = 'SN', invers = 0, ncount = 1 )
  cal_uv = get_aips_file( uv.disk, uv.name, 'CAL', - 1, 'UV' )
  call_aips_task( 'SPLIT', indata = uv, docalib = 100, gainuse = 0,
      outdisk = cal_uv.disk, outclass = cal_uv.klass, outseq = cal_uv.seq, douvcomp = 0 )

  # phase calibration against source model on shortest time scale
  calibrate_pb( cal_uv_1ch, pbm_facets, sigma = 0., conversion_method = 'DFT',
      reference_antenna = reference_antenna, calib_params = calib_params )

  # adjust facet definitions to include just the source
  facet_size = [ image_size, image_size ]
  facet_file_name = restore_parameter( cal_uv_1ch, 'pb_facet_file_name' )
  facet_file_name = path.expandvars( facet_file_name )
  temp_facet_file_name = facet_file_name + '.TEMP'
  facet_id = 1
  extract_facet_definitions( facet_file_name, [ facet_id ], temp_facet_file_name,
      include_clean_boxes = False )
  replace_facet( temp_facet_file_name, facet_id, facet_size, get_radec( cal_uv_1ch ),
      keep_boxes = False )
#  add_circular_clean_box( temp_facet_file_name, facet_id, [ 5, 5 ], 1 )
  add_rectangular_clean_box( temp_facet_file_name, facet_id, [ 0, 0 ], [ 0, 0 ] )
  move_file( temp_facet_file_name, facet_file_name )
  store_parameter( cal_uv_1ch, 'cpb_facet_count', 1 )
  store_parameter( cal_uv_1ch, 'cpb_image_size', facet_size[ 0 ] )
  store_parameter( cal_uv_1ch, 'pb_facet_count', 1 )
  store_parameter( cal_uv_1ch, 'pb_facet_size', facet_size[ 0 ] )

  # image and selfcal
  pb_facets = image_clean_pb_facets( cal_uv_1ch, sigma = 2., improvement_limit = 0.05,
      apply_solutions = True, conversion_method = 'DFT', do_sdi_clean = False,
      restore_components = True, add_boxes = True, box_sigma = 8., frequency_correction = False,
      imagr_params = imagr_params )
  selfcal_image_clean_pb( cal_uv_1ch, pb_facets, sigma = 2., try_final_amplitude = False,
      do_sdi_clean = False, reference_antenna = reference_antenna, improvement_limit = 0.05,
      restore_components = True, center_facets = False, add_boxes = True, box_sigma = 8.,
      frequency_correction = False, imagr_params = imagr_params, calib_params = calib_params )

  # perform bandpass calibration
  call_aips_task( 'TACOP', indata = cal_uv_1ch, outdata = cal_uv, inext = 'SN', invers = 0, ncount = 1 )
  call_aips_task( 'BPASS', indata = cal_uv, in2data = pb_facets, nmaps = 1, invers = 0, 
      cmethod = 'DFT', docalib = 100, gainuse = 0, solint = - 1, soltype = 'L1R',
      refant = reference_antenna, outvers = 0, ichansel = [ [ used_channel, used_channel, 1, 0 ] ],
      bpassprm = [ 0, 0, 0, 0, - 1, 0, 0, 0, 0, 0, 0 ] )
  call_aips_task( 'TACOP', indata = cal_uv, outdata = uv, inext = 'BP', invers = 0, ncount = 1 )

  # cleanup
  uv_1ch.zap()
  cal_uv_1ch.zap()
  cal_uv.zap()
  remove_facets( pbm_facets )
  remove_facets( pb_facets )
  remove_file( facet_file_name )
  remove_file( temp_facet_file_name )

  return

###############################################################################

def calibrate_instrumental_phases( full_uv, max_period = 1000., min_period = 5., max_gap = 5.,
    iterations = 10, rejection_limit = 2.0, height = 200.e3, equal_weights = True,
    reference_antenna = 0, print_info = False, image_size = 256, double = True, keep_scans = 5,
    polarization = '', imagr_params = {}, calib_params = {} ):
# uv contains single calibrator, all auto-polarizations
# note that cross-polarizations (LR,RL,XY,YX) are not calibrated
# max_gap in minutes

  # split off first polarization
  stokes_index = full_uv.header.ctype.index( 'STOKES' )
  stokes_list = [ full_uv.header.crval[ stokes_index ] + j * full_uv.header.cdelt[ stokes_index ]
      for j in range( full_uv.header.naxis[ stokes_index ] ) ]
  stokes = []
  if ( 1 in stokes_list ):
    stokes.append( 'I' )
  else:
    if ( - 1 in stokes_list ):
      stokes.append( 'RR' )
    if ( - 2 in stokes_list ):
      stokes.append( 'LL' )
    if ( - 5 in stokes_list ):
      stokes.append( 'XX' )
    if ( - 6 in stokes_list ):
      stokes.append( 'YY' )
  if ( len( stokes ) == 0 ):
    raise error( 'UV data does not contain suitable stokes/polarization' )
  elif ( polarization != '' ):
    if ( polarization in stokes ):
      # put selected polarization in front of list
      stokes.sort( cmp = lambda a, b: cmp( b == polarization, a == polarization ) )
    else:
      raise error( 'UV data does not contain stokes/polarization ' + polarization )
  for stoke in stokes:
    pol_uv = get_aips_file( full_uv.disk, full_uv.name, stoke, - 1, 'UV' )
    call_aips_task( 'SPLIT', indata = full_uv, outdisk = pol_uv.disk, outclass = pol_uv.klass,
        outseq = pol_uv.seq, stokes = stoke, douvcomp = 0 )
  uv = get_aips_file( pol_uv.disk, pol_uv.name, stokes[ 0 ], 0, 'UV' )

  # analyse UV data
  uv = prepare_uv_data( uv )

  # cover primary beam area with number of facets
  define_pb_facets( uv, max_smearing = 0.2 ) # primary beam

  # build source model
  source_radec = get_radec( uv )
  source_flux = 1.
  source_list = [ [ source_radec, source_flux ] ]
  pbm_facets = make_pb_model_from_source_list( uv, source_list )

  # first calibration against source model
  calibrate_pb( uv, pbm_facets, sigma = 0., conversion_method = 'DFT',
      reference_antenna = reference_antenna, calib_params = calib_params )

  # adjust facet definitions to include just the source
  facet_size = [ image_size, image_size ]
  facet_file_name = restore_parameter( uv, 'pb_facet_file_name' )
  facet_file_name = path.expandvars( facet_file_name )
  temp_facet_file_name = facet_file_name + '.TEMP'
  facet_id = 1
  extract_facet_definitions( facet_file_name, [ facet_id ], temp_facet_file_name,
      include_clean_boxes = False )
  replace_facet( temp_facet_file_name, facet_id, facet_size, source_radec, keep_boxes = False )
#  add_circular_clean_box( temp_facet_file_name, facet_id, [ 5, 5 ], 1 )
  add_rectangular_clean_box( temp_facet_file_name, facet_id, [ 0, 0 ], [ 0, 0 ] )
  move_file( temp_facet_file_name, facet_file_name )
  store_parameter( uv, 'cpb_facet_count', 1 )
  store_parameter( uv, 'cpb_image_size', facet_size[ 0 ] )
  store_parameter( uv, 'pb_facet_count', 1 )
  store_parameter( uv, 'pb_facet_size', facet_size[ 0 ] )

  # image and selfcal
  pb_facets = image_clean_pb_facets( uv, sigma = 2., improvement_limit = 0.05, apply_solutions = True,
      conversion_method = 'DFT', do_sdi_clean = False, restore_components = True, add_boxes = True,
      box_sigma = 5., frequency_correction = False, imagr_params = imagr_params )
  selfcal_image_clean_pb( uv, pb_facets, sigma = 2., try_final_amplitude = False, do_sdi_clean = False,
      reference_antenna = reference_antenna, improvement_limit = 0.05, restore_components = True,
      center_facets = False, add_boxes = True, box_sigma = 5., frequency_correction = False,
      imagr_params = imagr_params, calib_params = calib_params )

  # handle other polarization
  if ( len( stokes ) > 1 ):
    # calibrate complementary polarization product against same model
    # determine average phase difference
    pol_uv = get_aips_file( pol_uv.disk, pol_uv.name, stokes[ 1 ], 0, 'UV' )
    call_aips_task( 'TACOP', indata = uv, inext = 'PS', invers = 0, ncount = 1,
        outdata = pol_uv, outvers = 0 )
    call_aips_task( 'TACOP', indata = uv, inext = 'SN', invers = 0, ncount = 1,
        outdata = pol_uv, outvers = 0 )
    invert_solutions( pol_uv, in_version = 0, out_version = 0 )
    calibrate_pb( pol_uv, pbm_facets, sigma = 0., conversion_method = 'DFT',
        reference_antenna = reference_antenna, calib_params = calib_params )
    combine_solutions( pol_uv, in_version_1 = 0, in_version_2 = 0, out_version = 0 )
    time_window = ( restore_parameter( pol_uv, 'time_max' ) - restore_parameter( pol_uv, 'time_min' ) )
    time_window = time_window * 3. * 24. * 60. * 60. # window widths in seconds
    smooth_solutions_in_time( pol_uv, in_version = 0, out_version = 0, phase_window = time_window,
        amplitude_window = 0. )
  
  # filter instrumental phases
  call_aips_task( 'TACOP', indata = uv, inext = 'SN', invers = 0, ncount = 1,
      outdata = pb_facets, outvers = 0 )
  filter_instrumental_phases( uv, pb_facets, facet_list = [ 1 ], in_version = 0, 
    max_period = max_period, min_period = min_period, max_gap = max_gap,
    equal_weights = equal_weights, double_precision = double, iterations = iterations,
    height = height, out_version = 0, print_info = print_info, keep_scans = keep_scans )
  
  # handle other polarization
  if ( len( stokes ) > 1 ):

    # determine instrumental phases for other polarization
    pol_uv = get_aips_file( pol_uv.disk, pol_uv.name, stokes[ 1 ], 0, 'UV' )
    call_aips_task( 'TACOP', indata = uv, inext = 'SN', invers = 0, ncount = 1,
        outdata = pol_uv, outvers = 0 )
    combine_solutions( pol_uv, in_version_1 = 0, in_version_2 = 0, out_version = 0 )

    # copy polarization calibrations to target uv file
    solution_version = full_uv.table_highver( 'SN' )
    call_aips_task( 'SNDUP', indata = uv, invers = 0, outdata = full_uv,
        outvers = solution_version + 1 )
    call_aips_task( 'SNCOR', indata = full_uv, snver = solution_version + 1,
        stokes = stokes[ 1 ][ 0 ], opcode = 'CLPW', sncorprm = [ - 2., - 1. ] )
    call_aips_task( 'SNDUP', indata = pol_uv, invers = 0, outdata = full_uv,
        outvers = solution_version + 2 )
    call_aips_task( 'SNCOR', indata = full_uv, snver = solution_version + 2,
        stokes = stokes[ 0 ][ 0 ], opcode = 'CLPW', sncorprm = [ - 2., - 1. ] )
    time_window = ( restore_parameter( uv, 'time_max' ) - restore_parameter( uv, 'time_min' ) )
    time_window = time_window * 3 * 24. # window widths in hours
    call_aips_task( 'CLCAL', indata = full_uv, opcode = 'MERG', snver = solution_version + 1,
        invers = solution_version + 2, gainver = solution_version + 3 )
  else:
    call_aips_task( 'TACOP', indata = uv, inext = 'SN', invers = 0, ncount = 1,
        outdata = full_uv, outvers = 0 )

  # cleanup
  uv.zap()
  for stoke in stokes:
    pol_uv = get_aips_file( full_uv.disk, full_uv.name, stoke, 0, 'UV' )
    pol_uv.zap()
  remove_facets( pbm_facets )
  remove_facets( pb_facets )
  remove_file( facet_file_name )
  remove_file( temp_facet_file_name )

  return

###############################################################################

def fit_phi_poly_offset_model( P, dojac = None, antennas = None,
    X_table = None, za_table = None, phase_table = None,
    error_table = None, info_table = None, normalize_weights = None ):
  
  # handle input parameters
  if ( ( not dojac is None ) or ( antennas is None ) or
      ( X_table is None ) or ( za_table is None ) or 
      ( phase_table is None ) or ( error_table is None ) or
      ( info_table is None ) or ( normalize_weights is None ) ):
    return -1, None, None
  order = 5
  offsets = P[ 0 : antennas ]
  coefs = P[ antennas : ]
  times = len( coefs ) / order
  if ( len( coefs ) != order * times ):
    return -1, None, None
  coefs = coefs.reshape( times, order )
  
  # split pierce points per time
  time_list = []
  for time in info_table[ : , 0 ]:
    if ( not ( time in time_list ) ):
      time_list.append( time )
  if ( len( time_list ) != times ):
    return -1, None, None
  
  # fill model phase table with offsets
  phi_table = azeros( phase_table )
  for i in range( antennas ):
    sel = awhere( info_table[ : , 1 ] == i )
    if ( len( sel ) > 0 ):
      phi_table = aput( phi_table, sel, aget( phi_table, sel ) + offsets[ i ] )
  
  # calculate polynomial model phases per time stamp
  for n in range( times ):
    sel = awhere( info_table[ : , 0 ] == time_list[ n ] )
    X = aget( X_table, sel )
    za = aget( za_table, sel )
    phi_poly = phi_poly_model( X, coefs[ n ] ) / cos( aradians( za ) )
    phi_table = aput( phi_table, sel, aget( phi_table, sel ) + phi_poly )
  
  # calculate chi2 terms
  chi_list = []
  if normalize_weights:
    normalize_list = []
  for n in range( times ):
    sel = awhere( info_table[ : , 0 ] == time_list[ n ] )
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

def filter_instrumental_phases( uv, in_version = 0, max_period = 5.,
    min_period = 5., max_gap = 60., equal_weights = True, write_weights = False,
    height = 300.e3, out_version = 0, iterations = 1,  print_info = True,
    min_points = 5, time_range = [], drop_antennas = False, time_step = 1,
    density_scale = 1.e3, density_power = 1., max_range = 1.e6,
    phase_rms_criterion = None ):
# max_gap, max_period, min_period & max_range in minutes
# density scale in meters

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
  reference_array = array( calibration_data[ 4 ], dtype = int32 )
  ph_array = array( calibration_data[ 5 ], dtype = dtype )
  er_array = array( calibration_data[ 6 ], dtype = dtype )
  if ( len( time_range ) == 2 ):
    sel = awhere( ( time_array[ : , 0 ] >= time_range[ 0 ] ) &
        ( time_array[ : , 0 ] <= time_range[ 1 ] ) )
    time_array = aget( time_array, sel )
    reference_array = aget( reference_array, sel )
    ph_array = aget( ph_array, sel )
    er_array = aget( er_array, sel )
  time_count = len( time_array )
  
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
  
  # remove short scans
  short_scan_list = []
  for scan in scan_list:
    if ( scan[ 1 ] - scan[ 0 ] + 1 < min_points ):
      short_scan_list.append( scan )
  for scan in short_scan_list:
    scan_list.remove( scan )
  
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
  
  # remove entries within short, isolated scans
  time_list = []
  for scan in scan_list:
    time_list = time_list + range( scan[ 0 ], scan[ 1 ] + 1 )
  
  # convert to arrays
  antenna_count = len( antenna_table )
  scan_sel = array( scan_sel, dtype = int32 ).reshape( len( scan_sel ), 1 )
  
  # initial estimate instrumental phases and phase ambiguities
  uc_ph_array = azeros( ph_array )
  phi_instr_table = zeros( ( scan_count, antenna_count ), dtype = dtype )
  phi_ambig_table = zeros( ( time_count, antenna_count ), dtype = dtype )
  for i in range( antenna_count ):
    for s in range( scan_count ):
      scan = scan_list[ s ]
      sel = awhere( er_array[ scan[ 0 ] : scan[ 1 ] + 1, i ] > 0. )
      if ( len( sel ) > 0 ):
        p1 = aget( ph_array[ scan[ 0 ] : scan[ 1 ] + 1, i ], sel )
#        p2 = adegrees( unwrap_lms( aradians( p1 ), alpha = 1.e-4 ) )
        if ( len( p1 ) < 20 ):
          p2 = adegrees( unwrap( aradians( p1 ) ) )
        else:
          p2 = aunwrap_phase( p1, alpha = 0.01 )
          if sometrue( isnan( p2 ) ):
            p2 = aunwrap_phase( p1, alpha = 0.001 )
          if ( sometrue( isnan( p2 ) ) or ( max( fabs( p2 ) ) > 1.e4 ) ):
            p2 = adegrees( unwrap( aradians( p1 ) ) )
#        tck_phase = interpolate.splrep( aget( 
#            time_array[ scan[ 0 ] : scan[ 1 ] + 1, 0 ], sel ), p2 )
#        p2_interpolated = interpolate.splev( 
#            time_array[ scan[ 0 ] : scan[ 1 ] + 1, 0 ], tck_phase )
##        scan_center = int( ( scan[ 1 ] - scan[ 0 ] ) / 2 )
#        scan_center = 0 # scan[ 0 ]
#        p2 = p2 - p2_interpolated[ scan_center ]
        p2 = p2 - p2[ 0 ]
        uc_ph_array[ scan[ 0 ] : scan[ 1 ] + 1, i ] = aput(
            uc_ph_array[ scan[ 0 ] : scan[ 1 ] + 1, i ], sel, p2 )
        phi_instr_table[ s, i ] = amodulo( amean_phase( p1 - p2 ) + 180., 360. ) - 180.
        phi_ambig_table[ scan[ 0 ] : scan[ 1 ] + 1, i ] = aput( 
            phi_ambig_table[ scan[ 0 ] : scan[ 1 ] + 1, i ], sel, 
            around( ( p2 - p1 + phi_instr_table[ s, i ] ) / 360. ) * 360. )
  
  # solve per scan
  scan_rms_table = []
  scan_offset_table = []
  scan_time_table = []
  scan_coefs_table = []
  scan_phi_poly_table = []
  if print_info:
    print 'found %s scans' % ( repr( scan_count ) )
  for s in range( scan_count ):
    scan = scan_list[ s ]
    if print_info:
      print '... solving for scan %s' % ( repr( s + 1 ) )
    
    # containers for model input data
    X_table = []
    za_table = []
    ref_table = []
    phase_table = []
    error_table = []
    info_table = []
    
    # collect necessary data for fit
    for n in range( scan[ 0 ], scan[ 1 ] + 1, time_step ):
      
      # containers for model input data
      sX_table = []
      sza_table = []
      sref_table = []
      sphase_table = []
      serror_table = []
      sinfo_table = []
      
      # get pierce point coordinates
      pierce_table = calculate_pierce_coordinates( time_array[ n ],
          center_table, [ center_table ], array_table, antenna_table,
          height = height, iterations = iterations )
      
      # loop over pierce points
      ref_list = []
      j = 0
      for pierce_info in pierce_table:
        [ X, za, [ k, i ] ] = pierce_info
        if ( er_array[ n, i ] > 0. ):
          
          # store model fit input data
          sX_table.append( X )
          sza_table.append( za )
          sphase_table.append( uc_ph_array[ n, i ] )
          if equal_weights:
#            serror_table.append( 1. )
            serror_table.append( 1. / weight_table[ i ] )
          else:
#            serror_table.append( er_array[ n, i ] )
            serror_table.append( er_array[ n, i ] / weight_table[ i ] )
          if ( i == reference_array[ n ] ):
            ref_list.append( [ i, j ] )
          j = j + 1
          sinfo_table.append( [ n, i, reference_array[ n ] ] )
      
      # loop over pierce points
      for pierce_info in pierce_table:
        [ X, za, [ k, i ] ] = pierce_info
        if ( er_array[ n, i ] > 0. ):
          
          # store index to reference antenna
          sref_table.append( ref_list[ 0 ][ 1 ] )
      
      # if sufficient data available, add to big fit data
      if ( len( sphase_table ) <= 5 ):
        continue
      ref_offset = len( ref_table )
      for i in range( len( sphase_table ) ):
        X_table.append( [ x for x in sX_table[ i ] ] )
        za_table.append( sza_table[ i ] )
        ref_table.append( sref_table[ i ] + ref_offset )
        phase_table.append( sphase_table[ i ] )
        error_table.append( serror_table[ i ] )
        info_table.append( [ x for x in sinfo_table[ i ] ] )
    
    # check for empty result
    if ( len( info_table ) == 0 ):
      scan_rms_table.append( 180. )
      scan_offset_table.append( [] )
      scan_time_table.append( [] )
      scan_coefs_table.append( [] )
      scan_phi_poly_table.append( [] )
      continue
    
    # convert tables to arrays
    X_table = array( X_table, dtype = dtype )
    za_table = array( za_table, dtype = dtype )
    ref_table = array( ref_table, dtype = int32 )
    phase_table = array( phase_table, dtype = dtype )
    error_table = array( error_table, dtype = dtype )
    info_table = array( info_table, dtype = int32 )
    
    # check for varying reference antenna
    reference_antenna = info_table[ 0, 2 ]
    sel = awhere( info_table[ : , 2 ] != reference_antenna )
    if ( len( sel ) > 0 ):
      raise error( 'reference antenna change detected' )

    # check for dead antennas
    antenna_active = []
    for i in range( antenna_count ):
      sel = awhere( info_table[ : , 1 ] == i )
      if ( len( sel ) > 0 ):
        antenna_active.append( True )
      else:
        antenna_active.append( False )
    
    # split pierce points per time
    stime_list = []
    for time in info_table[ : , 0 ]:
      if ( not ( time in stime_list ) ):
        stime_list.append( time )
    stime_count = len( stime_list )
    stime_array = array( stime_list, dtype = int32 )
    
    # fit gradient + offset model to ionospheric phases
    # tie offset of reference antenna to zero
    # tie gradient of first time fit to zero
    function_keywords = { 'antennas' : antenna_count,
        'X_table' : X_table, 'za_table' : za_table,
        'phase_table' : phase_table, 'error_table' : error_table,
        'info_table' : info_table, 'normalize_weights' : normalize_weights }
    parameter_info = []
    for i in range( antenna_count ):
      if ( ( not antenna_active[ i ] ) or ( i == reference_antenna ) ):
        par_info = { 'parname' : 'P_%d' % ( i + 1 ), 'value' : 0.,
            'limits' : [ 0., 0. ] }
      else:
        par_info = { 'parname' : 'P_%d' % ( i + 1 ), 'value' : 0.,
            'limits' : [ None, None ] }
      parameter_info.append( par_info )
    for n in range( stime_count ):
      for o in range( 5 ):
        if ( ( o >= 2 ) or ( ( n == 0 ) and ( o < 2 ) ) ):
          par_info = { 'parname' : 'P_%d_%d' % ( n + 1, o + 1 ),
              'value' : 0., 'limits' : [ 0., 0. ] }
        else:
          par_info = { 'parname' : 'P_%d_%d' % ( n + 1, o + 1 ),
              'value' : 0., 'limits' : [ None, None ] }
        parameter_info.append( par_info )

    try:
      fit = mpfit( fit_phi_poly_offset_model, functkw = function_keywords,
          parinfo = parameter_info, quiet = True, autoderivative = True,
          debug = debug, fastnorm = True, nocovar = True, dblprec = double_precision )
    except:
      rms = 360.
      offsets = zeros( ( antenna_count ), dtype = dtype )
      coefs = zeros( ( stime_count, 5 ), dtype = dtype )
      if print_info:
        print '... gradient fit failed'
    else:
      rms = sqrt( fit.fnorm )
      offsets = fit.params[ 0 : antenna_count ]
      coefs = fit.params[ antenna_count : ].reshape( stime_count, 5 )
      if print_info:
        print '... gradient post-fit phase RMS = %s degrees' % repr( rms )
    
    if ( rms < 360. ):
      parameter_info = []
      for i in range( antenna_count ):
        if ( ( not antenna_active[ i ] ) or ( i == reference_antenna ) ):
          par_info = { 'parname' : 'P_%d' % ( i + 1 ), 'value' : 0.,
              'limits' : [ 0., 0. ] }
        else:
          par_info = { 'parname' : 'P_%d' % ( i + 1 ), 'value' : offsets[ i ],
              'limits' : [ None, None ] }
        parameter_info.append( par_info )
      for n in range( stime_count ):
        for o in range( 5 ):
#          if ( o < 2 ):
          if ( ( o < 2 ) and ( n > 0 ) ):
            coef = coefs[ n, o ]
            par_info = { 'parname' : 'P_%d_%d' % ( n + 1, o + 1 ),
#                'value' : coef, 'limits' : [ coef, coef ] }
                'value' : coef, 'limits' : [ None, None ] }
          elif ( n == 0 ):
            par_info = { 'parname' : 'P_%d_%d' % ( n + 1, o + 1 ),
                'value' : 0., 'limits' : [ 0., 0. ] }
          else:
            par_info = { 'parname' : 'P_%d_%d' % ( n + 1, o + 1 ),
                'value' : 0., 'limits' : [ None, None ] }
          parameter_info.append( par_info )
      try:
        fit = mpfit( fit_phi_poly_offset_model, functkw = function_keywords,
          parinfo = parameter_info, quiet = True, autoderivative = True,
          debug = debug, fastnorm = True, nocovar = True, dblprec = double_precision )
      except:
        rms = 360.
        if print_info:
          print '... polynomial fit failed'
      else:
        rms = sqrt( fit.fnorm )
        offsets = fit.params[ 0 : antenna_count ]
        coefs = fit.params[ antenna_count : ].reshape( stime_count, 5 )
        if print_info:
          print '... polynomial post-fit phase RMS = %s degrees' % repr( rms )
    
    # save fit results
    scan_rms_table.append( rms )
    scan_offset_table.append( offsets )
    scan_time_table.append( stime_array )
    scan_coefs_table.append( coefs )
    
    # calculate polynomial model phases per time stamp
    phi_table = azeros( phase_table )
#    for i in range( antenna_count ):
#      sel = awhere( info_table[ : , 1 ] == i )
#      phi_table = aput( phi_table, sel, aget( phi_table, sel ) + offsets[ i ] )
    for n in range( stime_count ):
      sel = awhere( info_table[ : , 0 ] == stime_list[ n ] )
      X = aget( X_table, sel )
      za = aget( za_table, sel )
      phi_poly = phi_poly_model( X, coefs[ n ] ) / cos( aradians( za ) )
      phi_table = aput( phi_table, sel, aget( phi_table, sel ) + phi_poly )
    sel = ref_table.reshape( len( ref_table ), 1 )
    phi_table = phi_table - aget( phi_table, sel )
    phi_poly_array = zeros( ( stime_count, antenna_count ), dtype = dtype )
    for x in range( len( info_table ) ):
      [ n, i, r ] = info_table[ x ]
      phi_poly_array[ stime_list.index( n ), i ] = phi_table[ x ]
    scan_phi_poly_table.append( phi_poly_array )
    
    if ( not phase_rms_criterion is None ):
      if ( rms < phase_rms_criterion ):
        if print_info:
          print '... phase RMS criterion achieved!'
        break

  s = scan_rms_table.index( min( scan_rms_table ) )
  scan = scan_list[ s ]
  rms = scan_rms_table[ s ]
  offset_table = scan_offset_table[ s ]
  time_table = scan_time_table[ s ]
#  coefs = scan_coefs_table[ s ]
  phi_poly_table = scan_phi_poly_table[ s ]
  
  if print_info:
    print 'selecting scan %s with phase RMS = %s degrees' % ( repr( s + 1 ),
        repr( rms ) )
  
  # update instrumental phases and phase ambiguities
  for i in range( antenna_count ):
    if ( offset_table[ i ] == 0. ):
      phi_instr_table[ s, i ] = 0.
    else:
      phi_instr_table[ s, i ] = phi_instr_table[ s, i ] + offset_table[ i ]
  sel = awhere( phi_instr_table[ s ] >= 180. )
  if ( len( sel ) > 0 ):
    for n in range( scan[ 0 ], scan[ 1 ] + 1 ):
      phi_ambig_table[ n ] = aput( phi_ambig_table[ n ], sel, 
          aget( phi_ambig_table[ n ], sel ) - 360. )
  sel = awhere( phi_instr_table[ s ] < -180. )
  if ( len( sel ) > 0 ):
    for n in range( scan[ 0 ], scan[ 1 ] + 1 ):
      phi_ambig_table[ n ] = aput( phi_ambig_table[ n ], sel,
          aget( phi_ambig_table[ n ], sel ) + 360. )
  phi_instr_table[ s ] = amodulo( phi_instr_table[ s ] + 180., 360. ) - 180.
  
  # calculate antenna-based RMS
  sel = time_table.reshape( len( time_table ), 1 )
  sub_time_array = aget( time_array, sel )
  sub_ph_array = aget( ph_array, sel )
  sub_er_array = aget( er_array, sel )
  sub_phi_ambig_table = aget( phi_ambig_table, sel )
  ant_res_table = zeros( ( len( sel ), antenna_count ), dtype = dtype ) - 360.
  for i in range( antenna_count ):
    sel = awhere( sub_er_array[ : , i ] > 0. )
    if ( len( sel ) > 0 ):
      p1 = ( aget( sub_ph_array[ : , i ], sel ) - phi_instr_table[ s, i ] +
          aget( sub_phi_ambig_table[ : , i ], sel ) )
      p2 = aget( phi_poly_table[ : , i ], sel )
      ant_res_table[ : , i ] = aput( ant_res_table[ : , i ], sel, p1 - p2 )
  ant_rms_list = []
  for i in range( antenna_count ):
    offset_list = []
    for n in range( len( time_table ) ):
#      sel = awhere( ant_res_table[ n ] > 0. )
      sel = awhere( sub_er_array[ n ] > 0. )
      if ( ( len( sel ) > 0 ) and ( ant_res_table[ n, i ] >= 0. ) ):
        offset_list.append( 
            ( aget( ant_res_table[ n ], sel ) - ant_res_table[ n, i ] ).mean() )
    if ( len( offset_list ) >= min_points ):
      rms = sqrt( ( array( offset_list, dtype = dtype )**2 ).mean() )
      ant_rms_list.append( rms )
      if print_info:
        print 'antenna %s phase RMS = %s degrees' % ( repr( i + 1 ), repr( rms ) )
    else:
      ant_rms_list.append( 0. )
      if print_info:
        print 'antenna %s phase RMS = undefined, too few data points' % (
            repr( i + 1 ) )
  
  # build solution table containing instrumental phase corrections
  solution_table = []
  time_array = array( calibration_data[ 0 ], dtype = dtype )
  reference_array = array( calibration_data[ 4 ], dtype = int32 )
  range_time = mean( sub_time_array[ : , 0 ] )
  for t in range( len( time_array ) ):
    solution_row = [ [ time_array[ t, 0 ], reference_array[ t ] + 1, 0., 0. ] ]
    for i in range( antenna_count ):
      if ( ant_rms_list[ i ] == 0. ):
        if drop_antennas:
          solution_row.append( [ 0., 0., 0., 0. ] )
        else:
          solution_row.append( [ 1., 0., 0., 1. ] )
      else:
        if ( abs( time_array[ t, 0 ] - range_time ) < max_range / ( 24. * 60 ) ):
          gain = r_phi_to_complex( [ 1., phi_instr_table[ s, i ] ] )
          if write_weights:
            weight = 1. / radians( ant_rms_list[ i ] )
          else:
            weight = 1.
          solution_row.append( [ gain.real, gain.imag, 0., weight ] )
        else:
          solution_row.append( [ 1., 0., 0., 1. ] )
      sel = awhere( array( solution_row )[ 1 : , 3 ] > 0. )
    if ( len( sel ) == 0 ):
      solution_row[ 0 ][ 1 ] = 0.
    solution_table.append( solution_row )
  write_solution_table( uv, solution_table, out_version = out_version )
  
  return

###############################################################################

###############################################################################

# import Python modules
from os import *
from math import *
from string import *
import datetime
import sys
import traceback
from numpy import *
import glob

# import user modules
from spam import *

###############################################################################

def process_target( uvfits_file_name, keep_log = True, minimize_storage = True,
    calib_params = { 'uvrang' : [ 1.,1.e6 ] }, imagr_params = { 'robust' : -1. },
    do_large_boxes = False, do_plots = True, do_cleanup = False, 
    model_source_list = [], model_resolution = 25., snapshot = False, 
    box_radius = 5, large_box_radius = 20, try_new = True ):
  
  # get full paths to data directories
  dat_path = path.expandvars( '${DAT}/' )
  fits_path = path.expandvars( '${FIT}/' )
  
  # check presence of input file(s)
  uvfits_file_names = glob.glob( uvfits_file_name )
  if ( len( uvfits_file_names ) == 0 ):
    raise error( 'UVFITS file(s) %s does not exist' % ( uvfits_file_name ) )
  pointing_name = uvfits_file_names[ 0 ].split( '/' )[ -1 ]
  pointing_name = pointing_name.split( '.' )[ 0 ]
  target_name = pointing_name.split( '_' )[ 0 ]
  
  # find free AIPS disk
  aips_disk = 0
  for i in range( 1, len( AIPS.disks ) ):
    if ( len( AIPSCat( i )[ i ] ) == 0 ):
      aips_disk = i
      break
  if ( aips_disk == 0 ):
    raise error( 'no free AIPS disk available' )
  
  # set up logging and benchmarking
  dt = datetime.datetime.now()
  dt_string = '%4d%02d%02d_%02d%02d%02d' % ( 
      dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second )
  print '%s: started processing of field %s' % ( dt_string, target_name )
  if keep_log:
    log_file_name = dat_path + 'spam_%s_%s.log' % ( pointing_name, dt_string )
    log_file = file( log_file_name, 'w' )
    old_stdout = sys.stdout
    old_stderr = sys.stderr
    sys.stdout = log_file
    sys.stderr = log_file
    print '%s: started processing of field %s' % ( dt_string, target_name )
    sys.stdout.flush()
  
  # process pointing
  try:
    t_iical_uv = get_aips_file( aips_disk, target_name, 'IICAL', -1, 'UV' )
    read_fits_uv( uvfits_file_names[ 0 ], t_iical_uv )
    if ( len(  uvfits_file_names ) > 1 ):
      change_uv_reference_date( t_iical_uv, t_iical_uv, offset = 0 )
      for i in range( 1, len( uvfits_file_names ) ):
        t_iical_uv1 = get_aips_file( aips_disk, target_name, 'IICAL', 0, 'UV' )
        t_iical_uv2 = get_aips_file( aips_disk, target_name, 'IICAL', -1, 'UV' )
        read_fits_uv( uvfits_file_names[ i ], t_iical_uv2 )
        change_uv_reference_date( t_iical_uv2, t_iical_uv1, offset = 2 * i )
        t_iical_uv = get_aips_file( aips_disk, target_name, 'IICAL', -1, 'UV' )
        call_aips_task( 'DBCON', indata = t_iical_uv1, in2data = t_iical_uv2,
            outdata = t_iical_uv, doarray = 1 )
        t_iical_uv1.zap()
        t_iical_uv2.zap()
    spam_target( t_iical_uv, try_new = try_new, minimize_storage = minimize_storage,
        calib_params = calib_params, imagr_params = imagr_params, do_plots = do_plots,
        do_large_boxes = do_large_boxes, model_source_list = model_source_list,
        snapshot = snapshot, model_resolution = model_resolution, 
        box_radius = box_radius, large_box_radius = large_box_radius )
  except KeyboardInterrupt:
    dt = datetime.datetime.now()
    dt_string = '%4d%02d%02d_%02d%02d%02d' % ( 
        dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second )
    if keep_log:
      print '%s: processing of field %s manually interrupted' % ( dt_string, target_name )
      sys.stdout.flush()
      sys.stdout = old_stdout
      sys.stderr = old_stderr
      log_file.close()
    print '%s: processing of field %s manually interrupted' % ( dt_string, target_name )
  except:
    print traceback.print_exc()
    dt = datetime.datetime.now()
    dt_string = '%4d%02d%02d_%02d%02d%02d' % ( 
        dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second )
    if keep_log:
      print '%s: processing of field %s failed' % ( dt_string, target_name )
      sys.stdout.flush()
      sys.stdout = old_stdout
      sys.stderr = old_stderr
      log_file.close()
    print '%s: processing of field %s failed' % ( dt_string, target_name )
  else:
    dt = datetime.datetime.now()
    dt_string = '%4d%02d%02d_%02d%02d%02d' % ( 
        dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second )
    if keep_log:
      print '%s: processing of field %s ended normally' % ( dt_string, target_name )
      sys.stdout.flush()
      sys.stdout = old_stdout
      sys.stderr = old_stderr
      log_file.close()
    print '%s: processing of field %s ended normally' % ( dt_string, target_name )
  
  if do_cleanup:
    clean_aips_disk( aips_disk )
  return

###############################################################################

def spam_target( t_iical_uv, minimize_storage = True, model_source_list = [],
    calib_params = { 'uvrang' : [ 1.,1.e6 ] }, imagr_params = { 'robust' : -1. },
    do_large_boxes = False, do_plots = True, snapshot = False, use_nvss = True,
    try_new = True, model_resolution = 25., box_radius = 5, large_box_radius = 20 ):
  
  # define target field processing parameters
  aips_disk = t_iical_uv.disk
  target_name = t_iical_uv.name
  uv = prepare_uv_data( t_iical_uv )
  if minimize_storage:
    t_iical_uv.zap()
  
  # get full paths to data directories
  dat_path = path.expandvars( '${DAT}/' )
  prt_path = path.expandvars( '${PRT}/' )
  fits_path = path.expandvars( '${FIT}/' )
  
  # determine instrument/frequency band
  if ( uv.header.telescop == 'GMRT' ):
    reference_antenna_list = [ 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,25 ]
    frequency = get_central_frequency( uv )
    if ( frequency > 140.e6 ) and ( frequency < 170.e6 ):
      max_smearing = 0.2
      outlier_fluxes = [ 3., 1.5 ]
      outlier_ratios = [ 3., 1.5 ]
      model_flux_min = 0.2
      outlier_flux_min = 0.2
      clip_levels = [ 50., 25. ]
      clip_uv_min = 0.5
      rfi_level = 10. # Jy
      amplitude_snr = 500
      use_catalog = [ True, True ]
      catalog_resolution = 50.
      plot_range = [ -0.005, 0.05 ]
    elif ( frequency > 200.e6 ) and ( frequency < 280.e6 ):
      max_smearing = 0.15
      outlier_fluxes = [ 1., 0.5 ]
      outlier_ratios = [ 3., 1.5 ]
      model_flux_min = 0.1
      outlier_flux_min = 0.1
      clip_levels = [ 10., 5. ]
      clip_uv_min = 1.
      rfi_level = 3. # Jy
      amplitude_snr = 700.
      use_catalog = [ True, True ]
      catalog_resolution = 35.
      plot_range = [ -0.0025, 0.025 ]
    elif ( frequency > 300.e6 ) and ( frequency < 350.e6 ):
      max_smearing = 0.15
      outlier_fluxes = [ 0.7, 0.2 ]
      outlier_ratios = [ 3., 1.5 ]
      model_flux_min = 0.05
      outlier_flux_min = 0.05
      clip_levels = [ 7., 3. ]
      clip_uv_min = 1.
      rfi_level = 0.3 # Jy
      amplitude_snr = 1200.
      use_catalog = [ True, False ]
      catalog_resolution = 25.
      plot_range = [ -0.00005, 0.0005 ]
    elif ( frequency > 580.e6 ) and ( frequency < 650.e6 ):
      max_smearing = 0.1
      outlier_fluxes = [ 0.3, 0.1 ] # Jy
      outlier_ratios = [ 5., 2. ]
      model_flux_min = 0.005 # Jy
      outlier_flux_min = 0.025 # Jy
      clip_levels = [ 1.5, 0.8 ] # Jy
      clip_uv_min = 2. # klambda
      rfi_level = 0.1 # Jy
      amplitude_snr = 1000.
      use_catalog = [ False, False ]
      catalog_resolution = 15.
      plot_range = [ -0.00025, 0.0025 ]
    elif ( frequency > 900.e6 ) and ( frequency < 1500.e6 ):
      max_smearing = 0.1
      outlier_fluxes = [ 0.3, 0.1 ] # Jy
      outlier_ratios = [ 5., 2. ]
      model_flux_min = 0.005 # Jy
      outlier_flux_min = 0.025 # Jy
      clip_levels = [ 1.5, 0.8 ] # Jy
      clip_uv_min = 4. # klambda
      rfi_level = 0.1 # Jy
      amplitude_snr = 1000.
      use_catalog = [ False, False ]
      catalog_resolution = 15.
      plot_range = [ -0.00025, 0.0025 ]
    else:
      raise error( 'unsupported frequency' )
  elif ( uv.header.telescop == 'VLA' ):
    frequency = get_central_frequency( uv )
    if ( frequency > 300.e6 ) and ( frequency < 400.e6 ):
      max_smearing = 0.15
      outlier_fluxes = [ 0.7, 0.2 ]
      outlier_ratios = [ 3., 1.5 ]
      model_flux_min = 0.1
      outlier_flux_min = 0.1
      clip_levels = [ 7., 3. ]
      clip_uv_min = 1.
      rfi_level = 1. # Jy
      amplitude_snr = 1000.
      use_catalog = [ True, False ]
      catalog_resolution = 20.
      plot_range = [ -0.00025, 0.0025 ]
  elif ( uv.header.telescop == 'EVLA' ):
    frequency = get_central_frequency( uv )
    if ( frequency > 40.e6 ) and ( frequency < 90.e6 ):
      max_smearing = 0.2
      outlier_fluxes = [ 3., 1.5 ]
      outlier_ratios = [ 3., 1.5 ]
      model_flux_min = 0.2
      outlier_flux_min = 0.2
      clip_levels = [ 50., 25. ]
      clip_uv_min = 0.5
      rfi_level = 10. # Jy
      amplitude_snr = 500
      use_catalog = [ True, True ]
      catalog_resolution = 50.
      plot_range = [ -0.005, 0.05 ]
  else:
    raise error( 'unsupported telescope' )
  if snapshot:
    time_step = 1
#    sidelobe_rejection = 1.25
    sidelobe_rejection = 1.5
    amplitude_snr = 0.5 * amplitude_snr
    plot_range = [ 2. * x for x in plot_range ]
    use_catalog[ 1 ] = False
  else:
    time_step = 4
#    sidelobe_rejection = 1.1
    sidelobe_rejection = 1.25
  [ uv_min, uv_max ] = [ 0., 1.e6 ]
  try:
    [ uv_min, uv_max ] = imagr_params[ 'uvrang' ]
  except KeyError:
    pass
  try:
    [ uv_min, uv_max ] = imagr_params[ 'uvrange' ]
  except KeyError:
    pass
  if ( len( model_source_list ) > 0 ):
    catalog_resolution = model_resolution
  if try_new:
    init_uv_min = 1.4 / get_central_wavelength( uv )
    init_uv_max = ( 4. * 60. / restore_parameter( uv, 'beam_size' ) ) / 2.
    init_uv_max = max( ( 4. * 60. / model_resolution ) * 2., init_uv_max )
    init_calib_params = { 'uvrange' : [ init_uv_min, init_uv_max ] }
  else:
    init_calib_params = { 'uvrange' : [ 0.5, 25. ] }
  
  # cover primary beam area with number of facets
  # add facets at places of interfering outlier sources
  define_pb_facets( uv, max_smearing = max_smearing )
  define_o_facets( uv, outlier_fluxes[ 0 ], field_ratio = outlier_ratios[ 0 ],
      use_nvss = use_nvss, use_vlss = True )
  define_o_facets( uv, outlier_fluxes[ 1 ], field_ratio = outlier_ratios[ 1 ],
      append = True, use_nvss = use_nvss, use_vlss = True )
  define_a_s_facets( uv )
  
  # check if target field is near anything bright
  a_facet_file_name_e = path.expandvars( restore_parameter( uv, 'a_facet_file_name' ) )
  s_facet_file_name_e = path.expandvars( restore_parameter( uv, 's_facet_file_name' ) )
  o_facet_file_name_e = path.expandvars( restore_parameter( uv, 'o_facet_file_name' ) )
  o_facet_size = restore_parameter( uv, 'o_facet_size' )
  as_facet_list = get_facet_list( a_facet_file_name_e )
  as_facet_list = as_facet_list + get_facet_list( s_facet_file_name_e )
  near_bright_radius = outlier_ratios[ 0 ] * restore_parameter( uv, 'field_size' ) / 2.
  near_bright = False
  target_radec = get_radec( uv )
  for [ i, radec ] in as_facet_list:
    [ r, p ] = calculate_angular_separation( target_radec, radec )
    if ( r < near_bright_radius ):
      print 'WARNING: target field center is %4.1f degrees from a very bright source' % ( r )
      near_bright = True
      add_facet( o_facet_file_name_e, radec, [ o_facet_size, o_facet_size ] )
      o_facet_count = restore_parameter( uv, 'o_facet_count' )
      store_parameter( uv, 'o_facet_count', o_facet_count + 1 )
  
  # make primary beam sky model
  pbm_facets = make_pb_model( uv, model_flux_min, imagr_params = imagr_params,
      use_nvss = use_nvss, use_vlss = True )
  model_radec_list = []
  model_flux = 0.
  for i in range( 1, 1 + restore_parameter( pbm_facets, 'facet_count' ) ):
    facet = get_facet( pbm_facets, i )
    model_radec_list.append( get_radec( facet ) )
    model_flux = model_flux + get_model_flux( facet )
  
  # clip visibility amplitudes
  flagver = max( 1, uv.table_highver( 'FG' ) )
  call_aips_task( 'CLIP', indata = uv, flagver = flagver, stokes = 'I', 
      aparm = [ 3. * model_flux, 0,0.001,0,0,0,0,0,0,0 ], outfgver = flagver )
  flag_uv = apply_flag_table( uv )
  uv.zap()
  flag_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  
  # make alternative sky model
  if ( len( model_source_list ) > 0 ):
    remove_facets( pbm_facets )
    pbm_facets = make_pb_model_from_source_list( uv, model_source_list,
      imagr_params = imagr_params )
    model_radec_list = []
    for i in range( 1, 1 + restore_parameter( pbm_facets, 'facet_count' ) ):
      facet = get_facet( pbm_facets, i )
      model_radec_list.append( get_radec( facet ) )
  
  # save model image
  if ( not minimize_storage ):
    determine_facet_overlap( pbm_facets )
    rpbm_facets = restore_model_components( pbm_facets, imagr_params = imagr_params )
    store_parameter( uv, 'cpb_noise', 1. )
    rpbm_image = combine_facets( uv, rpbm_facets, save_tables = False, save_info = True )
    write_fits_image( rpbm_image, fits_path + target_name + '.PBM.FITS',
        include_tables = True )
    remove_facets( rpbm_facets )
    rpbm_image.zap()
  
  # do full primary beam calibration
  integration_time = restore_parameter( uv, 'integration_time' )
  init_phase_interval_min = 0.1 * integration_time / 60.
  init_calib_params = { 'uvrange' : [ 0.5, 25. ] }
  reference_antenna = 0
  calibrate_model( uv, pbm_facets, sigma = 0., conversion_method = 'GRID',
      calib_params = init_calib_params, reference_antenna = reference_antenna, 
      sidelobe_rejection = 0., phase_interval = init_phase_interval_min,
      apply_solutions = False, snr_limit = 0.1 )
  remove_facets( pbm_facets )
  reference_antenna = select_reference_antenna( uv, count_power = 2.,
      antenna_list = reference_antenna_list )
  re_reference_solutions( uv, reference_antenna )
  
  # make gain phase plots
  if do_plots:
    day_min = int( floor( restore_parameter( uv, 'time_min' ) ) )
    day_max = int( ceil( restore_parameter( uv, 'time_max' ) ) )
    for day in range( day_min, day_max, 1 ):
      pl_version_low = uv.table_highver( 'PL' ) + 1
      try:
        call_aips_task( 'SNPLT', indata = uv, nplots = 8, dotv = -1, factor = 0.3,
            optype = 'PHAS', pixrange = [ -180,180 ], inext = 'SN', invers = 0,
            timerang = [ day,0,0,0, day + 1,0,0,0 ] )
      except RuntimeError:
        continue
      pl_version_high = uv.table_highver( 'PL' )
      plot_file_name = prt_path + aips_file_name_to_string( uv ) + '_DAY%02d.PS' % ( day )
      if file_exists( plot_file_name ):
        remove_file( plot_file_name )
      call_aips_task( 'LWPLA', indata = uv, plver = pl_version_low, 
          invers = pl_version_high, lpen = 2, outfile = plot_file_name )
      for pl_version in range( pl_version_low, pl_version_high + 1 ):
        uv.zap_table( 'PL', pl_version )
  
  # measure initial noise
  cpb_facets = image_cpb_facets( uv, imagr_params = imagr_params, 
     apply_solutions = True )
  cpb_noise = measure_cpb_noise( uv, cpb_facets, keep_image = ( not minimize_storage ) )
  remove_facets( cpb_facets )
  if ( not minimize_storage ):
    cpb_image = get_aips_file( aips_disk, 'CPB', 'FLATN', 0, 'MA' )
    write_fits_image( cpb_image, fits_path + target_name + '.CPB.FITS',
        include_tables = True )
  
  # image primary beam and nearby outliers
  [ pb_facets, o_facets ] = image_clean_pbo_facets( uv, sigma = 3.,
      conversion_method = 'GRID', imagr_params = imagr_params,
      sidelobe_rejection = sidelobe_rejection, clean_box_radius = box_radius )
  mc1_image = combine_facets( uv, pb_facets, save_tables = False, save_info = True )
  mc1_image.rename( klass = 'MC1' )
  write_fits_image( mc1_image, fits_path + target_name + '.MC1.FITS',
      include_tables = True )
  write_fits_uv( uv, fits_path + target_name + '.MC1.UVFITS' )
  
  # calibrate against true sky model
  integration_time = restore_parameter( uv, 'integration_time' )
  phase_interval_min = 0.1 * integration_time / 60.
  calibrate_model( uv, pb_facets, reference_antenna = reference_antenna,
      calib_params = calib_params, phase_interval = phase_interval_min,
      apply_solutions = False )
  sc_version = uv.table_highver( 'SN' )
  
  # remove instrumental phase
  day_min = int( floor( restore_parameter( uv, 'time_min' ) ) )
  day_max = int( ceil( restore_parameter( uv, 'time_max' ) ) )
  for day in range( day_min, day_max, 2 ):
    filter_instrumental_phases( uv, max_period = 30., max_range = 12. * 60.,
        time_range = dhmsdhms_to_timetime( [ day,0,0,0, day + 2,0,0,0 ] ),
        in_version = sc_version, time_step = time_step, min_points = 3, 
        min_period = 2., phase_rms_criterion = 7. )
    re_reference_solutions( uv, reference_antenna )
    if ( day > day_min ):
      combine_solutions( uv )
  cal_uv = apply_solution_table( uv )
  combine_solutions( uv, in_version_1 = sc_version, invert_2 = True )
  call_aips_task( 'TACOP', indata = uv, outdata = cal_uv, inext = 'SN', ncount = 1 )
  uv.zap()
  cal_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  
  # add some nearby outlier sources to PB area
  # subtract other nearby outliers
  add_o_facet_list = []
  for i in range( 1, 1 + restore_parameter( o_facets, 'facet_count' ) ):
    facet = get_facet( o_facets, i )
    print i, get_model_flux( facet )
    if ( get_model_flux( facet ) > outlier_flux_min ):
      print '... added to primary beam facets'
      add_o_facet_list.append( i )
  if ( len( add_o_facet_list ) > 0 ):
    add_centered_pb_facets( uv, pb_facets, 1.e6, imagr_params = imagr_params,
        center_facets = o_facets, facet_list = add_o_facet_list,
        keep_model_components = True )
    store_parameter( uv, 'added_pb_facet_count', 0 )
    store_parameter( pb_facets, 'added_facet_count', 0 )
  o_facet_list = range( 1, 1 + restore_parameter( o_facets, 'facet_count' ) )
  for x in add_o_facet_list:
    o_facet_list.remove( x )
  sub_o_uv = subtract_model( uv, o_facets, facet_list = o_facet_list )
  uv.zap()
  sub_o_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  remove_facets( o_facets )
  
  # subtract pb sources for RFI mitigation
  sub_pb_uv = subtract_model( uv, pb_facets )
  flag_uv = apply_flag_table( uv )
  uv.zap()
  flag_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  
  # flagging
  flag_image_baselines( sub_pb_uv, keep_images = ( not minimize_storage ),
      scaling_factor = 32., imagr_params = imagr_params, flag_version = -1 )
  flag_image_baselines( sub_pb_uv, keep_images = ( not minimize_storage ),
      scaling_factor = 16., imagr_params = imagr_params )
  flag_image_baselines( sub_pb_uv, keep_images = ( not minimize_storage ),
      scaling_factor = 8., imagr_params = imagr_params )
  flag_image_baselines( sub_pb_uv, keep_images = ( not minimize_storage ),
      scaling_factor = 4., imagr_params = imagr_params )
  flag_image_baselines( sub_pb_uv, keep_images = ( not minimize_storage ),
      scaling_factor = 2., imagr_params = imagr_params )
  flag_image_baselines( sub_pb_uv, keep_images = ( not minimize_storage ),
      scaling_factor = 1., imagr_params = imagr_params )
  call_aips_task( 'TACOP', indata = sub_pb_uv, outdata = uv, inext = 'FG',
      ncount = 1 )
  extend_flags( uv )
  flag_uv = apply_flag_table( uv )
  uv.zap()
  flag_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  sub_pb_uv.zap()
  
  # re-image primary beam
  re_image_clean_pb_facets( uv, pb_facets, sigma = 3., imagr_params = imagr_params,
      conversion_method = 'GRID' )
  sc1_image = combine_facets( uv, pb_facets, save_tables = False, save_info = True )
  sc1_image.rename( klass = 'SC1' )
  write_fits_image( sc1_image, fits_path + uv.name + '.SC1.FITS',
      include_tables = True )
  write_fits_uv( uv, fits_path + target_name + '.SC1.UVFITS' )
  if minimize_storage:
    mc1_image.zap()
    remove_file( fits_path + target_name + '.MC1.FITS' )
    remove_file( fits_path + target_name + '.MC1.UVFITS' )
  
  # determine amplitude calibration
  calibrate_model( uv, pb_facets, reference_antenna = reference_antenna,
      calib_params = calib_params, phase_interval = phase_interval_min,
      apply_solutions = False )
  flag_uv = apply_flag_table( uv )
  uv.zap()
  flag_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  reference_antenna = select_reference_antenna( uv, count_power = 2.,
      antenna_list = reference_antenna_list )
  re_reference_solutions( uv, reference_antenna )
  calibrate_model( uv, pb_facets, reference_antenna = reference_antenna,
      calib_params = calib_params, do_amplitude = True, normalize_gains = False,
      amplitude_interval = phase_interval_min )
  flag_uv = apply_flag_table( uv )
  uv.zap()
  flag_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  filter_solutions( uv, amplitude_limits = [ 0.,3. ], phase_limits = [ -50.,50. ],
      amplitude_sigma = None, phase_sigma = None )
  filter_solutions( uv, amplitude_sigma = 5., phase_sigma = 5.,
      amplitude_window = 12. * 60. * 60., phase_window = 12. * 60. * 60. )
  filter_solutions( uv, amplitude_sigma = 5., phase_sigma = 5.,
      amplitude_window = 60. * 60., phase_window = 60. * 60. )
  reference_antenna = select_reference_antenna( uv, count_power = 2.,
      antenna_list = reference_antenna_list )
  re_reference_solutions( uv, reference_antenna )
  
  # apply amplitude calibration
  smooth_solutions_in_time( uv, phase_window = 0., amplitude_window = 5. * 60., order = 1 )
  normalize_solutions( uv )
  sn_version = uv.table_highver( 'SN' )
  call_aips_task( 'SNCOR', indata = uv, snver = sn_version, opcode = 'ZPHS' )
  cal_uv = apply_solution_table( uv )
  uv.zap()
  cal_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  
  # calibrate and re-image
  calibrate_model( uv, pb_facets, reference_antenna = reference_antenna,
      calib_params = calib_params, phase_interval = phase_interval_min,
      apply_solutions = False )
  add_clean_boxes( uv, pb_facets, sidelobe_rejection = sidelobe_rejection,
      clean_box_radius = box_radius )
  if do_large_boxes:
    add_clean_boxes( uv, pb_facets, convolve_factor = 3., keep_boxes = True,
        clean_box_radius = large_box_radius, sidelobe_rejection = sidelobe_rejection )
  re_image_clean_pb_facets( uv, pb_facets, sigma = 3., imagr_params = imagr_params,
      conversion_method = 'GRID' )
  sc2_image = combine_facets( uv, pb_facets, save_tables = False, save_info = True )
  sc2_image.rename( klass = 'SC2' )
  write_fits_image( sc2_image, fits_path + target_name + '.SC2.FITS',
      include_tables = True )
  write_fits_uv( uv, fits_path + target_name + '.SC2.UVFITS' )
  if minimize_storage:
    sc1_image.zap()
    remove_file( fits_path + target_name + '.SC1.FITS' )
    remove_file( fits_path + target_name + '.SC1.UVFITS' )
  
  # subtract pb sources for RFI mitigation
  sub_pb_uv = subtract_model( uv, pb_facets )
  flag_uv = apply_flag_table( uv )
  uv.zap()
  flag_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  
  # clip visibility amplitudes
  flagver = max( 1, sub_pb_uv.table_highver( 'FG' ) )
  call_aips_task( 'CLIP', indata = sub_pb_uv, flagver = flagver, outfgver = flagver,
      aparm = [ clip_levels[ 0 ], 0,0.001,0,0,0,0,0,0,0 ], stokes = 'I',
      uvrange = [ uv_min, 1.e6 ] )
  call_aips_task( 'CLIP', indata = sub_pb_uv, flagver = flagver, outfgver = flagver,
      aparm = [ clip_levels[ 1 ], 0,0.001,0,0,0,0,0,0,0 ], stokes = 'I',
      uvrange = [ max( clip_uv_min, uv_min ), 1.e6 ] )
  call_aips_task( 'TACOP', indata = sub_pb_uv, outdata = uv, inext = 'FG',
      invers = 0, outvers = 0, ncount = 1 )
  flag_uv = apply_flag_table( uv )
  uv.zap()
  flag_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  flag_uv = apply_flag_table( sub_pb_uv )
  sub_pb_uv.zap()
  flag_uv.rename( name = sub_pb_uv.name, klass = sub_pb_uv.klass, seq = sub_pb_uv.seq )
  
  # flagging
  flag_image_baselines( sub_pb_uv, keep_images = ( not minimize_storage ),
      scaling_factor = 32., imagr_params = imagr_params, flag_version = -1 )
  flag_image_baselines( sub_pb_uv, keep_images = ( not minimize_storage ),
      scaling_factor = 16., imagr_params = imagr_params )
  flag_image_baselines( sub_pb_uv, keep_images = ( not minimize_storage ),
      scaling_factor = 8., imagr_params = imagr_params )
  flag_image_baselines( sub_pb_uv, keep_images = ( not minimize_storage ),
      scaling_factor = 4., imagr_params = imagr_params )
  flag_image_baselines( sub_pb_uv, keep_images = ( not minimize_storage ),
      scaling_factor = 2., imagr_params = imagr_params )
  flag_image_baselines( sub_pb_uv, keep_images = ( not minimize_storage ),
      scaling_factor = 1., imagr_params = imagr_params )
  extend_flags( sub_pb_uv )
  call_aips_task( 'TACOP', indata = sub_pb_uv, outdata = uv, inext = 'FG', ncount = 1 )
  flag_uv = apply_flag_table( uv )
  uv.zap()
  flag_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  flag_uv = apply_flag_table( sub_pb_uv )
  sub_pb_uv.zap()
  flag_uv.rename( name = sub_pb_uv.name, klass = sub_pb_uv.klass, seq = sub_pb_uv.seq )
  
  # run Obit task LowFRFI() to subtract lower level RFI
  rfi_uv = get_aips_file( aips_disk, target_name, 'LFRFI', -1, 'UV' )
  rfi_avg = 10. # min
  obit_script_name = 'python/obit_lowfrfi.py'
  if file_exists( obit_script_name ):
    remove_file( obit_script_name )
  obit_file = file( obit_script_name, 'w' )
  obit_file.write( "from ObitTask import *\n" )
  obit_file.write( "lr = ObitTask( 'LowFRFI' )\n" )
  obit_file.write( "lr.DataType = 'AIPS'\n" )
  obit_file.write( "lr.userno = %s\n" % ( repr( get_aips_userid() ) ) )
  obit_file.write( "lr.inDisk = %s\n" % ( repr( uv.disk ) ) )
  obit_file.write( "lr.inName = %s\n" % ( "'" + uv.name + "'" ) )
  obit_file.write( "lr.inClass = %s\n" % ( "'" + uv.klass + "'" ) )
  obit_file.write( "lr.inSeq = %s\n" % ( repr( uv.seq ) ) )
  obit_file.write( "lr.in2Disk = %s\n" % ( repr( sub_pb_uv.disk ) ) )
  obit_file.write( "lr.in2Name = %s\n" % ( "'" + sub_pb_uv.name + "'" ) )
  obit_file.write( "lr.in2Class = %s\n" % ( "'" + sub_pb_uv.klass + "'" ) )
  obit_file.write( "lr.in2Seq = %s\n" % ( repr( sub_pb_uv.seq ) ) )
  obit_file.write( "lr.outDisk = %s\n" % ( repr( rfi_uv.disk ) ) )
  obit_file.write( "lr.outName = %s\n" % ( "'" + rfi_uv.name + "'" ) )
  obit_file.write( "lr.outClass = %s\n" % ( "'" + rfi_uv.klass + "'" ) )
  obit_file.write( "lr.outSeq = %s\n" % ( repr( rfi_uv.seq ) ) )
  obit_file.write( "lr.flagVer = -1\n" )
  obit_file.write( "lr.flagRes = -1\n" )
  obit_file.write( "lr.timeInt = %s\n" % ( repr( integration_time ) ) )
  obit_file.write( "lr.doInvert = True\n" )
  obit_file.write( "lr.timeAvg = %s\n" % ( repr( rfi_avg ) ) )
  obit_file.write( "lr.minRFI = %s\n" % ( repr( rfi_level ) ) )
  obit_file.write( "lr.go()\n" )
  obit_file.close()
  if ( sys.stdout.name == '<stdout>' ):
    result = system( 'start_obittalk.sh %s %d %s' % ( getenv( 'PWD' ),
        get_aips_userid(), obit_script_name ) )
  else:
    result = system( 'start_obittalk.sh %s %d %s >> %s 2>&1' % ( getenv( 'PWD' ),
        get_aips_userid(), obit_script_name, sys.stdout.name ) )
    sys.stdout.flush()
    sys.stderr.flush()
  if ( result != 0 ):
#    raise error( 'Obit RFI subtraction failed' )
    print 'WARNING: Obit RFI subtraction failed'
  
  # copy essential tables and replace
  if rfi_uv.exists():
    call_aips_task( 'TACOP', indata = uv, outdata = rfi_uv, inext = 'SN', ncount = 1 )
    call_aips_task( 'TACOP', indata = uv, outdata = rfi_uv, inext = 'PS', ncount = 1 )
    uv.zap()
    rfi_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  sub_pb_uv.zap()
  
  # calibrate bandpass on target field
  model_uv = make_model_uv( uv, pb_facets, apply_solutions = True )
  div_uv = divide_uv( uv, model_uv )
  call_aips_task( 'INDXR', indata = div_uv, cparm = [ 12. * 60.,20. * 60.,-1,0,0,0 ] )
  call_aips_task( 'BPASS', indata = div_uv, flagver = 0, dopol = -1,
      blver = -1, docalib = -1, gainuse = 0, doband = -1, solint = 0.,
      soltype = 'L1', outvers = 0, bpassprm = [ 0,0,1,0,2,0,0,1,1,3,1 ],
      specindx = 1.e-3, refant = reference_antenna, cmethod = 'DFT',
      **calib_params )
  call_aips_task( 'TACOP', indata = div_uv, outdata = uv, inext = 'BP', ncount = 1 )
  div_uv.zap()
  model_uv.zap()
  temp_uv = get_aips_file( aips_disk, target_name, 'TEMP', -1, 'UV' )
  call_aips_task( 'SPLAT', indata = uv, smooth = [ 0 ], douvcomp = 0,
      aparm = [ 3,0,1,0,0,0,1 ], channel = 1, chinc = 1, docalib = -1,
      doband = 4, bpver = 0, flagver = 0, outdata = temp_uv )
  uv.zap()
  temp_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  
  # make global astrometric correction
  source_flux_min = 10. * restore_parameter( uv, 'cpb_noise' )
  source_list = get_source_list_from_facets( pb_facets, source_flux_min )
  radec_list = []
  for x in source_list:
    facet = get_facet( pb_facets, x[ 0 ] )
    radec_list.append( calculate_source_radec( facet, x[ 1 ] ) )
  print 'detected %d sources for astrometry check' % ( len( radec_list ) )
  print 'using %d catalog sources for astrometry check' % ( len( model_radec_list ) )
  global_offset = None
  for beam_factor in [ 2., 6., 18., 54. ]:
    try:
      global_offset = check_astrometry( sc2_image, radec_list, beam_factor = beam_factor,
        ref_radec_list = model_radec_list, plot_offsets = False, min_points = 8,
        use_nvss = use_nvss )
    except RuntimeError:
      continue
    except:
      raise error( traceback.print_exc() )
    else:
      [ dra_median, ddec_median, dr_std ] = global_offset
      if ( dr_std > catalog_resolution / 2. ):
        continue
      break
  if ( not global_offset is None ):
    offset = [ -dra_median / 3600., -ddec_median / 3600. ]
    pb_facet_file_name = restore_parameter( uv, 'pb_facet_file_name' )
    pb_facet_file_name_e = path.expandvars( pb_facet_file_name )
    for i in range( 1, 1 + restore_parameter( pb_facets, 'facet_count' ) ):
      facet = get_facet( pb_facets, i )
      cor_radec = calculate_offset_position_dradec( get_radec( facet ), offset )
      replace_facet( pb_facet_file_name_e, i, get_image_size( facet ), cor_radec,
          keep_boxes = True )
      set_radec( facet, cor_radec )
  
  # calibrate and re-image
  calibrate_model( uv, pb_facets, reference_antenna = reference_antenna,
      calib_params = calib_params, phase_interval = phase_interval_min,
      apply_solutions = False )
  add_clean_boxes( uv, pb_facets, sidelobe_rejection = sidelobe_rejection,
      clean_box_radius = box_radius )
  if do_large_boxes:
    add_clean_boxes( uv, pb_facets, convolve_factor = 3., keep_boxes = True,
        clean_box_radius = large_box_radius, sidelobe_rejection = sidelobe_rejection )
  re_image_clean_pb_facets( uv, pb_facets, sigma = 3., imagr_params = imagr_params,
      conversion_method = 'GRID' )
  sc3_image = combine_facets( uv, pb_facets, save_tables = False, save_info = True )
  sc3_image.rename( klass = 'SC3' )
  write_fits_image( sc3_image, fits_path + target_name + '.SC3.FITS',
      include_tables = True )
  write_fits_uv( uv, fits_path + target_name + '.SC3.UVFITS' )
  if minimize_storage:
    sc2_image.zap()
    remove_file( fits_path + target_name + '.SC2.FITS' )
    remove_file( fits_path + target_name + '.SC2.UVFITS' )
  
  # update calibration and reference antenna before peeling
  calibrate_model( uv, pb_facets, reference_antenna = reference_antenna,
      calib_params = calib_params, phase_interval = phase_interval_min,
      apply_solutions = False )
  reference_antenna = select_reference_antenna( uv, count_power = 2.,
      antenna_list = reference_antenna_list )
  re_reference_solutions( uv, reference_antenna )
  
  # peel sources in PB area
  sub_pb_uv = subtract_model( uv, pb_facets )
  flag_uv = apply_flag_table( uv )
  uv.zap()
  flag_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  if use_catalog[ 0 ]:
    convl_facets = convolve_facets( pb_facets, resolution = catalog_resolution )
    beam_size = get_beam_size( pb_facets )
    search_ratio = ( 10. * sqrt( beam_size[ 0 ]**2 + beam_size[ 1 ]**2 ) / 
        ( 2. * catalog_resolution ) )
  else:
    convl_facets = pb_facets
    search_ratio = 10.
  [ ppb_uv, ppb_facets ] = peel_pbo_facets( sub_pb_uv, convl_facets,
      phase_interval_min = phase_interval_min, reference_antenna = reference_antenna, 
      imagr_params = imagr_params, calib_params = calib_params, max_relativity = 1,
      sidelobe_rejection = sidelobe_rejection, conversion_method = 'GRID',
      amplitude_snr = 2. * amplitude_snr, double_search_ratio = search_ratio,
      use_catalog = use_catalog[ 0 ], catalog_resolution = catalog_resolution,
      use_nvss = use_nvss, input_catalog = model_source_list,
      clean_box_radius = box_radius )
  if use_catalog[ 0 ]:
    remove_facets( convl_facets )
  sub_pb_uv.zap()
  if ( ppb_facets is None ):
    raise error( 'no peeled sources' )
  copy_solution_tables( ppb_facets, ppb_uv )
  determine_facet_overlap( ppb_facets )
  
  # check peeled sources
  if use_catalog[ 0 ]:
    convl_facets = convolve_facets( ppb_facets, resolution = catalog_resolution )
  else:
    convl_facets = ppb_facets
  match_types = []
  offsets = []
  for i in range( 1, 1 + restore_parameter( convl_facets, 'facet_count' ) ):
    facet = get_facet( convl_facets, i )
    facet_ref_pos = get_pixel_reference( facet )
    [ facet_max, facet_max_pos ] = get_image_maximum( facet )
    offset = [ facet_max_pos[ 0 ] - facet_ref_pos[ 0 ], 
        facet_max_pos[ 1 ] - facet_ref_pos[ 1 ] ]
    offsets.append( offset )
    match_type = restore_parameter( facet, 'match_type' )
    match_types.append( match_type.split( ', ' ) )
    print '%2d: %6.4f Jy, %4.2f min, %s, %4d, %s' % ( i, get_model_flux( facet ),
        restore_parameter( facet, 'solution_interval' ), match_type, 
        sqrt( offset[ 0 ]**2 + offset[ 1 ]**2 ), repr( get_facet_overlap( facet ) ) )
#    if do_plots:
#      facet = get_facet( ppb_facets, i )
#      plot_image( facet, z_range = plot_range,
#          plot_title = '%s source %d' % ( target_name, i ) )
#      plot_source( pb_facets, get_radec( facet ), z_range = plot_range )
  if use_catalog[ 0 ]:
    remove_facets( convl_facets )
  
  # select peeled sources for ionospheric modeling
  if use_catalog[ 0 ]:
    offset_limits = [ 1., 2. ]
  else:
    offset_limits = [ 1.5, 3. ]
  facet_list = range( 1, 1 + restore_parameter( ppb_facets, 'facet_count' ) )
  if ( len( facet_list ) == 0 ):
    raise error( 'no peeled sources to fit ionosphere model' )
  if ( len( facet_list ) >= 4 ):
    for i in range( 1, 1 + restore_parameter( ppb_facets, 'facet_count' ) ):
      facet = get_facet( ppb_facets, i )
      if ( match_types[ i - 1 ][ 1 ] == 'no catalog match' ):
        facet_list.remove( i )
  if ( len( facet_list ) >= 4 ):
    sel = ( array( facet_list ) - 1 ).reshape( ( len( facet_list ), 1 ) )
    offsets = array( offsets )
    median_offset = median( aget( offsets, sel ), axis = 0 )
    offsets = sqrt( add.reduce( ( offsets - median_offset )**2, axis = 1 ) )
    for i in range( 1, 1 + restore_parameter( ppb_facets, 'facet_count' ) ):
      if ( not i in facet_list ):
        continue
      if ( match_types[ i - 1 ][ 0 ] == 'single source' ):
        if ( offsets[ i - 1 ] > offset_limits[ 0 ] ):
          facet_list.remove( i )
      elif ( match_types[ i - 1 ][ 0 ] == 'double source' ):
        if ( offsets[ i - 1 ] > offset_limits[ 1 ] ):
          facet_list.remove( i )
  if ( len( facet_list ) >= 4 ):
    for i in range( 1, 1 + restore_parameter( ppb_facets, 'facet_count' ) ):
      if ( not i in facet_list ):
        continue
      facet = get_facet( ppb_facets, i )
      source_facet_list = [ x[ 0 ] for x in 
          find_source_facets( ppb_facets, get_radec( facet ) ) ]
      while ( i != source_facet_list[ 0 ] ):
        if ( source_facet_list[ 0 ] in facet_list ):
          break
        source_facet_list = source_facet_list[ 1 : ]
      if ( i != source_facet_list[ 0 ] ):
        facet_list.remove( i )
  print 'selected facets = ', facet_list
  
  # what to do if there are too few peeled sources
  if ( len( facet_list ) < 4 ):
    
    # re-image using peeling solutions and selfcal solutions
    add_clean_boxes( uv, pb_facets, sidelobe_rejection = sidelobe_rejection,
        clean_box_radius = box_radius )
    if do_large_boxes:
      add_clean_boxes( uv, pb_facets, convolve_factor = 3., keep_boxes = True,
          clean_box_radius = large_box_radius, sidelobe_rejection = sidelobe_rejection )
    add_centered_pb_facets( uv, pb_facets, 1.e6, print_info = True,
        imagr_params = imagr_params, center_facets = ppb_facets, min_separation = 10. )
    cpb_noise = restore_parameter( uv, 'cpb_noise' )
    snr_max = 0.
    i_max = 1
    for i in range( 1, 1 + restore_parameter( ppb_facets, 'facet_count' ) ):
      facet = get_facet( ppb_facets, i )
      flux = max( get_model_flux( facet ), 0. )
      [ facet_max, facet_max_pos ] = get_image_maximum( facet )
      snr = facet_max * sqrt( flux / facet_max ) / cpb_noise
      if ( snr > snr_max ):
        snr_max = snr
        i_max = i
    if ( snr_max > amplitude_snr ):
      facet = get_facet( ppb_facets, i_max )
      call_aips_task( 'TACOP', indata = facet, outdata = uv, inext = 'SN', ncount = 1 )
      sn_version = uv.table_highver( 'SN' )
      call_aips_task( 'SNCOR', indata = uv, snver = sn_version, opcode = 'NORM' )
      if ( not i_max in facet_list ):
        facet_list = [ i_max ]
    for i in range( 1, 1 + restore_parameter( pb_facets, 'facet_count' ) ):
      facet = get_facet( pb_facets, i )
      call_aips_task( 'TACOP', indata = uv, outdata = facet, inext = 'SN', ncount = 1 )
    replace_model_solutions_with_peel_solutions( uv, pb_facets, ppb_facets,
        facet_list = facet_list )
    re_image_clean_pb_facets( uv, pb_facets, imagr_params = imagr_params,
        conversion_method = 'GRID' )
    sp0_image = combine_facets( uv, pb_facets, save_tables = False, save_info = True )
    sp0_image.rename( klass = 'SP0' )
    write_fits_image( sp0_image, fits_path + target_name + '.SP0.FITS',
        include_tables = True )
    write_fits_uv( uv, fits_path + target_name + '.SP0.UVFITS' )
    if minimize_storage:
      sc3_image.zap()
      remove_file( fits_path + target_name + '.SC3.FITS' )
      remove_file( fits_path + target_name + '.SC3.UVFITS' )
    
    # make global astrometric correction
    source_flux_min = 10. * restore_parameter( uv, 'cpb_noise' )
    source_list = get_source_list_from_facets( pb_facets, source_flux_min )
    radec_list = []
    for x in source_list:
      facet = get_facet( pb_facets, x[ 0 ] )
      radec_list.append( calculate_source_radec( facet, x[ 1 ] ) )
    print 'detected %d sources for astrometry check' % ( len( radec_list ) )
    print 'using %d catalog sources for astrometry check' % ( len( model_radec_list ) )
    global_offset = None
    for beam_factor in [ 2., 6., 18., 54. ]:
      try:
        global_offset = check_astrometry( sp0_image, radec_list, beam_factor = beam_factor,
          ref_radec_list = model_radec_list, plot_offsets = False, min_points = 5,
          use_nvss = use_nvss )
      except RuntimeError:
        continue
      except:
        raise error( traceback.print_exc() )
      else:
        [ dra_median, ddec_median, dr_std ] = global_offset
        if ( dr_std > catalog_resolution / 2. ):
          continue
        break
    if ( not global_offset is None ):
      offset = [ -dra_median / 3600., -ddec_median / 3600. ]
      pb_facet_file_name = restore_parameter( uv, 'pb_facet_file_name' )
      pb_facet_file_name_e = path.expandvars( pb_facet_file_name )
      for i in range( 1, 1 + restore_parameter( pb_facets, 'facet_count' ) ):
        facet = get_facet( pb_facets, i )
        cor_radec = calculate_offset_position_dradec( get_radec( facet ), offset )
        replace_facet( pb_facet_file_name_e, i, get_image_size( facet ), cor_radec,
            keep_boxes = True )
        set_radec( facet, cor_radec )
      calibrate_model( uv, pb_facets, reference_antenna = reference_antenna,
          calib_params = calib_params, phase_interval = phase_interval_min,
          apply_solutions = True )
      flag_uv = apply_flag_table( uv )
      uv.zap()
      flag_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
      smooth_solutions_in_time( uv, max_gap = 12. * 60., phase_window = 10. * 60. )
      re_sample_solutions( uv, gap_time = 12. * 60. * 60, weight_multiplier = 0.,
          interpolation_method = 'nearest' )
      re_reference_solutions( uv, reference_antenna )
      cor_uv = apply_solution_table( uv )
      uv.zap()
      cor_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
      calibrate_model( uv, pb_facets, reference_antenna = reference_antenna,
          calib_params = calib_params, phase_interval = phase_interval_min,
          apply_solutions = False, keep_flags = False )
    use_catalog[ 0 ] = False
    
    # peel sources in PB area
    sub_pb_uv = subtract_model( uv, pb_facets )
    flag_uv = apply_flag_table( uv )
    uv.zap()
    flag_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
    if use_catalog[ 0 ]:
      convl_facets = convolve_facets( pb_facets, resolution = catalog_resolution )
      beam_size = get_beam_size( pb_facets )
      search_ratio = ( 10. * sqrt( beam_size[ 0 ]**2 + beam_size[ 1 ]**2 ) / 
          ( 2. * catalog_resolution ) )
    else:
      convl_facets = pb_facets
      search_ratio = 10.
    [ ppb_uv, ppb_facets ] = peel_pbo_facets( sub_pb_uv, convl_facets, 
        phase_interval_min = phase_interval_min, reference_antenna = reference_antenna, 
        imagr_params = imagr_params, calib_params = calib_params, max_relativity = 3,
        sidelobe_rejection = sidelobe_rejection, conversion_method = 'GRID',
        amplitude_snr = amplitude_snr, double_search_ratio = search_ratio,
        use_catalog = use_catalog[ 0 ], catalog_resolution = catalog_resolution,
        phase_interval_max = 4., use_nvss = use_nvss, input_catalog = model_source_list,
        clean_box_radius = box_radius )
    if use_catalog[ 0 ]:
      remove_facets( convl_facets )
    sub_pb_uv.zap()
    if ( ppb_facets is None ):
      raise error( 'no peeled sources' )
    copy_solution_tables( ppb_facets, ppb_uv )
    determine_facet_overlap( ppb_facets )
    
    # check peeled sources
    if use_catalog[ 0 ]:
      convl_facets = convolve_facets( ppb_facets, resolution = catalog_resolution )
    else:
      convl_facets = ppb_facets
    match_types = []
    offsets = []
    for i in range( 1, 1 + restore_parameter( convl_facets, 'facet_count' ) ):
      facet = get_facet( convl_facets, i )
      facet_ref_pos = get_pixel_reference( facet )
      [ facet_max, facet_max_pos ] = get_image_maximum( facet )
      offset = [ facet_max_pos[ 0 ] - facet_ref_pos[ 0 ],
          facet_max_pos[ 1 ] - facet_ref_pos[ 1 ] ]
      offsets.append( offset )
      match_type = restore_parameter( facet, 'match_type' )
      match_types.append( match_type.split( ', ' ) )
      print '%2d: %6.4f Jy, %4.2f min, %s, %4d, %s' % ( i, get_model_flux( facet ),
          restore_parameter( facet, 'solution_interval' ), match_type, 
          sqrt( offset[ 0 ]**2 + offset[ 1 ]**2 ), repr( get_facet_overlap( facet ) ) )
#      if do_plots:
#        facet = get_facet( ppb_facets, i )
#        plot_image( facet, z_range = plot_range, 
#            plot_title = '%s source %d' % ( target_name, i ) )
#        plot_source( pb_facets, get_radec( facet ), z_range = plot_range )
    if use_catalog[ 0 ]:
      remove_facets( convl_facets )
    
    # select peeled sources for ionospheric modeling
    facet_list = range( 1, 1 + restore_parameter( ppb_facets, 'facet_count' ) )
    if ( len( facet_list ) < 3 ):
      raise error( 'too few peeled sources to fit ionosphere model' )
    for i in range( 1, 1 + restore_parameter( ppb_facets, 'facet_count' ) ):
      facet = get_facet( ppb_facets, i )
      if ( match_types[ i - 1 ][ 1 ] == 'no catalog match' ):
        facet_list.remove( i )
    if ( len( facet_list ) < 3 ):
      raise error( 'too few peeled sources to fit ionosphere model' )
    sel = ( array( facet_list ) - 1 ).reshape( ( len( facet_list ), 1 ) )
    offsets = array( offsets )
    median_offset = median( aget( offsets, sel ), axis = 0 )
    offsets = sqrt( add.reduce( ( offsets - median_offset )**2, axis = 1 ) )
    for i in range( 1, 1 + restore_parameter( ppb_facets, 'facet_count' ) ):
      if ( not i in facet_list ):
        continue
      if ( match_types[ i - 1 ][ 0 ] == 'single source' ):
        if ( offsets[ i - 1 ] > offset_limits[ 0 ] ):
          facet_list.remove( i )
      elif ( match_types[ i - 1 ][ 0 ] == 'double source' ):
        if ( offsets[ i - 1 ] > offset_limits[ 1 ] ):
          facet_list.remove( i )
    if ( len( facet_list ) < 3 ):
      raise error( 'too few peeled sources to fit ionosphere model' )
    for i in range( 1, 1 + restore_parameter( ppb_facets, 'facet_count' ) ):
      if ( not i in facet_list ):
        continue
      facet = get_facet( ppb_facets, i )
      source_facet_list = [ x[ 0 ] for x in 
          find_source_facets( ppb_facets, get_radec( facet ) ) ]
      while ( i != source_facet_list[ 0 ] ):
        if ( source_facet_list[ 0 ] in facet_list ):
          break
        source_facet_list = source_facet_list[ 1 : ]
      if ( i != source_facet_list[ 0 ] ):
        facet_list.remove( i )
    print 'selected facets = ', facet_list
    if ( len( facet_list ) < 3 ):
      raise error( 'too few peeled sources to fit ionosphere model' )
  
  # fit ionospheric model
  order = min( int( round( 1.5 * float( len( facet_list ) ) ) ), 20 )
#  order = min( int( round( 2. * float( len( facet_list ) ) - 5. ) ), 20 )
  sc_version = uv.table_highver( 'SN' )
  fit_ionospheric_pmkl_model( uv, ppb_facets, order = order, print_info = True,
      facet_list = facet_list )
  call_aips_task( 'TACOP', indata = uv, outdata = ppb_uv, inext = 'SN', ncount = 1 )
  uv.zap_table( 'SN', 0 )
#  if do_plots:
#    plot_fit_errors( uv )
  
  # remove instrumental phase effects
  if snapshot:
    window_list = [ 12. * 60., 12. * 60., 12. * 60., 60., 10., 10. ]
  else:
    window_list = [ 12. * 60., 60., 10. ]
  for window in window_list:
    
    # filter out systematic phase effects
    smooth_solutions_in_time( ppb_uv, max_gap = 12. * 60., phase_window = window * 60. )
    re_sample_solutions( ppb_uv, gap_time = 12. * 60. * 60, weight_multiplier = 0.,
        interpolation_method = 'nearest' )
    re_reference_solutions( ppb_uv, reference_antenna )
    for i in range( 1, 1 + restore_parameter( ppb_facets, 'facet_count' ) ):
      facet = get_facet( ppb_facets, i )
      sn_version = facet.table_highver( 'SN' )
      call_aips_task( 'TACOP', indata = ppb_uv, outdata = facet, inext = 'SN', ncount = 1 )
      combine_solutions( facet, invert_2 = True )
    call_aips_task( 'TACOP', indata = ppb_uv, outdata = uv, inext = 'SN', ncount = 1 )
    cor_uv = apply_solution_table( uv )
    combine_solutions( uv, in_version_1 = sc_version, invert_2 = True )
    call_aips_task( 'TACOP', indata = uv, outdata = cor_uv, inext = 'SN', ncount = 1 )
    uv.zap()
    cor_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
    copy_solution_tables( ppb_facets, ppb_uv )
    
    # re-fit ionospheric model
    uv.zap_table( 'NI', 0 )
    uv.zap_table( 'OB', 0 )
    sc_version = uv.table_highver( 'SN' )
    fit_ionospheric_pmkl_model( uv, ppb_facets, order = order, print_info = True,
        facet_list = facet_list )
    call_aips_task( 'TACOP', indata = uv, outdata = ppb_uv, inext = 'SN', ncount = 1 )
    uv.zap_table( 'SN', 0 )
  
  # filter solutions on outliers, and flag bad solutions
#  if do_plots:
#    plot_fit_errors( uv )
  filter_solutions( ppb_uv, phase_sigma = 5., amplitude_sigma = None )
  re_reference_solutions( ppb_uv, reference_antenna )
  flag_bad_solutions( ppb_uv )
  call_aips_task( 'TACOP', indata = ppb_uv, outdata = uv, inext = 'FG', ncount = 1 )
  flag_uv = apply_flag_table( uv )
  uv.zap()
  flag_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  
  # combine model phase solutions and peeling amplitude solutions and delays
  delay_facet_list = []
  model_facet_list = []
  cpb_noise = restore_parameter( uv, 'cpb_noise' )
  for i in range( 1, 1 + restore_parameter( ppb_facets, 'facet_count' ) ):
    facet = get_facet( ppb_facets, i )
    flux = max( get_model_flux( facet ), 0. )
    [ facet_max, facet_max_pos ] = get_image_maximum( facet )
    offset = offsets[ i - 1 ]
    snr =  facet_max * sqrt( flux / facet_max ) / cpb_noise
    print i, flux, facet_max, offset, flux / cpb_noise, facet_max / cpb_noise, snr
    if ( ( i in facet_list ) and ( snr > amplitude_snr ) ):
      delay_facet_list.append( i )
    else:
      model_facet_list.append( i )
  for i in range( 1, 1 + restore_parameter( ppb_facets, 'facet_count' ) ):
    facet = get_facet( ppb_facets, i )
    call_aips_task( 'TACOP', indata = facet, outdata = ppb_uv, inext = 'SN', ncount = 1 )
    snver = ppb_uv.table_highver( 'SN' )
    if ( i in model_facet_list ):
      call_aips_task( 'SNCOR', indata = ppb_uv, opcode = 'ZPHS', snver = snver )
    call_aips_task( 'SNCOR', indata = ppb_uv, opcode = 'ZDEL', snver = snver )
    call_aips_task( 'TACOP', indata = ppb_uv, outdata = facet, inext = 'SN',
        ncount = 1 )
    ppb_uv.zap_table( 'SN', snver )
  if ( len( delay_facet_list ) > 0 ):
    generate_solution_tables_from_pmkl_fit_table( uv, ppb_facets, rms_limit = 45.,
        reference_antenna = reference_antenna, facet_list = delay_facet_list,
        include_phases = False )
  if ( len( model_facet_list ) > 0 ):
    generate_solution_tables_from_pmkl_fit_table( uv, ppb_facets, rms_limit = 45.,
        reference_antenna = reference_antenna, facet_list = model_facet_list )
  for i in range( 1, 1 + restore_parameter( ppb_facets, 'facet_count' ) ):
    facet = get_facet( ppb_facets, i )
    combine_solutions( facet )
  copy_solution_tables( ppb_facets, ppb_uv )
  
  # generate solution tables and re-image
  add_clean_boxes( uv, pb_facets, sidelobe_rejection = sidelobe_rejection,
      clean_box_radius = box_radius )
  if do_large_boxes:
    add_clean_boxes( uv, pb_facets, convolve_factor = 3., keep_boxes = True,
        clean_box_radius = large_box_radius, sidelobe_rejection = sidelobe_rejection )
  add_centered_pb_facets( uv, pb_facets, 1.e6, print_info = True,
      imagr_params = imagr_params, center_facets = ppb_facets, min_separation = 10. )
  generate_solution_tables_from_pmkl_fit_table( uv, pb_facets, rms_limit = 45.,
      reference_antenna = reference_antenna )
  replace_model_solutions_with_peel_solutions( uv, pb_facets, ppb_facets )
  re_image_clean_pb_facets( uv, pb_facets, imagr_params = imagr_params,
      conversion_method = 'GRID' )
  sp1_image = combine_facets( uv, pb_facets, save_tables = False, save_info = True )
  sp1_image.rename( klass = 'SP1' )
  write_fits_image( sp1_image, fits_path + target_name + '.SP1.FITS',
      include_tables = True )
  write_fits_uv( uv, fits_path + target_name + '.SP1.UVFITS' )
  if minimize_storage:
    if file_exists( fits_path + target_name + '.SC3.FITS' ):
      sc3_image.zap()
      remove_file( fits_path + target_name + '.SC3.FITS' )
      remove_file( fits_path + target_name + '.SC3.UVFITS' )
    if file_exists( fits_path + target_name + '.SP0.FITS' ):
      sp0_image.zap()
      remove_file( fits_path + target_name + '.SP0.FITS' )
      remove_file( fits_path + target_name + '.SP0.UVFITS' )
    ppb_uv.zap()
    remove_facets( ppb_facets )
  
  # calibrate bandpass phases on target field
  model_uv = make_model_uv( uv, pb_facets, apply_solutions = True )
  div_uv = divide_uv( uv, model_uv )
  call_aips_task( 'INDXR', indata = div_uv, cparm = [ 12. * 60.,20. * 60.,-1,0,0,0 ] )
  call_aips_task( 'BPASS', indata = div_uv, flagver = 0, dopol = -1,
      blver = -1, docalib = -1, gainuse = 0, doband = -1, solint = 0.,
      soltype = 'L1', outvers = 0, bpassprm = [ 0,0,1,0,2,0,0,1,1,3,1 ],
      specindx = 1.e-3, refant = reference_antenna, cmethod = 'DFT',
      **calib_params )
  call_aips_task( 'TACOP', indata = div_uv, outdata = uv, inext = 'BP', ncount = 1 )
  div_uv.zap()
  temp_uv = get_aips_file( aips_disk, target_name, 'TEMP', -1, 'UV' )
  call_aips_task( 'SPLAT', indata = uv, smooth = [ 0 ], douvcomp = 0,
      aparm = [ 3,0,1,0,0,0,1 ], channel = 1, chinc = 1, docalib = -1,
      doband = 4, bpver = 0, flagver = 0, outdata = temp_uv )
  uv.zap()
  temp_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  
  # determine amplitude calibration
  if try_new:
    # exclude facets with amplitude corrections
    no_amp_facet_list = range( 1, 1 + restore_parameter( pb_facets, 'facet_count' ) )
    amp_facet_list = []
    for i in range( 1, 1 + restore_parameter( pb_facets, 'facet_count' ) ):
      facet = get_facet( pb_facets, i )
      if has_amplitude_solutions( facet ):
        no_amp_facet_list.remove( i )
        amp_facet_list.append( i )
    if ( len( amp_facet_list ) > 0 ):
      model_uv.zap()
      model_uv = make_model_uv( uv, pb_facets, facet_list = no_amp_facet_list,
          apply_solutions = True, keep_flags = False )
      sub_uv = subtract_model( uv, pb_facets, facet_list = amp_facet_list,
          keep_flags = False )
    else:
      sub_uv = uv
    calibrate_uv( sub_uv, model_uv, reference_antenna = reference_antenna,
        calib_params = calib_params, do_amplitude = True, apply_flags = False,
        amplitude_interval = phase_interval_min )
    model_uv.zap()
    if ( len( amp_facet_list ) > 0 ):
      call_aips_task( 'TACOP', indata = sub_uv, outdata = uv, inext = 'SN', ncount = 1 )
      sub_uv.zap()
  else:
    calibrate_uv( uv, model_uv, reference_antenna = reference_antenna,
        calib_params = calib_params, do_amplitude = True, apply_flags = False,
        amplitude_interval = phase_interval_min )
    model_uv.zap()
  
  # filter amplitude solutions on outliers
  filter_solutions( uv, amplitude_sigma = 5., phase_sigma = 5.,
      amplitude_window = 12. * 60. * 60., phase_window = 12. * 60. * 60. )
  filter_solutions( uv, amplitude_sigma = 5., phase_sigma = 5.,
      amplitude_window = 60. * 60., phase_window = 12. * 60. * 60. )
  smooth_solutions_in_time( uv, phase_window = 0., amplitude_window = 60., order = 1 )
  sn_version = uv.table_highver( 'SN' )
  call_aips_task( 'SNCOR', indata = uv, opcode = 'ZPHS', snver = sn_version )
  if try_new:
    # remove amplitude corrections from factes with previous amplitude corrections
    if ( len( amp_facet_list ) > 0 ):
      for i in amp_facet_list:
        facet = get_facet( pb_facets, i )
        call_aips_task( 'TACOP', indata = uv, outdata = facet, inext = 'SN', ncount = 1 )
        combine_solutions( facet, invert_2 = True )
  reference_antenna = select_reference_antenna( uv, count_power = 2.,
      antenna_list = reference_antenna_list )
  re_reference_solutions( uv, reference_antenna )
  cor_uv = apply_solution_table( uv )
  uv.zap()
  cor_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  
  # calibrate and re-image
#  calibrate_model( uv, pb_facets, reference_antenna = reference_antenna,
#      calib_params = calib_params, phase_interval = phase_interval_min,
#      apply_solutions = False )
  add_clean_boxes( uv, pb_facets, sidelobe_rejection = sidelobe_rejection,
      clean_box_radius = box_radius )
  if do_large_boxes:
    add_clean_boxes( uv, pb_facets, convolve_factor = 3., keep_boxes = True,
        clean_box_radius = large_box_radius, sidelobe_rejection = sidelobe_rejection )
  re_image_clean_pb_facets( uv, pb_facets, imagr_params = imagr_params,
      conversion_method = 'GRID' )
  sp1a_image = combine_facets( uv, pb_facets, save_tables = False, save_info = True )
  sp1a_image.rename( klass = 'SP1A' )
  write_fits_image( sp1a_image, fits_path + target_name + '.SP1A.FITS',
      include_tables = True )
  write_fits_uv( uv, fits_path + target_name + '.SP1A.UVFITS' )
  if minimize_storage:
    sp1_image.zap()
    remove_file( fits_path + target_name + '.SP1.FITS' )
    remove_file( fits_path + target_name + '.SP1.UVFITS' )
  
  # flag image undulations
  flag_image_undulations( uv, sp1a_image, imagr_params,
      keep_images = ( not minimize_storage ) )
  cpb_noise = restore_parameter( uv, 'cpb_noise' )
  measure_cpb_noise( uv, pb_facets, keep_image = True )
  cpb_image = get_aips_file( aips_disk, 'CPB', 'FLATN', 0, 'MA' )
  flag_image_undulations( uv, cpb_image, imagr_params,
      keep_images = ( not minimize_storage ) )
  if minimize_storage:
    cpb_image.zap()
  store_parameter( uv, 'cpb_noise', cpb_noise )
  for i in range( 1, 1 + restore_parameter( pb_facets, 'facet_count' ) ):
    facet = get_facet( pb_facets, i )
    flux = max( get_model_flux( facet ), 0. )
    [ facet_max, facet_max_pos ] = get_image_maximum( facet )
    facet_ref_pos = get_pixel_reference( facet )
    offset = int( round( sqrt( ( facet_max_pos[ 0 ] - facet_ref_pos[ 0 ] )**2 + 
        ( facet_max_pos[ 1 ] - facet_ref_pos[ 1 ] )**2 ) ) )
    snr =  facet_max * sqrt( flux / facet_max ) / cpb_noise
    if ( ( i == 1 ) or ( offset < 10 ) ):
      print i, flux, facet_max, offset, flux / cpb_noise, facet_max / cpb_noise, snr
      if ( ( i == 1 ) or ( snr > amplitude_snr ) ):
        flag_image_undulations( uv, facet, imagr_params,
            keep_images = ( not minimize_storage ) )
  flag_uv = apply_flag_table( uv )
  uv.zap()
  flag_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  
  # subtract pb sources and flag
  sub_pb_uv = subtract_model( uv, pb_facets )
  flag_uv = apply_flag_table( uv )
  uv.zap()
  flag_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  flag_image_baselines( sub_pb_uv, keep_images = ( not minimize_storage ),
      scaling_factor = 32., imagr_params = imagr_params, flag_version = -1 )
  flag_image_baselines( sub_pb_uv, keep_images = ( not minimize_storage ),
      scaling_factor = 16., imagr_params = imagr_params )
  flag_image_baselines( sub_pb_uv, keep_images = ( not minimize_storage ),
      scaling_factor = 8., imagr_params = imagr_params )
  flag_image_baselines( sub_pb_uv, keep_images = ( not minimize_storage ),
      scaling_factor = 4., imagr_params = imagr_params )
  flag_image_baselines( sub_pb_uv, keep_images = ( not minimize_storage ),
      scaling_factor = 2., imagr_params = imagr_params )
  flag_image_baselines( sub_pb_uv, keep_images = ( not minimize_storage ),
      scaling_factor = 1., imagr_params = imagr_params )
  extend_flags( sub_pb_uv )
  call_aips_task( 'TACOP', indata = sub_pb_uv, outdata = uv, inext = 'FG',
      ncount = 1 )
  flag_uv = apply_flag_table( uv )
  uv.zap()
  flag_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  flag_uv = apply_flag_table( sub_pb_uv )
  sub_pb_uv.zap()
  flag_uv.rename( name = sub_pb_uv.name, klass = sub_pb_uv.klass, seq = sub_pb_uv.seq )
  
  # flag based on residual amplitudes
#  total_clean_flux = restore_parameter( sp1a_image, 'total_clean_flux' )
#  if ( 
  call_aips_task( 'FLGIT', indata = sub_pb_uv, outdata = sub_pb_uv,
    flagver = 1, outfgver = 1, order = 0, ichansel = [ [ 1,0,1,1 ] ], opcode = '',
    aparm = [ clip_levels[ 1 ], 0,0,2.,2.,0, get_channel_count( sub_pb_uv ), 0,0,0 ] )
  call_aips_task( 'FLGIT', indata = sub_pb_uv, outdata = sub_pb_uv,
    flagver = 1, outfgver = 1, order = 0, ichansel = [ [ 1,0,1,1 ] ], opcode = '',
    aparm = [ clip_levels[ 1 ] ,0,0,2.,2.,0, get_channel_count( sub_pb_uv ), 0,0,0 ] )
  call_aips_task( 'TACOP', indata = sub_pb_uv, outdata = uv, inext = 'FG',
      ncount = 1 )
  flag_uv = apply_flag_table( uv )
  uv.zap()
  flag_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  sub_pb_uv.zap()
  
  # calibrate and re-image
#  calibrate_model( uv, pb_facets, reference_antenna = reference_antenna,
#      calib_params = calib_params, phase_interval = phase_interval_min,
#      apply_solutions = False )
  add_clean_boxes( uv, pb_facets, sidelobe_rejection = sidelobe_rejection,
      clean_box_radius = box_radius )
  if do_large_boxes:
    add_clean_boxes( uv, pb_facets, convolve_factor = 3., keep_boxes = True,
        clean_box_radius = large_box_radius, sidelobe_rejection = sidelobe_rejection )
  re_image_clean_pb_facets( uv, pb_facets, imagr_params = imagr_params,
      conversion_method = 'GRID' )
  sp1b_image = combine_facets( uv, pb_facets, save_tables = False, save_info = True )
  sp1b_image.rename( klass = 'SP1B' )
  write_fits_image( sp1b_image, fits_path + target_name + '.SP1B.FITS',
      include_tables = True )
  write_fits_uv( uv, fits_path + target_name + '.SP1B.UVFITS' )
  if minimize_storage:
    sp1a_image.zap()
    remove_file( fits_path + target_name + '.SP1A.FITS' )
    remove_file( fits_path + target_name + '.SP1A.UVFITS' )
  
  # peel sources in PB area
  sub_pb_uv = subtract_model( uv, pb_facets )
  flag_uv = apply_flag_table( uv )
  uv.zap()
  flag_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  if use_catalog[ 1 ]:
    convl_facets = convolve_facets( pb_facets, resolution = catalog_resolution )
    beam_size = get_beam_size( pb_facets )
    search_ratio = ( 10. * sqrt( beam_size[ 0 ]**2 + beam_size[ 1 ]**2 ) / 
        ( 2. * catalog_resolution ) )
  else:
    convl_facets = pb_facets
    search_ratio = 10.
  [ ppb_uv, ppb_facets ] = peel_pbo_facets( sub_pb_uv, convl_facets, 
      phase_interval_min = phase_interval_min, reference_antenna = reference_antenna, 
      imagr_params = imagr_params, calib_params = calib_params, max_relativity = 3,
      sidelobe_rejection = sidelobe_rejection, conversion_method = 'GRID',
      amplitude_snr = amplitude_snr, double_search_ratio = search_ratio,
      use_catalog = use_catalog[ 1 ], catalog_resolution = catalog_resolution,
      use_nvss = use_nvss, input_catalog = model_source_list,
      clean_box_radius = box_radius )
  if use_catalog[ 1 ]:
    remove_facets( convl_facets )
  sub_pb_uv.zap()
  copy_solution_tables( ppb_facets, ppb_uv )
  determine_facet_overlap( ppb_facets )
  
  # check peeled sources
  if use_catalog[ 1 ]:
    convl_facets = convolve_facets( ppb_facets, resolution = catalog_resolution )
  else:
    convl_facets = ppb_facets
  match_types = []
  offsets = []
  for i in range( 1, 1 + restore_parameter( convl_facets, 'facet_count' ) ):
    facet = get_facet( convl_facets, i )
    facet_ref_pos = get_pixel_reference( facet )
    [ facet_max, facet_max_pos ] = get_image_maximum( facet )
    offset = [ facet_max_pos[ 0 ] - facet_ref_pos[ 0 ],
        facet_max_pos[ 1 ] - facet_ref_pos[ 1 ] ]
    offsets.append( offset )
    match_type = restore_parameter( facet, 'match_type' )
    match_types.append( match_type.split( ', ' ) )
    print '%2d: %6.4f Jy, %4.2f min, %s, %4d, %s' % ( i, get_model_flux( facet ),
        restore_parameter( facet, 'solution_interval' ), match_type, 
        sqrt( offset[ 0 ]**2 + offset[ 1 ]**2 ), repr( get_facet_overlap( facet ) ) )
#    if do_plots:
#      facet = get_facet( ppb_facets, i )
#      plot_image( facet, z_range = plot_range, 
#          plot_title = '%s source %d' % ( target_name, i ) )
#      plot_source( pb_facets, get_radec( facet ), z_range = plot_range )
  if use_catalog[ 1 ]:
    remove_facets( convl_facets )
  
  # select peeled sources for ionospheric modeling
  if use_catalog[ 1 ]:
    offset_limits = [ 1., 2. ]
  else:
    offset_limits = [ 1.5, 3. ]
  facet_list = range( 1, 1 + restore_parameter( ppb_facets, 'facet_count' ) )
  if ( len( facet_list ) < 3 ):
    raise error( 'too few peeled sources to fit ionosphere model' )
  for i in range( 1, 1 + restore_parameter( ppb_facets, 'facet_count' ) ):
    facet = get_facet( ppb_facets, i )
    if ( match_types[ i - 1 ][ 1 ] == 'no catalog match' ):
      facet_list.remove( i )
  if ( len( facet_list ) < 3 ):
    raise error( 'too few peeled sources to fit ionosphere model' )
  sel = ( array( facet_list ) - 1 ).reshape( ( len( facet_list ), 1 ) )
  offsets = array( offsets )
  median_offset = median( aget( offsets, sel ), axis = 0 )
  offsets = sqrt( add.reduce( ( offsets - median_offset )**2, axis = 1 ) )
  for i in range( 1, 1 + restore_parameter( ppb_facets, 'facet_count' ) ):
    if ( not i in facet_list ):
      continue
    if ( match_types[ i - 1 ][ 0 ] == 'single source' ):
      if ( offsets[ i - 1 ] > offset_limits[ 0 ] ):
        facet_list.remove( i )
    elif ( match_types[ i - 1 ][ 0 ] == 'double source' ):
      if ( offsets[ i - 1 ] > offset_limits[ 1 ] ):
        facet_list.remove( i )
  if ( len( facet_list ) < 3 ):
    raise error( 'too few peeled sources to fit ionosphere model' )
  for i in range( 1, 1 + restore_parameter( ppb_facets, 'facet_count' ) ):
    if ( not i in facet_list ):
      continue
    facet = get_facet( ppb_facets, i )
    source_facet_list = [ x[ 0 ] for x in 
        find_source_facets( ppb_facets, get_radec( facet ) ) ]
    while ( i != source_facet_list[ 0 ] ):
      if ( source_facet_list[ 0 ] in facet_list ):
        break
      source_facet_list = source_facet_list[ 1 : ]
    if ( i != source_facet_list[ 0 ] ):
      facet_list.remove( i )
  print 'selected facets = ', facet_list
  if ( len( facet_list ) < 3 ):
    raise error( 'too few peeled sources to fit ionosphere model' )
  
  # fit ionospheric model
  order = min( int( round( 1.5 * float( len( facet_list ) ) ) ), 20 )
#  order = min( int( round( 2. * float( len( facet_list ) ) - 5. ) ), 20 )
  sc_version = uv.table_highver( 'SN' )
  fit_ionospheric_pmkl_model( uv, ppb_facets, order = order, print_info = True,
      facet_list = facet_list )
  call_aips_task( 'TACOP', indata = uv, outdata = ppb_uv, inext = 'SN', ncount = 1 )
  uv.zap_table( 'SN', 0 )
  
  # remove instrumental phase effects
  window = 10.
  smooth_solutions_in_time( ppb_uv, max_gap = 12. * 60., phase_window = window * 60. )
  re_sample_solutions( ppb_uv, gap_time = 12. * 60. * 60, weight_multiplier = 0.,
      interpolation_method = 'nearest' )
  re_reference_solutions( ppb_uv, reference_antenna )
  for i in range( 1, 1 + restore_parameter( ppb_facets, 'facet_count' ) ):
    facet = get_facet( ppb_facets, i )
    sn_version = facet.table_highver( 'SN' )
    call_aips_task( 'TACOP', indata = ppb_uv, outdata = facet, inext = 'SN', ncount = 1 )
    combine_solutions( facet, invert_2 = True )
  call_aips_task( 'TACOP', indata = ppb_uv, outdata = uv, inext = 'SN', ncount = 1 )
  cor_uv = apply_solution_table( uv )
#  combine_solutions( uv, in_version_1 = sc_version, invert_2 = True )
#  call_aips_task( 'TACOP', indata = uv, outdata = cor_uv, inext = 'SN', ncount = 1 )
  uv.zap()
  cor_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  copy_solution_tables( ppb_facets, ppb_uv )
  
  # re-fit ionospheric model
  uv.zap_table( 'NI', 0 )
  uv.zap_table( 'OB', 0 )
  sc_version = uv.table_highver( 'SN' )
  fit_ionospheric_pmkl_model( uv, ppb_facets, order = order, print_info = True,
      facet_list = facet_list )
  call_aips_task( 'TACOP', indata = uv, outdata = ppb_uv, inext = 'SN', ncount = 1 )
  uv.zap_table( 'SN', 0 )
  
  # filter solutions on outliers, and apply
#  if do_plots:
#    plot_fit_errors( uv )
  filter_solutions( ppb_uv, phase_sigma = 4., amplitude_sigma = None )
  re_reference_solutions( ppb_uv, reference_antenna )
  for i in range( 1, 1 + restore_parameter( ppb_facets, 'facet_count' ) ):
    facet = get_facet( ppb_facets, i )
    call_aips_task( 'TACOP', indata = ppb_uv, outdata = facet, inext = 'SN', ncount = 1 )
    combine_solutions( facet, invert_2 = True )
  call_aips_task( 'TACOP', indata = ppb_uv, outdata = uv, inext = 'SN', ncount = 1 )
  cal_uv = apply_solution_table( uv )
#  combine_solutions( uv, invert_2 = True )
#  call_aips_task( 'TACOP', indata = uv, outdata = cal_uv, inext = 'SN', ncount = 1 )
  uv.zap()
  cal_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  
  # combine model phase solutions and peeling amplitude solutions and delays
  delay_facet_list = []
  model_facet_list = []
  cpb_noise = restore_parameter( uv, 'cpb_noise' )
  for i in range( 1, 1 + restore_parameter( ppb_facets, 'facet_count' ) ):
    facet = get_facet( ppb_facets, i )
    flux = get_model_flux( facet )
    [ facet_max, facet_max_pos ] = get_image_maximum( facet )
    offset = offsets[ i - 1 ]
    snr =  facet_max * sqrt( flux / facet_max ) / cpb_noise
    print i, flux, facet_max, offset, flux / cpb_noise, facet_max / cpb_noise, snr
    if ( ( i in facet_list ) and ( snr > amplitude_snr ) ):
      delay_facet_list.append( i )
    else:
      model_facet_list.append( i )
  for i in range( 1, 1 + restore_parameter( ppb_facets, 'facet_count' ) ):
    facet = get_facet( ppb_facets, i )
    call_aips_task( 'TACOP', indata = facet, outdata = ppb_uv, inext = 'SN', ncount = 1 )
    snver = ppb_uv.table_highver( 'SN' )
    if ( i in model_facet_list ):
      call_aips_task( 'SNCOR', indata = ppb_uv, opcode = 'ZPHS', snver = snver )
    call_aips_task( 'SNCOR', indata = ppb_uv, opcode = 'ZDEL', snver = snver )
    call_aips_task( 'TACOP', indata = ppb_uv, outdata = facet, inext = 'SN',
        ncount = 1 )
    ppb_uv.zap_table( 'SN', snver )
  if ( len( delay_facet_list ) > 0 ):
    # Note: the delay phase offset from center to ref frequency is not added again
    generate_solution_tables_from_pmkl_fit_table( uv, ppb_facets,
        reference_antenna = reference_antenna, facet_list = delay_facet_list,
        include_phases = False, include_delay_phases = False )
  if ( len( model_facet_list ) > 0 ):
    generate_solution_tables_from_pmkl_fit_table( uv, ppb_facets,
        reference_antenna = reference_antenna, facet_list = model_facet_list )
  for i in range( 1, 1 + restore_parameter( ppb_facets, 'facet_count' ) ):
    facet = get_facet( ppb_facets, i )
    combine_solutions( facet )
  copy_solution_tables( ppb_facets, ppb_uv )
  
  # generate solution tables and re-image
  add_clean_boxes( uv, pb_facets, sidelobe_rejection = sidelobe_rejection,
      clean_box_radius = box_radius )
  if do_large_boxes:
    add_clean_boxes( uv, pb_facets, convolve_factor = 3., keep_boxes = True,
        clean_box_radius = large_box_radius, sidelobe_rejection = sidelobe_rejection )
  add_centered_pb_facets( uv, pb_facets, 1.e6, center_facets = ppb_facets,
       imagr_params = imagr_params, min_separation = 10. )
  generate_solution_tables_from_pmkl_fit_table( uv, pb_facets, 
      reference_antenna = reference_antenna )
  replace_model_solutions_with_peel_solutions( uv, pb_facets, ppb_facets )
  re_image_clean_pb_facets( uv, pb_facets, imagr_params = imagr_params,
      conversion_method = 'GRID' )
  sp2_image = combine_facets( uv, pb_facets, save_tables = False, save_info = True )
  sp2_image.rename( klass = 'SP2' )
  write_fits_image( sp2_image, fits_path + target_name + '.SP2.FITS',
      include_tables = True )
  write_fits_uv( uv, fits_path + target_name + '.SP2.UVFITS' )
  if minimize_storage:
    sp1b_image.zap()
    remove_file( fits_path + target_name + '.SP1B.FITS' )
    remove_file( fits_path + target_name + '.SP1B.UVFITS' )
    ppb_uv.zap()
    remove_facets( ppb_facets )
  
  # calibrate bandpass phases on target field
  model_uv = make_model_uv( uv, pb_facets, apply_solutions = True )
  div_uv = divide_uv( uv, model_uv )
  call_aips_task( 'INDXR', indata = div_uv, cparm = [ 12. * 60.,20. * 60.,-1,0,0,0 ] )
  call_aips_task( 'BPASS', indata = div_uv, flagver = 0, dopol = -1,
      blver = -1, docalib = -1, gainuse = 0, doband = -1, solint = 0.,
      soltype = 'L1', outvers = 0, bpassprm = [ 0,0,1,0,2,0,0,1,1,3,1 ],
      specindx = 1.e-3, refant = reference_antenna, cmethod = 'DFT',
      **calib_params )
  call_aips_task( 'TACOP', indata = div_uv, outdata = uv, inext = 'BP', ncount = 1 )
  div_uv.zap()
  temp_uv = get_aips_file( aips_disk, target_name, 'TEMP', -1, 'UV' )
  call_aips_task( 'SPLAT', indata = uv, smooth = [ 0 ], douvcomp = 0,
      aparm = [ 3,0,1,0,0,0,1 ], channel = 1, chinc = 1, docalib = -1,
      doband = 4, bpver = 0, flagver = 0, outdata = temp_uv )
  uv.zap()
  temp_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  
  # calibrate baselines on target field
  div_uv = divide_uv( uv, model_uv )
  call_aips_task( 'INDXR', indata = div_uv, cparm = [ 12. * 60.,20. * 60.,-1,0,0,0 ] )
  call_aips_task( 'BLCAL', indata = div_uv, flagver = 0, dopol = -1, bparm = [ 1,0 ],
      blver = 0, docalib = -1, gainuse = 0, doband = -1,  solint = 0,
      cmethod = 'DFT', ichansel = [ [ 1, get_channel_count( uv ), 1,1 ] ] )
  call_aips_task( 'TACOP', indata = div_uv, outdata = uv, inext = 'BL', ncount = 1 )
  div_uv.zap()
  bl_version = uv.table_highver( 'BL' )
  temp_uv = get_aips_file( aips_disk, target_name, 'TEMP', -1, 'UV' )
  call_aips_task( 'SPLAT', indata = uv, smooth = [ 0 ], douvcomp = 0,
      aparm = [ 3,0,1,0,0,0,1 ], channel = 1, chinc = 1, docalib = -1,
      doband = -1, flagver = 0, blver = bl_version, outdata = temp_uv )
  call_aips_task( 'TACOP', indata = uv, outdata = temp_uv, inext = 'SN', ncount = 1 )
  uv.zap()
  temp_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  
  # determine amplitude calibration
  # filter solutions on outliers
  calibrate_uv( uv, model_uv, reference_antenna = reference_antenna,
      calib_params = calib_params, do_amplitude = True, apply_flags = False,
      amplitude_interval = phase_interval_min )
  model_uv.zap()
  filter_solutions( uv, amplitude_sigma = 5., phase_sigma = 5.,
      amplitude_window = 12. * 60. * 60., phase_window = 12. * 60. * 60. )
  filter_solutions( uv, amplitude_sigma = 4., phase_sigma = 4.,
      amplitude_window = 60. * 60., phase_window = 60. * 60. )
  reference_antenna = select_reference_antenna( uv, count_power = 2.,
      antenna_list = reference_antenna_list )
  re_reference_solutions( uv, reference_antenna )
  smooth_solutions_in_time( uv, phase_window = 0., amplitude_window = 60., order = 1 )
  sn_version = uv.table_highver( 'SN' )
  call_aips_task( 'SNCOR', indata = uv, opcode = 'ZPHS', snver = sn_version )
  cor_uv = apply_solution_table( uv )
  uv.zap()
  cor_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  
  # calibrate and re-image
#  calibrate_model( uv, pb_facets, reference_antenna = reference_antenna,
#      calib_params = calib_params, phase_interval = phase_interval_min,
#      apply_solutions = False )
  add_clean_boxes( uv, pb_facets, sidelobe_rejection = sidelobe_rejection,
      clean_box_radius = box_radius )
  if do_large_boxes:
    add_clean_boxes( uv, pb_facets, convolve_factor = 3., keep_boxes = True,
        clean_box_radius = large_box_radius, sidelobe_rejection = sidelobe_rejection )
  re_image_clean_pb_facets( uv, pb_facets, imagr_params = imagr_params,
      conversion_method = 'GRID' )
  sp2a_image = combine_facets( uv, pb_facets, save_tables = False, save_info = True )
  sp2a_image.rename( klass = 'SP2A' )
  write_fits_image( sp2a_image, fits_path + target_name + '.SP2A.FITS',
      include_tables = True )
  write_fits_uv( uv, fits_path + target_name + '.SP2A.UVFITS' )
  if minimize_storage:
    sp2_image.zap()
    remove_file( fits_path + target_name + '.SP2.FITS' )
    remove_file( fits_path + target_name + '.SP2.UVFITS' )
  
  # flag image undulations
  flag_image_undulations( uv, sp2a_image, imagr_params,
      keep_images = ( not minimize_storage ) )
  cpb_noise = restore_parameter( uv, 'cpb_noise' )
  measure_cpb_noise( uv, pb_facets, keep_image = True )
  cpb_image = get_aips_file( aips_disk, 'CPB', 'FLATN', 0, 'MA' )
  flag_image_undulations( uv, cpb_image, imagr_params,
      keep_images = ( not minimize_storage ) )
  if minimize_storage:
    cpb_image.zap()
  store_parameter( uv, 'cpb_noise', cpb_noise )
  for i in range( 1, 1 + restore_parameter( pb_facets, 'facet_count' ) ):
    facet = get_facet( pb_facets, i )
    flux = get_model_flux( facet )
    [ facet_max, facet_max_pos ] = get_image_maximum( facet )
    facet_ref_pos = get_pixel_reference( facet )
    offset = int( round( sqrt( ( facet_max_pos[ 0 ] - facet_ref_pos[ 0 ] )**2 + 
        ( facet_max_pos[ 1 ] - facet_ref_pos[ 1 ] )**2 ) ) )
    snr =  facet_max * sqrt( flux / facet_max ) / cpb_noise
    if ( ( i == 1 ) or ( offset < 10 ) ):
      print i, flux, facet_max, offset, flux / cpb_noise, facet_max / cpb_noise, snr
      if ( ( i == 1 ) or ( snr > amplitude_snr ) ):
        flag_image_undulations( uv, facet, imagr_params,
            keep_images = ( not minimize_storage ) )
  flag_uv = apply_flag_table( uv )
  uv.zap()
  flag_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  
  # subtract pb sources and flag
  sub_pb_uv = subtract_model( uv, pb_facets )
  flag_uv = apply_flag_table( uv )
  uv.zap()
  flag_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  flag_image_baselines( sub_pb_uv, keep_images = ( not minimize_storage ),
      scaling_factor = 32., imagr_params = imagr_params, flag_version = -1 )
  flag_image_baselines( sub_pb_uv, keep_images = ( not minimize_storage ),
      scaling_factor = 16., imagr_params = imagr_params )
  flag_image_baselines( sub_pb_uv, keep_images = ( not minimize_storage ),
      scaling_factor = 8., imagr_params = imagr_params )
  flag_image_baselines( sub_pb_uv, keep_images = ( not minimize_storage ),
      scaling_factor = 4., imagr_params = imagr_params )
  flag_image_baselines( sub_pb_uv, keep_images = ( not minimize_storage ),
      scaling_factor = 2., imagr_params = imagr_params )
  flag_image_baselines( sub_pb_uv, keep_images = ( not minimize_storage ),
      scaling_factor = 1., imagr_params = imagr_params )
  extend_flags( sub_pb_uv )
  call_aips_task( 'TACOP', indata = sub_pb_uv, outdata = uv, inext = 'FG',
      ncount = 1 )
  flag_uv = apply_flag_table( uv )
  uv.zap()
  flag_uv.rename( name = uv.name, klass = uv.klass, seq = uv.seq )
  sub_pb_uv.zap()
  
  # calibrate and re-image
#  calibrate_model( uv, pb_facets, reference_antenna = reference_antenna,
#      calib_params = calib_params, phase_interval = phase_interval_min,
#      apply_solutions = False )
  add_clean_boxes( uv, pb_facets, sidelobe_rejection = sidelobe_rejection,
      clean_box_radius = box_radius )
  if do_large_boxes:
    add_clean_boxes( uv, pb_facets, convolve_factor = 3., keep_boxes = True,
        clean_box_radius = large_box_radius, sidelobe_rejection = sidelobe_rejection )
  re_image_clean_pb_facets( uv, pb_facets, imagr_params = imagr_params,
      conversion_method = 'GRID' )
  sp2b_image = combine_facets( uv, pb_facets, save_tables = False, save_info = True )
  sp2b_image.rename( klass = 'SP2B' )
  write_fits_image( sp2b_image, fits_path + target_name + '.SP2B.FITS',
      include_tables = True )
  write_fits_uv( uv, fits_path + target_name + '.SP2B.UVFITS' )
  if minimize_storage:
    sp2a_image.zap()
    remove_file( fits_path + target_name + '.SP2A.FITS' )
    remove_file( fits_path + target_name + '.SP2A.UVFITS' )
  
  ### do post-processing
  
  # make global astrometric correction
  source_flux_min = 10. * restore_parameter( uv, 'cpb_noise' )
  source_list = get_source_list_from_facets( pb_facets, source_flux_min )
  radec_list = []
  for x in source_list:
    facet = get_facet( pb_facets, x[ 0 ] )
    radec_list.append( calculate_source_radec( facet, x[ 1 ] ) )
  print 'detected %d sources for astrometry check' % ( len( radec_list ) )
  print 'using %d catalog sources for astrometry check' % ( len( model_radec_list ) )
  global_offset = None
  for beam_factor in [ 2., 6., 18., 54. ]:
    try:
      global_offset = check_astrometry( sp2b_image, radec_list, beam_factor = beam_factor,
        ref_radec_list = model_radec_list, plot_offsets = False, min_points = 8,
        use_nvss = use_nvss )
    except RuntimeError:
      continue
    except:
      raise error( traceback.print_exc() )
    else:
      [ dra_median, ddec_median, dr_std ] = global_offset
      if ( dr_std > catalog_resolution / 2. ):
        continue
      break
  if ( not global_offset is None ):
    offset = [ -dra_median / 3600., -ddec_median / 3600. ]
    pb_facet_file_name = restore_parameter( uv, 'pb_facet_file_name' )
    pb_facet_file_name_e = path.expandvars( pb_facet_file_name )
    for i in range( 1, 1 + restore_parameter( pb_facets, 'facet_count' ) ):
      facet = get_facet( pb_facets, i )
      cor_radec = calculate_offset_position_dradec( get_radec( facet ), offset )
      replace_facet( pb_facet_file_name_e, i, get_image_size( facet ), cor_radec,
          keep_boxes = True )
      set_radec( facet, cor_radec )
    cor_radec = calculate_offset_position_dradec( get_radec( uv ), offset )
    set_radec( uv, cor_radec )
    cor_radec = calculate_offset_position_dradec( get_radec( sp2b_image ), offset )
    set_radec( sp2b_image, cor_radec )
    write_fits_image( sp2b_image, fits_path + target_name + '.SP2B.FITS',
        include_tables = True )
    write_fits_uv( uv, fits_path + target_name + '.SP2B.UVFITS' )
  
  # create primary-beam corrected image
  sp2b_pbcor_image = apply_pb_correction( uv, sp2b_image )
  write_fits_image( sp2b_pbcor_image, fits_path + target_name + '.SP2B.PBCOR.FITS',
      include_tables = False )
  if minimize_storage:
    sp2b_pbcor_image.zap()
  
  # create flagged & calibrated residual UV data
  sub_pb_uv = subtract_model( uv, pb_facets )
  while table_exists( sub_pb_uv, 'FG', 0 ):
    flag_uv = apply_flag_table( sub_pb_uv )
    sub_pb_uv.zap()
    flag_uv.rename( name = sub_pb_uv.name, klass = sub_pb_uv.klass, seq = sub_pb_uv.seq )
  if table_exists( pb_facets, 'SN', 0 ):
    call_aips_task( 'TACOP', indata = pb_facets, outdata = sub_pb_uv, inext = 'SN',
        ncount = 1 )
  else:
    if ( not table_exists( uv, 'SN', 0 ) ):
      calibrate_model( uv, pb_facets, reference_antenna = reference_antenna,
          calib_params = calib_params, phase_interval = phase_interval_min,
          apply_solutions = False )
    call_aips_task( 'TACOP', indata = uv, outdata = sub_pb_uv, inext = 'SN',
        ncount = 1 )
  if table_exists( sub_pb_uv, 'SN', 0 ):
    cal_uv = apply_solution_table( sub_pb_uv )
    sub_pb_uv.zap()
    cal_uv.rename( name = sub_pb_uv.name, klass = sub_pb_uv.klass, seq = sub_pb_uv.seq )
  for table in [ 'SN', 'NI', 'OB' ]:
    while table_exists( sub_pb_uv, table, 0 ):
      sub_pb_uv.zap_table( table, 0 )
  write_fits_uv( sub_pb_uv, fits_path + target_name + '.SP2B.RES.UVFITS' )
  
  # create flagged & calibrated UV data with outliers removed
  radec = get_radec( uv )
  radius = restore_parameter( uv, 'field_size' ) * 60. * 0.6
  add_uv = add_model_in_circle( sub_pb_uv, pb_facets, radec, radius,
      apply_solutions = False )
  write_fits_uv( add_uv, fits_path + target_name + '.SP2B.CAL.UVFITS' )
  if minimize_storage:
    sub_pb_uv.zap()
    add_uv.zap()
  
  # create residual image
  rst_facets = restore_model_components( pb_facets, imagr_params = imagr_params,
      subtract = True )
  sp2br_image = combine_facets( uv, rst_facets, save_tables = False, save_info = True )
  sp2br_image.rename( name = pb_facets.name, klass = 'SP2BR' )
  while table_exists( sp2br_image, 'CC', 0 ):
    sp2br_image.zap_table( 'CC', 0 )
  call_aips_task( 'TACOP', indata = sp2b_image, outdata = sp2br_image, inext = 'CC',
      ncount = 1 )
  write_fits_image( sp2br_image, fits_path + target_name + '.SP2B.RES.FITS',
      include_tables = True )
  remove_facets( rst_facets )
  if minimize_storage:
    sp2br_image.zap()
  
  # print info
  cpb_noise =  restore_parameter( sp2b_image, 'background_rms' )
  beam = get_beam_size( sp2b_image )
  print '%.3f mJy/beam noise for a %.1F" x %.1F" (PA %d deg) beam' % (
      cpb_noise * 1.e3, beam[ 0 ], beam[ 1 ], beam[ 2 ] )
  
  return

###############################################################################

def summarize_spam_log( log_file_name ):
  
  # read info from log file
  host = '???'
  status = '???'
  flux = 0.
  dflux = 0.
  cc_count = 0
  dcc_count = 0
  first = True
  log_file_names = glob.glob( log_file_name )
  if ( len( log_file_names ) == 0 ):
    raise error( 'log file %s not found' % ( log_file_name ) )
  log_file_names.sort()
  try:
    image_data = []
    log_file = file( log_file_names[ -1 ], mode = 'r' )
    line_number = 0
    for line in log_file:
      line_number = line_number + 1
      if ( line_number == 1 ):
        start_dt = line[ 0 : 8 ] + line[ 9 : 15 ]
        start_dt = datetime.datetime( int( start_dt[ 0 : 4 ] ), int( start_dt[ 4 : 6 ] ),
            int( start_dt[ 6 : 8 ] ), int( start_dt[ 8 : 10 ] ), 
            int( start_dt[ 10 : 12 ] ), int( start_dt[ 12 : 14 ] ) )
        words = line.split()
        pointing_name = words[ 5 ]
        continue
      if ( line[ 0 ] == '*' ):
        continue
      words = line.split()
      if ( line[ 0 : 5 ] == 'UVCOP' ):
        if ( 'vis records' in line ):
          if ( not 'flagged' in line ):
            vis = int( float( words[ 2 ] ) )
      elif ( line[ 0 : 5 ] == 'SPLIT' ):
        if ( 'Visibilities written' in line ):
          vis = int( words[ 1 ] )
        elif ( 'Partially' in line ):
          vis = int( words[ 4 ] )
        elif ( 'Fully' in line ):
          vis = vis + int( words[ 4 ] )
      elif ( line[ 0 : 5 ] == 'IMEAN' ):
        if ( 'from histogram' in line ):
          if ( words[ 2 ] == 'Rms=' ):
            rms= float( words[ 3 ] )
          else:
            rms= float( words[ 4 ] )
      elif ( line[ 0 : 5 ] == 'IMAGR' ):
        if ( 'Total Cleaned flux density' in line ):
          dflux = float( words[ 6 ] )
          if ( words[ 7 ][ -2 : ] == 'Jy' ):
            flux_units = words[ 7 ].lower()
            dcc_count = int( words[ 8 ] )
          else:
            flux_units = words[ 7 ].lower() + words[ 8 ].lower()
            dcc_count = int( words[ 9 ] )
          if ( flux_units == 'jy' ):
            dflux = dflux * 1.e0
          elif ( flux_units == 'millijy' ):
            dflux = dflux * 1.e-3
          elif ( flux_units == 'microjy' ):
            dflux = dflux * 1.e-6
          elif ( flux_units == 'nanojy' ):
            dflux = dflux * 1.e-9
          elif ( flux_units == 'kilojy' ):
            dflux = dflux * 1.e3
          elif ( flux_units == 'megajy' ):
            dflux = dflux * 1.e6
          else:
            raise ValueError
        elif ( 'Merging the Clean components files' in line ):
          if first:
            flux = flux + dflux
            cc_count = cc_count + dcc_count
          else:
            flux = dflux
            cc_count = dcc_count
      elif ( line[ 0 : 3 ] == '...' ):
        if ( 'total # components' in line ):
          cc_count = int( words[ 5 ] )
        elif ( 'total cleaned flux' in line ):
          flux = float( words[ 5 ] )
          flux_units = 'Jy'
      elif ( line[ 0 : 5 ] == 'FITTP' ):
        if ( '.FITS' in line ):
          image = words[ -1 ].split( '.' )[ -2 ]
          image_data.append( [ image, flux, cc_count, rms, vis ] )
          first = False
      elif ( line[ 0 : 5 ] == 'SETFC' ):
        if ( 'Cpu=' in line ):
          host = words[ 1 ].upper()
    end_dt = line[ 0 : 8 ] + line[ 9 : 15 ]
    end_dt = datetime.datetime( int( end_dt[ 0 : 4 ] ), int( end_dt[ 4 : 6 ] ),
        int( end_dt[ 6 : 8 ] ), int( end_dt[ 8 : 10 ] ), 
        int( end_dt[ 10 : 12 ] ), int( end_dt[ 12 : 14 ] ) )
    if ( 'failed' in line ):
      status = 'FAILED'
    elif ( 'normally' in line ):
      status = 'DONE'
    elif ( 'interrupted' in line ):
      status = 'STOPPED'
    log_file.close()
  except KeyboardInterrupt:
    log_file.close()
  except IOError:
    raise error( 'error reading log file %s' % ( log_file_name ) )
  except:
    log_file.close()
    raise error( 'interpretation error on line %d' % ( line_number ) )
  
  print 'field %s is processed with status = %s' % ( 
      pointing_name, status )
  print '... processing took %s hours on host %s' % ( str( end_dt - start_dt ), host )
  for [ image, flux, cc_count, rms, vis ] in image_data:
    print( '...... %4s image: %6d visibilities, ' % ( image, vis ) + 
        '%7.3f Jy in %5d CLEAN components, ' % ( flux, cc_count ) +
        '%6.3f mJy/beam RMS noise ' % ( rms * 1.e3 ) )
  
  return True

###############################################################################

def create_central_uv_data( uv, pb_facets, radius, radec = None ):
# radius in arcminutes
  sub_pb_uv = subtract_model( uv, pb_facets )
  while table_exists( sub_pb_uv, 'FG', 0 ):
    flag_uv = apply_flag_table( sub_pb_uv )
    sub_pb_uv.zap()
    flag_uv.rename( name = sub_pb_uv.name, klass = sub_pb_uv.klass, seq = sub_pb_uv.seq )
  if table_exists( pb_facets, 'SN', 0 ):
    call_aips_task( 'TACOP', indata = pb_facets, outdata = sub_pb_uv, inext = 'SN',
        ncount = 1 )
  else:
    if ( not table_exists( uv, 'SN', 0 ) ):
      calibrate_model( uv, pb_facets, reference_antenna = reference_antenna,
          calib_params = calib_params, phase_interval = phase_interval_min,
          apply_solutions = False )
    call_aips_task( 'TACOP', indata = uv, outdata = sub_pb_uv, inext = 'SN',
        ncount = 1 )
  if table_exists( sub_pb_uv, 'SN', 0 ):
    cal_uv = apply_solution_table( sub_pb_uv )
    sub_pb_uv.zap()
    cal_uv.rename( name = sub_pb_uv.name, klass = sub_pb_uv.klass, seq = sub_pb_uv.seq )
  for table in [ 'SN', 'NI', 'OB' ]:
    while table_exists( sub_pb_uv, table, 0 ):
      sub_pb_uv.zap_table( table, 0 )
  if ( radec is None ):
    radec = get_radec( uv )
  add_uv = add_model_in_circle( sub_pb_uv, pb_facets, radec, radius,
      apply_solutions = False )
  sub_pb_uv.zap()
  return add_uv

###############################################################################

def add_model_in_circle( uv, pb_facets, radec, radius, apply_solutions = True ):
# radius in arcmin
  add_uv = get_aips_file( uv.disk, uv.name, 'ADD', -1, 'UV' )
  call_aips_task( 'MOVE', indata = uv, outdata = add_uv,
      userid = get_aips_userid() )
  add_facets = get_aips_file( pb_facets.disk, 'ADD', 'ICL001', -1, 'MA' )
  copy_facets( pb_facets, add_facets, include_beams = False )
  facet_file_name = restore_parameter( add_facets, 'facet_file_name' )
  facet_file_name_e = path.expandvars( facet_file_name )
  add_facet_file_name_e = facet_file_name_e + '.ADD'
  facet_count = restore_parameter( add_facets, 'facet_count' )
  extract_facet_definitions( facet_file_name_e, range( 1, 1 + facet_count ),
      add_facet_file_name_e, include_clean_boxes = False )
  cell_size = restore_parameter( add_uv, 'cell_size' )
  radius_pix = int( ceil( radius * 60. / cell_size ) )
  for i in range( 1, 1 + facet_count ):
    facet = get_facet( add_facets, i )
    pos = calculate_source_position( facet, radec, to_grid = True )
    add_circular_clean_box( add_facet_file_name_e, i, pos, radius_pix )
  for i in range( 1, 1 + facet_count ):
    extract_model_components( add_facets, i, add_facet_file_name_e )
  new_add_uv = add_model( add_uv, add_facets, keep_flags = False,
      apply_solutions = apply_solutions )
  if ( new_add_uv != None ):
    add_uv.zap()
    new_add_uv.rename( name = add_uv.name, klass = add_uv.klass,
        seq = add_uv.seq )
  remove_facets( add_facets )
  return add_uv

###############################################################################

def subtract_model_in_circle( uv, pb_facets, radec, radius, apply_solutions = True ):
# radius in arcmin
  sub_uv = get_aips_file( uv.disk, uv.name, 'SUB', -1, 'UV' )
  call_aips_task( 'MOVE', indata = uv, outdata = sub_uv,
      userid = get_aips_userid() )
  sub_facets = get_aips_file( pb_facets.disk, 'SUB', 'ICL001', -1, 'MA' )
  copy_facets( pb_facets, sub_facets, include_beams = False )
  facet_file_name = restore_parameter( sub_facets, 'facet_file_name' )
  facet_file_name_e = path.expandvars( facet_file_name )
  sub_facet_file_name_e = facet_file_name_e + '.SUB'
  facet_count = restore_parameter( sub_facets, 'facet_count' )
  extract_facet_definitions( facet_file_name_e, range( 1, 1 + facet_count ),
      sub_facet_file_name_e, include_clean_boxes = False )
  cell_size = restore_parameter( sub_uv, 'cell_size' )
  radius_pix = int( ceil( radius * 60. / cell_size ) )
  for i in range( 1, 1 + facet_count ):
    facet = get_facet( sub_facets, i )
    pos = calculate_source_position( facet, radec, to_grid = True )
    sub_circular_clean_box( sub_facet_file_name_e, i, pos, radius_pix )
  for i in range( 1, 1 + facet_count ):
    extract_model_components( sub_facets, i, sub_facet_file_name_e )
  new_sub_uv = subtract_model( sub_uv, sub_facets, keep_flags = False,
      apply_solutions = apply_solutions )
  if ( new_sub_uv != None ):
    sub_uv.zap()
    new_sub_uv.rename( name = sub_uv.name, klass = sub_uv.klass,
        seq = sub_uv.seq )
  remove_facets( sub_facets )
  return sub_uv

###############################################################################


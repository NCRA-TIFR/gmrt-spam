# import Python modules
from os import *
from math import *
from string import *
from datetime import *
import sys, os, glob, traceback
from numpy import *

# import user modules
from spam import *

###############################################################################

def process_lta( lta_file_name, pointings = [], keep_log = True, 
    primary_cal = None, correction_list = [] ):
  
  # get full paths to data directories
  dat_path = path.expandvars( '${DAT}/' )
  fits_path = path.expandvars( '${FIT}/' )
  
  # check presence of input file
  pointing_name = lta_file_name.split( '/' )[ -1 ]
  pointing_name = pointing_name.split( '.' )[ 0 ]
  
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
  print '%s: started processing of observation %s' % ( dt_string, pointing_name )
  if keep_log:
    log_file_name = dat_path + 'spam_%s_%s.log' % ( pointing_name, dt_string )
    log_file = file( log_file_name, mode = 'w' )
    old_stdout = sys.stdout
    old_stderr = sys.stderr
    sys.stdout = log_file
    sys.stderr = log_file
    print '%s: started processing of observation %s' % ( dt_string, pointing_name )
    sys.stdout.flush()
  
  # process pointing
  try:
    spam_lta( lta_file_name, pointings = pointings, primary_cal = primary_cal,
        correction_list = correction_list )
  except KeyboardInterrupt:
    dt = datetime.datetime.now()
    dt_string = '%4d%02d%02d_%02d%02d%02d' % ( 
        dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second )
    if keep_log:
      print '%s: processing of observation %s manually interrupted' % ( dt_string, pointing_name )
      sys.stdout.flush()
      sys.stdout = old_stdout
      sys.stderr = old_stderr
      log_file.close()
    print '%s: processing of observation %s manually interrupted' % ( dt_string, pointing_name )
  except:
    print traceback.print_exc()
    dt = datetime.datetime.now()
    dt_string = '%4d%02d%02d_%02d%02d%02d' % ( 
        dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second )
    if keep_log:
      print '%s: processing of observation %s failed' % ( dt_string, pointing_name )
      sys.stdout.flush()
      sys.stdout = old_stdout
      sys.stderr = old_stderr
      log_file.close()
    print '%s: processing of observation %s failed' % ( dt_string, pointing_name )
  else:
    dt = datetime.datetime.now()
    dt_string = '%4d%02d%02d_%02d%02d%02d' % ( 
        dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second )
    if keep_log:
      print '%s: processing of observation %s ended normally' % ( dt_string, pointing_name )
      sys.stdout.flush()
      sys.stdout = old_stdout
      sys.stderr = old_stderr
      log_file.close()
    print '%s: processing of observation %s ended normally' % ( dt_string, pointing_name )
  
  return

###############################################################################

def spam_lta( lta_file_name, pointings = [], do_tb_sort = False, max_scans = 128,
    primary_cal = None, correction_list = [] ):
  # correction_list = [ [ lta#, scan#, old_source, new_source ] ]
  
  # get full paths to data directories
  dat_path = path.expandvars( '${DAT}/' )
  prt_path = path.expandvars( '${PRT}/' )
  fits_path = path.expandvars( '${FIT}/' )
  
  # processing parameters
  project_name = 'TGSS'
  bpcal_names = [ '3C48','3C147','3C286' ]
  if ( not primary_cal is None ):
    if ( primary_cal in bpcal_names ):
      bpcal_names = [ primary_cal ]
    else:
      raise error( 'specified primary calibrator needs to be 3C48/147/286' )
  
  # check presence of LTA file(s)
  is_hw_correlator = False
  lta_file_names = []
  try:
    lta_index = lta_file_name.index( '.lt' )
  except:
    raise error( 'LTA file name %s has unexpected format (not .lta or .ltb)' % 
        ( lta_file_name ) )
  file_name = lta_file_name[ : lta_index + 4 ]
  if ( not file_name[ -3 : ] in [ 'lta', 'ltb' ] ):
    raise error( 'LTA file name %s has unexpected format (not .lta or .ltb)' % 
        ( lta_file_name ) )
  file_name = file_name[ : -4 ]
  flags_file_name = file_name + '.lta*##*.FLAGS*'
  for i in range( 10 ):
    if ( i == 0 ):
      extension = ''
    else:
      extension = '.%d' % ( i )
    if file_exists( file_name + '.lta' + extension ):
      lta_file_names.append( file_name + '.lta' + extension )
      if file_exists( file_name + '.ltb' + extension ):
        lta_file_names.append( file_name + '.ltb' + extension )
        is_hw_correlator = True
      elif is_hw_correlator:
        raise error( 'LTA file %s does not exist' % 
            ( file_name + '.ltb' + extension ) )
    elif file_exists( file_name + '.ltb' + extension ):
      raise error( 'LTA file %s does not exist' % 
          ( file_name + '.lta' + extension ) )
  if ( len( lta_file_names ) == 0 ):
    raise error( 'LTA file(s) %s does not exist' % ( lta_file_name ) )
  
  # convert LTA files into UVFITS files
  uvfits_file_names = []
  source_list = []
  for lta_file_name in lta_file_names:
    file_name = lta_file_name.split( '/' )[ -1 ]
    extension = '.' + file_name.split( '.' )[ -1 ]
    file_name = file_name[ : - len( extension ) ]
    if ( extension == '.lta' ):
      extcode = '.0A'
    elif ( extension == '.ltb' ):
      extcode = '.0B'
    else:
      extcode = extension
      extension = '.' + file_name.split( '.' )[ -1 ]
      file_name = file_name[ : - len( extension ) ]
      if ( extension == '.lta' ):
        extcode = extcode + 'A'
      else:
        extcode = extcode + 'B'
    lta_correction_list = []
    if ( len( correction_list ) > 0 ):
      lta_number = int( extcode[ 1 ] )
      for c in correction_list:
        if ( c[ 0 ] == lta_number ):
          lta_correction_list.append( c[ 1 : ] )
    for i in range( 1000 ):
      uvfits_file_name = fits_path + file_name.upper() + extcode + '.UVFITS.%03d' % ( i )
      if file_exists( uvfits_file_name ):
        remove_file( uvfits_file_name )
      new_source_list = convert_lta_to_uvfits( lta_file_name, uvfits_file_name, 
          max_scans = max_scans, scan_offset = i * max_scans, 
          correction_list = lta_correction_list, target_list = pointings )
      if ( new_source_list is None ):
        break
      source_list = source_list + new_source_list
      uvfits_file_names.append( uvfits_file_name )
  uvfits_file_names.sort()
  
  # filter source names
  target_names = []
  for source_name in source_list:
    if ( source_name in bpcal_names ):
      continue
    if ( source_name in target_names ):
      continue
    if ( len( pointings ) > 0 ):
      if ( source_name in pointings ):
        target_names.append( source_name )
      continue
    elif ( len( source_name ) == 6 ):
      if ( ( source_name[ 0 ] == 'R' ) and ( source_name[ 3 ] == 'D' ) ):
        target_names.append( source_name )
  target_names.sort()
  
  # pre-calibrate target fields
  if is_hw_correlator:
    half = len( uvfits_file_names ) / 2
    target_uvfits_file_names_1 = pre_calibrate_targets( uvfits_file_names[ : half ],
        flags_file_name = flags_file_name, bw_smearing = 0.5, keep_channel_one = True,
        target_names = target_names, do_tb_sort = do_tb_sort, bpcal_names = bpcal_names )
    target_uvfits_file_names_2 = pre_calibrate_targets( uvfits_file_names[ half : ],
        flags_file_name = flags_file_name, bw_smearing = 0.5, keep_channel_one = True,
        target_names = target_names, do_tb_sort = do_tb_sort, bpcal_names = bpcal_names )
    if ( len( target_uvfits_file_names_1 ) != len( target_uvfits_file_names_2 ) ):
      raise error( 'HW correlator LTA files do not have the same number of targets' )
    target_uvfits_file_names_1.sort()
    target_uvfits_file_names_2.sort()
    target_uvfits_file_names = []
    for i in range( len( target_uvfits_file_names_1 ) ):
      if ( '_LL_' in target_uvfits_file_names_1[ i ] ):
        uvfits_file_name = target_uvfits_file_names_1[ i ].replace( '_LL_', '_RRLL_' )
      else:
        uvfits_file_name = target_uvfits_file_names_2[ i ].replace( '_LL_', '_RRLL_' )
      combine_stokes( target_uvfits_file_names_1[ i ], target_uvfits_file_names_2[ i ],
          uvfits_file_name, time_offset = 0 )
      remove_file( target_uvfits_file_names_1[ i ] )
      remove_file( target_uvfits_file_names_2[ i ] )
      target_uvfits_file_names.append( uvfits_file_name )
  else:
    target_uvfits_file_names = pre_calibrate_targets( uvfits_file_names,
        flags_file_name = flags_file_name, bw_smearing = 0.5, do_tb_sort = do_tb_sort,
        target_names = target_names, bpcal_names = bpcal_names )
  
  # cleanup and exit
  for uvfits_file_name in uvfits_file_names:
    remove_file( fits_path + uvfits_file_name )
  
  return target_uvfits_file_names

###############################################################################

def process_tgss_pointing( uvfits_file_name, keep_log = True, minimize_storage = True,
    do_cleanup = False, use_nvss = True ):
  
  # get full paths to data directories
  dat_path = path.expandvars( '${DAT}/' )
  prt_path = path.expandvars( '${PRT}/' )
  fits_path = path.expandvars( '${FIT}/' )
  
  # check presence of input file(s)
  uvfits_file_names = glob.glob( uvfits_file_name )
  if ( len( uvfits_file_names ) == 0 ):
    raise error( 'UVFITS file(s) %s does not exist' % ( uvfits_file_name ) )
  pointing_name = uvfits_file_names[ 0 ].split( '/' )[ -1 ]
  pointing_name = pointing_name.split( '.' )[ 0 ]
  target_name = pointing_name.split( '_' )[ 0 ]
  
  # create symbolic links with short file names
  for i in range( len( uvfits_file_names ) ):
    uvfits_file_name = uvfits_file_names[ i ]
    new_uvfits_file_name = '%s%s_%d.UVFITS' % ( fits_path, target_name, i + 1 )
    if file_exists( new_uvfits_file_name ):
      remove_file( new_uvfits_file_name )
    symlink( uvfits_file_name, new_uvfits_file_name )
    uvfits_file_names[ i ] = new_uvfits_file_name
  
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
    spam_pointing( t_iical_uv, minimize_storage = minimize_storage,
        imagr_params = { 'robust' : -1., 'uvrang' : [ 0.2, 1.e6 ] },
        calib_params = { 'uvrang' : [ 1., 1.e6 ] }, use_nvss = use_nvss )
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

def process_pointing( uvfits_file_name, keep_log = True, minimize_storage = True,
    do_cleanup = False, model_file_name = None, snapshot = True, use_nvss = True,
    calib_params = { 'uvrang' : [ 1.,1.e6 ] }, imagr_params = { 'robust' : -1. },
    do_large_boxes = False, do_plots = False ):
  
  # get full paths to data directories
  dat_path = path.expandvars( '${DAT}/' )
  prt_path = path.expandvars( '${PRT}/' )
  fits_path = path.expandvars( '${FIT}/' )
  
  # check presence of input file(s)
  uvfits_file_names = glob.glob( uvfits_file_name )
  if ( len( uvfits_file_names ) == 0 ):
    raise error( 'UVFITS file(s) %s does not exist' % ( uvfits_file_name ) )
  pointing_name = uvfits_file_names[ 0 ].split( '/' )[ -1 ]
  pointing_name = pointing_name.split( '.' )[ 0 ]
  target_name = pointing_name.split( '_' )[ 0 ]
  model_source_list = []
  if ( not model_file_name is None ):
    if ( model_file_name != '' ):
      if ( not file_exists( model_file_name ) ):
        raise error( 'model file %s does not exist' % ( model_file_name ) )
      try:
        model_cat = read_starbase_catalog( model_file_name )
        ra_index = model_cat[ 0 ].index( 'RA' )
        dec_index = model_cat[ 0 ].index( 'DEC' )
        flux_index = model_cat[ 0 ].index( 'FLUX' )
        for source in model_cat[ 1 : ]:
          model_source_list.append( [ [ source[ ra_index ], source[ dec_index ] ], 
              source[ flux_index ] ] )
      except:
        raise error( 'model file %s could not be read' % ( model_file_name ) )
  
  # create symbolic links with short file names
  for i in range( len( uvfits_file_names ) ):
    uvfits_file_name = uvfits_file_names[ i ]
    new_uvfits_file_name = '%s%s_%d.UVFITS' % ( fits_path, target_name, i + 1 )
    if file_exists( new_uvfits_file_name ):
      remove_file( new_uvfits_file_name )
    symlink( uvfits_file_name, new_uvfits_file_name )
    uvfits_file_names[ i ] = new_uvfits_file_name
  
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
    spam_pointing( t_iical_uv, minimize_storage = minimize_storage,
        imagr_params = imagr_params, calib_params = calib_params,
        use_nvss = use_nvss, model_source_list = model_source_list,
        snapshot = snapshot, do_large_boxes = do_large_boxes, do_plots = do_plots )
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

def spam_pointing( t_iical_uv, minimize_storage = True, model_source_list = [],
    calib_params = { 'uvrang' : [ 1.,1.e6 ] }, imagr_params = { 'robust' : -1. },
    do_large_boxes = False, do_plots = False, snapshot = True, use_nvss = True ):
  
  # get full paths to data directories
  dat_path = path.expandvars( '${DAT}/' )
  prt_path = path.expandvars( '${PRT}/' )
  fits_path = path.expandvars( '${FIT}/' )
  
  # run pipeline
  spam_target( t_iical_uv, minimize_storage = minimize_storage, 
      model_source_list = model_source_list, calib_params = calib_params,
      imagr_params = imagr_params, do_large_boxes = do_large_boxes,
      do_plots = do_plots, snapshot = snapshot, use_nvss = use_nvss )
  
  return

###############################################################################

def create_tgss_mosaic( mosaic_name, mosaic_center, mosaic_size, do_galactic = False,
    tgss_path = '/lustre/hintema/tgss/adr1/', allow_missing_pointings = True,
    restoring_beam = [ 25.,0.,0. ], make_beam_plot = False, beam_radius = 2.,
    max_image_size = 16000, convolve_background = True, pixels_per_beam = 4.,
    get_pointings = False, pointing_model = [ 4.05, 24.6 ], weight_power = -2.,
    pbparms = [ 0.36, 1., -2.460, 10.76, -6.853, 35.73, -25.82 ], try_new = True ):
  gmrt_lonlat = dmsdms_to_degdeg( [ 74,2,59.07, 19,5,47.46 ] )
  if ( restoring_beam[ 1 ] > 0. ):
    mosaic_beam_size = restoring_beam
  else:
    bmin = restoring_beam[ 0 ]
    bmaj = around( bmin / cos( radians( min( mosaic_center[ 1 ] - gmrt_lonlat[ 1 ], 0. ) ) ), 2 )
    mosaic_beam_size = [ bmaj, bmin, 0. ]
  if try_new:
    result = create_mosaic_new( mosaic_name, mosaic_center, mosaic_size, 
        pointing_fits_path = tgss_path + 'fits/', do_galactic = do_galactic,
        pointing_file_name = tgss_path + 'pointing_grid.rdb',
        restoring_beam = mosaic_beam_size, max_image_size = max_image_size,
        convolve_background = convolve_background, pointing_model = pointing_model,
        allow_missing_pointings = allow_missing_pointings, beam_radius = beam_radius,
        make_beam_plot = make_beam_plot, pixels_per_beam = pixels_per_beam,
        get_pointings = get_pointings, pbparms = pbparms, weight_power = weight_power )
  else:
    result = create_mosaic( mosaic_name, mosaic_center, mosaic_size, 
        pointing_fits_path = tgss_path + 'fits/', do_galactic = do_galactic,
        pointing_file_name = tgss_path + 'pointing_grid.rdb',
        restoring_beam = mosaic_beam_size, max_image_size = max_image_size,
        convolve_background = convolve_background, pointing_model = pointing_model,
        allow_missing_pointings = allow_missing_pointings, beam_radius = beam_radius,
        make_beam_plot = make_beam_plot, pixels_per_beam = pixels_per_beam,
        get_pointings = get_pointings, pbparms = pbparms )
  
  return result

###############################################################################

def create_mosaic( mosaic_name, mosaic_center, mosaic_size, do_galactic = False,
    pointing_fits_path = '/lustre/hintema/tgss/fits/', get_pointings = False,
    pointing_file_name = '/lustre/hintema/tgss/pointing_grid.rdb',
    restoring_beam = [ 20.,0.,0. ], make_beam_plot = False, beam_radius = 1.8,
    max_image_size = 16000, convolve_background = False, pixels_per_beam = 4.,
    allow_missing_pointings = True, pbparms = None, pointing_model = None ):
# if do_galactic is False, mosaic_center is [ RA,DEC ]
# if do_galactic is True, mosaic_center is [ GLON,GLAT ]
# mosaic_size in degrees - box: hor x ver; circle: diam x -1
# pointing_model: [ arcsec/deg DEC, DEC_0 ]
  
  # get full paths to data directories
  dat_path = os.path.expandvars( '${DAT}/' )
  prt_path = os.path.expandvars( '${PRT}/' )
  fits_path = os.path.expandvars( '${FIT}/' )
  
  # identify pointings in mosaic
  if ( mosaic_center == [ 0.,0. ] ):
    mosaic_center = [ 0., 1. / 360000. ]
  if do_galactic:
    mosaic_radec = galactic_to_equatorial( mosaic_center )
    lb = equatorial_to_galactic( mosaic_radec )
    dlb = equatorial_to_galactic( [ mosaic_radec[ 0 ], mosaic_radec[ 1 ] + 1. ] )
    [ r, p ] = calculate_angular_separation( lb, dlb )
    rotation = -p
  else:
    mosaic_radec = mosaic_center
    rotation = 0.
  pointing_grid = read_starbase_catalog( pointing_file_name )
  mosaic_pointings = pointing_grid[ 0 : 1 ]
  for pointing in pointing_grid[ 1 : ]:
    [ pointing_name, pointing_ra, pointing_dec ] = pointing
    pointing_center = [ pointing_ra, pointing_dec ]
    if do_galactic:
      pointing_center = equatorial_to_galactic( pointing_center )
    [ r, p ] = calculate_angular_separation( mosaic_center, pointing_center )
    if ( mosaic_size[ 1 ] < 0. ):
      radius = mosaic_size[ 0 ] / 2. + beam_radius
      if ( r < radius ):
        mosaic_pointings.append( pointing )
    else:
      center = calculate_offset_position( [ 0.,0. ], r, p )
      if ( ( abs( center[ 1 ] ) < mosaic_size[ 1 ] / 2. + beam_radius ) and
          ( abs( amodulo( center[ 0 ] + 180., 360 ) - 180. ) < 
          mosaic_size[ 0 ] / 2. + beam_radius ) ):
        mosaic_pointings.append( pointing )
  if ( len( mosaic_pointings ) <= 1 ):
    raise error( 'no pointings found' )
  write_starbase_catalog( mosaic_pointings, '%s/%s.MOSAIC.rdb' % ( fits_path, mosaic_name ) )
  
  # check images
  i = 0
  for [ pointing_name, ra, dec ] in mosaic_pointings[ 1 : ]:
    pointing_name = pointing_name.split( '_' )[ 0 ]
    pointing_file_names = glob.glob( '%s/%s*.RES.FITS' % 
        ( pointing_fits_path, pointing_name ) )
    if ( len( pointing_file_names ) == 0 ):
      print 'pointing %s has no image' % pointing_name
  if get_pointings:
    return [ x[ 0 ] for x in mosaic_pointings[ 1 : ] ]
  
  # find free AIPS disk
  aips_disk = 0
  for i in range( 1, len( AIPS.disks ) ):
    if ( len( AIPSCat( i )[ i ] ) == 0 ):
      aips_disk = i
      break
  if ( aips_disk == 0 ):
    raise error( 'no free AIPS disk available' )
  
  # read in FITS residual images
  pointings = get_aips_file( aips_disk, mosaic_name[ 0 :12 ], 'IMG001', -1, 'MA' )
  i = 0
  for [ pointing_name, ra, dec ] in mosaic_pointings[ 1 : ]:
    pointing_name = pointing_name.split( '_' )[ 0 ]
    pointing_file_names = glob.glob( '%s/%s*.RES.FITS' % 
        ( pointing_fits_path, pointing_name ) )
    if ( len( pointing_file_names ) == 0 ):
      if ( not allow_missing_pointings ):
        raise error( 'pointing %s has no image' % pointing_name )
        if pointings.exists():
          remove_facets( pointings )
    for pointing_file_name in pointing_file_names:
      i = i + 1
      pointing = get_facet( pointings, i )
      read_fits_image( pointing_file_name, pointing )
  pointing_count = i
  if ( pointing_count == 0 ):
    raise error( 'no pointings found' )
  store_parameter( pointings, 'facet_count', pointing_count )
  
  # get some image statistics
  radecs = []
  beam_sizes = []
  noises = []
  for i in range( 1, 1 + pointing_count ):
    pointing = get_facet( pointings, i )
    radec = get_radec( pointing )
    radecs.append( radec )
    beam_size = get_beam_size( pointing )
    beam_sizes.append( beam_size )
    [ avg, noise ] = call_aips_task( 'IMEAN', indata = pointing,
        pixavg = 0., pixstd = 5.e-3, pixrange = [ -50.e-3, 50.e-3  ], 
        outputs = [ 'pixavg', 'pixstd' ] )
    noises.append( noise )
  radecs = array( radecs )
  beam_sizes = array( beam_sizes )
  noises = array( noises )
  
  # define restoring beam (and make beam plot)
  gmrt_lonlat = dmsdms_to_degdeg( [ 74,2,59.07, 19,5,47.46 ] )
  if ( restoring_beam[ 1 ] > 0. ):
    [ bmaj, bmin, bpa ] = restoring_beam
    mosaic_beam_size = restoring_beam
  else:
    bmin = restoring_beam[ 0 ]
    bmaj = around( bmin / cos( radians( mosaic_radec[ 1 ] - gmrt_lonlat[ 1 ] ) ), 1 )
    mosaic_beam_size = [ bmaj, bmin, 0. ]
  if make_beam_plot:
    dec_min = radecs[ : , 1 ].min()
    dec_max = radecs[ : , 1 ].max()
    decs = array( [ dec_min + ( dec_max - dec_min ) * ( float( x ) / 1000. )
        for x in range( 1001 ) ] )
    plot( decs, 0.9 * bmin / cos( radians( decs - gmrt_lonlat[ 1 ] ) ), 'b' )
    plot( decs, 1.0 * bmin / cos( radians( decs - gmrt_lonlat[ 1 ] ) ), 'g' )
    plot( decs, 1.1 * bmin / cos( radians( decs - gmrt_lonlat[ 1 ] ) ), 'r' )
    plot( radecs[ : , 1 ], beam_sizes[ : , 0 ], 'k.' )
    show()
    clear_aips_disk( aips_disk )
    return
  
  # define mosaic parameters
  mosaic_pixel_size = around( mosaic_beam_size[ 1 ] / pixels_per_beam, 1 )
  if ( mosaic_size[ 1 ] <= 0. ):
    mosaic_image_size = [ int( around( mosaic_size[ 0 ] * 3600. / mosaic_pixel_size ) ),
        int( around( mosaic_size[ 0 ] * 3600. / mosaic_pixel_size ) ) ]
  else:
    mosaic_image_size = [ int( around( mosaic_size[ 0 ] * 3600. / mosaic_pixel_size ) ),
        int( around( mosaic_size[ 1 ] * 3600. / mosaic_pixel_size ) ) ]
  if ( max( mosaic_image_size ) > max_image_size ):
    raise error( 'image size %s exceeds maximum' % ( repr( mosaic_image_size ) ) )
    remove_facets( pointings )
  
  # calculate restoring beam position angle per pointing
  pa_list = []
  for i in range( 1, 1 + pointing_count ):
    pointing = get_facet( pointings, i )
    pointing_radec = get_radec( pointing )
    p = amodulo( pointing_radec[ 0 ] - mosaic_radec[ 0 ] + 180., 360. ) - 180.
    p = asign( pointing_radec[ 1 ] ) * ( 1. - cos( aradians( pointing_radec[ 1 ] ) ) ) * p
    pa = amodulo( mosaic_beam_size[ 2 ] + p + 90., 180. ) - 90.
    pa_list.append( pa )
  
  # convolve pointings to mosaic resolution
  if convolve_background:
    for i in range( 1, 1 + pointing_count ):
      pointing = get_facet( pointings, i )
      convl_pointing = get_aips_file( aips_disk, mosaic_name[ 0 : 12 ], 'CONVL', -1, 'MA' )
      beam_size = get_beam_size( pointing )
      factor = ( ( mosaic_beam_size[ 0 ] * mosaic_beam_size[ 1 ] ) / 
          ( beam_size[ 0 ] * beam_size[ 1 ] ) )
      if ( ( beam_size[ 0 ] < mosaic_beam_size[ 0 ] ) and
          ( beam_size[ 1 ] < mosaic_beam_size[ 1 ] ) ):
        try:
          call_aips_task( 'CONVL', indata = pointing, outdata = convl_pointing,
              bmaj = mosaic_beam_size[ 0 ], bmin = mosaic_beam_size[ 1 ],
              bpa = pa_list[ i - 1 ], opcode = 'GAUS', doblank = 1, factor = factor )
          pixels = get_image_pixels( convl_pointing )
          if alltrue( pixels == get_aips_magic_value() ):
            raise error( 'convolution failed' )
        except:
          beam_size[ 0 ] = mosaic_beam_size[ 0 ]
      if ( ( beam_size[ 0 ] >= mosaic_beam_size[ 0 ] ) and
          ( beam_size[ 1 ] < mosaic_beam_size[ 1 ] ) ):
        if convl_pointing.exists():
          convl_pointing.zap()
        try:
          call_aips_task( 'CONVL', indata = pointing, outdata = convl_pointing,
              bmaj = 1.001 * beam_size[ 0 ], bmin = mosaic_beam_size[ 1 ],
              bpa = pa_list[ i - 1 ], opcode = 'GAUS', doblank = 1, factor = factor )
          pixels = get_image_pixels( convl_pointing )
          if alltrue( pixels == get_aips_magic_value() ):
            raise error( 'convolution failed' )
        except:
          beam_size[ 1 ] = mosaic_beam_size[ 1 ]
      if ( ( beam_size[ 0 ] >= mosaic_beam_size[ 0 ] ) and
          ( beam_size[ 1 ] >= mosaic_beam_size[ 1 ] ) ):
        if convl_pointing.exists():
          convl_pointing.zap()
        # should we still scale?
        pixels = get_image_pixels( pointing )
        sel = awhere( pixels != get_aips_magic_value() )
        pixels = aput( pixels, sel, aget( pixels, sel ) * factor )
        set_image_pixels( pointing, pixels )
      if convl_pointing.exists():
        pointing.zap()
        convl_pointing.rename( name = pointing.name, klass = pointing.klass,
            seq = pointing.seq )
      set_beam_size( pointing, mosaic_beam_size[ 0 : 2 ] + [ pa_list[ i - 1 ] ] )
  
  # restore CLEAN components
  for i in range( 1, 1 + pointing_count ):
    pointing = get_facet( pointings, i )
    pixels = get_image_pixels( pointing )
    sel = awhere( pixels == get_aips_magic_value() )
    rst_pointing = restore_model_components( pointings, facet_list = [ i ],
        cross_restore = False, imagr_params = { 'bmaj' : mosaic_beam_size[ 0 ], 
        'bmin' : mosaic_beam_size[ 1 ], 'bpa' : pa_list[ i - 1 ] } )
    pixels = get_image_pixels( rst_pointing )
    pixels = aput( pixels, sel, get_aips_magic_value() )
    set_image_pixels( rst_pointing, pixels )
    pointing = get_facet( pointings, i )
    pointing.zap()
    rst_pointing.rename( name = pointing.name, klass = pointing.klass,
        seq = pointing.seq )
  
  # resample to final resolution
  # convert to galactic, if requested
  for i in range( 1, 1 + pointing_count ):
    pointing = get_facet( pointings, i )
    if do_galactic:
      convert_to_galactic( pointing )
    scale_pointing = scale_image( pointing, rotate = True,
        new_pixel_size = [ mosaic_pixel_size, mosaic_pixel_size ] )
    pointing.zap()
    scale_pointing.rename( name = pointing.name, klass = pointing.klass,
        seq = pointing.seq )
    if do_galactic:
      for ctype in pointing.header.ctype:
        if ( ctype.find( 'GLON' ) != - 1 ):
          ra_index = pointing.header.ctype.index( ctype )
        if ( ctype.find( 'GLAT' ) != - 1 ):
          dec_index = pointing.header.ctype.index( ctype )
      wiz_im = wizardry( pointing )
      wiz_im.header.ctype[ ra_index ] = 'RA--' + wiz_im.header.ctype[ ra_index ][ 4 : ]
      wiz_im.header.ctype[ dec_index ] = 'DEC-' + wiz_im.header.ctype[ dec_index ][ 4 : ]
      wiz_im.header.update()
      wiz_im = wizardry( pointing )
      del wiz_im
      set_epoch( pointing, 2000. )
  
  # select primary beam model parameters
  frequency = get_frequency( pointings )
  if ( pbparms is None ):
    if ( frequency > 140.e6 ) and ( frequency < 170.e6 ):
      pbparms = [ 0., 1., -4.04, 76.2, -68.8, 22.03, 0. ]
    elif ( frequency > 200.e6 ) and ( frequency < 280.e6 ):
      pbparms = [ 0., 1., -3.366, 46.159, -29.963, 7.529, 0. ]
    elif ( frequency > 300.e6 ) and ( frequency < 350.e6 ):
      pbparms = [ 0., 1., -3.397, 47.192, -30.931, 7.803, 0. ] 
    elif ( frequency > 580.e6 ) and ( frequency < 650.e6 ):
      pbparms = [ 0., 1., -3.486, 47.749, -35.203, 10.399, 0. ] 
    elif ( frequency > 1.e9 ) and ( frequency < 1.8e9 ):
      pbparms = [ 0., 1., -2.27961, 21.4611, -9.7929, 1.80153, 0. ] 
    else:
      raise error( 'beam shape is unknown' )
  else:
    pbparms[ 0 ] = 0.
    pbparms[ 1 ] = 1.
  pb_cutoff = calculate_pbparm_attenuation( frequency, beam_radius, pbparms )
  pbparms[ 0 ] = pb_cutoff
  
  # make RMS maps
  rms_pointings = get_aips_file( aips_disk, 'RMS', pointings.klass, -1, 'MA' )
  for i in range( 1, 1 + pointing_count ):
    pointing = get_facet( pointings, i )
    rms_pointing = get_facet( rms_pointings, i )
    call_aips_task( 'RMSD', indata = pointing, outdata = rms_pointing,
        imsize = [ 91, -1 ], xinc = 20, yinc = 20, optype = 'ROBS' )
  
  # make PB-corrected RMS maps
  pbrms_pointings = get_aips_file( aips_disk, 'PBRMS', pointings.klass, -1, 'MA' )
  for i in range( 1, 1 + pointing_count ):
    rms_pointing = get_facet( rms_pointings, i )
    pbrms_pointing = get_facet( pbrms_pointings, i )
    center = get_radec( rms_pointing )
    coordinates = degdeg_to_hmsdms( center )
    if ( center[ 1 ] < 0. ):
      coordinates[ 3 ] = -abs( coordinates[ 3 ] )
      coordinates[ 4 ] = -abs( coordinates[ 4 ] )
      coordinates[ 5 ] = -abs( coordinates[ 5 ] )
    call_aips_task( 'PBCOR', indata = rms_pointing, outdata = pbrms_pointing,
        pbparm = pbparms, coordina = coordinates )
  
  # make mask maps
  ones_pointings = get_aips_file( aips_disk, 'ONES', pointings.klass, -1, 'MA' )
  for i in range( 1, 1 + pointing_count ):
    pbrms_pointing = get_facet( pbrms_pointings, i )
    ones_pointing = get_facet( ones_pointings, i )
    units = pbrms_pointing.header.bunit
    call_aips_task( 'MATHS', indata = pbrms_pointing, outdata = ones_pointing,
        opcode = 'POLY', cparm = [ 1,0,0,0 ] )
    set_header_keyword( ones_pointing, 'bunit', units )
  
  # apply pointing model
  if ( not pointing_model is None ):
    
    # make masked RMS maps
    masked_rms_pointings = get_aips_file( aips_disk, 'MRMS', pointings.klass, -1, 'MA' )
    for i in range( 1, 1 + pointing_count ):
      rms_pointing = get_facet( rms_pointings, i )
      ones_pointing = get_facet( ones_pointings, i )
      masked_rms_pointing = get_facet( masked_rms_pointings, i )
      units = rms_pointing.header.bunit
      call_aips_task( 'COMB', indata = ones_pointing, in2data = rms_pointing,
          outdata = masked_rms_pointing, opcode = 'MULT', aparm = [ 1,0,0,0,0,0,0,0,0,0 ] )
      set_header_keyword( masked_rms_pointing, 'bunit', units )
    remove_facets( rms_pointings )
    rms_pointings = masked_rms_pointings
    
    # make shifted PB-corrected RMS maps
    # pointing_model: [ arcsec/deg DEC, DEC_0 ]
    remove_facets( pbrms_pointings )
    for i in range( 1, 1 + pointing_count ):
      rms_pointing = get_facet( rms_pointings, i )
      pbrms_pointing = get_facet( pbrms_pointings, i )
      center = get_radec( rms_pointing )
      dec_shift = pointing_model[ 0 ] * ( center[ 1 ] - pointing_model[ 1 ] ) / 3600.
      center[ 1 ] = center[ 1 ] + dec_shift
      pb_cutoff = calculate_pbparm_attenuation( frequency, 
          beam_radius + abs( dec_shift ), [ 0. ] + pbparms[ 1 : ] )
      coordinates = degdeg_to_hmsdms( center )
      if ( center[ 1 ] < 0. ):
        coordinates[ 3 ] = -abs( coordinates[ 3 ] )
        coordinates[ 4 ] = -abs( coordinates[ 4 ] )
        coordinates[ 5 ] = -abs( coordinates[ 5 ] )
      call_aips_task( 'PBCOR', indata = rms_pointing, outdata = pbrms_pointing,
          pbparm = [ pb_cutoff ] + pbparms[ 1 : ], coordina = coordinates )
  
  # multiply RMS maps with PB-corrected RMS maps
  pbrms2_pointings = get_aips_file( aips_disk, 'PBRMS2', pointings.klass, -1, 'MA' )
  for i in range( 1, 1 + pointing_count ):
    rms_pointing = get_facet( rms_pointings, i )
    pbrms_pointing = get_facet( pbrms_pointings, i )
    pbrms2_pointing = get_facet( pbrms2_pointings, i )
    units = pbrms_pointing.header.bunit
    call_aips_task( 'COMB', indata = rms_pointing, in2data = pbrms_pointing,
        outdata = pbrms2_pointing, opcode = 'MULT', aparm = [ 1,0,0,0,0,0,0,0,0,0 ] )
    set_header_keyword( pbrms2_pointing, 'bunit', units )
  remove_facets( rms_pointings )
  
  # make inverted PB-corrected RMS map
  ipbrms_pointings = get_aips_file( aips_disk, 'IPBRMS', pointings.klass, -1, 'MA' )
  for i in range( 1, 1 + pointing_count ):
    pbrms_pointing = get_facet( pbrms_pointings, i )
    ones_pointing = get_facet( ones_pointings, i )
    ipbrms_pointing = get_facet( ipbrms_pointings, i )
    units = pbrms_pointing.header.bunit
    call_aips_task( 'COMB', indata = ones_pointing, in2data = pbrms_pointing,
        outdata = ipbrms_pointing, opcode = 'DIV', aparm = [ 1,0,0,0,0,0,0,0,0,0 ] )
    set_header_keyword( ipbrms_pointing, 'bunit', units )
  remove_facets( ones_pointings )
  remove_facets( pbrms_pointings )
  
  # make weight maps
  weight_pointings = get_aips_file( aips_disk, 'WEIGHT', pointings.klass, -1, 'MA' )
  for i in range( 1, 1 + pointing_count ):
    ipbrms_pointing = get_facet( ipbrms_pointings, i )
    weight_pointing = get_facet( weight_pointings, i )
    units = ipbrms_pointing.header.bunit
    call_aips_task( 'COMB', indata = ipbrms_pointing, in2data = ipbrms_pointing,
        outdata = weight_pointing, opcode = 'MULT', aparm = [ 1,0,0,0,0,0,0,0,0,0 ] )
    set_header_keyword( weight_pointing, 'bunit', units )
  remove_facets( ipbrms_pointings )
  
  # make field maps
  field_pointings = get_aips_file( aips_disk, 'FIELD', pointings.klass, -1, 'MA' )
  for i in range( 1, 1 + pointing_count ):
    pointing = get_facet( pointings, i )
    pbrms2_pointing = get_facet( pbrms2_pointings, i )
    field_pointing = get_facet( field_pointings, i )
    units = pbrms2_pointing.header.bunit
    call_aips_task( 'COMB', indata = pointing, in2data = pbrms2_pointing,
        outdata = field_pointing, opcode = 'DIV', aparm = [ 1,0,0,0,0,0,0,0,0,0 ] )
    set_header_keyword( field_pointing, 'bunit', units )
  remove_facets( pointings )
  remove_facets( pbrms2_pointings )
  
  # combine field images and weight images
  coordinates = degdeg_to_hmsdms( mosaic_center )
  if ( mosaic_center[ 1 ] < 0. ):
    coordinates[ 3 ] = -abs( coordinates[ 3 ] )
    coordinates[ 4 ] = -abs( coordinates[ 4 ] )
    coordinates[ 5 ] = -abs( coordinates[ 5 ] )
  field_mosaic = get_aips_file( aips_disk, mosaic_name[ 0 : 12 ], 'FLDMOS', -1, 'MA' )
  call_aips_task( 'FLATN', indata = field_pointings, nfield = pointing_count,
      outdata = field_mosaic, imsize = mosaic_image_size, coordina = coordinates )
#      reweight = [ 3, 0.5 ] )
  remove_facets( field_pointings )
  weight_mosaic = get_aips_file( aips_disk, mosaic_name[ 0 : 12 ], 'WGTMOS', -1, 'MA' )
  call_aips_task( 'FLATN', indata = weight_pointings, nfield = pointing_count,
      outdata = weight_mosaic, imsize = mosaic_image_size, coordina = coordinates )
#      reweight = [ 3, 0.5 ] )
  remove_facets( weight_pointings )
  
  # make mosaic image
  mosaic = get_aips_file( aips_disk, mosaic_name[ 0 : 12 ], 'MOSAIC', -1, 'MA' )
  units = field_mosaic.header.bunit
  call_aips_task( 'COMB', indata = field_mosaic, in2data = weight_mosaic,
      outdata = mosaic, opcode = 'DIV', aparm = [ 1,0,0,0,0,0,0,0,0,0 ] )
  set_header_keyword( mosaic, 'bunit', units )
  field_mosaic.zap()
  weight_mosaic.zap()
  
  # write mosaic image to disk
  dummy = get_source_names( mosaic )[ 0 ]
  change_source_name( mosaic, dummy, mosaic_name )
  mosaic_beam_size[ 2 ] = amodulo( mosaic_beam_size[ 2 ] + rotation + 90., 180. ) - 90.
  set_beam_size( mosaic, mosaic_beam_size )
  set_observer( mosaic, 'TGSS' )
  if do_galactic:
    for ctype in mosaic.header.ctype:
      if ( ctype.find( 'RA--' ) != - 1 ):
        ra_index = mosaic.header.ctype.index( ctype )
      if ( ctype.find( 'DEC-' ) != - 1 ):
        dec_index = mosaic.header.ctype.index( ctype )
    wiz_im = wizardry( mosaic )
    wiz_im.header.ctype[ ra_index ] = 'GLON' + wiz_im.header.ctype[ ra_index ][ 4 : ]
    wiz_im.header.ctype[ dec_index ] = 'GLAT' + wiz_im.header.ctype[ dec_index ][ 4 : ]
    wiz_im.header.update()
    wiz_im = wizardry( mosaic )
    del wiz_im
  clear_history( mosaic )
  write_fits_image( mosaic, fits_path + mosaic_name + '.MOSAIC.FITS' )
  mosaic.zap()
  
  clear_aips_disk( aips_disk )
  return

###############################################################################

def create_mosaic_new( mosaic_name, mosaic_center, mosaic_size, do_galactic = False,
    pointing_fits_path = '/lustre/hintema/tgss/fits/', get_pointings = False,
    pointing_file_name = '/lustre/hintema/tgss/pointing_grid.rdb', weight_power = -2.,
    restoring_beam = [ 20.,0.,0. ], make_beam_plot = False, beam_radius = 1.8,
    max_image_size = 16000, convolve_background = False, pixels_per_beam = 4.,
    allow_missing_pointings = True, pbparms = None, pointing_model = None ):
# if do_galactic is False, mosaic_center is [ RA,DEC ]
# if do_galactic is True, mosaic_center is [ GLON,GLAT ]
# mosaic_size in degrees - box: hor x ver; circle: diam x -1
# pointing_model: [ arcsec/deg DEC, DEC_0 ]
  
  # get full paths to data directories
  dat_path = os.path.expandvars( '${DAT}/' )
  prt_path = os.path.expandvars( '${PRT}/' )
  fits_path = os.path.expandvars( '${FIT}/' )
  
  # identify pointings in mosaic
  if ( mosaic_center == [ 0.,0. ] ):
    mosaic_center = [ 0., 1. / 360000. ]
  if do_galactic:
    mosaic_radec = galactic_to_equatorial( mosaic_center )
    lb = equatorial_to_galactic( mosaic_radec )
    dlb = equatorial_to_galactic( [ mosaic_radec[ 0 ], mosaic_radec[ 1 ] + 1. ] )
    [ r, p ] = calculate_angular_separation( lb, dlb )
    rotation = -p
  else:
    mosaic_radec = mosaic_center
    rotation = 0.
  pointing_grid = read_starbase_catalog( pointing_file_name )
  mosaic_pointings = pointing_grid[ 0 : 1 ]
  for pointing in pointing_grid[ 1 : ]:
    [ pointing_name, pointing_ra, pointing_dec ] = pointing
    pointing_center = [ pointing_ra, pointing_dec ]
    if do_galactic:
      pointing_center = equatorial_to_galactic( pointing_center )
    [ r, p ] = calculate_angular_separation( mosaic_center, pointing_center )
    if ( mosaic_size[ 1 ] < 0. ):
      radius = mosaic_size[ 0 ] / 2. + beam_radius
      if ( r < radius ):
        mosaic_pointings.append( pointing )
    else:
      center = calculate_offset_position( [ 0.,0. ], r, p )
      if ( ( abs( center[ 1 ] ) < mosaic_size[ 1 ] / 2. + beam_radius ) and
          ( abs( amodulo( center[ 0 ] + 180., 360 ) - 180. ) < 
          mosaic_size[ 0 ] / 2. + beam_radius ) ):
        mosaic_pointings.append( pointing )
  if ( len( mosaic_pointings ) <= 1 ):
    raise error( 'no pointings found' )
  write_starbase_catalog( mosaic_pointings, '%s/%s.MOSAIC.rdb' % ( fits_path, mosaic_name ) )
  
  # check images
  i = 0
  for [ pointing_name, ra, dec ] in mosaic_pointings[ 1 : ]:
    pointing_name = pointing_name.split( '_' )[ 0 ]
    pointing_file_names = glob.glob( '%s/%s*.RES.FITS' % 
        ( pointing_fits_path, pointing_name ) )
    if ( len( pointing_file_names ) == 0 ):
      print 'pointing %s has no image' % pointing_name
  if get_pointings:
    return [ x[ 0 ] for x in mosaic_pointings[ 1 : ] ]
  
  # find free AIPS disk
  aips_disk = 0
  for i in range( 1, len( AIPS.disks ) ):
    if ( len( AIPSCat( i )[ i ] ) == 0 ):
      aips_disk = i
      break
  if ( aips_disk == 0 ):
    raise error( 'no free AIPS disk available' )
  
  # read in FITS residual images
  pointings = get_aips_file( aips_disk, mosaic_name[ 0 :12 ], 'IMG001', -1, 'MA' )
  i = 0
  for [ pointing_name, ra, dec ] in mosaic_pointings[ 1 : ]:
    pointing_name = pointing_name.split( '_' )[ 0 ]
    pointing_file_names = glob.glob( '%s/%s*.RES.FITS' % 
        ( pointing_fits_path, pointing_name ) )
    if ( len( pointing_file_names ) == 0 ):
      if ( not allow_missing_pointings ):
        raise error( 'pointing %s has no image' % pointing_name )
        if pointings.exists():
          remove_facets( pointings )
    for pointing_file_name in pointing_file_names:
      i = i + 1
      pointing = get_facet( pointings, i )
      read_fits_image( pointing_file_name, pointing )
  pointing_count = i
  if ( pointing_count == 0 ):
    raise error( 'no pointings found' )
  store_parameter( pointings, 'facet_count', pointing_count )
  
  # get some image statistics
  radecs = []
  beam_sizes = []
  noises = []
  for i in range( 1, 1 + pointing_count ):
    pointing = get_facet( pointings, i )
    radec = get_radec( pointing )
    radecs.append( radec )
    beam_size = get_beam_size( pointing )
    beam_sizes.append( beam_size )
    [ avg, noise ] = call_aips_task( 'IMEAN', indata = pointing,
        pixavg = 0., pixstd = 5.e-3, pixrange = [ -50.e-3, 50.e-3  ], 
        outputs = [ 'pixavg', 'pixstd' ] )
    noises.append( noise )
  radecs = array( radecs )
  beam_sizes = array( beam_sizes )
  noises = array( noises )
  
  # define restoring beam (and make beam plot)
  gmrt_lonlat = dmsdms_to_degdeg( [ 74,2,59.07, 19,5,47.46 ] )
  if ( restoring_beam[ 1 ] > 0. ):
    [ bmaj, bmin, bpa ] = restoring_beam
    mosaic_beam_size = restoring_beam
  else:
    bmin = restoring_beam[ 0 ]
    bmaj = around( bmin / cos( radians( mosaic_radec[ 1 ] - gmrt_lonlat[ 1 ] ) ), 1 )
    mosaic_beam_size = [ bmaj, bmin, 0. ]
  if make_beam_plot:
    dec_min = radecs[ : , 1 ].min()
    dec_max = radecs[ : , 1 ].max()
    decs = array( [ dec_min + ( dec_max - dec_min ) * ( float( x ) / 1000. )
        for x in range( 1001 ) ] )
    plot( decs, 0.9 * bmin / cos( radians( decs - gmrt_lonlat[ 1 ] ) ), 'b' )
    plot( decs, 1.0 * bmin / cos( radians( decs - gmrt_lonlat[ 1 ] ) ), 'g' )
    plot( decs, 1.1 * bmin / cos( radians( decs - gmrt_lonlat[ 1 ] ) ), 'r' )
    plot( radecs[ : , 1 ], beam_sizes[ : , 0 ], 'k.' )
    show()
    clear_aips_disk( aips_disk )
    return
  
  # define mosaic parameters
  mosaic_pixel_size = around( mosaic_beam_size[ 1 ] / pixels_per_beam, 1 )
  if ( mosaic_size[ 1 ] <= 0. ):
    mosaic_image_size = [ int( around( mosaic_size[ 0 ] * 3600. / mosaic_pixel_size ) ),
        int( around( mosaic_size[ 0 ] * 3600. / mosaic_pixel_size ) ) ]
  else:
    mosaic_image_size = [ int( around( mosaic_size[ 0 ] * 3600. / mosaic_pixel_size ) ),
        int( around( mosaic_size[ 1 ] * 3600. / mosaic_pixel_size ) ) ]
  if ( max( mosaic_image_size ) > max_image_size ):
    raise error( 'image size %s exceeds maximum' % ( repr( mosaic_image_size ) ) )
    remove_facets( pointings )
  
  # calculate restoring beam position angle per pointing
  pa_list = []
  for i in range( 1, 1 + pointing_count ):
    pointing = get_facet( pointings, i )
    pointing_radec = get_radec( pointing )
    p = amodulo( pointing_radec[ 0 ] - mosaic_radec[ 0 ] + 180., 360. ) - 180.
    p = asign( pointing_radec[ 1 ] ) * ( 1. - cos( aradians( pointing_radec[ 1 ] ) ) ) * p
    pa = amodulo( mosaic_beam_size[ 2 ] + p + 90., 180. ) - 90.
    pa_list.append( pa )
  
  # convolve pointings to mosaic resolution
  if convolve_background:
    for i in range( 1, 1 + pointing_count ):
      pointing = get_facet( pointings, i )
      convl_pointing = get_aips_file( aips_disk, mosaic_name[ 0 : 12 ], 'CONVL', -1, 'MA' )
      beam_size = get_beam_size( pointing )
      factor = ( ( mosaic_beam_size[ 0 ] * mosaic_beam_size[ 1 ] ) / 
          ( beam_size[ 0 ] * beam_size[ 1 ] ) )
      if ( ( beam_size[ 0 ] < mosaic_beam_size[ 0 ] ) and
          ( beam_size[ 1 ] < mosaic_beam_size[ 1 ] ) ):
        try:
          call_aips_task( 'CONVL', indata = pointing, outdata = convl_pointing,
              bmaj = mosaic_beam_size[ 0 ], bmin = mosaic_beam_size[ 1 ],
              bpa = pa_list[ i - 1 ], opcode = 'GAUS', doblank = 1, factor = factor )
          pixels = get_image_pixels( convl_pointing )
          if alltrue( pixels == get_aips_magic_value() ):
            raise error( 'convolution failed' )
        except:
          beam_size[ 0 ] = mosaic_beam_size[ 0 ]
      if ( ( beam_size[ 0 ] >= mosaic_beam_size[ 0 ] ) and
          ( beam_size[ 1 ] < mosaic_beam_size[ 1 ] ) ):
        if convl_pointing.exists():
          convl_pointing.zap()
        try:
          call_aips_task( 'CONVL', indata = pointing, outdata = convl_pointing,
              bmaj = 1.001 * beam_size[ 0 ], bmin = mosaic_beam_size[ 1 ],
              bpa = pa_list[ i - 1 ], opcode = 'GAUS', doblank = 1, factor = factor )
          pixels = get_image_pixels( convl_pointing )
          if alltrue( pixels == get_aips_magic_value() ):
            raise error( 'convolution failed' )
        except:
          beam_size[ 1 ] = mosaic_beam_size[ 1 ]
      if ( ( beam_size[ 0 ] >= mosaic_beam_size[ 0 ] ) and
          ( beam_size[ 1 ] >= mosaic_beam_size[ 1 ] ) ):
        if convl_pointing.exists():
          convl_pointing.zap()
        # should we still scale?
        pixels = get_image_pixels( pointing )
        sel = awhere( pixels != get_aips_magic_value() )
        pixels = aput( pixels, sel, aget( pixels, sel ) * factor )
        set_image_pixels( pointing, pixels )
      if convl_pointing.exists():
        pointing.zap()
        convl_pointing.rename( name = pointing.name, klass = pointing.klass,
            seq = pointing.seq )
      set_beam_size( pointing, mosaic_beam_size[ 0 : 2 ] + [ pa_list[ i - 1 ] ] )
  
  # restore CLEAN components
  for i in range( 1, 1 + pointing_count ):
    pointing = get_facet( pointings, i )
    pixels = get_image_pixels( pointing )
    sel = awhere( pixels == get_aips_magic_value() )
    rst_pointing = restore_model_components( pointings, facet_list = [ i ],
        cross_restore = False, imagr_params = { 'bmaj' : mosaic_beam_size[ 0 ], 
        'bmin' : mosaic_beam_size[ 1 ], 'bpa' : pa_list[ i - 1 ] } )
    pixels = get_image_pixels( rst_pointing )
    pixels = aput( pixels, sel, get_aips_magic_value() )
    set_image_pixels( rst_pointing, pixels )
    pointing = get_facet( pointings, i )
    pointing.zap()
    rst_pointing.rename( name = pointing.name, klass = pointing.klass,
        seq = pointing.seq )
  
  # resample to final resolution
  # convert to galactic, if requested
  for i in range( 1, 1 + pointing_count ):
    pointing = get_facet( pointings, i )
    if do_galactic:
      convert_to_galactic( pointing )
    scale_pointing = scale_image( pointing, rotate = True,
        new_pixel_size = [ mosaic_pixel_size, mosaic_pixel_size ] )
    pointing.zap()
    scale_pointing.rename( name = pointing.name, klass = pointing.klass,
        seq = pointing.seq )
    if do_galactic:
      for ctype in pointing.header.ctype:
        if ( ctype.find( 'GLON' ) != - 1 ):
          ra_index = pointing.header.ctype.index( ctype )
        if ( ctype.find( 'GLAT' ) != - 1 ):
          dec_index = pointing.header.ctype.index( ctype )
      wiz_im = wizardry( pointing )
      wiz_im.header.ctype[ ra_index ] = 'RA--' + wiz_im.header.ctype[ ra_index ][ 4 : ]
      wiz_im.header.ctype[ dec_index ] = 'DEC-' + wiz_im.header.ctype[ dec_index ][ 4 : ]
      wiz_im.header.update()
      wiz_im = wizardry( pointing )
      del wiz_im
      set_epoch( pointing, 2000. )
  
  # select primary beam model parameters
  frequency = get_frequency( pointings )
  if ( pbparms is None ):
    if ( frequency > 140.e6 ) and ( frequency < 170.e6 ):
      pbparms = [ 0., 1., -4.04, 76.2, -68.8, 22.03, 0. ]
    elif ( frequency > 200.e6 ) and ( frequency < 280.e6 ):
      pbparms = [ 0., 1., -3.366, 46.159, -29.963, 7.529, 0. ]
    elif ( frequency > 300.e6 ) and ( frequency < 350.e6 ):
      pbparms = [ 0., 1., -3.397, 47.192, -30.931, 7.803, 0. ] 
    elif ( frequency > 580.e6 ) and ( frequency < 650.e6 ):
      pbparms = [ 0., 1., -3.486, 47.749, -35.203, 10.399, 0. ] 
    elif ( frequency > 1.e9 ) and ( frequency < 1.8e9 ):
      pbparms = [ 0., 1., -2.27961, 21.4611, -9.7929, 1.80153, 0. ] 
    else:
      raise error( 'beam shape is unknown' )
  else:
    pbparms[ 0 ] = 0.
    pbparms[ 1 ] = 1.
  pb_cutoff = calculate_pbparm_attenuation( frequency, beam_radius, pbparms )
  pbparms[ 0 ] = pb_cutoff
  
  # make RMS maps
  rms_pointings = get_aips_file( aips_disk, 'RMS', pointings.klass, -1, 'MA' )
  for i in range( 1, 1 + pointing_count ):
    pointing = get_facet( pointings, i )
    rms_pointing = get_facet( rms_pointings, i )
    call_aips_task( 'RMSD', indata = pointing, outdata = rms_pointing,
        imsize = [ 91, -1 ], xinc = 20, yinc = 20, optype = 'ROBS' )
  
  # replace (near) zeros with blanks
  for i in range( 1, 1 + pointing_count ):
    pointing = get_facet( pointings, i )
    rms_pointing = get_facet( rms_pointings, i )
    pixels = get_image_pixels( pointing )
    rms_pixels = get_image_pixels( rms_pointing )
    sel = awhere( ( pixels == 0. ) | 
        ( rms_pixels < 0.25 * median( rms_pixels ) ) )
    pixels = aput( pixels, sel, get_aips_magic_value() )
    rms_pixels = aput( rms_pixels, sel, get_aips_magic_value() )
    set_image_pixels( pointing, pixels )
    set_image_pixels( rms_pointing, rms_pixels )
  del pixels
  del rms_pixels
  
  # make mask maps
  pbcor_pointings = get_aips_file( aips_disk, 'PBCOR', pointings.klass, -1, 'MA' )
  mask_pointings = get_aips_file( aips_disk, 'ONES', pointings.klass, -1, 'MA' )
  for i in range( 1, 1 + pointing_count ):
    pointing = get_facet( pointings, i )
    pbcor_pointing = get_facet( pbcor_pointings, i )
    center = get_radec( pointing )
    coordinates = degdeg_to_hmsdms( center )
    if ( center[ 1 ] < 0. ):
      coordinates[ 3 ] = -abs( coordinates[ 3 ] )
      coordinates[ 4 ] = -abs( coordinates[ 4 ] )
      coordinates[ 5 ] = -abs( coordinates[ 5 ] )
    call_aips_task( 'PBCOR', indata = pointing, outdata = pbcor_pointing,
        pbparm = pbparms, coordina = coordinates )
  for i in range( 1, 1 + pointing_count ):
    pbcor_pointing = get_facet( pbcor_pointings, i )
    mask_pointing = get_facet( mask_pointings, i )
    units = pbcor_pointing.header.bunit
    call_aips_task( 'MATHS', indata = pbcor_pointing, outdata = mask_pointing,
        opcode = 'POLY', cparm = [ 1,0,0,0 ] )
    set_header_keyword( mask_pointing, 'bunit', units )
  
  if ( pointing_model is None ):
    
    # make PB-corrected RMS maps
    pbrms_pointings = get_aips_file( aips_disk, 'PBRMS', pointings.klass, -1, 'MA' )
    for i in range( 1, 1 + pointing_count ):
      rms_pointing = get_facet( rms_pointings, i )
      pbrms_pointing = get_facet( pbrms_pointings, i )
      center = get_radec( rms_pointing )
      coordinates = degdeg_to_hmsdms( center )
      if ( center[ 1 ] < 0. ):
        coordinates[ 3 ] = -abs( coordinates[ 3 ] )
        coordinates[ 4 ] = -abs( coordinates[ 4 ] )
        coordinates[ 5 ] = -abs( coordinates[ 5 ] )
      call_aips_task( 'PBCOR', indata = rms_pointing, outdata = pbrms_pointing,
          pbparm = pbparms, coordina = coordinates )
  
  else:
    # make shifted, PB-corrected image and RMS maps
    # pointing_model: [ arcsec/deg DEC, DEC_0 ]
    remove_facets( pbcor_pointings )
    pbrms_pointings = get_aips_file( aips_disk, 'PBRMS', pointings.klass, -1, 'MA' )
    for i in range( 1, 1 + pointing_count ):
      rms_pointing = get_facet( rms_pointings, i )
      pbrms_pointing = get_facet( pbrms_pointings, i )
      center = get_radec( rms_pointing )
      dec_shift = pointing_model[ 0 ] * ( center[ 1 ] - pointing_model[ 1 ] ) / 3600.
      center[ 1 ] = center[ 1 ] + dec_shift
      pb_cutoff = calculate_pbparm_attenuation( frequency, 
          beam_radius + abs( dec_shift ), [ 0. ] + pbparms[ 1 : ] )
      coordinates = degdeg_to_hmsdms( center )
      if ( center[ 1 ] < 0. ):
        coordinates[ 3 ] = -abs( coordinates[ 3 ] )
        coordinates[ 4 ] = -abs( coordinates[ 4 ] )
        coordinates[ 5 ] = -abs( coordinates[ 5 ] )
      call_aips_task( 'PBCOR', indata = rms_pointing, outdata = pbrms_pointing,
          pbparm = [ pb_cutoff ] + pbparms[ 1 : ], coordina = coordinates )
      pointing = get_facet( pointings, i )
      pbcor_pointing = get_facet( pbcor_pointings, i )
      call_aips_task( 'PBCOR', indata = pointing, outdata = pbcor_pointing,
          pbparm = [ pb_cutoff ] + pbparms[ 1 : ], coordina = coordinates )

    # mask outer regions
    masked_pbrms_pointings = get_aips_file( aips_disk, 'MPBRMS', pointings.klass, -1, 'MA' )
    masked_pbcor_pointings = get_aips_file( aips_disk, 'MPBCOR', pointings.klass, -1, 'MA' )
    for i in range( 1, 1 + pointing_count ):
      mask_pointing = get_facet( mask_pointings, i )
      pbrms_pointing = get_facet( pbrms_pointings, i )
      masked_pbrms_pointing = get_facet( masked_pbrms_pointings, i )
      units = pbrms_pointing.header.bunit
      call_aips_task( 'COMB', indata = mask_pointing, in2data = pbrms_pointing,
          outdata = masked_pbrms_pointing, opcode = 'MULT', aparm = [ 1,0,0,0,0,0,0,0,0,0 ] )
      set_header_keyword( masked_pbrms_pointing, 'bunit', units )
      pbcor_pointing = get_facet( pbcor_pointings, i )
      masked_pbcor_pointing = get_facet( masked_pbcor_pointings, i )
      units = pbcor_pointing.header.bunit
      call_aips_task( 'COMB', indata = mask_pointing, in2data = pbcor_pointing,
          outdata = masked_pbcor_pointing, opcode = 'MULT', aparm = [ 1,0,0,0,0,0,0,0,0,0 ] )
      set_header_keyword( masked_pbcor_pointing, 'bunit', units )
    remove_facets( pbrms_pointings )
    pbrms_pointings = masked_pbrms_pointings
    remove_facets( pbcor_pointings )
    pbcor_pointings = masked_pbcor_pointings
  
  # make weight maps
  weight_pointings = get_aips_file( aips_disk, 'WEIGHT', pointings.klass, -1, 'MA' )
  for i in range( 1, 1 + pointing_count ):
    pbrms_pointing = get_facet( pbrms_pointings, i )
    weight_pointing = get_facet( weight_pointings, i )
    units = pbrms_pointing.header.bunit
    call_aips_task( 'MATHS', indata = pbrms_pointing, outdata = weight_pointing,
        opcode = 'POWR', cparm = [ 0,1,1,0, weight_power ] )
    set_header_keyword( weight_pointing, 'bunit', units )
  remove_facets( pbrms_pointings )
  
  # make field maps
  field_pointings = get_aips_file( aips_disk, 'FIELD', pointings.klass, -1, 'MA' )
  for i in range( 1, 1 + pointing_count ):
    pbcor_pointing = get_facet( pbcor_pointings, i )
    weight_pointing = get_facet( weight_pointings, i )
    field_pointing = get_facet( field_pointings, i )
    units = pbcor_pointing.header.bunit
    call_aips_task( 'COMB', indata = pbcor_pointing, in2data = weight_pointing,
        outdata = field_pointing, opcode = 'MULT', aparm = [ 1,0,0,0,0,0,0,0,0,0 ] )
    set_header_keyword( field_pointing, 'bunit', units )
  remove_facets( pbcor_pointings )
  
  # combine field images and weight images
  coordinates = degdeg_to_hmsdms( mosaic_center )
  if ( mosaic_center[ 1 ] < 0. ):
    coordinates[ 3 ] = -abs( coordinates[ 3 ] )
    coordinates[ 4 ] = -abs( coordinates[ 4 ] )
    coordinates[ 5 ] = -abs( coordinates[ 5 ] )
  reference_x = int( floor( float( mosaic_image_size[ 0 ] + 1 ) / 2. ) )
  reference_y = int( ceil( float( mosaic_image_size[ 1 ] + 1 ) / 2. ) )
  field_mosaic = get_aips_file( aips_disk, mosaic_name[ 0 : 12 ], 'FLDMOS', -1, 'MA' )
  call_aips_task( 'FLATN', indata = field_pointings, nfield = pointing_count,
      outdata = field_mosaic, imsize = mosaic_image_size, coordina = coordinates,
      cooref = [ reference_x, reference_y ] )
  remove_facets( field_pointings )
  weight_mosaic = get_aips_file( aips_disk, mosaic_name[ 0 : 12 ], 'WGTMOS', -1, 'MA' )
  call_aips_task( 'FLATN', indata = weight_pointings, nfield = pointing_count,
      outdata = weight_mosaic, imsize = mosaic_image_size, coordina = coordinates,
      cooref = [ reference_x, reference_y ] )
  remove_facets( weight_pointings )
  
  # make mosaic image
  mosaic = get_aips_file( aips_disk, mosaic_name[ 0 : 12 ], 'MOSAIC', -1, 'MA' )
  units = field_mosaic.header.bunit
  call_aips_task( 'COMB', indata = field_mosaic, in2data = weight_mosaic,
      outdata = mosaic, opcode = 'DIV', aparm = [ 1,0,0,0,0,0,0,0,0,0 ] )
  set_header_keyword( mosaic, 'bunit', units )
  field_mosaic.zap()
  weight_mosaic.zap()
  
  # write mosaic image to disk
  dummy = get_source_names( mosaic )[ 0 ]
  change_source_name( mosaic, dummy, mosaic_name )
  mosaic_beam_size[ 2 ] = amodulo( mosaic_beam_size[ 2 ] + rotation + 90., 180. ) - 90.
  set_beam_size( mosaic, mosaic_beam_size )
  set_observer( mosaic, 'TGSS' )
  if do_galactic:
    for ctype in mosaic.header.ctype:
      if ( ctype.find( 'RA--' ) != - 1 ):
        ra_index = mosaic.header.ctype.index( ctype )
      if ( ctype.find( 'DEC-' ) != - 1 ):
        dec_index = mosaic.header.ctype.index( ctype )
    wiz_im = wizardry( mosaic )
    wiz_im.header.ctype[ ra_index ] = 'GLON' + wiz_im.header.ctype[ ra_index ][ 4 : ]
    wiz_im.header.ctype[ dec_index ] = 'GLAT' + wiz_im.header.ctype[ dec_index ][ 4 : ]
    wiz_im.header.update()
    wiz_im = wizardry( mosaic )
    del wiz_im
  clear_history( mosaic )
  write_fits_image( mosaic, fits_path + mosaic_name + '.MOSAIC.FITS' )
  mosaic.zap()
  
  clear_aips_disk( aips_disk )
  return

###############################################################################

###############################################################################

# import Python modules
from os import *
from sys import *
from math import *
from string import *
from random import *

# import 3rd party modules
from numpy import *

# import user modules
from spam import *

###############################################################################

def pre_calibrate_targets_vla( uvfits_file_name,
    target_names = [], bw_smearing = 0.7, ta_smearing = 0.35, 
    keep_channel_one = False, do_tb_sort = False, force_date = False,
    force_stokes = True, force_frequency = True, make_plots = True,
    bpcal_names = [ '3C48','3C147','3C286' ], bpcal_snr_power = 1.,
    bpcal_per_scan = True, initial_reference_antenna = 0, amplitude_limit = 20.,
    channel_range = None ):
  
  # get full paths to data directories
  dat_path = path.expandvars( '${DAT}/' )
  prt_path = path.expandvars( '${PRT}/' )
  fits_path = path.expandvars( '${FIT}/' )
  
  # get list of UVFITS files
  if ( type( uvfits_file_name ) == type( [] ) ):
    uvfits_file_names = []
    for file_name in uvfits_file_name:
      uvfits_file_names = uvfits_file_names + glob.glob( file_name )
  else:
    uvfits_file_names = glob.glob( uvfits_file_name )
  if ( len( uvfits_file_names ) == 0 ):
    raise error( 'UVFITS file(s) %s does not exist' % ( repr( uvfits_file_name ) ) )
  uvfits_file_names.sort()
  
  # find free AIPS disk
  aips_disk = 0
  for i in range( 1, len( AIPS.disks ) ):
    if ( len( AIPSCat( i )[ i ] ) == 0 ):
      aips_disk = i
      break
  if ( aips_disk == 0 ):
    raise error( 'no free AIPS disk available' )
  
  # read UVFITS files into AIPS
  multi_uv_list = []
  date_list = []
  stokes_list = []
  freqs_list = []
#  project_name = 'VLA'
  for uvfits_file_name in uvfits_file_names:
    
    # read UV data
    multi_uv = get_aips_file( aips_disk, project_name, 'UV', -1, 'UV' )
    read_fits_uv( uvfits_file_name, multi_uv )
    if ( not multi_uv.header.telescop in [ 'VLA' ] ):
      raise error( 'UVFITS file %s is not a VLA observation' % ( uvfits_file_name ) )
    multi_uv_list.append( multi_uv )
    date_list.append( get_reference_date( multi_uv ).replace( '-', '' ) )
    freqs_list.append( get_frequency_list( multi_uv ) )
    stokes = get_stokes( multi_uv )
    stokes_string = ''
    for stoke in stokes:
      if ( not stoke in [ 'RR', 'LL', 'RL', 'LR' ] ):
        raise error( 'UVFITS file %s contains unsupported stokes' % ( uvfits_file_name ) )
      stokes_string = stokes_string + stoke
    stokes_list.append( stokes_string )
  
  # check meta data
  if ( len( multi_uv_list ) > 1 ):
    for i in range( 1, len( multi_uv_list ) ):
      if ( force_date and ( date_list[ i ] != date_list[ 0 ] ) ):
        raise error( 'UVFITS files are from different observing days' )
      if ( force_stokes and ( stokes_list[ i ] != stokes_list[ 0 ] ) ):
        raise error( 'UVFITS files have different stokes definitions' )
      if ( force_frequency and ( freqs_list[ i ] != freqs_list[ 0 ] ) ):
        raise error( 'UVFITS files have different frequency definitions' )
  extension = '_' + date_list[ 0 ]
  extension = extension + '_' + stokes_list[ 0 ]
  if ( freqs_list[ 0 ][ -1 ] > freqs_list[ 0 ][ 0 ] ):
    extension = extension + '_USB'
  else:
    extension = extension + '_LSB'
  
  # concatenate data
  multi_uv = multi_uv_list[ 0 ]
  if ( len( multi_uv_list ) > 1 ):
    for i in range( 1, len( multi_uv_list ) ):
      multi_uv1 = multi_uv
      multi_uv2 = multi_uv_list[ i ]
      multi_uv = get_aips_file( aips_disk, project_name, 'UV', -1, 'UV' )
      call_aips_task( 'DBCON', indata = multi_uv1, in2data = multi_uv2,
          outdata = multi_uv, doarray = 1 )
      multi_uv1.zap()
      multi_uv2.zap()
  fix_su_table( multi_uv )
  
  # fix calibrator names
  change_source_name( multi_uv, '0137+331', '3C48' )
  change_source_name( multi_uv, '0542+498', '3C147' )
  change_source_name( multi_uv, '1331+305', '3C286' )
  
  # determine processing parameters
  integration_time = around( find_integration_time( multi_uv ), decimals = 2 )
  dfrequency = get_channel_width( multi_uv )
  frequency = mean( freqs_list[ 0 ] )
  channel_count = len( freqs_list[ 0 ] )
  if keep_channel_one:
    channel_increment = int( round( ( freqs_list[ 0 ][ 0 ] / dfrequency ) * 
        bw_smearing * ( ( 2. * 25. ) / 35.e3 ) ) )
  else:
    channel_increment = int( round( ( frequency / dfrequency ) * 
        bw_smearing * ( ( 2. * 25. ) / 35.e3 ) ) )
  channel_increment = max( channel_increment, 1 )
  new_channel_count = channel_count / channel_increment
  dchannel_count = channel_count - new_channel_count * channel_increment
  if keep_channel_one:
    channel_low = 1
    channel_high = new_channel_count * channel_increment
  else:
    if ( not channel_range is None ):
      channel_count = channel_range[ 1 ] - channel_range[ 0 ] + 1
      new_channel_count = channel_count / channel_increment
      dchannel_count = channel_count - new_channel_count * channel_increment
      channel_high = channel_range[ 1 ] - dchannel_count / 2
    else:
      channel_high = channel_count - dchannel_count / 2
    channel_low = channel_high - new_channel_count * channel_increment + 1
  central_channel = ( channel_high + channel_low ) / 2
  central_channels = [ max( central_channel - 4, channel_low ), 
      min( central_channel + 4, channel_high ) ]
  time_increment = max( 1, int( round( ( 1. / integration_time ) *
      ta_smearing * ( ( 2. * 25. ) / 35.e3 ) * ( 86400. / ( 2. * pi ) ) ) ) )
  if ( frequency > 300.e6 ) and ( frequency < 400.e6 ):
    uv_range = [ 1.,1.e6 ]
    bp_dfluxes = [ 12.,7.,5. ]
    clip_level = 40.
    extension = '_VLA325' + extension
  else:
    raise error( 'unsupported frequency' )
  
  # put data in TB order
  if do_tb_sort:
    sort_uv = get_aips_file( aips_disk, project_name, 'TBSORT', -1, 'UV' )
    call_aips_task( 'UVSRT', indata = multi_uv, outdata = sort_uv, sort = 'TB' )
    multi_uv.zap()
    sort_uv.rename( name = multi_uv.name, klass = multi_uv.klass, seq = multi_uv.seq )
  
  # recalculate UV coordinates
  an_table = wizardry( multi_uv ).table( 'AN', 0 )
  an_table.keywords[ 'DATUTC' ] = an_table.keywords[ 'IATUTC' ]
  an_table.close()
#  fix_uv = get_aips_file( aips_disk, project_name, 'UVCOR', -1, 'UV' )
#  call_aips_task( 'UVFIX', indata = multi_uv, outdata = fix_uv )
#  multi_uv.zap()
#  fix_uv.rename( name = multi_uv.name, klass = multi_uv.klass, seq = multi_uv.seq )
  
  # get overview of raw UV content
  while table_exists( multi_uv, 'CL', 0 ):
    multi_uv.zap_table( 'CL', 0 )
  while table_exists( multi_uv, 'NX', 0 ):
    multi_uv.zap_table( 'NX', 0 )
  call_aips_task( 'INDXR', indata = multi_uv, bparm = [ -1,-1 ],
      cparm = [ 0,9999, integration_time / 60., 0,0,0 ] )
  listr_file_name = '%s%s.LISTR' % ( prt_path, aips_file_name_to_string( multi_uv ) )
  if file_exists( listr_file_name ):
    remove_file( listr_file_name )
  call_aips_task( 'LISTR', indata = multi_uv, optype = 'SCAN', docrt = -1,
      outprint = listr_file_name )
  prtan_file_name = '%s%s.PRTAN' % ( prt_path, aips_file_name_to_string( multi_uv ) )
  if file_exists( prtan_file_name ):
    remove_file( prtan_file_name )
  call_aips_task( 'PRTAN', indata = multi_uv, docrt = -1, outprint = prtan_file_name )
  
  # create separate UV files for different sources
  source_names = get_source_names( multi_uv )
  if ( len( target_names ) > 0 ):
    new_source_names = []
    for source in source_names:
      if ( source in source_names + bpcal_names ):
        new_source_names.append( source )
    source_names = new_source_names
  for i in range( len( source_names ) )[ : : 30 ]:
    try:
      call_aips_task( 'SPLAT', indata = multi_uv, smooth = [ 0 ], douvcomp = 0,
          aparm = [ 3,0,1,0,0,0,1 ], channel = 1, chinc = 1, flagver = 0,
          sources = source_names[ i : min( i + 30, len( source_names ) ) ], 
          doband = -1, stokes = '', docalib = -1,
          outdisk = multi_uv.disk, outclass = multi_uv.klass, outseq = multi_uv.seq )
    except RuntimeError:
      print 'WARNING: sources not found in multi-UV file: %s' % ( 
          repr( source_names[ i : min( i + 30, len( source_names ) ) ] ) )
  
  # process bandpass calibrators
  bpcal_snrs = []
  for bpcal_name in bpcal_names:
    try:
      main_bp_uv = get_aips_file( aips_disk, bpcal_name, 'UV', 0, 'UV' )
    except:
      continue
    
    # initial flagging
    flags_bpcal = [
        { 'bchan' : 1, 'echan' : 1 },
        { 'antennas' : [ 31,32 ] }]
    add_flags( main_bp_uv, flags = flags_bpcal, flag_version = 1 )
    
    # process per scan
    scan_list = find_scans( main_bp_uv )
    if ( ( not bpcal_per_scan ) and ( len( scan_list ) > 1 ) ):
      scan_list = [ [ scan_list[ 0 ][ 0 ], scan_list[ -1 ][ 1 ] ] ]
    for scan_id, scan_times in enumerate( scan_list ):
      
      try:
        bp_uv = get_aips_file( aips_disk, bpcal_name, 'UV', -1, 'UV' )
        time_range = timetime_to_dhmsdhms( scan_times )
        call_aips_task( 'SPLAT', indata = main_bp_uv, docalib = -1, doband = -1,
            flagver = 0, blver = -1, douvcomp = 0, aparm = [ 3,0,1,0,0,0,0 ],
            stokes = '', outdata = bp_uv, timerang = time_range )
        
        # calculate calibrator flux
        frequency = mean( get_frequency_list( bp_uv )[
            central_channels[ 0 ] - 1 : central_channels[ 1 ] ] )
        [ bp_flux, bp_spix ] = get_SH_model_flux_spix( bpcal_name, frequency )
        
        # do short interval calibration against model
        dummy_uv = get_aips_file( bp_uv.disk, bp_uv.name, 'DUMMY', -1, 'UV' )
        call_aips_task( 'CALIB', indata = bp_uv, flagver = 0, smooth = [ 0 ],
            ichansel = [ central_channels + [ 1,1 ] ], weightit = 1,
            cmethod = 'DFT', smodel = [ bp_flux, 0,0,0 ], uvrang = uv_range,
            solint = 1. / 60., soltype = 'L1', solmode = 'A&P', outdata = dummy_uv,
            snver = 0, aparm = [ 0,0,0,0,0,0,1.,0,0,0 ], docalib = -1,
            refant = initial_reference_antenna )
        dummy_uv.zap()
        
        # filter gain amplitudes and phases
        filter_solutions( bp_uv, amplitude_limits = [ 0., amplitude_limit ],
            amplitude_sigma = 10., phase_sigma = 10. )
        filter_solutions( bp_uv, amplitude_sigma = 3., phase_sigma = None,
            amplitude_window = 60. * 60., phase_window = None )
        filter_solutions( bp_uv, amplitude_sigma = None, phase_sigma = 3.,
            amplitude_window = None, phase_window = 5. * 60. )
        filter_dead_antennas( bp_uv )
        reference_antenna = select_reference_antenna( bp_uv, 
            antenna_list = [ 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,25 ] )  
        re_reference_solutions( bp_uv, reference_antenna )
        
        # calibrate bandpass
        call_aips_task( 'BPASS', indata = bp_uv, flagver = 0, dopol = -1, blver = -1,
            docalib = 100, gainuse = 0, smodel = [ bp_flux,0,0,0 ], doband = -1,
            solint = -1, soltype = 'L1', outvers = 0, refant = reference_antenna,
            bpassprm = [ 0,0,0,0,2,0,0,1,1,3,1 ], cmethod = 'DFT', specindx = 1.e-3,
            ichansel = [ central_channels + [ 1,1 ] ], uvrang = uv_range )
        
        # generate calibrated, flagged, collapsed data for phase calibrator
        # apply amplitude and bandpass calibration
        smooth_solutions_in_time( bp_uv, phase_window = 20. * 60. * 60., 
            amplitude_window = 20. * 60. * 60., full_polarization = True,
            order = 0, max_gap = 10. )
        bp_cal_uv = get_aips_file( aips_disk, bpcal_name, 'CAL', -1, 'UV' )
        call_aips_task( 'SPLAT', indata = bp_uv, docalib = 100, gainuse = 0,
            doband = 5, bpver = 0, flagver = 0, blver = -1, douvcomp = 0,
            aparm = [ 3,0,1,0,0,0,0 ], stokes = '', outdata = bp_cal_uv,
            bchan = channel_low, echan = channel_high )
        
        # recalculate calibrator flux
        frequency = mean( get_frequency_list( bp_cal_uv ) )
        [ bp_flux, bp_spix ] = get_SH_model_flux_spix( bpcal_name, frequency )
        
        # clip outliers and apply flags
        call_aips_task( 'CLIP', indata = bp_cal_uv, docalib = -1, gainuse = 0,
            doband = -1, blver = 0, flagver = 1, outfgver = 1, stokes = '',
            aparm = [ 3. * bp_flux, 0,0,0,0,0,0,0,0,0 ] )
        temp_uv = get_aips_file( bp_uv.disk, bp_uv.name, 'TEMP', -1, 'UV' )
        call_aips_task( 'UVCOP', indata = bp_cal_uv, outdata = temp_uv, flagver = 1,
            uvcopprm = [ 0,0,0,0,0,3,0 ] )
        
        # calibration and clipping
        for bp_dflux in bp_dfluxes:
          
          # do short interval phase calibration against model
          call_aips_task( 'CALIB', indata = temp_uv, flagver = 0, smooth = [ 0 ],
              refant = reference_antenna, smodel = [ bp_flux, 0,0,0 ], blver = -1,
              cmethod = 'DFT', solint = 1. / 60., soltype = 'L1', solmode = 'A&P',
              outdata = dummy_uv, snver = 0, doband = -1, bpver = 0, uvrang = uv_range,
              aparm = [ 4,0,0,0,0,0,1.,0,0,0 ], docalib = -1, gainuse = 0,
              ichansel = [ [ 1, get_channel_count( temp_uv ), 1,1 ] ], weightit = 1 )
          dummy_uv.zap()
          
          # clip outliers
          bp_flux_max = bp_flux + bp_dflux
          bp_flux_min = max( 0., bp_flux - bp_dflux )
          call_aips_task( 'CLIP', indata = temp_uv, docalib = 100, gainuse = 0,
              doband = -1, blver = -1, flagver = 1, outfgver = 1, stokes = '',
              aparm = [ bp_flux_max, 0, bp_flux_min, 0,0,0,0,0,0,0 ] )
          
          # apply flags
          temp_uv.zap_table( 'SN', 0 )
          temp_uv2 = get_aips_file( bp_uv.disk, bp_uv.name, 'TEMP', -1, 'UV' )
          call_aips_task( 'UVCOP', indata = temp_uv, outdata = temp_uv2, flagver = 1,
              uvcopprm = [ 0,0,0,0,0,3,0 ] )
          temp_uv.zap()
          temp_uv = temp_uv2
        
        bp_cal_uv.zap()
        temp_uv.rename( name = bp_cal_uv.name, klass = bp_cal_uv.klass,
            seq = bp_cal_uv.seq )
        
        # recalculate calibrator flux
        frequency = mean( get_frequency_list( bp_cal_uv ) )
        [ bp_flux, bp_spix ] = get_SH_model_flux_spix( bpcal_name, frequency )
        
        # do short interval phase calibration against model
        call_aips_task( 'CALIB', indata = bp_cal_uv, flagver = -1, smooth = [ 0 ],
            refant = reference_antenna, smodel = [ bp_flux, 0,0,0 ], blver = -1,
            cmethod = 'DFT', solint = 1. / 60., soltype = 'L1', solmode = 'A&P',
            outdata = dummy_uv, snver = 0, doband = -1, bpver = 0, uvrang = uv_range,
            aparm = [ 4,0,0,0,0,0,2.,0,0,0 ], docalib = -1, gainuse = 0,
            ichansel = [ [ 1, get_channel_count( bp_cal_uv ), 1,1 ] ], weightit = 1 )
        dummy_uv.zap()
        
        # calibrate bandpass
        call_aips_task( 'BPASS', indata = bp_cal_uv, flagver = 0, dopol = -1,
            docalib = 100, gainuse = 0, smodel = [ bp_flux, 0,0,0 ], doband = -1,
            solint = -1, soltype = 'L1', outvers = 0, refant = reference_antenna,
            bpassprm = [ 0,0,0,0,2,0,0,1,1,3,1 ], cmethod = 'DFT',
            ichansel = [ central_channels + [ 1,1 ] ], specindx = bp_spix + 0.7, # 0.,
            uvrang = uv_range, blver = -1 )
        
        # determine polarization offsets
        # apply solutions
        # average in time and frequency
        try:
          determine_polarization_solution_offsets( bp_cal_uv, normalize_gains = False )
        except:
          pass
        smooth_solutions_in_time( bp_cal_uv, phase_window = 20. * 60. * 60.,
            amplitude_window = 20. * 60. * 60., full_polarization = True,
            order = 0, max_gap = 20. * 60. )
        re_sample_solutions( bp_cal_uv, interpolation_method = 'nearest',
            gap_time = 20. * 60, force_reference = False, full_polarization = True,
            add_mirror_points = False, weight_multiplier = 0. )
        bp_ical_uv = get_aips_file( aips_disk, bpcal_name, 'ICAL', -1, 'UV' )
        call_aips_task( 'SPLAT', indata = bp_cal_uv, docalib = 100, gainuse = 0,
            doband = 5, bpver = 0, flagver = 0, blver = -1, douvcomp = 0,
            aparm = [ 3,0,1,0,0,0,0 ], stokes = 'I', outdata = bp_ical_uv,
            channel = channel_increment, chinc = channel_increment )
        tavg_uv = time_average_uv_data( bp_ical_uv, time_increment )
        bp_ical_uv.zap()
        tavg_uv.rename( name = bp_ical_uv.name, klass = bp_ical_uv.klass, 
            seq = bp_ical_uv.seq )
        
        # do short interval phase calibration against model
        call_aips_task( 'CALIB', indata = bp_ical_uv, flagver = -1, smooth = [ 0 ],
            refant = reference_antenna, smodel = [ bp_flux, 0,0,0 ], blver = -1,
            cmethod = 'DFT', solint = 1. / 60., soltype = 'L1', solmode = 'A&P',
            outdata = dummy_uv, snver = 0, doband = -1, bpver = 0, uvrang = uv_range,
            aparm = [ 4,0,0,0,0,0,2.,0,0,0 ], docalib = -1, gainuse = 0,
            ichansel = [ [ 1, get_channel_count( bp_ical_uv ), 1,1 ] ], weightit = 1 )
        dummy_uv.zap()
        
        # determine instrumental phase offsets
        # combine amplitudes and phases
        smooth_solutions_in_time( bp_ical_uv, phase_window = 0.,
            amplitude_window = 20. * 60. * 60., full_polarization = True,
            order = 0, max_gap = 30. )
        amp_version = bp_ical_uv.table_highver( 'SN' )
#        time_range = []
        try:
          filter_instrumental_phases( bp_ical_uv, max_period = 10., min_period = 5.,
              max_gap = 5., iterations = 1, height = 300.e3, print_info = True,
              time_step = 1, density_scale = 1.e3, density_power = 1., in_version = 0,
              drop_antennas = False, time_range = [], phase_rms_criterion = None )
        except:
          continue
        re_reference_solutions( bp_ical_uv, reference_antenna )
        phase_version = bp_ical_uv.table_highver( 'SN' )
        call_aips_task( 'SNCOR', indata = bp_ical_uv, snver = amp_version,
            opcode = 'ZPHS' )
        combine_solutions( bp_ical_uv, in_version_1 = amp_version, 
            in_version_2 = phase_version, force_match = True )
        
        # assess quality of calibration
        amplitude_stats = calculate_solution_table_amplitude_median_rms( bp_ical_uv,
            version = 1 )
        total_snr = 0.
        for stat in amplitude_stats:
          if ( stat[ 0 ] > 0. ):
            stat2 = get_robust_mean_deviations( 1. / 
                draw_from_gaussian_distribution( stat[ 0 ], stat[ 1 ] ) )
            amp = stat2[ 0 ]
            damp = max( stat2[ 1 ], -stat2[ 2 ] )
            total_snr = total_snr + sqrt( amp / damp )
          if ( len( stat ) > 2 ):
            if ( stat[ 2 ] > 0. ):
              stat2 = get_robust_mean_deviations( 1. / 
                  draw_from_gaussian_distribution( stat[ 2 ], stat[ 3 ] ) )
              amp = stat2[ 0 ]
              damp = max( stat2[ 1 ], -stat2[ 2 ] )
              total_snr = total_snr + ( amp / damp )**( bpcal_snr_power )
        bpcal_snrs.append( [ bpcal_name, scan_id, total_snr, bp_uv, bp_cal_uv, bp_ical_uv ] )
        
      # skip calibrator scan when an error occurs
      except:
        continue
  
  # select best calibrator scan
  if ( len( bpcal_snrs ) == 0 ):
    raise error( 'no suitable calibrator scans found' )
  for [ bpcal_name, scan_id, total_snr, bp_uv, bp_cal_uv, bp_ical_uv ] in bpcal_snrs:
    print 'calibrator %s scan %d has an SNR of %d' % (
        bpcal_name, scan_id + 1, int( total_snr ) )
  total_snrs = [ x[ 2 ] for x in bpcal_snrs ]
  bpcal_index = total_snrs.index( max( total_snrs ) )
  [ bpcal_name, scan_id, total_snr, bp_uv, bp_cal_uv, bp_ical_uv ] = bpcal_snrs[ bpcal_index ]
  print '... selecting calibrator %s scan %d' % ( bpcal_name, scan_id + 1 )
  
  # make gain amplitude plots
  if make_plots:
    pl_version_low = bp_cal_uv.table_highver( 'PL' ) + 1
    call_aips_task( 'SNPLT', indata = bp_cal_uv, nplots = 8, dotv = -1, factor = 0.3,
        optype = 'AMP', pixrange = [ 0, 2 ], inext = 'SN', invers = 1 )
    pl_version_high = bp_cal_uv.table_highver( 'PL' )
    plot_file_name = prt_path + aips_file_name_to_string( bp_cal_uv ) + '_GAIN_AMP.PS'
    if file_exists( plot_file_name ):
      remove_file( plot_file_name )
    call_aips_task( 'LWPLA', indata = bp_cal_uv, plver = pl_version_low, 
        invers = pl_version_high, lpen = 2, outfile = plot_file_name )
    for pl_version in range( pl_version_low, pl_version_high + 1 ):
      bp_cal_uv.zap_table( 'PL', pl_version )
  
  # make gain phase plots
  if make_plots:
    pl_version_low = bp_cal_uv.table_highver( 'PL' ) + 1
    call_aips_task( 'SNPLT', indata = bp_cal_uv, nplots = 8, dotv = -1, factor = 0.3,
        optype = 'PHAS', pixrange = [ -180,180 ], inext = 'SN', invers = 1 )
    pl_version_high = bp_cal_uv.table_highver( 'PL' )
    plot_file_name = prt_path + aips_file_name_to_string( bp_cal_uv ) + '_GAIN_PHS.PS'
    if file_exists( plot_file_name ):
      remove_file( plot_file_name )
    call_aips_task( 'LWPLA', indata = bp_cal_uv, plver = pl_version_low, 
        invers = pl_version_high, lpen = 2, outfile = plot_file_name )
    for pl_version in range( pl_version_low, pl_version_high + 1 ):
      bp_cal_uv.zap_table( 'PL', pl_version )
  
  # make bandpass plot
  if make_plots:
    pl_version_low = bp_uv.table_highver( 'PL' ) + 1
    call_aips_task( 'POSSM', indata = bp_uv, nplots = 4, dotv = -1, factor = 0.3,
        aparm = [ 0,1,0,2,-180,180,0,2,0 ] )
    pl_version_high = bp_uv.table_highver( 'PL' )
    plot_file_name = prt_path + aips_file_name_to_string( bp_uv ) + '_BANDPASS.PS'
    if file_exists( plot_file_name ):
      remove_file( plot_file_name )
    call_aips_task( 'LWPLA', indata = bp_uv, plver = pl_version_low, 
        invers = pl_version_high, lpen = 2, outfile = plot_file_name )
    for pl_version in range( pl_version_low, pl_version_high + 1 ):
      bp_uv.zap_table( 'PL', pl_version )
  
  # process target fields
  target_names.sort()
  target_uvfits_file_names = []
  for target_name in target_names:
    try:
      t_uv = get_aips_file( aips_disk, target_name, 'UV', 0, 'UV' )
    except:
      continue
    
    # apply amplitude & bandpass calibration
    call_aips_task( 'TACOP', indata = bp_uv, inext = 'BP', ncount = 1, 
        outdata = t_uv )
    call_aips_task( 'TACOP', indata = bp_uv, inext = 'SN', ncount = 1,
        outdata = t_uv )
    re_sample_solutions( t_uv, interpolation_method = 'linear',
        gap_time = 20. * 60., force_reference = False,
        full_polarization = True, add_mirror_points = False )
    add_flags( t_uv, flags = flags_bpcal, flag_version = 1 )
    t_cal_uv = get_aips_file( aips_disk, target_name, 'CAL', -1, 'UV' )
    call_aips_task( 'SPLAT', indata = t_uv, docalib = 100, gainuse = 0,
        doband = 5, bpver = 0, flagver = 0, blver = -1, douvcomp = 0,
        aparm = [ 3,0,1,0,0,0,0 ], stokes = '', outdata = t_cal_uv,
        bchan = channel_low, echan = channel_high )
    
    # apply instrumental calibration
    call_aips_task( 'TACOP', indata = bp_cal_uv, inext = 'BP', ncount = 1, 
        outdata = t_cal_uv )
    call_aips_task( 'TACOP', indata = bp_cal_uv, inext = 'SN', ncount = 1,
        outdata = t_cal_uv )
    re_sample_solutions( t_cal_uv, interpolation_method = 'linear',
        gap_time = 20. * 60, force_reference = False,
        full_polarization = True, add_mirror_points = False )
    t_ical_uv = get_aips_file( aips_disk, target_name, 'ICAL', -1, 'UV' )
    call_aips_task( 'SPLAT', indata = t_cal_uv, docalib = 100, gainuse = 0,
        doband = 5, bpver = 0, flagver = 0, blver = -1, douvcomp = 0,
        aparm = [ 3,0,1,0,0,0,0 ], stokes = '', outdata = t_ical_uv )
    
    # clip visibility amplitudes
    flagver = max( 1, t_ical_uv.table_highver( 'FG' ) )
    call_aips_task( 'CLIP', indata = t_ical_uv, flagver = flagver, stokes = '',
        aparm = [ clip_level, 0,0.001,0,0,0,0,0,0,0 ], outfgver = flagver )
    flag_uv = apply_flag_table( t_ical_uv )
    t_ical_uv.zap()
    flag_uv.rename( name = t_ical_uv.name, klass = t_ical_uv.klass,
        seq = t_ical_uv.seq )
    
    # create time- and frequency-collapsed data
    favg_uv = get_aips_file( aips_disk, target_name, 'FAVG', -1, 'UV' )
    call_aips_task( 'SPLAT', indata = t_ical_uv, docalib = -1, gainuse = 0,
        doband = -1, bpver = 0, flagver = 0, blver = 0, douvcomp = 0,
        aparm = [ 3,0,1,0,0,0,0 ], stokes = 'I', outdata = favg_uv,
        channel = channel_increment, chinc = channel_increment )
    t_ical_uv.zap()
    favg_uv.rename( name = t_ical_uv.name, klass = t_ical_uv.klass,
        seq = t_ical_uv.seq )
    tavg_uv = time_average_uv_data( t_ical_uv, time_increment )
    t_ical_uv.zap()
    tavg_uv.rename( name = t_ical_uv.name, klass = t_ical_uv.klass,
        seq = t_ical_uv.seq )
    
    # apply instrumental phase calibration
    call_aips_task( 'TACOP', indata = bp_ical_uv, inext = 'SN', 
        invers = 0, ncount = 1, outdata = t_ical_uv, outvers = 0 )
    re_sample_solutions( t_ical_uv, interpolation_method = 'linear',
        gap_time = 20. * 60, force_reference = False,
        full_polarization = True, add_mirror_points = False )
    t_iical_uv = get_aips_file( aips_disk, target_name, 'IICAL', -1, 'UV' )
    call_aips_task( 'SPLAT', indata = t_ical_uv, docalib = 100, gainuse = 0,
        doband = -1, bpver = 0, flagver = 0, blver = 0, douvcomp = 0,
        aparm = [ 3,0,1,0,0,0,0 ], stokes = 'I', outdata = t_iical_uv )
    
    # determine and apply Tsys correction
    cor_uv = calculate_haslam_tsys_correction( t_iical_uv, bpcal_name,
        apply_correction = True )
    t_iical_uv.zap()
    cor_uv.rename( name = t_iical_uv.name, klass = t_iical_uv.klass,
        seq = t_iical_uv.seq )
    
    # save target UVFITS
    target_uvfits_file_name = '%s%s%s.UVFITS' % ( fits_path, target_name, extension )
    write_fits_uv( t_iical_uv, target_uvfits_file_name )
    target_uvfits_file_names.append( target_uvfits_file_name )
  
  # cleanup
  clear_aips_disk( aips_disk )
  
  return target_uvfits_file_names

###############################################################################

project_name = 'AF295'
aips_disk = 1

multi_uv = get_aips_file( aips_disk, project_name, 'P BAND', -1, 'UV' )
call_aips_task( 'FILLM', datain = './fits/%s_' % ( project_name ), band = 'P', outdisk = multi_uv.disk, outname = multi_uv.name, outseq = multi_uv.seq, ncount = 2, cparm = [ 0,1,16,0,0,0,0,10,10 ], bparm = [ 0 ] )
multi_uv.rename( klass = 'VLAPA' )

temp_uv = get_aips_file( aips_disk, project, 'TEMP', -1, 'UV' )
call_aips_task( 'NOIFS', indata = multi_uv, outdata = temp_uv, flagver = 0,
    bif = 1, eif = 2 )
multi_uv.zap()
temp_uv.rename( name = multi_uv.name, klass = multi_uv.klass, seq = multi_uv.seq )






###############################################################################

# import Python modules
from sys import *
from os import *
from datetime import *
from math import *

# import user modules
from files import *
from aips import *
from sphere import *
from parameter import *
from skymodel import *
from image import *
from calibrate import *
from solutions import *
from error import *

###############################################################################

def get_source_list_from_facets( facets, flux_min, peak_flux_ratio_max = None,
    area_ratio_max = None, size_ratio_min = 0.75, blank_factor = 2.0, edge_size = 8 ):

  facet_count = restore_parameter( facets, 'facet_count' )
  facet_list = range( 1, facet_count + 1 )
#  facet_size = get_image_size( facets )
  source_list = []

  # create copies of facets to blank
  peak_facets = get_aips_file( facets.disk, 'PEAK', facets.klass, -1, 'MA' )
  for i in facet_list:
    facet_i = get_facet( facets, i )
    peak_facet_i = get_facet( peak_facets, i )
    call_aips_task( 'MOVE', indata = facet_i, outdata = peak_facet_i, 
        userid = get_aips_userid(), opcode = '' )
    fill_image( peak_facet_i, do_edge_circle = True, edge_size = edge_size )

  # determine minimum flux for sources
  if ( not peak_flux_ratio_max is None ):
    peak_flux_max = 0.
    for i in facet_list:
      peak_facet_i = get_facet( peak_facets, i )
      max_facet = get_image_extremum( peak_facet_i, force_positive = True )
      if ( max_facet is None ):
        max_facet = [ 0., [ 0, 0 ] ]
      [ max_flux, max_pos ] = max_facet
      if ( max_flux > peak_flux_max ):
        peak_flux_max = max_flux
    peak_flux_min = max( [ flux_min, peak_flux_max / peak_flux_ratio_max ] )
  else:
    peak_flux_min = flux_min

  # find and fit sources
  res_facets = get_aips_file( facets.disk, 'RES', facets.klass, -1, 'MA' )
  for i in facet_list:
    facet_i = get_facet( facets, i )
    [ beam_bmaj, beam_bmin, beam_bpa ] = get_beam_size( facet_i )
    res_facet_i = get_facet( res_facets, i )
    call_aips_task( 'MOVE', indata = facet_i, outdata = res_facet_i, userid = get_aips_userid(),
        opcode = '' )
    peak_facet_i = get_facet( peak_facets, i )
    max_facet = get_image_extremum( peak_facet_i, force_positive = True )
    if ( max_facet is None ):
      max_facet = [ 0., [ 0, 0 ] ]
    [ max_flux, max_pos ] = max_facet
    peak_facet_i.zap()

    while ( max_flux >= peak_flux_min ):

      fit_results = fit_gaussian_to_peak( res_facet_i, pos = max_pos, return_double_fit = True ) 
      blank_list = []

      if ( len( fit_results ) == 0 ):

        [ peak_x, peak_y ] = [ float( p ) for p in max_pos ]
        peak_flux = 0.
        if ( beam_bmaj != 0. ):
          [ int_bmaj, int_bmin, int_bpa ] = [ beam_bmaj, beam_bmin, beam_bpa ]
        else:
          pixel_size = max( get_pixel_size( facet_i ) )
          [ int_bmaj, int_bmin, int_bpa ] = [ pixel_size * 4., pixel_size * 4., 0. ]
        blank_list.append( [ [ peak_x, peak_y ], peak_flux, [ int_bmaj, int_bmin, int_bpa ] ] )

      else:

        for fit in fit_results:
          [ [ peak_x, peak_y ], peak_flux, [ int_bmaj, int_bmin, int_bpa ] ] = fit
          
          # check fit against clean beam size
          if ( beam_bmaj != 0. ):
            if ( ( int_bmaj < size_ratio_min * beam_bmaj ) or 
                ( int_bmin < size_ratio_min * beam_bmin ) ):
              [ peak_x, peak_y ] = max_pos
              peak_flux = 0.
              [ int_bmaj, int_bmin, int_bpa ] = [ beam_bmaj, beam_bmin, beam_bpa ]
              blank_list.append( [ [ peak_x, peak_y ], peak_flux, [ int_bmaj, int_bmin, int_bpa ] ] )
              continue

          # check fit against selection criteria
          if ( beam_bmaj != 0. ):
            area_ratio = ( int_bmaj * int_bmin ) / ( beam_bmaj * beam_bmin )
            int_flux = peak_flux * area_ratio
            if ( not area_ratio_max is None ):
              if ( area_ratio > area_ratio_max ):
                blank_list.append( [ [ peak_x, peak_y ], peak_flux, [ int_bmaj, int_bmin, int_bpa ] ] )
                continue
          else: # special case when fitting to dirty beam facets
            int_flux = peak_flux

          # check if source is nearest to this facet center, and if so, add to source list
          blank_list.append( [ [ peak_x, peak_y ], peak_flux, [ int_bmaj, int_bmin, int_bpa ] ] )
          peak_radec = calculate_source_radec( res_facet_i, [ peak_x, peak_y ] )
          overlapping_facet_list = [ i ] + get_facet_overlap( facet_i )
          [ [ main_i, main_pos ] ] = find_source_facets( facets, peak_radec, primary_facet_only = True,
              facet_list = overlapping_facet_list )
          if ( i == main_i ):
            source_list.append( [ i, [ peak_x, peak_y ], peak_flux, int_flux, 
                [ int_bmaj, int_bmin, int_bpa ] ] )

      for blank in blank_list:
        [ [ peak_x, peak_y ], peak_flux, [ int_bmaj, int_bmin, int_bpa ] ] = blank

        # subtract fitted source from residual facet
        # put central source area to zero to prevent fitting on residuals
        [ bmaj_pix, bmin_pix, bpa_pix ] = convert_beam_size( facet_i, 
            beam = [ int_bmaj, int_bmin, int_bpa ], to_pixel = True )
        old_res_facet_i = get_aips_file( res_facet_i.disk, res_facet_i.name, 
            res_facet_i.klass, res_facet_i.seq, 'MA' )
        old_res_facet_i.rename( name = 'OLDRES', seq = 0 )
        call_aips_task( 'IMMOD', indata = old_res_facet_i, outdata = res_facet_i,
            opcode = 'GAUS', ngaus = 1, fmax = [ - peak_flux, 0,0,0 ],
            fpos = [ [ peak_x, peak_y ], [ 0,0 ], [ 0,0 ], [ 0,0 ] ],
            fwidth = [ [ bmaj_pix, bmin_pix, bpa_pix ], [ 0,0,0 ], [ 0,0,0 ], [ 0,0,0 ] ],
            factor = 1 )
        old_res_facet_i.zap()
        fill_source( res_facet_i, [ peak_x, peak_y ], value = 0.,
            beam = [ blank_factor * int_bmaj, blank_factor * int_bmin, int_bpa ] )

      # blank facet boundaries and get next extremum
      call_aips_task( 'MOVE', indata = res_facet_i, outdata = peak_facet_i,
          userid = get_aips_userid() )
      fill_image( peak_facet_i, do_edge_circle = True, edge_size = edge_size )
      max_facet = get_image_extremum( peak_facet_i, force_positive = True )
      if ( max_facet is None ):
        max_facet = [ 0., [ 0, 0 ] ]
      [ max_flux, max_pos ] = max_facet
      peak_facet_i.zap()

    res_facet_i.zap()

  # sort source by integrated flux
  if ( len( source_list ) > 0 ):
    source_list.sort( cmp = lambda a, b: cmp( b[ 3 ], a[ 3 ] ) )

  return source_list

###############################################################################

def get_source_list( image, flux_min, peak_flux_ratio_max = None, area_ratio_max = None,
    size_ratio_min = 0.75, blank_factor = 2.0, edge_size = 8 ):
  
  # create copy of image
  facets = get_aips_file( image.disk, image.name, 'ICL001', -1, 'MA' )
  call_aips_task( 'MOVE', indata = image, outdata = facets, userid = get_aips_userid() )
  store_parameter( facets, 'facet_count', 1 )
  determine_facet_overlap( facets )
  
  # extract sources
  source_list = get_source_list_from_facets( facets, flux_min, 
      peak_flux_ratio_max = peak_flux_ratio_max, area_ratio_max = area_ratio_max,
      size_ratio_min = size_ratio_min, blank_factor = blank_factor,
      edge_size = edge_size )
  facets.zap()
  
  return source_list

###############################################################################

def peel_s_facets_old( uv, facets, sigma_min = 5., signal_to_noise_min = 5., 
    solution_interval_max = 60., improvement_limit = 0.02, apply_solutions = True,
    reference_antenna = 0, imagr_params = {}, calib_params = {}, snr_limit = 2.5,
    sidelobe_rejection = 0. ):

  # calculate required signal per integration interval
#  time_count = restore_parameter( uv, 'time_count' )
  time_count = len( get_time_list( uv ) )
  cpb_noise = restore_parameter( uv, 'cpb_noise' )
  integration_time = restore_parameter( uv, 'integration_time' )
  facet_file_name = restore_parameter( facets, 'facet_file_name' )
  noise_per_interval_min = cpb_noise * sqrt( float( time_count ) )
  signal_per_interval_min = signal_to_noise_min * noise_per_interval_min
  signal_min = signal_per_interval_min / sqrt( floor( solution_interval_max / integration_time ) )

  # get list of bright enough, compact enough sources to peel
  source_list = get_source_list_from_facets( facets, signal_min )

  # make a copy of input UV
  peel_uv = get_aips_file( uv.disk, uv.name, 'PEELS', -1, 'UV' )
  call_aips_task( 'MOVE', indata = uv, outdata = peel_uv, userid = get_aips_userid(), opcode = '' )

  # if model is bright enough, try peeling it
  if ( len( source_list ) != 0 ):

    # split UV to get timeranges where source is above horizon
    rts_times = calculate_rise_transit_set_times( uv, radec = get_radec( facets ) )
    [ up_uv, down_uv ] = split_uv_on_time_range( uv, rts_times )

    if ( not up_uv is None ):
      peel_uv.zap()

      # add model back to UV
      add_up_uv = add_model( up_uv, facets, sigma = 0., apply_solutions = apply_solutions,
          keep_solutions = True, flag_solutions = False )
      up_uv.zap()

      # use an initial point source model for peeling
      total_model_flux = 0.
      facet_count = restore_parameter( facets, 'facet_count' )
      for i in range( 1, facet_count + 1 ):
        facet_i = get_facet( facets, i )
        total_model_flux = total_model_flux + get_model_flux( facet_i )
      for i in range( 1, facet_count + 1 ):
        facet_i = get_facet( facets, i )
        facet_pixel_ref = get_pixel_reference( facet_i )
        if ( i == 1 ):
          call_aips_task( 'CCMOD', indata = facets, invers = -1, opcode = 'POIN', 
              flux = total_model_flux, pixxy = facet_pixel_ref )
        else:
          call_aips_task( 'CCMOD', indata = facets, invers = -1, opcode = 'POIN',
              flux = 0., pixxy = facet_pixel_ref )

      # selfcal and image facets
      selfcal_image_clean_facets_old( add_up_uv, facets, sigma_min = sigma_min,
          reference_antenna = reference_antenna, do_sdi_clean = False,
          signal_to_noise_min = signal_to_noise_min, snr_limit = snr_limit,
          improvement_limit = improvement_limit, imagr_params = imagr_params,
          calib_params = calib_params, sidelobe_rejection = sidelobe_rejection )

      # subtract updated model from UV
      peel_up_uv = subtract_model( add_up_uv, facets, sigma = sigma_min,
          apply_solutions = True, keep_solutions = True, flag_solutions = True )
      add_up_uv.zap()

      # combine UV
      dummy_uv = merge_uv( peel_up_uv, down_uv )
      peel_up_uv.zap()
      dummy_uv.rename( name = peel_uv.name, klass = peel_uv.klass, seq = peel_uv.seq )
      peel_uv = get_aips_file( peel_uv.disk, peel_uv.name, peel_uv.klass, peel_uv.seq, 'UV' )

    if ( not down_uv is None ):
      down_uv.zap()

  return peel_uv

###############################################################################

def peel_a_facets_old( uv, facets, sigma_min = 5., signal_to_noise_min = 10., 
    solution_interval_max = 60., improvement_limit = 0.05, apply_solutions = True,
    reference_antenna = 0, imagr_params = {}, calib_params = {}, snr_limit = 2.,
    signal_multiplier = 10., print_info = True, try_final_amplitude = False,
    amplitude_interval = 5., convergence_limit = -0.5, conversion_method = 'DFT',
    phase_interval_min = 0., resolve_power = 0.5, sidelobe_rejection = 0. ):
  
  facet_file_name = restore_parameter( facets, 'facet_file_name' )
  
  # calculate required signal per integration interval
#  time_count = restore_parameter( uv, 'time_count' )
  time_count = len( get_time_list( uv ) )
  cpb_noise = restore_parameter( uv, 'cpb_noise' )
  integration_time = restore_parameter( uv, 'integration_time' )
  facet_file_name = restore_parameter( facets, 'facet_file_name' )
  noise_per_interval_min = cpb_noise * sqrt( float( time_count ) )
  signal_per_interval_min = signal_to_noise_min * noise_per_interval_min
  signal_min = signal_per_interval_min / sqrt( floor( solution_interval_max / 
      integration_time ) )
  
  # allow for a broader initial selection of sources, as the initial flux measurement
  # might be affected by ionospheric phase errors
  # fainter sources which need too long a calibration solution interval are rejected later on
  signal_min = signal_min / signal_multiplier
  if print_info:
    print '... minimum peeling peak flux = %s Jy' % ( repr( signal_min ) )
  
  # get list of bright enough, compact enough sources to peel
  source_list = get_source_list_from_facets( facets, signal_min, area_ratio_max = None )

  
  # make a copy of input UV
  peel_uv = get_aips_file( uv.disk, uv.name, 'PEELA', - 1, 'UV' )
  call_aips_task( 'MOVE', indata = uv, outdata = peel_uv, 
      userid = get_aips_userid(), opcode = '' )
  
  # if some sources are bright enough, try to peel them
  peel_facets = get_aips_file( facets.disk, 'PEELA', facets.klass, - 1, 'MA' )
  peel_facet_file_name = facet_file_name + '.PEELA'
  j = 0
  facet_list = []
  facets_done = []
  for source in source_list:
    # split UV to get timeranges where source is above horizon
    [ i, pos, peak_flux, int_flux, shape ] = source
    if ( i in facets_done ):
      continue
    facet_i = get_facet( facets, i )
    rts_times = calculate_rise_transit_set_times( peel_uv, radec = get_radec( facet_i ) )
    [ up_uv, down_uv ] = split_uv_on_time_range( peel_uv, rts_times )
    
    if ( ( not up_uv is None ) and ( not i in facet_list ) ):
      
      # copy facet and beam
      j = j + 1
      facet_list.append( i )
      beam_i = get_facet_beam( facet_i )
      facet_j = get_facet( peel_facets, j )
      beam_j = get_facet_beam( facet_j )
      call_aips_task( 'MOVE', indata = facet_i, outdata = facet_j, 
          userid = get_aips_userid() )
      call_aips_task( 'MOVE', indata = beam_i, outdata = beam_j, 
          userid = get_aips_userid() )
      if ( j == 1 ):
        call_aips_task( 'TACOP', indata = get_facet( facets, 1 ), inext = 'PS', 
            invers = 0, ncount = 1, outdata = facet_j, outvers = 0 )
      extract_facet_definitions( facet_file_name, facet_list, peel_facet_file_name )
#      store_parameter( peel_facets, 'facet_count', j )
      
      # add model back to UV
      add_up_uv = add_model( up_uv, peel_facets, facet_list = [ j ], sigma = 0.,
          apply_solutions = apply_solutions, keep_solutions = True,
          flag_solutions = False )
      up_uv.zap()
      
      # use an initial central point source model for peeling
      facet_pixel_ref = get_pixel_reference( facet_j )
      call_aips_task( 'CCMOD', indata = facet_j, invers = - 1, opcode = 'POIN', 
          flux = int_flux, pixxy = facet_pixel_ref )
      
      # selfcal and image facets
      converge = selfcal_image_clean_facet_old( add_up_uv, facet_j, 
          facet_file_name = peel_facet_file_name, snr_limit = snr_limit,
          sigma_min = sigma_min, do_sdi_clean = False, re_center_model = False, 
          signal_to_noise_min = signal_to_noise_min, print_info = print_info, 
          try_final_amplitude = try_final_amplitude, restore_components = True,
          amplitude_interval = amplitude_interval, flux_rejection_ratio = 0.2,
          improvement_limit = improvement_limit, resolve_power = resolve_power,
          convergence_limit = convergence_limit, frequency_correction = False,
          conversion_method = conversion_method, imagr_params = imagr_params,
          calib_params = calib_params, reference_antenna = reference_antenna,
          phase_interval_min = phase_interval_min, 
          sidelobe_rejection = sidelobe_rejection )
      
      if converge:
        facets_done.append( i )
        solution_interval = 60. * restore_parameter( facet_j, 'solution_interval' )
        if ( solution_interval > solution_interval_max ):
          converge = False
      
      if converge:
        peel_uv.zap()
        
        # subtract updated model from UV
        peel_up_uv = subtract_model( add_up_uv, peel_facets, facet_list = [ j ],
          sigma = sigma_min, apply_solutions = True, keep_solutions = True,
          flag_solutions = True )
        add_up_uv.zap()
        
        # combine UV
        if ( not down_uv is None ):
          dummy_uv = merge_uv( peel_up_uv, down_uv )
          peel_up_uv.zap()
          down_uv.zap()
          dummy_uv.rename( name = peel_uv.name, klass = peel_uv.klass, seq = peel_uv.seq )
        else:
          peel_up_uv.rename( name = peel_uv.name, klass = peel_uv.klass, seq = peel_uv.seq )
      
      else:
        add_up_uv.zap()
        if ( not down_uv is None ):
          down_uv.zap()
        facet_j.zap()
        beam_j.zap()
        remove_facet( peel_facet_file_name, j )
        j = j - 1
        facet_list = facet_list[ : - 1 ]
    
    else:
      if ( not up_uv is None ):
        up_uv.zap()
      if ( not down_uv is None ):
        down_uv.zap()
  
  if ( j == 0 ):
    peel_facets = None
  else:
    store_parameter( peel_facets, 'facet_count', j )
  
  return [ peel_uv, peel_facets ]

###############################################################################

def peel_o_facets_old( uv, facets, sigma_min = 5., signal_to_noise_min = 8.,
    solution_interval_max = 60., try_final_amplitude = False, snr_limit = 2.5,
    improvement_limit = 0.02, apply_solutions = True, reference_antenna = 0,
    imagr_params = {}, calib_params = {}, sidelobe_rejection = 0. ):

  # calculate required signal per integration interval
#  time_count = restore_parameter( uv, 'time_count' )
  time_count = len( get_time_list( uv ) )
  cpb_noise = restore_parameter( uv, 'cpb_noise' )
  integration_time = restore_parameter( uv, 'integration_time' )
  facet_file_name = restore_parameter( facets, 'facet_file_name' )
  noise_per_interval_min = cpb_noise * sqrt( float( time_count ) )
  signal_per_interval_min = signal_to_noise_min * noise_per_interval_min
  signal_min = signal_per_interval_min / sqrt( floor( solution_interval_max / 
      integration_time ) )

  # get list of bright enough, compact enough sources to peel
  source_list = get_source_list_from_facets( facets, signal_min, area_ratio_max = 3.**2 )

  # make a copy of input UV
  peel_uv = get_aips_file( uv.disk, uv.name, 'PEELO', -1, 'UV' )
  call_aips_task( 'MOVE', indata = uv, outdata = peel_uv,
      userid = get_aips_userid(), opcode = '' )

  # if some sources are bright enough, try to peel them
  for source in source_list:
    [ i, pos, peak_flux, int_flux, shape ] = source
    facet_i = get_facet( facets, i )

    # add model back to UV
    add_uv = add_model( peel_uv, facets, facet_list = [ i ], sigma = 0., 
        apply_solutions = apply_solutions, keep_solutions = True,
        flag_solutions = False )
    peel_uv.zap()

    # use an initial central point source model for peeling
    facet_pixel_ref = get_pixel_reference( facet_i )
    call_aips_task( 'CCMOD', indata = facet_i, invers = -1, opcode = 'POIN', 
        flux = int_flux, pixxy = facet_pixel_ref )

    # selfcal and image facet
    selfcal_image_clean_facets_old( add_uv, facets, facet_list = [ i ], 
        sigma_min = sigma_min, signal_to_noise_min = signal_to_noise_min,
        try_final_amplitude = try_final_amplitude, do_sdi_clean = False,
        improvement_limit = improvement_limit, imagr_params = imagr_params,
        reference_antenna = reference_antenna, calib_params = calib_params,
        snr_limit = snr_limit, sidelobe_rejection = sidelobe_rejection )

    # subtract updated model from UV
    dummy_uv = subtract_model( add_uv, facets, facet_list = [ i ], sigma = sigma_min,
        apply_solutions = True, keep_solutions = True, flag_solutions = True )
    add_uv.zap()
    dummy_uv.rename( name = peel_uv.name, klass = peel_uv.klass, seq = peel_uv.seq )
    peel_uv = get_aips_file( peel_uv.disk, peel_uv.name, peel_uv.klass, peel_uv.seq, 'UV' )

  return peel_uv

###############################################################################

def replace_model_solutions_with_peel_solutions( uv, facets, peel_facets,
    facet_list = [], version = 0, re_grid_solutions = True, max_separation = None,
    normalize_amplitudes = False ):
  # max_separation in arcsec
  
  if ( len( facet_list ) == 0 ):
    peel_facet_count = restore_parameter( peel_facets, 'facet_count' )
    peel_facet_list = range( 1, 1 + peel_facet_count )
  else:
    peel_facet_count = len( facet_list )
    peel_facet_list = [ i for i in facet_list ]
  if re_grid_solutions:
    time_list = get_time_list( uv )
    gap_time = restore_parameter( uv, 'integration_time' ) / ( 4. * 60. )
  if ( not max_separation is None ):
    max_radius = max_separation / 3600.
  else:
    max_radius = 0.45 * min( array( get_pixel_size( facets ) ) * 
        array( get_image_size( facets ) ) ) / 3600.
  used_facet_list = []
  for i in peel_facet_list:
    peel_facet_i = get_facet( peel_facets, i )
    peel_radec = get_radec( peel_facet_i )
    source_facet_list = find_source_facets( facets, peel_radec,
        primary_facet_only = False )
#    [ bmaj, bmin, bpa ] = get_beam_size( peel_facet_i )
#    if ( len( source_facet_list ) > 0 ):
#      [ j, pos ] = source_facet_list[ 0 ]
    for [ j, pos ] in source_facet_list:
      if ( j in used_facet_list ):
        continue
      facet_j = get_facet( facets, j )
      radec = get_radec( facet_j )
      [ radius, angle ] = calculate_angular_separation( peel_radec, radec )
#      if ( 3600. * radius < max_beam_offset * bmaj ):
      if ( radius > max_radius ):
        continue
      if normalize_amplitudes:
        snver = uv.table_highver( 'SN' ) + 1
        call_aips_task( 'TACOP', indata = peel_facet_i, inext = 'SN',
            invers = version, ncount = 1, outdata = uv, outvers = snver )
        call_aips_task( 'SNCOR', indata = uv, snver = snver, opcode = 'NORM' )
        call_aips_task( 'TACOP', indata = uv, inext = 'SN',
            invers = 0, ncount = 1, outdata = facet_j, outvers = 0 )
        for ver in range( snver, 1 + uv.table_highver( 'SN' ) ):
          uv.zap_table( 'SN', ver )
      else:
        call_aips_task( 'TACOP', indata = peel_facet_i, inext = 'SN',
            invers = version, ncount = 1, outdata = facet_j, outvers = 0 )
      if re_grid_solutions:
        try:
          re_sample_solutions( facet_j, gap_time = gap_time, time_list = time_list,
              interpolation_method = 'nearest', force_reference = True )
        except RuntimeError:
          reference_antenna = get_reference_antenna( facet_j )
          re_reference_solutions( facet_j, reference_antenna )
          re_sample_solutions( facet_j, gap_time = gap_time, time_list = time_list,
              interpolation_method = 'nearest', force_reference = True )
        except:
          raise error( 'unable to resample solutions' )
        else:
          pass
      used_facet_list.append( j )
  return

###############################################################################
###############################################################################

def calibrate_facet( uv, facets, facet_id, sigma_min = 0., snr_limit = 2.5,
    signal_to_noise_min = 15., reference_antenna = 0, conversion_method = 'DFT',
    do_amplitude = False, amplitude_interval = 5., phase_interval_min = 0.,
    calib_params = {}, sidelobe_rejection = 0., frequency_correction = False,
    normalize_gains = True ):
  
  # determine total flux in model
  facet_i = get_facet( facets, facet_id )
  total_model_flux = get_model_flux( facet_i )
  
  cpb_noise = restore_parameter( uv, 'cpb_noise' )
  if ( do_amplitude ):
    solution_interval = amplitude_interval
    interval_count = 1.
  else:
    # calculate noise per integration interval
    time_count = len( get_time_list( uv ) )
    noise_per_interval = cpb_noise * sqrt( float( time_count ) )
    
    # calculate solution interval
    integration_time = restore_parameter( uv, 'integration_time' ) / 60.
    sn_per_interval = total_model_flux / noise_per_interval
    interval_count = ceil( ( signal_to_noise_min / sn_per_interval )**2 )
    solution_interval = integration_time * interval_count
    solution_interval = max( [ solution_interval, phase_interval_min ] )
  
  # calibrate UV data 
  calibrate_model( uv, facets, reference_antenna = reference_antenna,
      phase_interval = solution_interval, do_amplitude = do_amplitude,
      amplitude_interval = amplitude_interval, apply_solutions = False,
      keep_flags = False, sigma = sigma_min, calib_params = calib_params,
      sidelobe_rejection = sidelobe_rejection, facet_list = [ facet_id ],
      frequency_correction = frequency_correction, snr_limit = snr_limit,
      conversion_method = conversion_method, normalize_gains = normalize_gains )
  
  return solution_interval

###############################################################################

def image_clean_facet( uv, facets, facet_id, clean_flux_min, facet_file_name = '',  
    conversion_method = 'DFT', model_version = 0, apply_solutions = True,
    do_sdi_clean = False, restore_components = False, frequency_correction = False,
    gain = 0.1, factor = 0., imagr_params = {} ):
  
  # process input params
  i_params = imagr_params.copy()
  if ( not restore_components ):
    i_params[ 'bmaj' ] = -1
  if ( facet_file_name == '' ):
    used_facet_file_name = restore_parameter( facets, 'facet_file_name' )
  else:
    used_facet_file_name = facet_file_name
  used_facet_file_name = path.expandvars( used_facet_file_name )
  
  # get selected facet
  sel_facet_file_name = used_facet_file_name + '.SEL'
  extract_facet_definitions( used_facet_file_name, [ facet_id ], sel_facet_file_name )
  facet = get_facet( facets, facet_id )
  beam = get_facet_beam( facet )
  sel_facet = get_aips_file( facet.disk, 'SEL', 'ICL001', - 1, 'MA' )
  sel_beam = get_facet_beam( sel_facet )
  call_aips_task( 'MOVE', indata = facet, outdata = sel_facet, userid = get_aips_userid() )
  call_aips_task( 'MOVE', indata = beam, outdata = sel_beam, userid = get_aips_userid() )
  
  # set parameters
  if frequency_correction:
    dish_diameter = restore_parameter( uv, 'dish_diameter' )
  else:
    dish_diameter = 0.
  cell_size = get_pixel_size( sel_facet, make_absolute = True )
  facet_size = get_image_size( sel_facet )
  channel_count = get_channel_count( uv )
  uv_size = restore_parameter( uv, 'pb_image_size' )
  if apply_solutions:
    docalib = 100
    gainuse = 0
  else:
    docalib = -1
    gainuse = -1
  if do_sdi_clean:
    sdi_param = 0.1
  else:
    sdi_param = 0.

  # remove old model
  if ( ( model_version != 0 ) and ( table_exists( sel_facet, 'CC', model_version ) ) ):
    sel_facet.zap_table( 'CC', model_version )
  
  # image facet
  call_aips_task( 'IMAGR', indata = uv, nchav = channel_count, docalib = docalib,
      gainuse = gainuse, outver = model_version, niter = 100000, flagver = -1,
      outdisk = sel_facet.disk, outname = sel_facet.name, outseq = sel_facet.seq,
      in2disk = uv.disk, cellsize = cell_size, imsize = facet_size, do3dimag = 1,
      flux = 0.95 * clean_flux_min, boxfile = sel_facet_file_name, dotv = 0,
      cmethod = conversion_method, minpatch = facet_size[ 0 ] - 1, overlap = 2,
      gain = gain, nfield = 1, bcomp = [ 0 for i in range( 64 ) ], allokay = 1,
      imagrprm = [ dish_diameter, 0,0, sdi_param, 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0 ],
      maxpixel = 0, factor = factor, uvsize = [ uv_size, uv_size ], **i_params )
  
  # replace original facet & cleanup
  facet.zap()
  beam.zap()
  sel_facet.rename( name = facet.name, klass = facet.klass, seq = facet.seq )
  sel_beam.rename( name = beam.name, klass = beam.klass, seq = beam.seq )
  remove_file( sel_facet_file_name )
  
  return

###############################################################################

def selfcal_image_clean_facet( uv, facets, facet_id, facet_file_name = '',
    sigma_min = 3., signal_to_noise_min = 15., improvement_limit = 0.02, 
    selfcal_cycle_min = 4, reference_antenna = 0, conversion_method = 'DFT',
    add_boxes = True, box_sigma = 5., clean_box_radius = 5, do_sdi_clean = False,
    re_center_model = False, imagr_params = {}, calib_params = {}, snr_limit = 2.,
    convergence_limit = -0.5, frequency_correction = False, print_info = True, 
    amplitude_interval = 5., offset_ratio_max = 2., allow_pixel_peak = True,
    flux_rejection_ratio = 0.5, phase_interval_min = 0.1, resolve_power = 0.5,
    amplitude_noise_factor = 0., amplitude_snr = 300., sidelobe_rejection = 0.,
    edge_size = 8, normalize_gains = True, convolve_size = None ):
  
  # initialise some parameters
  channel_count = get_channel_count( uv )
  cpb_noise = restore_parameter( uv, 'cpb_noise' )
  if ( facet_file_name == '' ):
    sel_facet_file_name = restore_parameter( facets, 'facet_file_name' )
  else:
    sel_facet_file_name = facet_file_name
  sel_facet_file_name = path.expandvars( sel_facet_file_name )
  sel_facet = get_facet( facets, facet_id )
  facet_size = get_image_size( sel_facet )
  model_version = sel_facet.table_highver( 'CC' ) + 1
  facet_radec = get_radec( sel_facet )
  if re_center_model:
    cmethod = 'DFT'
  else:
    cmethod = conversion_method
  
  # determine signal/noise level in facet
  signal = get_model_flux( sel_facet )
  initial_signal = signal
  rms = get_image_rms( sel_facet )
  [ avg, noise ] = call_aips_task( 'IMEAN', indata = sel_facet, pixavg = 0.,
      pixstd = rms, pixrange = [ -5. * rms, 5. * rms  ], 
      outputs = [ 'pixavg', 'pixstd' ] )
  if ( ( noise <= 0. ) or ( noise > 2. * rms ) ):
    if print_info:
      print '... WARNING: histogram noise fit failed, using image RMS instead'
    noise = rms
  initial_noise = noise
  sn = signal / noise
  initial_sn = sn
  if print_info:
    print '... signal = %s, noise = %s, S/N = %s' % ( repr( signal ), repr( noise ), 
        repr( sn ) )
  
  # enter selfcal loop
  sn_improvement = 1.
  selfcal_i = 0
  while ( ( selfcal_i < selfcal_cycle_min ) or ( sn_improvement > improvement_limit ) ):
    selfcal_i = selfcal_i + 1
    last_sn = sn
    if print_info:
      print '... selfcal loop %s' % repr( selfcal_i )
    
    # determine total / peak flux ratio
    total_flux = get_model_flux( sel_facet )
    fit_results = fit_gaussian_to_peak( sel_facet, return_double_fit = True,
        offset_ratio_max = offset_ratio_max )
    if ( ( len( fit_results ) == 0 ) and allow_pixel_peak ):
      # use maximum pixel position
      if print_info:
        print 'WARNING: using maximum pixel position instead'
      fill_image( sel_facet, do_edge_circle = True, edge_size = edge_size )
      max_results = get_image_maximum( sel_facet )
      if ( max_results is None ):
        max_results = [ 0., [ 0, 0 ] ]
      fit_results = [ [ [ float( max_results[ 1 ][ 0 ] ), 
          float( max_results[ 1 ][ 1 ] ) ], max_results[ 0 ], [ 0.,0.,0. ] ] ]
    if ( len( fit_results ) == 0 ):
      if re_center_model:
        if print_info:
          print ( 'WARNING: gaussian fit to peeling source failed while ' + 
              're-centering model' )
        sn = 1.
        sn_improvement = - 100.
        break
      else:
        if print_info:
          print ( 'WARNING: gaussian fit to peeling source failed while ' + 
              'determining total / peak ratio' )
        total_peak_ratio = 1.
    else:
      [ fit_pos, fit_peak, fit_beam ] = fit_results[ 0 ]
      total_peak_ratio = total_flux / fit_peak
      if print_info:
        print '...... total / peak ratio = %s' % ( repr( total_peak_ratio ) )
    
    # center facet on flux peak
    if re_center_model:
      if ( not convolve_size is None ):
        [ bmaj, bmin, bpa ] = get_beam_size( sel_facet )
        if ( bmin < convolve_size ):
          bmin = convolve_size
          bmaj = max( bmaj * 1.01, convolve_size )
          temp_facet = get_aips_file( sel_facet.disk, 'TEMP', sel_facet.klass, -1, 'MA' )
          call_aips_task( 'CONVL', indata = sel_facet, outdata = temp_facet,
              opcode = 'GAUS', bmaj = bmaj, bmin = bmin, bpa = bpa, factor = 0. )
          new_fit_results = fit_gaussian_to_peak( temp_facet, return_double_fit = True,
              offset_ratio_max = offset_ratio_max )
          temp_facet.zap()
          if ( len( new_fit_results ) > 0 ):
            fit_results = new_fit_results
            [ fit_pos, fit_peak, fit_beam ] = fit_results[ 0 ]
      if ( selfcal_i == 1 ):
        fit_count = len( fit_results )
      elif ( len( fit_results ) != fit_count ):
        if print_info:
          print ( 'WARNING: number of gaussians in fit to peeling source ' + 
              'changed while re-centering source model' )
        sn = 1.
        sn_improvement = - 100.
        break
      fit_radec = calculate_source_radec( sel_facet, fit_pos )
      if ( len( fit_results ) > 1 ):
        [ fit_pos_2, fit_peak_2, fit_beam_2 ] = fit_results[ 1 ]
        fit_radec_2 = calculate_source_radec( sel_facet, fit_pos_2 )
        if ( fit_peak_2 > fit_peak ):
          fit_radec = fit_radec_2
      set_radec( sel_facet, fit_radec, shift_model = True )
      set_radec( sel_facet, facet_radec, shift_model = False )
    
    # self-calibrate UV on model
    # minimum SNR is increased to reflect extendedness of source
    resolve_factor = max( 1., total_peak_ratio )**resolve_power
    solution_interval = calibrate_facet( uv, facets, facet_id, sigma_min = 0.,
        signal_to_noise_min = signal_to_noise_min * resolve_factor,
        calib_params = calib_params, reference_antenna = reference_antenna, 
        conversion_method = cmethod, phase_interval_min = phase_interval_min,
        frequency_correction = frequency_correction, sidelobe_rejection = 0.,
        snr_limit = snr_limit )
    solution_version = uv.table_highver( 'SN' )
    
    # image and clean
    if ( add_boxes  ): # and ( selfcal_i > 1 ) ):
      add_clean_boxes( uv, facets, facet_list = [ facet_id ], keep_boxes = False,
          box_sigma = box_sigma, clean_box_radius = clean_box_radius,
          facet_based_boxes = True, sidelobe_rejection = sidelobe_rejection,
          print_info = print_info, facet_file_name = sel_facet_file_name )
    image_clean_facet( uv, facets, facet_id, model_version = model_version, 
        clean_flux_min = sigma_min * cpb_noise, restore_components = True,
        imagr_params = imagr_params, frequency_correction = frequency_correction,
        do_sdi_clean = do_sdi_clean, conversion_method = conversion_method,
        facet_file_name = sel_facet_file_name )
    
    # determine signal/noise level and S/N improvement
    rms = get_image_rms( sel_facet )
    [ avg, noise ] = call_aips_task( 'IMEAN', indata = sel_facet, pixavg = 0.,
        pixstd = rms, pixrange = [ -5. * rms, 5. * rms  ], 
        outputs = [ 'pixavg', 'pixstd' ] )
    if ( ( noise <= 0. ) or ( noise > 2. * rms ) ):
      if print_info:
        print '... WARNING: histogram noise fit failed, using image RMS instead'
      noise = rms
    
    signal = get_model_flux( sel_facet )
    flux_ratio = signal / initial_signal
    if ( flux_ratio < flux_rejection_ratio ):
      if print_info:
        print 'WARNING: source flux dropped below threshold, aborting selfcal'
      sn = 1.
      sn_improvement = -100.
      break
    sn = signal / noise
    sn_improvement = ( sn / last_sn ) - 1.
    if print_info:
      print '...... signal = %s, noise = %s, S/N = %s' % ( repr( signal ),
          repr( noise ), repr( sn ) )
      print '...... S/N improvement = %s' % ( repr( sn_improvement ) )

    # only keep solutions when S/N improvement was detected
    if ( selfcal_i == 1 ):
      last_solution_version = solution_version
      last_solution_interval = solution_interval
    elif ( ( selfcal_i <= selfcal_cycle_min ) or ( sn_improvement >= 0. ) ):
      uv.zap_table( 'SN', last_solution_version - 1 )
      uv.zap_table( 'SN', last_solution_version )
      last_solution_version = solution_version
      last_solution_interval = solution_interval
    else:
      uv.zap_table( 'SN', solution_version - 1 )
      uv.zap_table( 'SN', solution_version )
      solution_version = last_solution_version
      solution_interval = last_solution_interval
      sn = last_sn
      image_clean_facet( uv, facets, facet_id, model_version = model_version,
          clean_flux_min = sigma_min * cpb_noise, restore_components = True,
          do_sdi_clean = do_sdi_clean, frequency_correction = frequency_correction,
          imagr_params = imagr_params, conversion_method = conversion_method,
          facet_file_name = sel_facet_file_name )
  
  # allow for single amplitude & phase selfcal
  if ( ( noise > amplitude_noise_factor * cpb_noise ) and 
      ( sn > amplitude_snr ) and ( sn_improvement > - 99. ) ):
    while True: # so we can use break command
      
      last_noise = noise
      convergence = True
      selfcal_i = selfcal_i + 1
      last_sn = sn
      if print_info:
        print '... (a&p) selfcal loop %s' % repr( selfcal_i )
      
      # center facet on flux peak
      if re_center_model:
        fit_results = fit_gaussian_to_peak( sel_facet, return_double_fit = True,
            offset_ratio_max = offset_ratio_max )
        if ( ( len( fit_results ) == 0 ) and allow_pixel_peak ):
          # use maximum pixel position
          if print_info:
            print 'WARNING: using maximum pixel position instead'
          fill_image( sel_facet, do_edge_circle = True, edge_size = edge_size )
          max_results = get_image_maximum( sel_facet )
          if ( max_results is None ):
            max_results = [ 0., [ 0, 0 ] ]
          fit_results = [ [ [ float( max_results[ 1 ][ 0 ] ),
              float( max_results[ 1 ][ 1 ] ) ], max_results[ 0 ], [ 0.,0.,0. ] ] ]
        if ( not convolve_size is None ):
          [ bmaj, bmin, bpa ] = get_beam_size( sel_facet )
          if ( bmin < convolve_size ):
            bmin = convolve_size
            bmaj = max( bmaj * 1.01, convolve_size )
            temp_facet = get_aips_file( sel_facet.disk, 'TEMP', sel_facet.klass, -1, 'MA' )
            call_aips_task( 'CONVL', indata = sel_facet, outdata = temp_facet,
                opcode = 'GAUS', bmaj = bmaj, bmin = bmin, bpa = bpa, factor = 0. )
            new_fit_results = fit_gaussian_to_peak( temp_facet, return_double_fit = True,
                offset_ratio_max = offset_ratio_max )
            temp_facet.zap()
            if ( len( new_fit_results ) > 0 ):
              fit_results = new_fit_results
        if ( len( fit_results ) != fit_count ):
          if print_info:
            print ( 'WARNING: number of gaussians in peeling source fit ' +
                'changed while re-centering source model' )
          sn = 1.
          sn_improvement = - 100.
          break
        [ fit_pos, fit_peak, fit_beam ] = fit_results[ 0 ]
        fit_radec = calculate_source_radec( sel_facet, fit_pos )
        if ( len( fit_results ) > 1 ):
          [ fit_pos_2, fit_peak_2, fit_beam_2 ] = fit_results[ 1 ]
          fit_radec_2 = calculate_source_radec( sel_facet, fit_pos_2 )
          if ( fit_peak_2 > fit_peak ):
            fit_radec = fit_radec_2
        set_radec( sel_facet, fit_radec, shift_model = True )
        set_radec( sel_facet, facet_radec, shift_model = False )
      
      # apply phase solutions to UV data
      temp_uv = apply_solution_table( uv )
      
      # self-calibrate UV on model
      ref_ant = get_reference_antenna( uv )
      temp_solution_interval = calibrate_facet( temp_uv, facets, facet_id,
          signal_to_noise_min = signal_to_noise_min, calib_params = calib_params,
          reference_antenna = ref_ant, conversion_method = cmethod,
          amplitude_interval = amplitude_interval, do_amplitude = True, 
          sidelobe_rejection = 0., frequency_correction = frequency_correction,  
          normalize_gains = normalize_gains, snr_limit = snr_limit, sigma_min = 0. )
      temp_solution_version = temp_uv.table_highver( 'SN' )
      
      # image and clean
      if ( add_boxes ):
        add_clean_boxes( uv, facets, facet_list = [ facet_id ], keep_boxes = False,
            box_sigma = box_sigma, clean_box_radius = clean_box_radius,
            facet_based_boxes = True, sidelobe_rejection = sidelobe_rejection,
            print_info = print_info, facet_file_name = sel_facet_file_name )
      image_clean_facet( temp_uv, facets, facet_id, model_version = model_version,
          clean_flux_min = sigma_min * cpb_noise, restore_components = True,
          do_sdi_clean = do_sdi_clean, frequency_correction = frequency_correction,
          imagr_params = imagr_params, conversion_method = conversion_method,
          facet_file_name = sel_facet_file_name )
      
      # determine noise level and noise improvement
      rms = get_image_rms( sel_facet )
      [ avg, noise ] = call_aips_task( 'IMEAN', indata = sel_facet, pixavg = 0.,
          pixstd = rms, pixrange = [ -5. * rms, 5. * rms  ], 
          outputs = [ 'pixavg', 'pixstd' ] )
      if ( ( noise <= 0. ) or ( noise > 2. * rms ) ):
        if print_info:
          print '... WARNING: histogram noise fit failed, using image RMS instead'
        noise = rms
      signal = get_model_flux( sel_facet )
      flux_ratio = signal / initial_signal
      if ( flux_ratio < flux_rejection_ratio ):
        if print_info:
          print 'WARNING: source flux dropped below threshold, aborting selfcal'
        sn = 1.
        sn_improvement = -100.
        break
      
      sn = signal / noise
      sn_improvement = ( sn / last_sn ) - 1.
      if print_info:
        print '...... signal = %s, noise = %s, S/N = %s' % ( repr( signal ),
            repr( noise ), repr( sn ) )
        print '...... S/N improvement = %s' % ( repr( sn_improvement ) )
      break
    
    # only keep solutions when noise improvement was detected
    if ( ( sn_improvement >= 0. ) or ( noise <= last_noise ) ):
      # combine tables
      call_aips_task( 'TACOP', indata = uv, inext = 'SN', ncount = 1,
          invers = solution_version, outdata = temp_uv, outvers = 0 )
      combine_solutions( temp_uv, force_match = True )
      call_aips_task( 'TACOP', indata = temp_uv, inext = 'SN', ncount = 1,
          invers = 0, outdata = uv, outvers = 0 )
      solution_version = uv.table_highver( 'SN' )
      uv.zap_table( 'SN', last_solution_version - 1 )
    else:
      sn = last_sn
      image_clean_facet( uv, facets, facet_id, model_version = model_version,
          clean_flux_min = sigma_min * cpb_noise, restore_components = True,
          do_sdi_clean = do_sdi_clean, frequency_correction = frequency_correction,
          imagr_params = imagr_params, conversion_method = conversion_method,
          facet_file_name = sel_facet_file_name )
    temp_uv.zap()
  
  # do final check
  overall_sn_improvement = ( sn / initial_sn ) - 1.
  if ( overall_sn_improvement < convergence_limit ):
    try:
      uv.zap_table( 'SN', solution_version - 1 )
    except:
      pass
    try:
      uv.zap_table( 'SN', solution_version )
    except:
      pass
    return False
  
  # save solutions to facet
  # store parameters
  # remove solutions from UV data
  call_aips_task( 'TACOP', indata = uv, inext = 'SN', ncount = 1,
      invers = solution_version - 1, outdata = sel_facet, outvers = 0 )
  call_aips_task( 'TACOP', indata = uv, inext = 'SN', invers = solution_version,
      ncount = 1, outdata = sel_facet, outvers = 0 )
  uv.zap_table( 'SN', solution_version - 1 )
  uv.zap_table( 'SN', solution_version )
  store_parameter( sel_facet, 'solution_interval', solution_interval )
  
  return True

###############################################################################

def peel_pbo_facets( uv, facets, sigma_min = 3., signal_to_noise_min = 15.,
    apply_solutions = True, print_info = True, use_catalog = True, gain = 0.1,
    relative_solutions = True, max_relativity = 1, keep_solutions = True,
    clean_box_ratio = 5., double_search_ratio = 10., improvement_limit = 0.05,
    signal_multiplier = sqrt( 2. ), input_catalog = [], catalog_resolution = 50.,
    assoc_ratio = 1.e9, frequency_correction = False, conversion_method = 'DFT',
    max_peel_count = 20, convergence_limit = -0.5, phase_interval_max = 2.,
    phase_interval_min = 0., reject_multiples = False, offset_ratio_max = 1.,
    reference_antenna = 0, resolve_power = 0.5, snr_limit = 0.1, add_boxes = True,
    amplitude_interval = 5., amplitude_noise_factor = 0., amplitude_snr = 500.,
    model_subtracted = True, sort = '', force_simple_match = False, edge_size = 8,
    imagr_params = {}, calib_params = {}, sidelobe_rejection = 1.05, 
    clean_box_radius = 5, normalize_gains = True, flux_rejection_ratio = 0.5,
    catalog_search_ratio = 2., use_nvss = True, normalize_relative_gains = False ):
  
  # calculate required signal per integration interval
#  time_count = restore_parameter( uv, 'time_count' )
  time_count = len( get_time_list( uv ) )
  cpb_noise = restore_parameter( uv, 'cpb_noise' )
  integration_time = restore_parameter( uv, 'integration_time' )
  facet_file_name = restore_parameter( facets, 'facet_file_name' )
  facet_file_name_e = path.expandvars( facet_file_name )
  noise_per_interval_min = cpb_noise * sqrt( float( time_count ) )
  signal_per_interval_min = signal_to_noise_min * noise_per_interval_min
  signal_min = signal_per_interval_min / sqrt( 60. * phase_interval_max / integration_time )
  
  # allow for a broader initial selection of sources, as the initial flux measurement
  # might be affected by ionospheric phase errors
  # fainter sources which need too long a calibration solution interval are rejected later on
  signal_min = signal_min / signal_multiplier
  if print_info:
    print '... minimum peeling peak flux = %s Jy' % ( repr( signal_min ) )
  
  # get list of bright enough, compact enough sources to peel
  facet_count = restore_parameter( facets, 'facet_count' )
  source_list = get_source_list_from_facets( facets, signal_min, area_ratio_max = 3.**2 )
  if print_info:
    print '... found %s possible peeling sources' % ( repr( len( source_list ) ) )
  
  # detect double component sources
  close_double_list = []
  for s in range( len( source_list ) ): # source_i in source_list:
    # get source info
    source_i = source_list[ s ]
    [ i, pos_i, peak_flux_i, int_flux_i, shape_i ] = source_i
    facet_i = get_facet( facets, i )
    source_radec_i = calculate_source_radec( facet_i, pos_i )
    [ bmaj_i, bmin_i, bpa_i ] = get_beam_size( facet_i )
    # check if source is close double with other source(s)
    for source_j in source_list[ s + 1 : ]:
      [ j, pos_j, peak_flux_j, int_flux_j, shape_j ] = source_j
      facet_j = get_facet( facets, j )
      source_radec_j = calculate_source_radec( facet_j, pos_j )
      [ radius, angle ] = calculate_angular_separation( source_radec_i, source_radec_j )
      if ( 3600. * radius < double_search_ratio * bmaj_i ):
        close_double_list.append( source_i )
        close_double_list.append( source_j )
        if print_info:
          print '... sources at RADECs %s and %s are part of multiple component source' % (
              radec_to_string( source_radec_i ), radec_to_string( source_radec_j ) )

  # detect multiple (>2) component sources
  multi_source = [ False for source in close_double_list ]
  for s in range( 0, len( close_double_list ), 2 ):
    source_i = close_double_list[ s ]
    source_j = close_double_list[ s + 1 ]
    for k in range( len( close_double_list ) ):
      source_k = close_double_list[ k ]
      if ( not ( k in [ s, s + 1 ] ) ):
        if ( ( source_k == source_i ) or ( source_k == source_j ) ):
          multi_source[ s ] = True
          multi_source[ s + 1 ] = True
  
  # if requested, keep first double association of multiples
  if ( not reject_multiples ):
    for s in range( 0, len( close_double_list ), 2 ):
      if ( multi_source[ s ] or multi_source[ s + 1 ] ):
        source_i = close_double_list[ s ]
        source_j = close_double_list[ s + 1 ]
        if ( ( not source_i in close_double_list[ 0 : s ] ) and 
             ( not source_j in close_double_list[ 0 : s ] ) ):
          multi_source[ s ] = False
          multi_source[ s + 1 ] = False
  
  # merge doubles only, add to temporary source list
  temp_double_list = []
  for s in range( 0, len( close_double_list ), 2 ):
    source_i = close_double_list[ s ]
    [ i, pos_i, peak_flux_i, int_flux_i, shape_i ] = source_i
    facet_i = get_facet( facets, i )
    source_radec_i = calculate_source_radec( facet_i, pos_i )
    if multi_source[ s ]:
      if print_info:
        print '... discarding multiple (>2) component source at RADEC %s' % ( 
            radec_to_string( source_radec_i ) )
    else:
      source_j = close_double_list[ s + 1 ]
      [ j, pos_j, peak_flux_j, int_flux_j, shape_j ] = source_j
      facet_j = get_facet( facets, j )
      source_radec_j = calculate_source_radec( facet_j, pos_j )
      if print_info:
        print '... combining double components at RADECs %s and %s' % ( 
            radec_to_string( source_radec_i ), radec_to_string( source_radec_j ) )
        if ( ( ( peak_flux_i > peak_flux_j ) and ( int_flux_i < int_flux_j ) ) or 
             ( ( peak_flux_i < peak_flux_j ) and ( int_flux_i > int_flux_j ) ) ): 
          print ( '...... WARNING: order of peak fluxes and total fluxes are ' + 
              'different, using peak fluxes' )
      if ( j == i ):
        k = i
        facet_k = facet_i
      else:
        if ( ( pos_i[ 0 ]**2 + pos_i[ 1 ]**2 ) < ( pos_j[ 0 ]**2 + pos_j[ 1 ]**2 ) ):
          k = i
          facet_k = facet_i
          pos_j = calculate_source_position( facet_k, source_radec_j )
        else:
          k = j
          facet_k = facet_j
          pos_i = calculate_source_position( facet_k, source_radec_i )
      # put clean box position in between sources
      [ radius, angle ] = calculate_angular_separation( source_radec_i, source_radec_j )
      r_ik = 0.5 * radius
      clean_radec_k = calculate_offset_position( source_radec_i, r_ik, angle )
      clean_pos_k = calculate_source_position( facet_k, clean_radec_k )
      bmajmin_k = ( 3600. * radius + clean_box_ratio * 
          ( shape_i[ 0 ] + shape_j[ 0 ] ) ) / clean_box_ratio
      shape_k = [ bmajmin_k, bmajmin_k, angle ]
      # put new facet at position of component with highest peak
      int_flux_k = int_flux_i + int_flux_j
      if ( peak_flux_i > peak_flux_j ):
        peak_flux_k = peak_flux_i
        pos_k = pos_i
      else:
        peak_flux_k = peak_flux_j
        pos_k = pos_j
      source_k = [ k, pos_k, peak_flux_k, int_flux_k, shape_k, clean_pos_k ]
      temp_double_list.append( source_k )
  
  # remove multiples from original source list
  for s in range( len( close_double_list ) ):
    source_i = close_double_list[ s ]
    try:
      source_list.remove( source_i )
    except:
      pass
  
  # add double and single sources to temporary source list
  temp_source_list = []
  for s in temp_double_list:
    temp_source_list.append( s )
  for s in source_list:
    [ i, pos, peak_flux, int_flux, shape ] = s
    temp_source_list.append( [ i, pos, peak_flux, int_flux, shape, pos ] )
  if print_info:
    print '... using %s peeling sources' % ( repr( len( temp_source_list ) ) )
  
  # sort list of candidate peeling sources by model component flux
  # TODO: what about double sources with both components in different facets?
  source_list = []
  double_list = []
  for s in temp_source_list:
    [ i, pos, peak_flux, int_flux, shape, clean_pos ] = s
    facet_i = get_facet( facets, i )
    [ bmaj, bmin, bpa ] = convert_beam_size( facet_i, beam = shape, to_pixel = True )
    clean_radius = int( ceil( 
        max( 2. * clean_box_radius, clean_box_ratio * bmaj ) / 2. ) ) + 1
    model_flux = get_model_flux_from_position_area( facet_i, clean_pos, clean_radius )
    source_list.append( [ i, pos, peak_flux, int_flux, model_flux, shape, clean_pos ] )
    if s in temp_double_list:
      double_list.append( [ i, pos, peak_flux, int_flux, model_flux, shape, clean_pos ] )
  if ( sort == 'model' ):
    source_list.sort( cmp = lambda a, b: cmp( b[ 4 ], a[ 4 ] ) )
  elif ( sort == 'total' ):
    source_list.sort( cmp = lambda a, b: cmp( b[ 3 ], a[ 3 ] ) )
  else: # peak
    source_list.sort( cmp = lambda a, b: cmp( b[ 2 ], a[ 2 ] ) )
  
  # associate candidate peeling sources with catalog sources
  assoc_list = []
  match_list = []
  complex_match_list = []
  if use_catalog:
    if ( len( input_catalog ) > 0 ):
      # copy input list
      model_list = []
      for s in input_catalog:
        model_list.append( [ s[ 0 ][ 0 : 2 ], s[ 1 ] ] )
    else:
      radec = get_radec( uv )
      field_size = restore_parameter( uv, 'field_size' )
      search_radius = 0.6 * field_size
      for source in source_list:
        [ i, pos_i ] = source[ 0 : 2 ]
        source_radec = calculate_source_radec( get_facet( facets, i ), pos_i )
        [ r, p ] = calculate_angular_separation( radec, source_radec )
        if ( r > search_radius ):
          search_radius = 1.05 * r
      frequency = get_central_frequency( uv )
      epoch = get_epoch( uv )
      model_list = generate_source_list( radec, search_radius, frequency, epoch = epoch,
          use_nvss = use_nvss, use_wenss = ( not use_nvss ), spectral_index = -0.8 )
      # TODO: use other catalog with more accurate positions ?
  for source in source_list:
    if ( source in double_list ):
      match_string = 'double source'
    else:
      match_string = 'single source'
    [ i, pos_i, peak_flux_i, int_flux_i, model_flux_i, shape_i, clean_box_i ] = source
    source_radec = calculate_source_radec( get_facet( facets, i ), pos_i )
    if use_catalog:
      min_flux = int_flux_i / assoc_ratio
      candidate_list = []
      for model in model_list:
        [ model_radec, model_flux ] = model
        if ( model_flux > min_flux ):
#        if True:
          [ radius, angle ] = calculate_angular_separation( source_radec, model_radec )
          if ( radius < catalog_search_ratio * catalog_resolution / ( 2. * 3600. ) ):
            candidate_list.append( [ model, radius ] )
      if ( len( candidate_list ) == 0 ):
        if print_info:
          print ( '... WARNING: found no catalog match at RADEC %s' % 
              ( radec_to_string( source_radec ) ) + ', using measured position' ) 
        assoc_list.append( [ source, False, [ source_radec, int_flux_i ], 0. ] )
        match_string = match_string + ', no catalog match, used measured position'
        complex_match_list.append( source )
      elif ( len( candidate_list ) > 1 ):
        # handle multiple candidates, sort by flux
        candidate_list.sort( cmp = lambda a, b: cmp( b[ 0 ][ 1 ], a[ 0 ][ 1 ] ) )
        candidate_1 = candidate_list[ 0 ]
        candidate_2 = candidate_list[ 1 ]
        if ( source in double_list ):
          # when original is double, pick the closest of the brightest two candidates
          if print_info:
            print ( '... WARNING: found multiple catalog matches at RADEC %s' % 
                ( radec_to_string( source_radec ) ) + ', using closest' )
          if ( candidate_1[ 1 ] < candidate_2[ 1 ] ):
            assoc_list.append( [ source, True ] + candidate_1 )
            model_list.remove( candidate_1[ 0 ] )
          else:
            assoc_list.append( [ source, True ] + candidate_1 )
            model_list.remove( candidate_2[ 0 ] )
          match_string = match_string + ', multiple catalog match, used closest'
          complex_match_list.append( source )
        # when original is single, merge brightest two candidate positions with flux as weight
        else:
          if print_info:
            print ( '... WARNING: found multiple catalog matches at RADEC %s' %
                ( radec_to_string( source_radec ) ) + ', using average' )
          [ radius, angle ] = calculate_angular_separation( candidate_1[ 0 ][ 0 ],
              candidate_2[ 0 ][ 0 ] )
          model_flux = candidate_1[ 0 ][ 1 ] + candidate_2[ 0 ][ 1 ]
          radius = radius * candidate_2[ 0 ][ 1 ] / model_flux
          model_radec = calculate_offset_position( candidate_1[ 0 ][ 0 ], radius, angle )
          [ radius, angle ] = calculate_angular_separation( source_radec, model_radec )
          assoc_list.append( [ source, True, [ model_radec, model_flux ], radius ] )
          model_list.remove( candidate_1[ 0 ] )
          model_list.remove( candidate_2[ 0 ] )
          match_string = match_string + ', multiple catalog match, used average'
          complex_match_list.append( source )
      else: # ( len( candidate_list ) == 1 )
        if print_info:
          print ( '... found single catalog match for source at RADEC %s' %
              ( radec_to_string( source_radec ) ) )
        assoc_list.append( [ source, True ] + candidate_list[ 0 ] )
        model_list.remove( candidate_list[ 0 ][ 0 ] )
        match_string = match_string + ', single catalog match'
    else:
      assoc_list.append( [ source, False, [ source_radec, int_flux_i ], 0. ] )
      match_string = match_string + ', used measured position'
    match_list.append( match_string )
  
  # if sources are bright enough, try to peel them
  peel_facet_file_name = facet_file_name + '.PEELPB'
  peel_facet_file_name_e = path.expandvars( peel_facet_file_name )
  if file_exists( peel_facet_file_name_e ):
    remove_file( peel_facet_file_name_e )
  peel_facets = get_aips_file( facets.disk, 'PEELPB', facets.klass, - 1, 'MA' )
  
  # make a copy of input UV
  peel_uv = get_aips_file( uv.disk, uv.name, 'PEELPB', - 1, 'UV' )
  call_aips_task( 'MOVE', indata = uv, outdata = peel_uv, 
      userid = get_aips_userid(), opcode = '' )
  while table_exists( peel_uv, 'NI', 0 ):
    peel_uv.zap_table( 'NI', 0 )
  while table_exists( peel_uv, 'OB', 0 ):
    peel_uv.zap_table( 'OB', 0 )
  
  peel_i = 0
  for source_i in source_list:
    if ( force_simple_match and ( source_i in complex_match_list ) ):
      if print_info:
        print '... skipping source %s/%s due to complex catalog match' % (
            repr( source_list.index( source_i ) + 1 ), repr( len( source_list ) ) )
      continue
    old_peel_uv = get_aips_file( peel_uv.disk, peel_uv.name, peel_uv.klass, 
        peel_uv.seq, 'UV' )
    peel_uv = get_aips_file( peel_uv.disk, peel_uv.name, peel_uv.klass, - 1, 'UV' )
    call_aips_task( 'MOVE', indata = old_peel_uv, outdata = peel_uv, 
        userid = get_aips_userid(), opcode = '' )
    
    peel_i = peel_i + 1
    peel_facet_i = get_facet( peel_facets, peel_i )
    peel_beam_i = get_facet_beam( peel_facet_i )
    
    # get source info
    [ i, pos, peak_flux, int_flux, model_flux, shape, clean_pos ] = source_i
    facet_i = get_facet( facets, i )
    source_radec = calculate_source_radec( facet_i, pos )
    clean_radec = calculate_source_radec( facet_i, clean_pos )
    [ bmaj, bmin, bpa ] = convert_beam_size( facet_i, beam = shape, to_pixel = True )
    if print_info:
      if ( pos == clean_pos ):
        print '... peeling single source %s (%s/%s) at RADEC %s' % ( repr( peel_i ),
            repr( source_list.index( source_i ) + 1 ), repr( len( source_list ) ),
            radec_to_string( source_radec ) )
      else:
        print '... peeling double source %s (%s/%s) at RADEC %s' % ( repr( peel_i ),
            repr( source_list.index( source_i ) + 1 ), repr( len( source_list ) ),
            radec_to_string( clean_radec ) )
    
    # define new facet with peel source in center
    # add facet definition to new facet file
    peel_facet_i_radec = source_radec
    if use_catalog:
      for assoc in assoc_list:
        [ assoc_source, assoc_found, [ assoc_radec, assoc_flux ], assoc_radius ] = assoc
        if ( assoc_source == source_i ):
          if assoc_found:
            # set facet center to catalog position
            peel_facet_i_radec = assoc_radec
            if print_info:
              print '... replacing source RADEC %s with catalog RADEC %s' % ( 
                  radec_to_string( source_radec ), radec_to_string( assoc_radec ) )
          break
    peel_facet_i_size = get_image_size( facet_i )
    peel_facet_i_pixel_ref = get_pixel_reference( facet_i )
    add_facet( peel_facet_file_name_e, peel_facet_i_radec, peel_facet_i_size,
        facet_id = peel_i )
    
    # make dirty image of residual facet
    peel_facet_i_pixel_size = get_pixel_size( facet_i, make_absolute = True )
    temp_facet = get_aips_file( peel_facets.disk, 'TEMP', 'ICL001', -1, 'MA' )
    temp_facet_file_name = peel_facet_file_name_e + '.TEMP'
    extract_facet_definitions( peel_facet_file_name_e, [ peel_i ], temp_facet_file_name )
    channel_count = get_channel_count( peel_uv )
    if ( apply_solutions and table_exists( peel_uv, 'SN', 0 ) ):
      sol_switch = 100
      sol_vers = peel_uv.table_highver( 'SN' )
    else:
      sol_switch = -1
      sol_vers = -1
    uv_size = restore_parameter( peel_uv, 'pb_image_size' )
    call_aips_task( 'IMAGR', indata = peel_uv, nchav = channel_count, flagver = -1,
        outdisk = temp_facet.disk, outname = temp_facet.name, outseq = temp_facet.seq,
        outver = 0, in2disk = peel_uv.disk, cellsize = peel_facet_i_pixel_size, 
        imsize = peel_facet_i_size, overlap = 2, cmethod = conversion_method,
        do3dimag = 1, niter = 100000, flux = 100000., boxfile = temp_facet_file_name, 
        nfield = 1, maxpixel = 0, minpatch = peel_facet_i_size[ 0 ] - 1,
        allokay = 0, docalib = sol_switch, gainuse = sol_vers, dotv = 0, gain = gain,
        imagrprm = [ 0 for idx in range( 20 ) ], uvsize = [ uv_size, uv_size ],
        **imagr_params )
    remove_file( temp_facet_file_name )
    temp_beam = get_facet_beam( temp_facet )
    temp_beam.rename( name = peel_beam_i.name, klass = peel_beam_i.klass,
        seq = peel_beam_i.seq )
    temp_facet.rename( name = peel_facet_i.name, klass = peel_facet_i.klass,
        seq = peel_facet_i.seq )
    
    # add clean box to new facet file
    if ( peel_i == 1 ):
      store_parameter( peel_facets, 'facet_file_name', peel_facet_file_name )
    store_parameter( peel_facets, 'facet_count', peel_i )
    clean_radius = int( around( ( min( peel_facet_i_size ) - ( 2 * edge_size ) + 1 ) / 2 ) )
#    peel_facet_i_clean_box_pos = calculate_source_position( peel_facet_i, clean_radec )
    peel_facet_i_clean_box_pos = calculate_source_position( peel_facet_i, peel_facet_i_radec )
    add_circular_clean_box( peel_facet_file_name_e, peel_i, peel_facet_i_clean_box_pos,
        clean_radius )
    
    # look for instances of source in all facets
    overlapping_facet_list = [ i ] + get_facet_overlap( facet_i )
    facet_source_list = find_source_facets( facets, peel_facet_i_radec,
        facet_list = overlapping_facet_list )
    for facet_source in facet_source_list:
      [ j, pos_j ] =  facet_source
      facet_j = get_facet( facets, j )
      pos_j = [ int( around( x ) ) for x in pos_j ]
      
      # extract source model and transfer model components to peel facet
      # add source model back to UV
      temp_facet_file_name = facet_file_name_e + '.TEMP'
      extract_facet_definitions( facet_file_name_e, [ j ], temp_facet_file_name, 
          include_clean_boxes = False, renumber_facets = False )
      add_circular_clean_box( temp_facet_file_name, j, pos_j, clean_radius )
      model_version = extract_model_components( facets, j, temp_facet_file_name )
      if ( not model_table_empty( facet_j ) ):
        transfer_model_components( facet_j, peel_facet_i, model_version = model_version )
        if model_subtracted:
          add_uv = add_model( peel_uv, facets, model_version = model_version,
              facet_list = [ j ], frequency_correction = frequency_correction,
              conversion_method = conversion_method, apply_solutions = apply_solutions )
          peel_uv.zap()
          add_uv.rename( name = peel_uv.name, klass = peel_uv.klass, seq = peel_uv.seq )
      remove_file( temp_facet_file_name )
      facet_j.zap_table( 'CC', model_version )
    
    # combine model component tables into one table
    model_version = combine_model_tables( peel_facet_i )
    model_versions = range( 1, 1 + peel_facet_i.table_highver( 'CC' ) )
    model_versions.remove( model_version )
    for i in model_versions:
      if table_exists( peel_facet_i, 'CC', i ):
        peel_facet_i.zap_table( 'CC', i )
    call_aips_task( 'TACOP', indata = peel_facet_i, outdata = peel_facet_i,
        inext = 'CC', invers = model_version, outvers = 1, ncount = 1 )
    peel_facet_i.zap_table( 'CC', model_version )
    
    # check for model flux
    if ( get_model_flux( peel_facet_i ) == 0. ):
      if print_info:
        print '... skipping source %s/%s which has no model flux' % (
            repr( source_list.index( source_i ) + 1 ), repr( len( source_list ) ) )
      peel_beam_i = get_facet_beam( peel_facet_i )
      peel_facet_i.zap()
      peel_beam_i.zap()
      remove_facet( peel_facet_file_name_e, peel_i )
      peel_i = peel_i - 1
      peel_uv.zap()
      peel_uv = get_aips_file( old_peel_uv.disk, old_peel_uv.name, old_peel_uv.klass,
          old_peel_uv.seq, 'UV' )
      continue
    
    # restore model components to image
    rst_facet = restore_model_components( peel_facets, facet_list = [ peel_i ],
        cross_restore = False, imagr_params = imagr_params )
    peel_facet_i.zap()
    rst_facet.rename( name = peel_facet_i.name, klass = peel_facet_i.klass,
        seq = peel_facet_i.seq )
    
    # determine peeling solutions relative to best available previous solutions
    if relative_solutions:
      cal_uv = get_aips_file( peel_uv.disk, peel_uv.name, 'CAL', - 1, 'UV' )
      # relative to previous facet solutions
      if ( ( max_relativity >= 3 ) and table_exists( facet_i, 'SN', 0 ) ):
        relativity = 3
        call_aips_task( 'TACOP', indata = facet_i, inext = 'SN', invers = 0, ncount = 1,
            outdata = peel_uv, outvers = 0 )
        ref_ant = get_reference_antenna( peel_uv )
        re_sample_solutions( peel_uv, interpolation_method = 'spline' )
        snver = peel_uv.table_highver( 'SN' )
        if normalize_relative_gains:
          call_aips_task( 'SNCOR', indata = peel_uv, opcode = 'NORM',
              snver = snver )
        call_aips_task( 'SPLIT', indata = peel_uv, docalib = 100, gainuse = 0, 
            douvcomp = 0, flagver = -1, outdisk = cal_uv.disk, 
            outclass = cal_uv.klass, outseq = cal_uv.seq )
        peel_uv.zap_table( 'SN', 0 )
        if normalize_relative_gains:
          peel_uv.zap_table( 'SN', 0 )
      # relative to peeling solutions brightest source
      elif ( ( max_relativity >= 2 ) and table_exists( peel_facets, 'SN', 0 ) ):
        relativity = 2
        call_aips_task( 'TACOP', indata = get_facet( peel_facets, 1 ),
            inext = 'SN', invers = 0, ncount = 1, outdata = peel_uv, outvers = 0 )
        ref_ant = get_reference_antenna( peel_uv )
        re_sample_solutions( peel_uv, interpolation_method = 'spline' )
        snver = peel_uv.table_highver( 'SN' )
        if normalize_relative_gains:
          call_aips_task( 'SNCOR', indata = peel_uv, opcode = 'NORM',
              snver = snver )
        call_aips_task( 'SPLIT', indata = peel_uv, docalib = 100, gainuse = 0,
            douvcomp = 0, flagver = -1, outdisk = cal_uv.disk,
            outclass = cal_uv.klass, outseq = cal_uv.seq )
        peel_uv.zap_table( 'SN', 0 )
        if normalize_relative_gains:
          peel_uv.zap_table( 'SN', 0 )
      # relative to selfcal solutions
      elif ( ( max_relativity >= 1 ) and table_exists( peel_uv, 'SN', 0 ) ):
        relativity = 1
        ref_ant = get_reference_antenna( peel_uv )
        re_sample_solutions( peel_uv, interpolation_method = 'spline' )
        snver = peel_uv.table_highver( 'SN' )
        if normalize_relative_gains:
          call_aips_task( 'SNCOR', indata = peel_uv, opcode = 'NORM',
              snver = snver )
        call_aips_task( 'SPLIT', indata = peel_uv, docalib = 100, gainuse = 0,
            douvcomp = 0, flagver = -1, outdisk = cal_uv.disk,
            outclass = cal_uv.klass, outseq = cal_uv.seq )
        if normalize_relative_gains:
          peel_uv.zap_table( 'SN', 0 )
      else:
        relativity = 0
        cal_uv = get_aips_file( peel_uv.disk, peel_uv.name, peel_uv.klass,
            peel_uv.seq, 'UV' )
        ref_ant = reference_antenna
    else:
      relativity = 0
      cal_uv = get_aips_file( peel_uv.disk, peel_uv.name, peel_uv.klass,
          peel_uv.seq, 'UV' )
      ref_ant = reference_antenna
    
    # selfcal and image peel facet
    # use the source model as starting model
    converge = selfcal_image_clean_facet( cal_uv, peel_facets, peel_i,
        facet_file_name = peel_facet_file_name_e, sigma_min = sigma_min,
        do_sdi_clean = False, signal_to_noise_min = signal_to_noise_min,
        snr_limit = snr_limit, re_center_model = use_catalog, print_info = print_info,
        amplitude_interval = amplitude_interval, improvement_limit = improvement_limit,
        convergence_limit = convergence_limit, reference_antenna = ref_ant,
        conversion_method = conversion_method, phase_interval_min = phase_interval_min,
        imagr_params = imagr_params, calib_params = calib_params, add_boxes = add_boxes,
        frequency_correction = frequency_correction, resolve_power = resolve_power,
        offset_ratio_max = offset_ratio_max, sidelobe_rejection = sidelobe_rejection,
        amplitude_noise_factor = amplitude_noise_factor, amplitude_snr = amplitude_snr,
        edge_size = edge_size, clean_box_radius = clean_box_radius,
        normalize_gains = normalize_gains, flux_rejection_ratio = flux_rejection_ratio,
        convolve_size = catalog_resolution )
    
    if ( relativity > 0 ):
      cal_uv.zap()
    
    # check final solution interval
    if ( not converge ):
      if print_info:
        print '... source rejected - peeling diverged'
    else:
      solution_interval = restore_parameter( peel_facet_i, 'solution_interval' )
      if ( solution_interval > phase_interval_max ):
        converge = False
        if print_info:
          print '... source rejected - peeling interval (%5.2f min) too high' % (
                solution_interval )
    if ( not converge ):
      peel_beam_i = get_facet_beam( peel_facet_i )
      peel_facet_i.zap()
      peel_beam_i.zap()
      remove_facet( peel_facet_file_name_e, peel_i )
      peel_i = peel_i - 1
      peel_uv.zap()
      peel_uv = get_aips_file( old_peel_uv.disk, old_peel_uv.name, old_peel_uv.klass,
          old_peel_uv.seq, 'UV' )
      continue
    
    # combine relative peeling solutions with previous solutions solutions
    if ( relativity > 0 ):
      selfcal_version = peel_uv.table_highver( 'SN' )
      if ( relativity == 3 ):
        call_aips_task( 'TACOP', indata = facet_i, inext = 'SN', invers = 0, 
            ncount = 1, outdata = peel_uv, outvers = 0 )
      elif ( relativity == 2 ):
        call_aips_task( 'TACOP', indata = get_facet( peel_facets, 1 ),
            inext = 'SN', invers = 0, ncount = 1, outdata = peel_uv, outvers = 0 )
      ref_ant_1 = get_reference_antenna( peel_uv )
      ref_ant_2 = get_reference_antenna( peel_facet_i )
      if ( ref_ant_1 != ref_ant_2 ):
        re_reference_solutions( peel_uv, ref_ant_2 )
      re_sample_solutions( peel_uv, interpolation_method = 'spline' )
      ref_version = peel_uv.table_highver( 'SN' )
      if normalize_relative_gains:
        call_aips_task( 'SNCOR', indata = peel_uv, opcode = 'NORM', snver = ref_version )
      call_aips_task( 'TACOP', indata = peel_facet_i, inext = 'SN', ncount = 1,
          outdata = peel_uv )
      solution_interval = restore_parameter( peel_facet_i, 'solution_interval' )
      re_sample_solutions( peel_uv, force_reference = True,
        gap_time = 1.5 * solution_interval, add_mirror_points = False )
      combine_solutions( peel_uv, force_match = True, in_version_1 = ref_version )
      if ( reference_antenna == 0 ):
        reference_antenna = ref_ant_2
      elif ( reference_antenna != ref_ant_2 ):
        re_reference_solutions( peel_uv, reference_antenna )
      call_aips_task( 'TACOP', indata = peel_uv, inext = 'SN', ncount = 1,
          outdata = peel_facet_i )
      while (  peel_uv.table_highver( 'SN' ) > selfcal_version ):
        peel_uv.zap_table( 'SN', 0 )
    
    # store info on catalog matching
    store_parameter( peel_facet_i, 'match_type', 
        match_list[ source_list.index( source_i ) ] )
    
    # subtract updated model from UV
    old_peel_uv.zap()
    sub_uv = subtract_model( peel_uv, peel_facets, facet_list = [ peel_i ],
        sigma = 0., apply_solutions = True, keep_solutions = True,
        frequency_correction = frequency_correction, flag_solutions = True,
        conversion_method = conversion_method )
    peel_uv.zap()
    sub_uv.rename( name = old_peel_uv.name, klass = old_peel_uv.klass,
        seq = old_peel_uv.seq )
    peel_uv = get_aips_file( old_peel_uv.disk, old_peel_uv.name, old_peel_uv.klass, 
        old_peel_uv.seq, 'UV' )
    
    # hard limit the number of peel sources
    if ( peel_i >= max_peel_count ):
      break
  
  # store parameters
  if ( peel_i > 0 ):
    store_parameter( peel_facets, 'facet_count', peel_i )
#    store_parameter( peel_facets, 'facet_file_name', peel_facet_file_name )
  else:
    peel_facets = None
  
  # remove solutions from peel uv
  while ( peel_uv.table_highver( 'SN' ) > 0 ):
    peel_uv.zap_table( 'SN', 0 )
  if keep_solutions:
    call_aips_task( 'TACOP', indata = uv, outdata = peel_uv, inext = 'SN', invers = 0,
        outvers = 0, ncount = 1 )
  
  # copy back NI and OB tables
  ni_version = uv.table_highver( 'NI' )
  for v in range( 1, 1 + ni_version ):
    if table_exists( uv, 'NI', v ):
      call_aips_task( 'TACOP', indata = uv, outdata = peel_uv, invers = v,
          outvers = v, inext = 'NI', ncount = 1 )
  ob_version = uv.table_highver( 'OB' )
  for v in range( 1, 1 + ob_version ):
    if table_exists( uv, 'OB', v ):
      call_aips_task( 'TACOP', indata = uv, outdata = peel_uv, invers = v,
          outvers = v, inext = 'OB', ncount = 1 )
  
  return [ peel_uv, peel_facets ]

###############################################################################

def convolve_facets( facets, resolution = 50. ):
# resolution in arcsec
  convl_facets = get_aips_file( facets.disk, 'CONVL', facets.klass, -1, 'MA' )
  for i in range( 1, 1 + restore_parameter( facets, 'facet_count' ) ):
    facet = get_facet( facets, i )
    convl_facet = get_facet( convl_facets, i )
    [ bmaj, bmin, bpa ] = get_beam_size( facet )
    if ( bmaj > resolution ):
      if ( bmin > resolution ):
        call_aips_task( 'MOVE', indata = facet, outdata = convl_facet,
            userid = get_aips_userid() )
        continue
      else:
        bmaj = bmaj * 1.001
        bmin = resolution
    else:
      bmaj = resolution
      bmin = resolution
      bpa = 0.
    call_aips_task( 'CONVL', indata = facet, outdata = convl_facet,
        bmaj = bmaj, bmin = bmin, bpa = bpa, opcode = 'GAUS', doblank = 1 )
    call_aips_task( 'TACOP', indata = facet, outdata = convl_facet,
      inext = 'PS', ncount = 1 )
  return convl_facets

###############################################################################


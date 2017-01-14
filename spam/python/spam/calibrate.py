###############################################################################

# import Python modules
from sys import *
from os import *
from datetime import *
from math import *

# import 3rd party modules
#import numarray as na
from numpy import *

# import user modules
from files import *
from aips import *
from sphere import *
from parameter import *
from skymodel import *
from image import *
from solutions import *
from error import *

###############################################################################

def divide_uv( uv, model_uv, apply_flags = True, flag_version = 0 ):
# flag table is expected to be attached to uv
  
  # apply flags
  if ( apply_flags and table_exists( uv, 'FG', flag_version ) ):
    flag_uv = apply_flag_table( uv )
  else:
    flag_uv = uv
  
  # divide UV by model UV
#  div_uv = get_aips_file( uv.disk, uv.name, 'DIV', -1, 'UV' )
#  call_aips_task( 'DIFUV', indata = flag_uv, in2data = model_uv,
#      outdata = div_uv, optype = 'DIV', solint = 0 )
  
  # adjust model weights to huge value
  # divide UV by model UV
  temp_uv = get_aips_file( uv.disk, uv.name, 'TEMP', - 1, 'UV' )
  call_aips_task( 'WTMOD', indata = model_uv, outdata = temp_uv, aparm = [ 1.e12,0,0 ] )
  div_uv = get_aips_file( uv.disk, uv.name, 'DIV', -1, 'UV' )
  call_aips_task( 'DIFUV', indata = flag_uv, in2data = temp_uv, solint = 0,
      outdata = div_uv, optype = 'DIV' )
  temp_uv.zap()
  
  # cleanup
  if ( flag_uv != uv ):
    flag_uv.zap()
  
  return div_uv

###############################################################################

def calibrate_model( uv, facets, reference_antenna = 0, phase_interval = 0.,
    do_amplitude = False, amplitude_interval = 5., interpolation_method = '',
    apply_solutions = True, keep_solutions = True, solution_version = 0,
    flag_solutions = True, fast_flag = False, keep_flags = True, sigma = 0., 
    model_version = 0, facet_list = [], conversion_method = 'DFT',
    calib_params = {}, frequency_correction = False, snr_limit = 2.,
    normalize_gains = True, sidelobe_rejection = 0. ):
  
  # work copies
  flag_solns = flag_solutions
  flag_uv = uv
  
  # determine solution mode
  if flag_solns:
    if apply_solutions:
      if ( len( facet_list ) > 0 ):
        facet = get_facet( facets, facet_list[ 0 ] )
      else:
        facet = facets
      if table_exists( facet, 'SN', solution_version ):
        sol_mode = 2
      elif table_exists( uv, 'SN', solution_version ):
        sol_mode = 1
      else:
        sol_mode = 0
        flag_solns = False
        
  # determine solution flagging mode
  if ( flag_solns and fast_flag and ( sol_mode == 2 ) ):
    # get selected facet list with non-empty models
    if ( len( facet_list ) > 0 ):
      temp_facet_count = len( facet_list )
      temp_facet_list = [ i for i in facet_list ]
    else:
      try:
        temp_facet_count = restore_parameter( facets, 'facet_count' )
      except:
        temp_facet_count = 1
      temp_facet_list = range( 1, temp_facet_count + 1 )
    for i in range( 1, 1 + max( temp_facet_list ) ):
      if ( i in temp_facet_list ):
        facet = get_facet( facets, i )
        if model_table_empty( facet, model_version = model_version ):
          temp_facet_list.remove( i )
    temp_facet_count = len( temp_facet_list )
    if ( temp_facet_count > 0 ):
      # try to do fast flagging
      try:
        facet_count = restore_parameter( facets, 'facet_count' )
        added_facet_count = restore_parameter( facets, 'added_facet_count' )
      except:
        pass
      else:
        flag_solns = False
        start_j = temp_facet_count - 1
        if ( added_facet_count > 0 ):
          orig_facet_count = facet_count - added_facet_count
          for j in range( temp_facet_count ):
            if ( temp_facet_list[ j ] > orig_facet_count ):
              start_j = max( 0, j - 1 )
              break
        for j in range( start_j, temp_facet_count ):
          if ( j == start_j ):
            flag_version = -1
          else:
            flag_version = 0
          if ( j == temp_facet_count - 1 ):
            apply_flags = True
          else:
            apply_flags = False
          facet = get_facet( facets, temp_facet_list[ j ] )
          flag_uv = flag_bad_solutions( uv, uvim = facet,
              apply_flags = apply_flags, flag_version = flag_version,
              solution_version = solution_version, keep_flags = keep_flags )
  
  # create model uv
  model_uv = make_model_uv( flag_uv, facets, facet_list = facet_list,
      apply_solutions = apply_solutions, keep_solutions = keep_solutions,
#      conversion_method = conversion_method, model_version = model_version,
      conversion_method = 'DFT', model_version = model_version,
      solution_version = solution_version, flag_solutions = flag_solns,
      frequency_correction = frequency_correction, keep_flags = keep_flags,
      flux_scale = 1., sigma = sigma, sidelobe_rejection = sidelobe_rejection )
  
  # calibrate uv data
  out_version = calibrate_uv( flag_uv, model_uv, snr_limit = snr_limit,
      phase_interval = phase_interval, do_amplitude = do_amplitude,
      amplitude_interval = amplitude_interval, calib_params = calib_params,
      interpolation_method = interpolation_method, apply_flags = True,
      reference_antenna = reference_antenna, normalize_gains = normalize_gains )
  
  # clean up
  model_uv.zap()
  if ( flag_uv != uv ):
    call_aips_task( 'TACOP', indata = flag_uv, inext = 'SN', invers = out_version,
         ncount = 1, outdata = uv )
    flag_uv.zap()
    out_version = uv.table_highver( 'SN' )
  
  return out_version

###############################################################################

def calibrate_uv( uv, model_uv, reference_antenna = 0, phase_interval = 0.,
    do_amplitude = False, amplitude_interval = 5., interpolation_method = '',
    apply_flags = True, calib_params = {}, snr_limit = 2., normalize_gains = True ):
  
  # divide model out
  div_uv = divide_uv( uv, model_uv, apply_flags = apply_flags )
  
  # calibrate UV data against 1 Jy point source model
  if normalize_gains:
    norm = 1
  else:
    norm = 0
  if do_amplitude:
    solmode = 'A&P'
    solution_interval = amplitude_interval
  else:
    solmode = 'P'
    solution_interval = phase_interval
  channel_count = get_channel_count( uv )
#  dummy_uv = get_aips_file( uv.disk, uv.name, 'DUMMY', -1, 'UV' )
#  call_aips_task( 'CALIB', indata = div_uv, smodel = [ 1,0,0,0,0,0,0 ],
#      outdata = dummy_uv, cmethod = 'DFT', snver = 0,
#      refant = reference_antenna, solint = solution_interval, solsub = 1,
#      solmin = 1, weightit = 1, soltype = 'L1R', solmode = solmode,
#      aparm = [ 4,0,0,0,0,0, snr_limit, 0,0 ], cparm = [ 0,0,0,0,0,0,0 ],
#      docalib = -1, gainuse = -1, ichansel = [ [ 1, channel_count, 1,1 ] ], 
#      cmodel = 'COMP', normaliz = norm, **calib_params )
#  dummy_uv.zap()
  call_aips_task( 'CALIB', indata = div_uv, cmethod = 'DFT', snver = 0,
      refant = reference_antenna, solint = solution_interval, weightit = 1,
      soltype = 'L1R', solmode = solmode, aparm = [ 4,1,0,0,0,0, snr_limit, 0,0 ],
      cparm = [ 0,0,0,0,0,0,0 ], ichansel = [ [ 1, channel_count, 1,1 ] ], 
      docalib = -1, gainuse = -1, normaliz = norm, **calib_params )
  
  # copy solutions to original UV data
  call_aips_task( 'TACOP', indata = div_uv, inext = 'SN', invers = 0, ncount = 1,
      outdata = uv, outvers = 0 )
  div_uv.zap()
  
  # resample solutions to UV database time grid
  if ( interpolation_method == '' ):
      im = 'spline'
  else:
    im = interpolation_method
  in_version = uv.table_highver( 'SN' )
  out_version = in_version + 1
  try:
    re_sample_solutions( uv, in_version = in_version, out_version = out_version,
        interpolation_method = im, gap_time = solution_interval,
        force_reference = True )
#    uv.zap_table( 'SN', in_version )
  except:
    re_sample_solutions( uv, in_version = in_version, out_version = out_version,
        interpolation_method = 'linear', gap_time = solution_interval,
        force_reference = True )
#    uv.zap_table( 'SN', in_version )
  if ( not do_amplitude ):
    sn_version = uv.table_highver( 'SN' )
    call_aips_task( 'SNCOR', indata = uv, snver = sn_version, opcode = 'NORM' )
  
  return out_version

###############################################################################

def make_model_uv( uv, facets, facet_list = [], sigma = 0., flux_scale = 1.,
    apply_solutions = False, keep_solutions = True, conversion_method = 'DFT',
    model_version = 0, solution_version = 0, frequency_correction = False,
    flag_solutions = True, keep_flags = True, merge_model = False,
    print_info = True, sidelobe_rejection = 0. ):
  # sidelobe_rejection:
  # > 0. rejects |components| < |largest negative * sidelobe_rejection|
  # < 0. rejects components < 0.
  
  # initialise some parameters
  if ( len( facet_list ) > 0 ):
    temp_facet_count = len( facet_list )
    temp_facet_list = [ i for i in facet_list ]
  else:
    try:
      temp_facet_count = restore_parameter( facets, 'facet_count' )
    except:
      temp_facet_count = 1
    temp_facet_list = range( 1, temp_facet_count + 1 )
  if ( sigma > 0. ):
    cpb_noise = restore_parameter( uv, 'cpb_noise' )
    flux_limit = sigma * cpb_noise
  else:
    flux_limit = 0.
  if frequency_correction:
    dish_diameter = restore_parameter( uv, 'dish_diameter' )
  
  # handle different solution possibilities
  sol_switch = 0
  if apply_solutions:
    if table_exists( get_facet( facets, temp_facet_list[ 0 ] ), 'SN', solution_version ):
      sol_switch = 2
    elif table_exists( uv, 'SN', solution_version ):
      sol_switch = 1
  
  # copy and rename selected facets
  if print_info:
    print '... selecting facets'
  temp_facets = get_aips_file( facets.disk, 'TEMP', 'ICL001', -1, 'MA' )
  j = 0
  for i in range( 1, 1 + max( temp_facet_list ) ):
    if ( i in temp_facet_list ):
      facet = get_facet( facets, i )
      if model_table_empty( facet, model_version = model_version ):
        temp_facet_list.remove( i )
        continue
      j = j + 1
      temp_facet = get_facet( temp_facets, j )
      call_aips_task( 'MOVE', indata = facet, outdata = temp_facet,
          userid = get_aips_userid() )
  temp_facet_count = len( temp_facet_list )
  if ( temp_facet_count == 0 ):
    # return empty model
    model_uv = get_aips_file( uv.disk, uv.name, 'MODEL', - 1, 'UV' )
    sub_uv = subtract_uv( uv, uv, apply_flags = False )
    sub_uv.rename( name = model_uv.name, klass = model_uv.klass,
        seq = model_uv.seq )
    return model_uv
  
  # merge, select and scale model components
  if print_info:
    print '... scaling model components'
  j = 1
  for i in range( 1, 1 + temp_facet_count ):
    temp_facet = get_facet( temp_facets, i )
    call_aips_task( 'CCMRG', indata = temp_facet, invers = model_version,
        outvers = 0 )
    threshold = abs( get_model_minimum( temp_facet ) ) * sidelobe_rejection
    threshold = max( threshold, flux_limit )
    flux_limit = 0.
    if ( ( sol_switch == 2 ) and ( j > 1 ) ):
      scale_model_flux( temp_facet, -flux_scale, threshold = threshold )
    else:
      scale_model_flux( temp_facet, flux_scale, threshold = threshold )
    if model_table_empty( temp_facet ):
      temp_facet.zap()
      temp_facet_count = temp_facet_count - 1
      continue
    if ( i != j ):
      new_facet = get_facet( temp_facets, j )
      temp_facet.rename( name = new_facet.name, klass = new_facet.klass,
          seq = new_facet.seq )
    j = j + 1
  if ( temp_facet_count == 0 ):
    # return empty model
    model_uv = get_aips_file( uv.disk, uv.name, 'MODEL', - 1, 'UV' )
    sub_uv = subtract_uv( uv, uv, apply_flags = False )
    sub_uv.rename( name = model_uv.name, klass = model_uv.klass,
        seq = model_uv.seq )
    return model_uv
  
  # apply no solution
  if print_info:
    print '... creating visibilities from model components'
  model_uv = get_aips_file( uv.disk, uv.name, 'MODEL', - 1, 'UV' )
  if ( sol_switch == 0 ):
    if frequency_correction:
      call_aips_task( 'OOSUB', indata = uv, in2data = temp_facets,
          invers = 0, nmaps = temp_facet_count, outdata = model_uv,
          flux = flux_limit, cmodel = 'COMP', cmethod = conversion_method,
          bparm = [ dish_diameter, 0, 0 ], opcode = 'MODL', factor = 1. )
    else:
      call_aips_task( 'UVSUB', indata = uv, in2data = temp_facets,
          invers = 0, nmaps = temp_facet_count, outdata = model_uv,
          flux = flux_limit, cmodel = 'COMP', cmethod = conversion_method,
          opcode = 'MODL', factor = 1. )
  
  # apply single solution for all facets
  elif ( sol_switch == 1 ):
    # flag bad solutions (important for DIFUV)
    if flag_solutions:
      flag_uv = flag_bad_solutions( uv, solution_version = solution_version,
          flag_version = -1, apply_flags = True, keep_flags = keep_flags )
    else:
      flag_uv = uv
    # generate model UV data set
    if frequency_correction:
      call_aips_task( 'OOSUB', indata = flag_uv, in2data = temp_facets,
          invers = 0, nmaps = temp_facet_count, outdata = model_uv,
          flux = flux_limit, cmodel = 'COMP', cmethod = conversion_method,
          bparm = [ dish_diameter, 0, 0 ], opcode = 'MODL', factor = 1. )
    else:
      call_aips_task( 'UVSUB', indata = flag_uv, in2data = temp_facets,
          invers = 0, nmaps = temp_facet_count, outdata = model_uv,
          flux = flux_limit, cmodel = 'COMP', cmethod = conversion_method,
          opcode = 'MODL', factor = 1. )
    if ( flag_uv != uv ):
      flag_uv.zap()
    # apply inverted SN to model UV
    call_aips_task( 'CLINV', indata = uv, inext = 'SN', invers = solution_version,
        outdata = model_uv, outvers = 0 )
#    snver = model_uv.table_highver( 'SN' ) + 1
#    call_aips_task( 'TACOP', indata = uv, inext = 'SN', invers = solution_version,
#        outdata = model_uv, outvers = snver )
#    call_aips_task( 'SNCOR', indata = model_uv, snver = snver, opcode = 'PHNEG' )
    cal_model_uv = get_aips_file( uv.disk, uv.name, 'CALMDL', - 1, 'UV' )
    call_aips_task( 'SPLIT', indata = model_uv, docalib = 100, gainuse = 0,
        douvcomp = 0, flagver = - 1, outdisk = cal_model_uv.disk,
        outclass = cal_model_uv.klass, outseq = cal_model_uv.seq )
    model_uv.zap()
    cal_model_uv.rename( name = model_uv.name, klass = model_uv.klass,
        seq = model_uv.seq )
  
  # apply solution per facet
  else: # ( sol_switch == 2 )
    # define some names
    temp_uv = get_aips_file( uv.disk, uv.name, 'TEMP', - 1, 'UV' )
    cal_model_uv = get_aips_file( uv.disk, uv.name, 'CALMDL', - 1, 'UV' )
    # flag bad solutions (important for DIFUV)
    if flag_solutions:
      for j in range( 1, 1 + temp_facet_count ):
        if ( j == 1 ):
          flag_version = -1
        else:
          flag_version = 0
        if ( j == temp_facet_count ):
          apply_flags = True
        else:
          apply_flags = False
        temp_facet = get_facet( temp_facets, j )
        flag_uv = flag_bad_solutions( uv, uvim = temp_facet,
            apply_flags = apply_flags, flag_version = flag_version,
            solution_version = solution_version, keep_flags = keep_flags )
    else:
      flag_uv = get_aips_file( uv.disk, uv.name, 'FLAG', -1, 'UV' )
      call_aips_task( 'MOVE', indata = uv, outdata = flag_uv,
          userid = get_aips_userid() )
    # remove larger tables to speed up things
    while table_exists( flag_uv, 'NI', 0 ):
      flag_uv.zap_table( 'NI', 0 )
    while table_exists( flag_uv, 'OB', 0 ):
      flag_uv.zap_table( 'OB', 0 )
    while table_exists( flag_uv, 'SN', 0 ):
      flag_uv.zap_table( 'SN', 0 )
    # cycle over facets
    for j in range( 1, 1 + temp_facet_count ):
      temp_facet = get_facet( temp_facets, j )
      if model_table_empty( temp_facet ):
        continue
      # generate model UV data set per facet
      if frequency_correction:
        call_aips_task( 'OOSUB', indata = flag_uv, in2data = temp_facet,
            invers = 0, nmaps = 1, outdata = temp_uv, flux = flux_limit,
            cmodel = 'COMP', cmethod = conversion_method, opcode = 'MODL',
            bparm = [ dish_diameter, 0, 0 ], factor = 1. )
      else:
        call_aips_task( 'UVSUB', indata = flag_uv, in2data = temp_facet,
            invers = 0, nmaps = 1, outdata = temp_uv, flux = flux_limit,
            cmodel = 'COMP', cmethod = conversion_method, opcode = 'MODL',
            factor = 1. )
      # apply inverted SN to model UV
      call_aips_task( 'CLINV', indata = temp_facet, inext = 'SN',
          invers = 0, outdata = temp_uv, outvers = 0 )
      call_aips_task( 'SPLIT', indata = temp_uv, docalib = 100, gainuse = 0,
          flagver = -1, douvcomp = 0, outdisk = cal_model_uv.disk,
          outclass = cal_model_uv.klass, outseq = cal_model_uv.seq )
      temp_uv.zap()
      # combine model UV from different facets
      if ( j == 1 ):
        cal_model_uv.rename( name = model_uv.name, klass = model_uv.klass,
            seq = model_uv.seq )
        cal_model_uv = get_aips_file( uv.disk, uv.name, 'CALMDL', -1, 'UV' )
      else:
        # subtract model UV from UV
        sub_uv = subtract_uv( model_uv, cal_model_uv, apply_flags = False )
        model_uv.zap()
        cal_model_uv.zap()
        sub_uv.rename( name = model_uv.name, klass = model_uv.klass,
            seq = model_uv.seq )
        if print_info:
          fraction_done = float( j ) / float( temp_facet_count )
          print '...... %d percent done' % ( int( round( 100. * fraction_done ) ) )
    # copy back large tables (SN tables are handled below)
    flag_uv.zap()
    for version in range( 1, 1 + uv.table_highver( 'NI' ) ):
      if table_exists( uv, 'NI', version ):
        call_aips_task( 'TACOP', indata = uv, inext = 'NI', ncount = 1,
            invers = version, outdata = model_uv, outvers = version )
    for version in range( 1, 1 + uv.table_highver( 'OB' ) ):
      if table_exists( uv, 'OB', version ):
        call_aips_task( 'TACOP', indata = uv, inext = 'OB', ncount = 1,
              invers = version, outdata = model_uv, outvers = version )
  
  # remove facet copies
  for j in range( 1, 1 + temp_facet_count ):
    temp_facet = get_facet( temp_facets, j )
    temp_facet.zap()
  
  # copy or remove SN tables when requested
  if ( keep_solutions and ( not table_exists( model_uv, 'SN', 0 ) ) ):
    for version in range( 1, 1 + uv.table_highver( 'SN' ) ):
      if table_exists( uv, 'SN', version ):
        call_aips_task( 'TACOP', indata = uv, inext = 'SN', ncount = 1,
            invers = version, outdata = model_uv, outvers = version )
  if ( ( not keep_solutions ) and table_exists( model_uv, 'SN', 0 ) ):
    while table_exists( model_uv, 'SN', 0 ):
      model_uv.zap_table( 'SN', 0 )
  
  return model_uv

###############################################################################

def subtract_uv( uv, model_uv, apply_flags = True, flag_version = 0 ):
# flag table is expected to be attached to uv
  
  # apply flags
  if ( apply_flags and table_exists( uv, 'FG', flag_version ) ):
    flag_uv = apply_flag_table( uv )
  else:
    flag_uv = uv
  
  # subtract model UV from UV
  # adjust weights to original value
#  temp_uv = get_aips_file( uv.disk, uv.name, 'TEMP', - 1, 'UV' )
#  sub_uv = get_aips_file( uv.disk, uv.name, 'SUB', - 1, 'UV' )
#  call_aips_task( 'DIFUV', indata = flag_uv, in2data = model_uv,
#      outdata = temp_uv, solint = 0 )
#  call_aips_task( 'WTMOD', indata = temp_uv, outdata = sub_uv, aparm = [ 2,0,0 ] )
#  temp_uv.zap()
  
  # adjust model weights to huge value
  # subtract model UV from UV
  temp_uv = get_aips_file( uv.disk, uv.name, 'TEMP', - 1, 'UV' )
  call_aips_task( 'WTMOD', indata = model_uv, outdata = temp_uv, aparm = [ 1.e12,0,0 ] )
  sub_uv = get_aips_file( uv.disk, uv.name, 'SUB', -1, 'UV' )
  call_aips_task( 'DIFUV', indata = flag_uv, in2data = temp_uv, solint = 0,
      outdata = sub_uv )
  temp_uv.zap()
  
  # cleanup
  if ( flag_uv != uv ):
    flag_uv.zap()
  
  return sub_uv

###############################################################################

def subtract_model( uv, facets, facet_list = [], sigma = 0., merge_model = False,
    apply_solutions = True, keep_solutions = True, conversion_method = 'DFT',
    model_version = 0, solution_version = 0, frequency_correction = False,
    flag_solutions = True, fast_flag = False, keep_flags = True,
    print_info = True, sidelobe_rejection = 0. ):
  
  # work copies
  flag_solns = flag_solutions
  flag_uv = uv
  
  # determine solution mode
  if apply_solutions:
    if ( len( facet_list ) > 0 ):
      facet = get_facet( facets, facet_list[ 0 ] )
    else:
      facet = facets
    if table_exists( facet, 'SN', solution_version ):
      sol_mode = 2
    elif table_exists( uv, 'SN', solution_version ):
      sol_mode = 1
    else:
      sol_mode = 0
      flag_solns = False
  else:
    sol_mode = 0
    flag_solns = False
        
  # determine solution flagging mode
  if ( flag_solns and fast_flag and ( sol_mode == 2 ) ):
    if print_info:
      print '... flagging bad solutions'
    # get selected facet list with non-empty models
    if ( len( facet_list ) > 0 ):
      temp_facet_count = len( facet_list )
      temp_facet_list = [ i for i in facet_list ]
    else:
      try:
        temp_facet_count = restore_parameter( facets, 'facet_count' )
      except:
        temp_facet_count = 1
      temp_facet_list = range( 1, temp_facet_count + 1 )
    for i in range( 1, 1 + max( temp_facet_list ) ):
      if ( i in temp_facet_list ):
        facet = get_facet( facets, i )
        if model_table_empty( facet, model_version = model_version ):
          temp_facet_list.remove( i )
    temp_facet_count = len( temp_facet_list )
    if ( temp_facet_count > 0 ):
      # try to do fast flagging
      try:
        facet_count = restore_parameter( facets, 'facet_count' )
        added_facet_count = restore_parameter( facets, 'added_facet_count' )
      except:
        pass
      else:
        flag_solns = False
        start_j = temp_facet_count - 1
        if ( added_facet_count > 0 ):
          orig_facet_count = facet_count - added_facet_count
          for j in range( temp_facet_count ):
            if ( temp_facet_list[ j ] > orig_facet_count ):
              start_j = max( 0, j - 1 )
              break
        for j in range( start_j, temp_facet_count ):
          if ( j == start_j ):
            flag_version = -1
          else:
            flag_version = 0
          if ( j == temp_facet_count - 1 ):
            apply_flags = True
          else:
            apply_flags = False
          facet = get_facet( facets, temp_facet_list[ j ] )
          flag_uv = flag_bad_solutions( uv, uvim = facet,
              apply_flags = apply_flags, flag_version = flag_version,
              solution_version = solution_version, keep_flags = keep_flags )
          flag_solns = False
  
  # create model uv
  if print_info:
    print '... creating source model data'
  model_uv = make_model_uv( flag_uv, facets, facet_list = facet_list,
    apply_solutions = apply_solutions, keep_solutions = keep_solutions,
#    conversion_method = conversion_method, model_version = model_version,
    conversion_method = 'DFT', model_version = model_version, # bug fix
    solution_version = solution_version, flag_solutions = flag_solns,
    frequency_correction = frequency_correction, keep_flags = flag_solns,
    flux_scale = 1., sigma = sigma, merge_model = merge_model,
    print_info = print_info, sidelobe_rejection = sidelobe_rejection )
  
  # subtract model uv from flagged uv
  if ( model_uv is None ):
#    sub_uv = None
    sub_uv = get_aips_file( uv.disk, uv.name, 'SUB', -1, 'UV' )
    call_aips_task( 'MOVE', indata = uv, outdata = sub_uv,
        userid = get_aips_userid() )
  else:
    if print_info:
      print '... subtracting model from real data'
    sub_uv = subtract_uv( flag_uv, model_uv, apply_flags = flag_solns )
    model_uv.zap()
    
    # when requested delete solutions
    if ( not keep_solutions ):
      while table_exists( sub_uv, 'SN', 0 ):
        sub_uv.zap_table( 'SN', 0 )
    
    # when requested, keep flags
    if flag_solns:
      if ( not keep_flags ):
        flag_uv.zap_table( 'FG', 0 )
      elif ( flag_uv != uv ):
        call_aips_task( 'TACOP', indata = flag_uv, outdata = uv, inext = 'FG',
            ncount = 1 )
  
  # cleanup
  if ( flag_uv != uv ):
    flag_uv.zap()
  
  return sub_uv

###############################################################################

def add_model( uv, facets, facet_list = [], sigma = 0., merge_model = False,
    apply_solutions = True, keep_solutions = True, conversion_method = 'DFT',
    model_version = 0, solution_version = 0, frequency_correction = False,
    flag_solutions = True, fast_flag = False, keep_flags = True,
    sidelobe_rejection = 0. ):
  
  # work copies
  flag_solns = flag_solutions
  flag_uv = uv
  
  # determine solution mode
  if flag_solns:
    if apply_solutions:
      if ( len( facet_list ) > 0 ):
        facet = get_facet( facets, facet_list[ 0 ] )
      else:
        facet = facets
      if table_exists( facet, 'SN', solution_version ):
        sol_mode = 2
      elif table_exists( uv, 'SN', solution_version ):
        sol_mode = 1
      else:
        sol_mode = 0
        flag_solns = False
        
  # determine solution flagging mode
  if ( flag_solns and fast_flag and ( sol_mode == 2 ) ):
    # get selected facet list with non-empty models
    if ( len( facet_list ) > 0 ):
      temp_facet_count = len( facet_list )
      temp_facet_list = [ i for i in facet_list ]
    else:
      try:
        temp_facet_count = restore_parameter( facets, 'facet_count' )
      except:
        temp_facet_count = 1
      temp_facet_list = range( 1, temp_facet_count + 1 )
    for i in range( 1, 1 + max( temp_facet_list ) ):
      if ( i in temp_facet_list ):
        facet = get_facet( facets, i )
        if model_table_empty( facet, model_version = model_version ):
          temp_facet_list.remove( i )
    temp_facet_count = len( temp_facet_list )
    if ( temp_facet_count > 0 ):
      # try to do fast flagging
      try:
        facet_count = restore_parameter( facets, 'facet_count' )
        added_facet_count = restore_parameter( facets, 'added_facet_count' )
      except:
        pass
      else:
        flag_solns = False
        start_j = temp_facet_count - 1
        if ( added_facet_count > 0 ):
          orig_facet_count = facet_count - added_facet_count
          for j in range( temp_facet_count ):
            if ( temp_facet_list[ j ] > orig_facet_count ):
              start_j = max( 0, j - 1 )
              break
        for j in range( start_j, temp_facet_count ):
          if ( j == start_j ):
            flag_version = -1
          else:
            flag_version = 0
          if ( j == temp_facet_count - 1 ):
            apply_flags = True
          else:
            apply_flags = False
          facet = get_facet( facets, temp_facet_list[ j ] )
          flag_uv = flag_bad_solutions( uv, uvim = facet,
              apply_flags = apply_flags, flag_version = flag_version,
              solution_version = solution_version, keep_flags = keep_flags )
  
  # create inverted model uv
  imodel_uv = make_model_uv( flag_uv, facets, facet_list = facet_list,
    apply_solutions = apply_solutions, keep_solutions = keep_solutions,
#    conversion_method = conversion_method, model_version = model_version,
    conversion_method = 'DFT', model_version = model_version, # bug fix
    solution_version = solution_version, flag_solutions = flag_solns,
    frequency_correction = frequency_correction, keep_flags = keep_flags,
    flux_scale = -1., sigma = sigma, merge_model = merge_model,
    sidelobe_rejection = sidelobe_rejection )
  
  # subtract inverted model uv from flagged uv ( = add model )
  add_uv = get_aips_file( uv.disk, uv.name, 'ADD', -1, 'UV' )
  if ( imodel_uv is None ):
    call_aips_task( 'MOVE', indata = uv, outdata = add_uv,
        userid = get_aips_userid() )
  else:
    sub_uv = subtract_uv( flag_uv, imodel_uv, apply_flags = True )
    sub_uv.rename( name = add_uv.name, klass = add_uv.klass, seq = add_uv.seq )
    imodel_uv.zap()
    
    # when requested delete solutions
    if ( not keep_solutions ):
      while table_exists( add_uv, 'SN', 0 ):
        add_uv.zap_table( 'SN', 0 )
  
  # cleanup
  if ( flag_uv != uv ):
    flag_uv.zap()
  
  return add_uv

###############################################################################

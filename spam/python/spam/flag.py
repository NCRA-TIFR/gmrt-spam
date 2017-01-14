###############################################################################

# import Python modules
from sys import *
from os import *
from datetime import *
from math import *
import pdb

# import user modules
from files import *
from aips import *
from acalc import *
from parameter import *
from solutions import *
from plot import *
from error import *

###############################################################################

def add_flags( uv, flags = [], flag_version = 0, use_uvflg = False,
    print_info = False ):
  if ( flag_version > 0 ):
    flagver = flag_version
  else:
    flagver = uv.table_highver( 'FG' ) + 1
  while ( uv.table_highver( 'HI' ) > 1 ):
    uv.zap_table( 'HI', 0 )
  if print_info:
    print '... adding %d flag commands to table version %d' % ( len( flags ), flagver )
  if use_uvflg:
    call_aips_task( 'TACOP', indata = uv, outdata = uv, inext = 'HI', ncount = 1,
        invers = 1, outvers = 2 )
    for uvflg_params in flags:
      call_aips_task( 'UVFLG', indata = uv, outfgver = flagver, opcode = 'FLAG',
           dohist = -1, **uvflg_params )
    uv.zap_table( 'HI', 1 )
    call_aips_task( 'TACOP', indata = uv, outdata = uv, inext = 'HI', ncount = 1,
        invers = 2, outvers = 1 )
    uv.zap_table( 'HI', 2 )
    return
  uv_sources = [ uv.header.object ]
  if table_exists( uv, 'SU', 0 ):
    uv_sources.append( uv.sources )
  stokes_table = [ 'YX','XY','YY','XX', 'LR','RL','LL','RR', 'HI', 'I','Q','U','V' ]
  uv_stokes = [ ( stokes_table.index( s ) - 8 ) for s in uv.stokes ]
  uv_stokes_min = min( uv_stokes )
  uv_stokes_max = max( uv_stokes )
  if ( ( uv_stokes_min >= 1 ) and ( uv_stokes_max <= 4 ) ):
    uv_stokes_offset = 1
  elif ( ( uv_stokes_min >= -4 ) and ( uv_stokes_max <= -1 ) ):
    uv_stokes_offset = -1
  elif ( ( uv_stokes_min >= -8 ) and ( uv_stokes_max <= -5 ) ):
    uv_stokes_offset = -5
  else:
    raise error( 'unknown UV stokes range "%s"' % ( repr( uv.stokes ) ) )
  if table_exists( uv, 'FG', flagver ):
    flag_table = wizardry( uv ).table( 'FG', flagver )
  else:
    flag_table = new_table( uv, 'FG', version = flagver )
  flag_row = new_table_row( flag_table )
  for keywords in flags:
    # put in row defaults
    flag_row.reason = ''
    flag_row.subarray = 0
    flag_row.chans = [ 0,0 ]
    flag_row.time_range = [ 0.,9999. ]
    flag_row.freq_id = -1
    flag_row.ifs = [ 1,1 ]
    flag_row.pflags = [ 1,1,1,1 ]
    flag_row.source = 0
    flag_row.ants = [ 0,0 ]
    source_list = []
    antennas = []
    baseline = []
    # expand sources and baseline values
    for key in keywords.keys():
      keyword = keywords[ key ]
      if ( key == 'sources' ):
        if ( len( keyword ) > 0 ):
          if ( len( keyword[ 0 ] ) > 0 ):
            for word in keyword:
              if ( len( word ) > 0 ):
                if ( not word in uv_sources ):
                  raise error( 'unknown source "%s"' % ( word ) )
                source_list.append( uv_sources.index( word ) )
              else:
                break
      elif ( key == 'antennas' ):
        if ( len( keyword ) > 0 ):
          if ( keyword[ 0 ] > 0 ):
            for word in keyword:
              if ( word > 0 ):
                antennas.append( word )
              else:
                break
      elif ( key == 'baseline' ):
        if ( len( keyword ) > 0 ):
          if ( keyword[ 0 ] > 0 ):
            for word in keyword:
              if ( word > 0 ):
                baseline.append( word )
              else:
                break
    if ( len( source_list ) == 0 ):
      source_list.append( 0 )
    if ( ( len( antennas ) == 0 ) and ( len( baseline ) == 0 ) ):
      ants_list = [ [ 0,0 ] ]
    elif ( ( len( antennas ) > 0 ) and ( len( baseline ) == 0 ) ):
      ants_list = [ [ ant, 0 ] for ant in antennas ]
    elif ( ( len( antennas ) == 0 ) and ( len( baseline ) > 0 ) ):
      ants_list = [ [ ant, 0 ] for ant in baseline ]
    else:
      ants_list = []
      for ant1 in antennas:
        for ant2 in baseline:
          ants_list.append( [ ant1, ant2 ] )
    # fill in column values and add flags to table
    for key in keywords.keys():
      keyword = keywords[ key ]
      if ( key in [ 'antennas', 'baseline', 'sources' ] ):
        continue
      elif ( key == 'reason' ):
        reason = keyword.strip()
        flag_row.reason = reason[ 0 : min( len( reason ), 24 ) ]
      elif ( key == 'subarray' ):
        flag_row.subarray = int( keyword )
      elif ( key == 'bchan' ):
        flag_row.chans[ 0 ] = int( keyword )
      elif ( key == 'echan' ):
        flag_row.chans[ 1 ] = int( keyword )
      elif ( key == 'timerang' ):
        flag_row.time_range[ 0 ] = float32( dhms_to_time( keyword[ 0 : 4 ] ) )
        flag_row.time_range[ 1 ] = float32( dhms_to_time( keyword[ 4 : 8 ] ) )
      elif ( key == 'freqid' ):
        flag_row.freq_id = int( keyword )
      elif ( key == 'bif' ):
        flag_row.ifs[ 0 ] = int( keyword )
      elif ( key == 'eif' ):
        flag_row.ifs[ 1 ] = int( keyword )
      elif ( key == 'stokes' ):
        # interpret different stokes options
        if ( uv_stokes_offset == 1 ):
          if ( len( keyword ) > 0 ):
            if ( keyword == 'I' ):
              flag_row.pflags = [ 1,0,0,0 ]
            elif ( keyword == 'Q' ):
              flag_row.pflags = [ 0,1,0,0 ]
            elif ( keyword == 'U' ):
              flag_row.pflags = [ 0,0,1,0 ]
            elif ( keyword == 'V' ):
              flag_row.pflags = [ 0,0,0,1 ]
            elif ( keyword == 'IV' ):
              flag_row.pflags = [ 1,0,0,1 ]
            elif ( keyword == 'IQU' ):
              flag_row.pflags = [ 1,1,1,0 ]
            elif ( ( keyword == 'IQUV' ) or ( keyword == 'FULL' ) ):
              pass
            else:
              raise error( 'cannot apply stokes flag "%s" to current UV data' %
                  ( keyword ) )
        elif ( uv_stokes_offset == -1 ):
          if ( len( keyword ) > 0 ):
            if ( keyword == 'RR' ):
              flag_row.pflags = [ 1,0,0,0 ]
            elif ( keyword == 'LL' ):
              flag_row.pflags = [ 0,1,0,0 ]
            elif ( keyword == 'RL' ):
              flag_row.pflags = [ 0,0,1,0 ]
            elif ( keyword == 'LR' ):
              flag_row.pflags = [ 0,0,0,1 ]
            elif ( ( keyword == 'RRLL' ) or ( keyword == 'HALF' ) ):
              flag_row.pflags = [ 1,1,0,0 ]
            elif ( ( keyword == 'RLLR' ) or ( keyword == 'CROS' ) ):
              flag_row.pflags = [ 0,0,1,1 ]
            elif ( ( keyword == 'RLRL' ) or ( keyword == 'FULL' ) ):
              pass
            else:
              raise error( 'cannot apply stokes flag "%s" to current UV data' %
                  ( keyword ) )
        else: # ( uv_stokes_offset == -5 ):
          if ( len( keyword ) > 0 ):
            if ( keyword == 'XX' ):
              flag_row.pflags = [ 1,0,0,0 ]
            elif ( keyword == 'YY' ):
              flag_row.pflags = [ 0,1,0,0 ]
            elif ( keyword == 'XY' ):
              flag_row.pflags = [ 0,0,1,0 ]
            elif ( keyword == 'YX' ):
              flag_row.pflags = [ 0,0,0,1 ]
            elif ( ( keyword == 'XXYY' ) or ( keyword == 'HALF' ) ):
              flag_row.pflags = [ 1,1,0,0 ]
            elif ( ( keyword == 'XYYX' ) or ( keyword == 'CROS' ) ):
              flag_row.pflags = [ 0,0,1,1 ]
            elif ( ( keyword == 'XYXY' ) or ( keyword == 'FULL' ) ):
              pass
            else:
              raise error( 'cannot apply stokes flag "%s" to current UV data' %
                  ( keyword ) )
      else:
        raise error( 'unknown flag table column "%s"' % ( key ) )
    for source in source_list:
      flag_row.source = source
      for ants in ants_list:
        flag_row.ants = ants
        flag_table.append( flag_row )
  flag_table.close()
  return

###############################################################################

def flag_uv_data( uv, flags = [], keep_solutions = True, flag_version = 0 ):
  add_flags( uv, flags = flags, flag_version = flag_version )
  flag_uv = apply_flag_table( uv, version = flag_version,
      keep_solutions = keep_solutions )
  return flag_uv

###############################################################################

def make_time_images( uv, image_size = 1024, apply_solutions = True, 
    solution_version = 0, imagr_params = {}, apply_flags = True, flag_version = 0,
    print_info = False, epsilon = 1.e-8, keep_images = False, time_step = 10. ):
# works best if source model is subtracted from uv
  
  if apply_solutions:
    if table_exists( uv, 'SN', solution_version ):
      docalib = 100
      gainuse = solution_version
    else:
      docalib = -1
      gainuse = -1
  else:
    docalib = -1
    gainuse = -1
  nchav = get_channel_count( uv )
  cell_size = restore_parameter( uv, 'cell_size' )
  uv_size = restore_parameter( uv, 'pb_image_size' )
  time_array = array( get_time_list( uv ) )
#  dtime = time_step / ( 60. * 24. )
  
  # apply previous flags
  if ( apply_flags and table_exists( uv, 'FG', flag_version ) ):
    flag_uv = apply_flag_table( uv, version = flag_version, keep_solutions = True )
  else:
    flag_uv = uv
  
  # get reference noise
  temp_image = get_aips_file( uv.disk, 'T0000', 'IIM001', -1, 'MA' )
  temp_beam = get_facet_beam( temp_image )
  call_aips_task( 'IMAGR', indata = flag_uv, nchav = nchav, nfield = 1, niter = 0,
      cellsize = [ cell_size, cell_size ], do3dimag = 1, outdisk = temp_image.disk,
      outname = temp_image.name, outseq = temp_image.seq, docalib = docalib,
      imsize = [ image_size, image_size ], gainuse = gainuse, flagver = -1,
      dotv = 0, allok = 0, uvsize = [ uv_size, uv_size ], **imagr_params )
  fill_image_edge( temp_image, do_edge_circle = True )
  rms = get_image_rms( temp_image )
  [ temp_avg, temp_noise ] = call_aips_task( 'IMEAN', indata = temp_image,
      pixavg = 0., pixstd = rms, pixrange = [ -5. * rms, 5. * rms  ], 
      outputs = [ 'pixavg', 'pixstd' ] )
  if ( ( temp_noise <= 0. ) or ( temp_noise > 2. * rms ) ):
    if print_info:
      print '... WARNING: histogram noise fit failed, using image RMS instead'
    temp_noise = rms
  if ( not keep_images ):
    temp_image.zap()
  temp_beam.zap()
  reference_noise = temp_noise
  if print_info:
    print '... reference noise = %s' % ( repr( reference_noise ) )
  
  # loop over time
  noise_list = []
  timerang_list = []
  t = 0
#  time = time_array[ 0 ]
  int_time = restore_parameter( uv, 'integration_time' )
  n = 0
  dn = int( ceil( time_step * 60. / int_time ) )
#  while ( time < time_array[ -1 ] ):
  while ( n < len( time_array ) ):
#    sel = awhere( ( time_array > time ) & ( time_array < time + dtime ) )
#    if ( len( sel ) == 0 ):
#      time = time + dtime
#      continue
#    sel = arange( n, n + dn ).reshape( dn, 1 )
    t = t + 1
    time_low = time_to_dhms( time_array[ n ] - 0.5 * int_time / ( 24. * 3600. ) )
    time_high = time_to_dhms( time_array[ min( n + dn, len( time_array ) ) - 1 ] + 
        0.5 * int_time / ( 24. * 3600. ) )
    timerang_list.append( time_low + time_high )
#    time = time + dtime
    n = n + dn
    if print_info:
      print '... time range %s = %s - %s' % ( repr( t ), repr( time_low ), 
          repr( time_high ) )
    call_aips_task( 'UVFLG', indata = flag_uv, outfgver = 0, opcode = 'FLAG', 
        reason = 'test', timerang = time_low + time_high, dohist = -1 )
    flagver = flag_uv.table_highver( 'FG' )
    temp_image = get_aips_file( uv.disk, 'T%04d' % ( t ), 'IIM001', -1, 'MA' )
    temp_beam = get_facet_beam( temp_image )
    call_aips_task( 'IMAGR', indata = flag_uv, nchav = nchav, nfield = 1, niter = 0,
        cellsize = [ cell_size, cell_size ], uvsize = [ uv_size, uv_size ],
        outdisk = temp_image.disk, outname = temp_image.name, outseq = temp_image.seq, 
        imsize = [ image_size, image_size ], do3dimag = 1, docalib = docalib,
        gainuse = gainuse, flagver = flagver, dotv = 0, allok = 0, **imagr_params )
    fill_image_edge( temp_image, do_edge_circle = True )
    rms = get_image_rms( temp_image )
    [ temp_avg, temp_noise ] = call_aips_task( 'IMEAN', indata = temp_image,
        pixavg = 0., pixstd = rms, pixrange = [ -5. * rms, 5. * rms  ], 
        outputs = [ 'pixavg', 'pixstd' ] )
    if ( ( temp_noise <= 0. ) or ( temp_noise > 2. * rms ) ):
      if print_info:
        print '... WARNING: histogram noise fit failed, using image RMS instead'
      temp_noise = rms
    if ( not keep_images ):
      temp_image.zap()
    temp_beam.zap()
    if ( temp_noise == reference_noise ):
      noise_list.append( 1.e9 )
    else:
      noise_list.append( temp_noise )
  min_t = 1 + noise_list.index( min( noise_list ) )
  min_timerang = timerang_list[ min_t - 1 ]
  min_noise = noise_list[ min_t - 1 ]
  if print_info:
    print '... minimum noise time range %s (%s)= %s' % ( repr( min_t ), repr( min_timerang ),
        repr( min_noise ) )
  
  if ( apply_flags and table_exists( uv, 'FG', flag_version ) ):
    flag_uv.zap()
  
  return noise_list

###############################################################################

def read_gmrt_flag_file( uv, flag_file_name, ignore = False, 
    day_offset = -5.5 / 24., skip_flags = [ [ 6,4 ], [ 7,4 ], [ 8,1 ] ] ):
  
  # prepare
  if ( not file_exists( flag_file_name ) ):
    raise error( 'flag file %s does not exist' % flag_file_name )
  if ( not 'GMRT' in uv.header.instrume ):
    raise error( 'UV data is not GMRT data' )
  antenna_names = uv.antennas
  reference_ymd = uv.header.date_obs.split( '-' )
  reference_day = date( int( reference_ymd[ 0 ] ), int( reference_ymd[ 1 ] ),
      int( reference_ymd[ 2 ] ) ).toordinal()
  
  # scan file
  code_lines = [ -1, -1 ]
  antenna_lines = [ -1, -1 ]
  flag_lines = [ -1, -1 ]
  i = 0
  ff = file( flag_file_name, 'r' )
  for line in ff:
    i = i + 1
    if ( '*{' in line ):
      if ( 'FlagCodes' in line ):
        code_lines[ 0 ] = i
      elif ( 'AntennaMap' in line ):
        antenna_lines[ 0 ] = i
      elif ( 'OnlineFlags' in line ):
        flag_lines[ 0 ] = i
    elif ( '*}' in line ):
      if ( 'FlagCodes' in line ):
        code_lines[ 1 ] = i
      elif ( 'AntennaMap' in line ):
        antenna_lines[ 1 ] = i
      elif ( 'OnlineFlags' in line ):
        flag_lines[ 1 ] = i
  ff.close()
  if ( -1 in code_lines ):
    raise error( 'Definition of flag codes is incomplete'  )
  if ( -1 in antenna_lines ):
    raise error( 'Definition of antenna map is incomplete'  )
  if ( flag_lines[ 0 ] == -1 ):
    raise error( 'Definition of online flags is incomplete'  )
  if ( flag_lines[ 1 ] == -1 ):
    print 'WARNING: definition of online flags is incomplete, adding end'
    flag_lines[ 1 ] = i
  
  # read flag codes
  flag_codes = []
  flag_comments = []
  i = 0
  ff = file( flag_file_name, 'r' )
  for line in ff:
    i = i + 1
    if ( ( i <= code_lines[ 0 ] ) or ( i >= code_lines[ 1 ] ) ): 
      continue
    line = line.strip( ' \r\n' ).strip()
    if ( '#' in line ):
      line = line[ : line.index( '#' ) ]
    if ( len( line ) == 0 ):
      continue
    parts = line.split( ':' )
    if ( len( parts ) < 2 ):
      if ignore:
        continue
      raise error( 'flag definition in line %d has unknown format' % ( i ) )
    codes = parts[ 0 ].split( ',' )
    if ( len( codes ) != 2 ):
      if ignore:
        continue
      raise error( 'flag definition in line %d has unknown format' % ( i ) )
    flag_codes.append( [ int( codes[ 0 ] ), int( codes[ 1 ] ) ] )
    comment = line[ line.index( ':' ) + 1 : ].strip()
    flag_comments.append( comment[ 0 : min( len( comment ), 24 ) ] )
  ff.close()
  flag_groups = []
  flag_states = []
  for code in flag_codes:
    if ( not code[ 0 ] in flag_groups ):
      flag_groups.append( code[ 0 ] )
      flag_states.append( [ code[ 1 ] ] )
    else:
      flag_states[ flag_groups.index( code[ 0 ] ) ].append( code[ 1 ] )
  
  # read antenna map
  antenna_codes = []
  antenna_match = []
  i = 0
  ff = file( flag_file_name, 'r' )
  for line in ff:
    i = i + 1
    if ( ( i <= antenna_lines[ 0 ] ) or ( i >= antenna_lines[ 1 ] ) ): 
      continue
    line = line.strip( ' \r\n' )
    if ( '#' in line ):
      line = line[ : line.index( '#' ) ].strip()
    if ( len( line ) == 0 ):
      continue
    while ( True ):
      if ( not '=>' in line ):
        break
      antenna_codes.append( int( line[ : line.index( '=>' ) ] ) )
      name = ( line[ line.index( '=>' ) + 2 : ].split() )[ 0 ]
      antenna_found = False
      for j in range( len( antenna_names ) ):
        if ( name in antenna_names[ j ] ):
          antenna_found = True
          break
      if ( not antenna_found ):
        j = -1
      antenna_match.append( j + 1 )
      line = line[ line.index( name ) + len( name ) : ]
  ff.close()
  
  # process flags
  flag_list = []
  months = [ 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
             'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec' ]
  i = 0
  ff = file( flag_file_name, 'r' )
  for line in ff:
    i = i + 1
    if ( ( i <= flag_lines[ 0 ] ) or ( i >= flag_lines[ 1 ] ) ): 
      continue
    line = line.strip( ' \r\n' )
    if ( len( line ) == 0 ):
      continue
    if ( line[ 0 ] != '$' ):
      continue
    words = line[ 1 : ].split()
    if ( len( words ) != 7 ):
      if ignore:
        continue
      raise error( 'line %d has unknown format' % ( i ) )
    ymd = [ int( words[ 4 ] ), months.index( words[ 1 ] ) + 1, int( words[ 2 ] ) ]
    day = date( ymd[ 0 ], ymd[ 1 ], ymd[ 2 ] ).toordinal() - reference_day
    if ( day < 0 ):
      if ignore:
        continue
      raise error( 'day offset in line %d is less than zero' % ( i ) )
    if ( day > 999 ):
      if ignore:
        continue
      raise error( 'day offset in line %d is larger than 999' % ( i ) )
    hms = words[ 3 ].split( ':' )
    t = dhms_to_time( [ day, int( hms[ 0 ] ), int( hms[ 1 ] ), int( hms[ 2 ] ) ] )
    t = t + day_offset
    mask = int( words[ 5 ], 16 )
    antenna_list = []
    for j in antenna_codes:
      if( mask & 2**j ):
        antenna_list.append( j )
    flags = words[ 6 ].split( '#' )[ 1 : ]
    if ( len( flags ) != len( antenna_list ) ):
      if ignore:
        continue
      raise error( 'no. of flags differs from no. of antennas in line %d' % i )
    for j in range( len( flags ) ):
      k = antenna_codes.index( antenna_list[ j ] )
      switches = flags[ j ].split( '@' )
      for switch in switches:
        codes = switch.split( ',' )
        if ( len( codes ) != 2 ):
          if ignore:
            continue
          raise error( 'flag in line %d has unknown format' % ( i ) )
        flag_group = flag_groups.index( int( codes[ 0 ] ) )
        flag_state = int( codes[ 1 ] )
        if ( not [ flag_group, flag_state ] in flag_codes ):
          if ignore:
            continue
          raise error( 'flag in line %d is not listed in flag codes' % ( i ) )
        flag_list.append( [ t, k, flag_group, flag_state ] )
  ff.close()
  flag_list.sort( cmp = lambda a, b: cmp( a[ 0 ], b[ 0 ] ) )
  
  # generate flag commands
  flags = []
  time_machine = zeros( ( len( antenna_codes ), len( flag_groups ) ),
      dtype = float64 )
  state_machine = zeros( ( len( antenna_codes ), len( flag_groups ) ),
      dtype = int64 )
  for flag in flag_list:
    [ time, antenna, group, state ] = flag
    if ( state != state_machine[ antenna, group ] ):
      if ( state_machine[ antenna, group ] == 0 ):
        # moving from non-error to error state
        state_machine[ antenna, group ] = state
        time_machine[ antenna, group ] = time
      else:
        # moving from error to (non-)error state
        dhms1 = [ int( x ) for x in round_dhms( time_to_dhms( 
            time_machine[ antenna, group ] ) ) ]
        dhms2 = [ int( x ) for x in round_dhms( time_to_dhms( time ) ) ]
        flag_code = [ flag_groups[ group ], state_machine[ antenna, group ] ]
        if ( not flag_code in skip_flags ):
          flags.append( { 'antennas' : [ antenna_match[ antenna ] ],
              'timerang' : dhms1 + dhms2,
              'reason' : flag_comments[ flag_codes.index( flag_code ) ] } )
        state_machine[ antenna, group ] = state
        time_machine[ antenna, group ] = time
#    else:
#      raise error( 'old and new state are the same' )
  
  # handle uncleared states
  sel = awhere( state_machine != 0 )
  if ( len( sel ) > 0 ):
#    raise error( 'some error states were never cleared' )
    for s in sel:
      [ antenna, group ] = s
      dhms1 = [ int( x ) for x in round_dhms( time_to_dhms( 
          time_machine[ antenna, group ] ) ) ]
      dhms2 = [ int( x ) for x in round_dhms( time_to_dhms( 999. ) ) ]
      flag_code = [ flag_groups[ group ], state_machine[ antenna, group ] ]
      if ( not flag_code in skip_flags ):
        flags.append( { 'antennas' : [ antenna_match[ antenna ] ],
            'timerang' : dhms1 + dhms2,
            'reason' : flag_comments[ flag_codes.index( flag_code ) ] } )
  
  return flags

###############################################################################

def write_flag_file( flags, flag_file_name, append = False ):
  
  if file_exists( flag_file_name ):
    if ( not append ):
      remove_file( flag_file_name )
      flag_file = file( flag_file_name, 'w' )
    else:
      flag_file = file( flag_file_name, 'a' )
  else:
    flag_file = file( flag_file_name, 'w' )
  
  for flag in flags:
    line = ""
    for key in flag:
      line = line + key.upper() + '='
      value = flag[ key ]
      if ( ( type( value ) == type( 1 ) ) or ( type( value ) == type( 1. ) ) ):
        line = line + repr( value )
      elif ( type( value ) == type( '1' ) ):
        line = line + '"' + value + '"'
      elif ( type( value ) == type( [ 1 ] ) ):
        for x in value:
          if ( ( type( x ) == type( 1 ) ) or ( type( x ) == type( 1. ) ) ):
            line = line + repr( x )
          elif ( type( x ) == type( '1' ) ):
            line = line + '"' + value + '"'
          elif ( type( x ) == type( [ 1 ] ) ):
            for y in x:
              if ( ( type( x ) == type( 1 ) ) or ( type( x ) == type( 1. ) ) ):
                line = line + repr( x )
              elif ( type( x ) == type( '1' ) ):
                line = line + '"' + value + '"'
              else:
                raise error( 'unknown data type in key ' + repr( key ) )
              line = line + ','
            line = line[ : -1 ]
          else:
            raise error( 'unknown data type in key ' + repr( key ) )
          line = line + ','
        line = line[ : -1 ]
      line = line + ' '
    line = line + '/\n'
    flag_file.write( line )
  flag_file.close()
  
  return

###############################################################################

def flag_bad_solutions( uv, uvim = None, solution_version = 0, flag_version = 0,
    apply_flags = False, keep_flags = False, int_factor = 0.6, 
    weight_cutoff = 0.001 ):
# flagver < 0: create new table; flagver = 0: append to highest table
  
  # if needed, copy solution table from other file
  if ( uvim is None ):
    sn_version = solution_version
  else:
    call_aips_task( 'TACOP', indata = uvim, inext = 'SN', outdata = uv, 
        invers = solution_version, ncount = 1, outvers = 0 )
    sn_version = uv.table_highver( 'SN' )
  
  # pick flag table version
  if ( flag_version > 0 ):
    if table_exists( uv, 'FG', flag_version ):
      fl_version = flag_version
    else:
      fl_version = -1
  elif ( flag_version == 0 ):
    fl_version = uv.table_highver( 'FG' )
    if ( fl_version == 0 ):
      fl_version = -1
  else: #( flag_version < 0 )
    fl_version = -1
  
  # make flag table based on failed solutions
  dt = restore_parameter( uv, 'integration_time' ) * int_factor
  call_aips_task( 'SNFLG', indata = uv, invers = sn_version, optype = 'JUMP',
      cutoff = weight_cutoff, dparm = [ 1000., dt ], flagver = fl_version )
  if ( not uvim is None ):
    uv.zap_table( 'SN', sn_version )
  new_fl_version = uv.table_highver( 'FG' )
  if ( fl_version > 0 ):
    uv.zap_table( 'FG', fl_version )
    call_aips_task( 'TACOP', indata = uv, outdata = uv, inext = 'FG',
        invers = new_fl_version, ncount = 1, outvers = fl_version )
    uv.zap_table( 'FG', new_fl_version )
  else:
    fl_version = new_fl_version
  
  # apply flag table
  flag_uv = None
  if apply_flags:
    flag_uv = apply_flag_table( uv, version = fl_version, keep_solutions = True )
    if ( not keep_flags ):
      uv.zap_table( 'FG', fl_version )
  
  return flag_uv

###############################################################################

def apply_flag_table( uv, version = 0, keep_solutions = True ):
  
  if ( version == 0 ):
    flag_version = uv.table_highver( 'FG' )
  else:
    flag_version = version
  
  # apply flag table
  flag_uv = get_aips_file( uv.disk, uv.name, 'FLAG', - 1, 'UV' )
  call_aips_task( 'UVCOP', indata = uv, outdata = flag_uv, flagver = flag_version )
  
  if ( not keep_solutions ):
    while table_exists( flag_uv, 'SN', 0 ):
      flag_uv.zap_table( 'SN', 0 )
  
  return flag_uv

###############################################################################

def image_cpb_facets( uv, apply_solutions = False, solution_version = 0,
    gain = 0.1, factor = 0., imagr_params = {} ):
  
  # make initial central primary beam facets (no cleaning)
  cell_size = restore_parameter( uv, 'cell_size' )
  uv_size = restore_parameter( uv, 'pb_image_size' )
  cpb_facet_count = restore_parameter( uv, 'cpb_facet_count' )
  cpb_facet_size = restore_parameter( uv, 'pb_facet_size' )
  cpb_facet_file_name = restore_parameter( uv, 'pb_facet_file_name' )
  cpb_facet_file_name_e = path.expandvars( cpb_facet_file_name )
  channel_count = get_channel_count( uv )
  cpb_facets = get_aips_file( uv.disk, 'CPB', 'ICL001', - 1, 'MA' )
  if ( apply_solutions and table_exists( uv, 'SN', solution_version ) ):
    docalib = 100
    gainuse = solution_version
  else:
    docalib = -1
    gainuse = -1
  call_aips_task( 'IMAGR', indata = uv, nchav = channel_count, in2disk = uv.disk,
      outdisk = cpb_facets.disk, overlap = 2, outver = 0, docalib = docalib,
      gainuse = gainuse, outname = cpb_facets.name, outseq = cpb_facets.seq,
      cellsize = [ cell_size, cell_size ], do3dimag = 1, niter = 100000,
      flux = 100000., imsize = [ cpb_facet_size, cpb_facet_size ], dotv = 0,
      boxfile = cpb_facet_file_name_e, cmethod = '', minpatch = cpb_facet_size - 1,
      gain = gain, nfield = cpb_facet_count, maxpixel = 0, factor = factor,
      imagrprm = [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0 ], flagver = -1,
      bcomp = [ 0 for j in range( 64 ) ], uvsize = [ uv_size, uv_size ],
      **imagr_params )
  
  # store parameters  
  store_parameter( cpb_facets, 'facet_count', cpb_facet_count )
  store_parameter( cpb_facets, 'facet_file_name', cpb_facet_file_name )
  
  return cpb_facets

###############################################################################

def combine_facets( uv, facets, facet_list = [], edge_size = 8, edge_weight = 0., 
    save_tables = False, save_info = False, projection = '', radec = None,
     image_size = None, print_info = True ):

  # process inputs
  # select facets to combine
  if ( len( facet_list ) == 0 ):
    used_facet_list = range( 1, 1 + restore_parameter( facets, 'facet_count' ) )
  else:
    used_facet_list = facet_list
  facet_count = len( used_facet_list )
  facet = get_facet( facets, used_facet_list[ 0 ] )
  if ( projection != '' ):
    if ( projection != get_image_projection( facet ) ):
      if ( not projection in [ 'SIN','TAN','ARC','NCP','AIT','GLS','STG','MER' ] ):
        raise error( 'new projection type %s unknown' % ( projection ) )
      cootype = '-' + projection
  else:
    cootype = ''
  
  # set new image center
  if ( radec is None ):
    radec = get_radec( uv )
  coordinates = degdeg_to_hmsdms( radec )
  if ( radec[ 1 ] < 0. ):
    coordinates[ 3 ] = -abs( coordinates[ 3 ] )
    coordinates[ 4 ] = -abs( coordinates[ 4 ] )
    coordinates[ 5 ] = -abs( coordinates[ 5 ] )
  
  # smooth overlapping facets
  if ( ( edge_weight is None ) or ( edge_weight >= 1. ) ):
    weightit = 0
  else:
    facet_size = get_image_size( facet )
    pixel_size = get_pixel_size( facet )
    xc = float( facet_size[ 0 ] - 1 ) / 2.
    yc = float( facet_size[ 1 ] - 1 ) / 2.
    rc = min( [ ( xc - float( edge_size ) ) * pixel_size[ 0 ],
                ( yc - float( edge_size ) ) * pixel_size[ 1 ] ] )
    weightit = ( 1. - max( [ 0.0001, edge_weight ] ) ) / rc
  
  # determine output image size
  if ( image_size is None ):
    used_image_size = restore_parameter( uv, 'pb_image_size' )
    used_image_size = [ used_image_size, used_image_size ]
  else:
    used_image_size = image_size
  
  # make work copies of facets
  # blank circular edge area of facets
  # sum total flux
  temp_facets = get_aips_file( facets.disk, 'TEMP', facets.klass, -1, 'MA' )
  j = 0
  total_clean_flux = 0.
  for i in used_facet_list:
    j = j + 1
    facet = get_facet( facets, i )
    temp_facet = get_facet( temp_facets, j )
    if ( not table_exists( facet, 'CC', 0 ) ):
      raise error( 'CC table missing on input facet %d' % ( i ) )
    total_clean_flux = total_clean_flux + get_model_flux( facet )
    call_aips_task( 'MOVE', indata = facet, outdata = temp_facet,
        userid = get_aips_userid() )
    fill_image( temp_facet, do_edge_circle = True, edge_size = edge_size )
  
  # combine facets into image
  image = get_aips_file( facets.disk, facets.name, 'FLATN', -1, 'MA' )
  call_aips_task( 'FLATN', indata = temp_facets, nfield = facet_count,
      outdata = image, imsize = used_image_size, coordina = coordinates,
      edgskp = edge_size, weightit = weightit, reweight = [ 3, 0.5 ],
       cootype = cootype )
  
  # save model tables and solution tables to image
  while table_exists( image, 'CC', 0 ):
    image.zap_table( 'CC', 0 )
  while table_exists( image, 'SN', 0 ):
    image.zap_table( 'SN', 0 )
  if save_tables:
    for i in used_facet_list:
      facet = get_facet( facets, i )
      call_aips_task( 'TACOP', indata = facet, outdata = image, inext = 'CC',
          invers = 0, ncount = 1, outvers = i )
      if table_exists( facet, 'SN', 0 ):
        call_aips_task( 'TACOP', indata = facet, outdata = image, inext = 'SN',
            invers = 0, ncount = 1, outvers = i )
  else:
    for i in range( 1, 1 + facet_count ):
      facet = get_facet( temp_facets, i )
      transfer_model_components( facet, image )
    combine_model_tables( image, model_versions = range( 1, 1 + facet_count ) )
    for i in range( 1, 1 + facet_count ):
      image.zap_table( 'CC', i )
    call_aips_task( 'TACOP', indata = image, outdata = image, inext = 'CC',
        invers = facet_count + 1, outvers = 1, ncount = 1 )
    image.zap_table( 'CC', facet_count + 1 )
  
  # clean up temporary facets
  for i in range( 1, 1 + facet_count ):
    facet = get_facet( temp_facets, i )
    facet.zap()
  
  # save info to image
  if save_info:
    try:
      clean_box_flux_min = restore_parameter( facets, 'clean_box_flux_min' )
    except:
      pass
    else:
      store_parameter( image, 'clean_box_flux_min', clean_box_flux_min )
    try:
      clean_flux_min = restore_parameter( facets, 'clean_flux_min' )
    except:
      pass
    else:
      store_parameter( image, 'clean_flux_min', clean_flux_min )
    cpb_noise = restore_parameter( uv, 'cpb_noise' )
    store_parameter( image, 'background_rms', cpb_noise )
    store_parameter( image, 'total_clean_flux', total_clean_flux )
    pixel_min = get_image_minimum( image )
    store_parameter( image, 'pixel_minimum', pixel_min[ 0 ] )
    pixel_max = get_image_maximum( image )
    store_parameter( image, 'pixel_maximum', pixel_max[ 0 ] )
  
  return image

###############################################################################

def measure_cpb_noise( uv, facets, facet_list = [], edge_size = 8, keep_image = False,
    print_info = True ):
  
  # process inputs
  if ( len( facet_list ) > 0 ):
    used_facet_list = facet_list
    if ( facet_list[ 0 ] != 1 ):
      facet = get_facet( facets, facet_list[ 0 ] )
      radec = get_radec( facet )
  else:
    radec = None
    cpb_facet_count = restore_parameter( uv, 'cpb_facet_count' )
    used_facet_list = range( 1, 1 + cpb_facet_count )
  cpb_image_size = restore_parameter( uv, 'cpb_image_size' )
  image_size = [ cpb_image_size, cpb_image_size ]

  # combine facets into image
  cpb_image = combine_facets( uv, facets, image_size = image_size, radec = radec,
      facet_list = used_facet_list, edge_size = edge_size, save_tables = False,
      save_info = False, print_info = print_info )
  if ( cpb_image.name != 'CPB' ):
    image = get_aips_file( facets.disk, 'CPB', 'FLATN', -1, 'MA' )
    cpb_image.rename( name = image.name, klass = image.klass, seq = image.seq )
  
  # measure noise in image
  rms = get_image_rms( cpb_image )
  [ cpb_avg, cpb_noise ] = call_aips_task( 'IMEAN', indata = cpb_image,
      pixavg = 0., pixstd = rms, pixrange = [ -5. * rms, 5. * rms  ], 
      outputs = [ 'pixavg', 'pixstd' ] )
  if ( ( cpb_noise <= 0. ) or ( cpb_noise > 2. * rms ) ):
    if print_info:
      print '... WARNING: histogram noise fit failed, using image RMS instead'
    cpb_noise = rms
  store_parameter( uv, 'cpb_noise', cpb_noise )
  
  # clean up
  if ( not keep_image ):
    cpb_image.zap()
  
  return cpb_noise

###############################################################################

def extend_flags( uv, flag_version = 0, apply_flags = True, print_info = True,
    cutoff = 0.1, do_channels = True, do_baselines = True, flag_factor = 0.6 ):
  
  # retrieve basic info
  channel_count = get_channel_count( uv )
  time_list = get_time_list( uv )
  baseline_count = get_baseline_count( uv )
  antenna_count = get_antenna_count( uv )  
  integration_time = restore_parameter( uv, 'integration_time' ) / 86400.
  if ( flag_version == -1 ):
    flagver = uv.table_highver( 'FG' ) + 1
  elif ( flag_version == 0 ):
    flagver = max( 1, uv.table_highver( 'FG' ) )
  else:
    flagver = flag_version
  if ( apply_flags and table_exists( uv, 'FG', 0 ) ):
    flag_uv = apply_flag_table( uv, keep_solutions = True )
  else:
    flag_uv = uv
  
  # flag times
  min_times = int( ceil( cutoff * float( len( time_list ) ) ) )
  time_counts = zeros( ( baseline_count ), dtype = int64 )
  min_channels = int( ceil( cutoff * float( channel_count ) ) )
  flag_list = []
  for group in wizardry( flag_uv ):
    bl = group.baseline
    sel = awhere( group.visibility[ 0, : , 0, 2 ] > 0. )
    if ( ( len( sel ) == 0 ) or ( do_channels and ( len( sel ) < min_channels ) ) ):
      time1 = time_to_dhms( group.time - flag_factor * integration_time )
      time2 = time_to_dhms( group.time + flag_factor * integration_time )
      flag_list.append( { 'antennas' : bl[ 0 : 1 ], 'baseline' : bl[ 1 : 2 ],
            'timerang' : time1 + time2 } )
    else:
      index = baseline_to_index( flag_uv, bl, antenna_count = antenna_count )
      time_counts[ index ] = time_counts[ index ] + 1
  if ( flag_uv != uv ):
    flag_uv.zap()
  
  # flag baselines
  if do_baselines:
    sel = awhere( time_counts < min_times )
  else:
    sel = awhere( time_counts == 0 )
  for s in sel:
    flag_bl = index_to_baseline( uv, s[ 0 ], antenna_count = antenna_count )
    new_flag_list = []
    for flag in flag_list:
      bl = [ flag[ 'antennas' ][ 0 ], flag[ 'baseline' ][ 0 ] ]
      if ( bl != flag_bl ):
        new_flag_list.append( flag )
    new_flag_list.append( { 'antennas' : flag_bl[ 0 : 1 ],
        'baseline' : flag_bl[ 1 : 2 ] } )
    flag_list = [ x for x in new_flag_list ]
  
  if print_info:
    print '... (re-)flagging %d baselines' % ( len( sel ) )
    print '... flagging %d visibility groups' % ( len( flag_list ) - len( sel ) )
  add_flags( uv, flags = flag_list, flag_version = flagver, 
      print_info = print_info )
  
  return

###############################################################################

def flag_image_undulations( uv, image, imagr_params, cutoff = 3., rotate = False,
    print_info = True, uv_range = [ 0., 1.e6 ], clip_sigma = 5., bin_power = 1., 
    keep_images = False, min_data = -1, flag_version = 0, fft_trunc = 0.7,
    power = 0., max_per_bin = 0.05, apply_flags = True, do_phase = False,
    clip_negative = False, count_power = 0.5, flag_factor = 0.6,
    apply_smoothing = False ):
  
  # check UV range limits in imagr_params
  [ uv_min, uv_max ] = uv_range
  try:
    uvrang = imagr_params[ 'uvrang' ]
    uv_min = max( uv_min, uvrang[ 0 ] * 1.e3 )
    uv_max = min( uv_max, uvrang[ 1 ] * 1.e3 )
  except KeyError:
    pass
  try:
    uvrange = imagr_params[ 'uvrange' ]
    uv_min = max( uv_min, uvrange[ 0 ] * 1.e3 )
    uv_max = min( uv_max, uvrange[ 1 ] * 1.e3 )
  except KeyError:
    pass
  try:
    uvtaper = imagr_params[ 'uvtaper' ]
    uv_max = min( uv_max, min( uvtaper ) * 1.e3 )
  except KeyError:
    pass
  min_uv2 = uv_min**2
  max_uv2 = uv_max**2
  if print_info:
    print '... flagging UV range %d to %d lambda' % ( int( uv_min ), int( uv_max ) )
  
  # cut image to nearest power of two
  if print_info:
    print '... Fourier-transforming image'
  image_size = get_image_size( image )
  new_image_size = [ int( 2.**ceil( log( image_size[ 0 ] ) / log( 2. ) ) ), 
      int( 2.**ceil( log( image_size[ 1 ] ) / log( 2. ) ) ) ]
  imsize = [ min( new_image_size[ 0 ], 8192 ), min( new_image_size[ 1 ], 8192 ) ]
  padded_image = get_aips_file( image.disk, image.name, 'PADDED', -1, 'MA' )
  call_aips_task( 'OGEOM', indata = image, outdata = padded_image, imsize = imsize,
      aparm = [ 0, 0, 0, 0, 0, 0, 0 ] )
  
  # remove tables to speed up things
  cc_max = padded_image.table_highver( 'CC' )
  for i in range( 1, 1 + cc_max ):
    if table_exists( padded_image, 'CC', i ):
      padded_image.zap_table( 'CC', i )
  sn_max = padded_image.table_highver( 'SN' )
  for i in range( 1, 1 + sn_max ):
    if table_exists( padded_image, 'SN', i ):
      padded_image.zap_table( 'SN', i )
  
  # clip flux above N sigma
  # replace image blanks by zeros
  pix = get_image_pixels( padded_image )
  new_pix = abs( pix )
  sel = awhere( ( pix != get_aips_magic_value() ) & ( new_pix > 0. ) )
  sel_pix = aget( new_pix, sel )
  sel_rms = 1.48 * median( sel_pix )
  while ( sel_pix.max() > clip_sigma * sel_rms ):
    sel2 = awhere( sel_pix < clip_sigma * sel_rms )
    sel = aget( sel, sel2 )
    sel_pix = aget( new_pix, sel )
    sel_rms = 1.48 * median( sel_pix )
  if ( not clip_negative ):
    sel = awhere( pix < clip_sigma * sel_rms )
  new_pix[ : , : ] = get_aips_magic_value()
  new_pix = aput( new_pix, sel, aget( pix, sel ) )
  set_image_pixels( padded_image, new_pix )
  del pix
  del new_pix
  
  # FFT image
  # combine FFT results into amplitude image
  fft_input = get_aips_file( image.disk, image.name, 'REMAG', -1, 'MA' )
  call_aips_task( 'REMAG', indata = padded_image, outdata = fft_input )
  padded_image.zap()
  if rotate:
    pix = get_image_pixels( fft_input, flip = False )
    set_image_pixels( fft_input, pix )
    del pix
  fft_real = get_aips_file( image.disk, image.name, 'UVREAL', -1, 'MA' )
  fft_imag = get_aips_file( image.disk, image.name, 'UVIMAG', -1, 'MA' )
  outseq = max( fft_real.seq, fft_imag.seq )
  call_aips_task( 'FFT', indata = fft_input, opcode = 'MARE',
      outdisk = fft_real.disk, outname = fft_real.name, outseq = outseq )
  fft_real = get_aips_file( image.disk, image.name, 'UVREAL', outseq, 'MA' )
  fft_imag = get_aips_file( image.disk, image.name, 'UVIMAG', outseq, 'MA' )
  fft_amp = get_aips_file( image.disk, image.name, 'FFTA', -1, 'MA' )
  if do_phase:
    call_aips_task( 'COMB', indata = fft_real, in2data = fft_imag,
        outdata = fft_amp, opcode = 'POLA', aparm = [ 1., 0., 0., 0. ] )
    pix = get_image_pixels( fft_amp, flip = rotate )
    pix = abs( amodulo( pix + pi, 2. * pi ) - pi )
  else:
    call_aips_task( 'COMB', indata = fft_real, in2data = fft_imag,
        outdata = fft_amp, opcode = 'POLI', aparm = [ 1., 1., 0., 0. ] )
    pix = get_image_pixels( fft_amp, flip = rotate ) - 1.
  fft_real.zap()
  fft_imag.zap()
  
  # get FFT image info
  for ctype in fft_amp.header.ctype:
    if ( ctype.find( 'UU' ) != - 1 ):
      u_index = fft_amp.header.ctype.index( ctype )
    if ( ctype.find( 'VV' ) != - 1 ):
      v_index = fft_amp.header.ctype.index( ctype )
  uv_ref = [ int( round( fft_amp.header.crpix[ u_index ] ) ) - 1,
      int( round( fft_amp.header.crpix[ v_index ] ) ) - 1 ]
  uv_size = [ fft_amp.header.cdelt[ u_index ], fft_amp.header.cdelt[ v_index ] ]
  if ( not keep_images ):
    fft_amp.zap()
  u_min = int( round( 0.5 * ( 1. - fft_trunc ) * pix.shape[ 0 ] ) )
  u_max = pix.shape[ 0 ] - 1 - u_min
  v_min = int( round( 0.5 * ( 1. - fft_trunc ) * pix.shape[ 1 ] ) )
  v_max = pix.shape[ 1 ] - 1 - v_min
  # TODO: implement boundary check
  
  # apply smoothing
  if apply_smoothing:
    if print_info:
      print '... smoothing Fourier transform'
    pix2 = pix.copy()
    pix2[ 1 : ][ : ] = pix2[ 1 : ][ : ] + 0.5 * pix[ : -1 ][ : ]
    pix2[ : ][ 1 : ] = pix2[ : ][ 1 : ] + 0.5 * pix[ : ][ : -1 ]
    pix2[ : -1 ][ : ] = pix2[ : -1 ][ : ] + 0.5 * pix[ 1 : ][ : ]
    pix2[ : ][ : -1 ] = pix2[ : ][ : -1 ] + 0.5 * pix[ : ][ 1 : ]
    pix2[ 1 : ][ 1 : ] = pix2[ 1 : ][ 1 : ] + 0.25 * pix[ : -1 ][ : -1 ]
    pix2[ 1 : ][ : -1 ] = pix2[ 1 : ][ : -1 ] + 0.25 * pix[ : -1 ][ 1 : ]
    pix2[ : -1 ][ 1 : ] = pix2[ : -1 ][ 1 : ] + 0.25 * pix[ 1 : ][ : -1 ]
    pix2[ : -1 ][ : -1 ] = pix2[ : -1 ][ : -1 ] + 0.25 * pix[ 1 : ][ 1 : ]
    pix2 = 0.25 * pix2
    del pix
    pix = pix2
  
  # apply UV track mask
  if print_info:
    print '... generating UV data mask'
  if ( flag_version == -1 ):
    flagver = uv.table_highver( 'FG' ) + 1
  elif ( flag_version == 0 ):
    flagver = max( 1, uv.table_highver( 'FG' ) )
  else:
    flagver = flag_version
  if ( apply_flags and table_exists( uv, 'FG', flagver ) ):
    flag_uv = apply_flag_table( uv, version = flagver )
  else:
    flag_uv = uv
  factors = get_frequency( uv ) / array( get_frequency_list( uv ), dtype = float32 )
  count = len( factors )
  counts = zeros( pix.shape, dtype = int64 )
  i = 0
  big_sel = []
  for group in wizardry( flag_uv ):
    i = i + 1
    uvw = group.uvw
    dus = array( ( uvw[ 0 ] / ( uv_size[ 0 ] * factors ) ).round(), dtype = int64 )
    dvs = array( ( uvw[ 1 ] / ( uv_size[ 1 ] * factors ) ).round(), dtype = int64 )
    sel = awhere( group.visibility[ 0, : , 0, 2 ] > 0. )
    count = len( sel )
    if ( count == 0 ):
      continue
    dus = aget( dus, sel )
    dvs = aget( dvs, sel )
    sel = zeros( ( 2 * count, 2 ), dtype = sel.dtype )
    sel[ 0 : count, 0 ] = uv_ref[ 0 ] + dus
    sel[ 0 : count, 1 ] = uv_ref[ 1 ] + dvs
    sel[ count : 2 * count, 0 ] = uv_ref[ 0 ] - dus
    sel[ count : 2 * count, 1 ] = uv_ref[ 1 ] - dvs
    big_sel = big_sel + sel.tolist()
    if ( i % 1000 == 0 ):
      big_sel = array( big_sel )
      big_counts = acount( big_sel )
      counts = aput( counts, big_sel, aget( counts, big_sel ) + big_counts )
      big_sel = []
    if ( print_info and ( i % 10000 == 0 ) ):
      print '... %7d visibility groups done' % ( i )
  if ( len( big_sel ) > 0 ):
    big_sel = array( big_sel )
    big_counts = acount( big_sel )
    counts = aput( counts, big_sel, aget( counts, big_sel ) + big_counts )
  mask = ( counts > 0 )
  counts = ( float32( counts )**count_power ) * float32( mask )
  counts = counts + ( float32( mask ) - 1. )
  pix = pix / counts
  pix = pix / pix.max()
  if keep_images:
    set_image_pixels( fft_amp, pix )
  pix = pix * float32( mask )
  if ( min_data < 0 ):
    min_data = int( around( sqrt( float( len( awhere( mask ) ) ) ) ) )
  
  # identify amplitude outliers per UV bin
  if print_info:
    print '... flagging per UV distance bin'
  flag_mask = mask.copy()
  sel = awhere( flag_mask )
  data_list = []
  for s in sel:
    u = ( s[ 0 ] - uv_ref[ 0 ] ) * uv_size[ 0 ]
    v = ( s[ 1 ] - uv_ref[ 1 ] ) * uv_size[ 1 ]
    uv2 = u**2 + v**2
    data = pix[ s[ 0 ], s[ 1 ] ]
    data_list.append( [ uv2, s, data ] )
  data_list.sort( cmp = lambda a, b: cmp( a[ 0 ], b[ 0 ] ) )
  sel = array( [ d[ 1 ] for d in data_list ], dtype = int64 )
  data = array( [ d[ 2 ] for d in data_list ], dtype = float32 )
  ndata = len( data )
  if ( ndata <= 2. * min_data ):
    nbins = 1
    dbin = ndata
  else:
#    nbins = int( floor( sqrt( ndata / min_data ) ) )
    nbins = int( ceil( sqrt( float( ndata ) ) ) )
    dbin = int( ceil( float( ndata ) / ( float( nbins )**bin_power ) ) )
    dbin = max( min_data, dbin )
    nbins = int( ceil( ( float( ndata ) / float( dbin ) )**( 1. / bin_power ) ) )
#  for i in range( nbins ):
  i = -1
  while True:
    i = i + 1
    lim1 = int( round( ( float( i )**bin_power ) * float( dbin ) ) )
    lim2 = int( round( ( float( i + 1 )**bin_power ) * float( dbin ) ) )
    if ( lim1 >= ndata ):
      break
    if ( lim2 > ndata ):
      lim2 = ndata
    bin_sel = sel[ lim1 : lim2 ]
    bin_data = aget( pix, bin_sel )
    sel2 = awhere( bin_data > 0. )
    bin_sel = aget( bin_sel, sel2 )
    start_count = len( bin_sel )
    min_data2 = int( round( ( 1. - max_per_bin ) * float( len( bin_sel ) ) ) )
    while ( True ):
      bin_data = aget( pix, bin_sel )
#      if ( len( bin_data ) <= min_data ):
#         break
      if ( len( bin_data ) <= min_data2 ):
         break
      rms = 1.48 * median( abs( bin_data ) )
      sel2 = awhere( bin_data <= cutoff * rms )
      while ( len( sel2 ) < min_data2 ):
        rms = 1.01 * rms
        sel2 = awhere( bin_data <= cutoff * rms )
      sel2 = awhere( bin_data > cutoff * rms )
      if ( len( sel2 ) == 0 ):
        break
      sel2 = aget( bin_sel, sel2 )
      flag_mask = aput( flag_mask, sel2, False )
      sel2 = awhere( bin_data <= cutoff * rms )
      bin_sel = aget( bin_sel, sel2 )
    end_count = len( bin_sel )
    if print_info:
      if ( start_count > 0 ):
        fraction = 1. - float( end_count ) / float( start_count )
      else:
        fraction = 0.
      if ( fraction > 0.4 ):
        print '...... flagging bin %d (%d points)' % ( i + 1, len( bin_sel ) )
        print '...... flagged %s percent of data' % ( repr( 100. * fraction ) )
        print '...... WARNING: THIS IS A RATHER HIGH PERCENTAGE'
        pdb.set_trace()
  if print_info:
    old_count = len( awhere( mask ) )
    new_count = len( awhere( flag_mask ) )
    fraction = 1. - float( new_count ) / float( old_count )
    print '... flagged %s percent of data' % ( repr( 100. * fraction ) )
  if keep_images:
    fft_amp_flagged = get_aips_file( image.disk, image.name, 'FFTAF', -1, 'MA' )
    call_aips_task( 'MOVE', indata = fft_amp, outdata = fft_amp_flagged,
        userid = get_aips_userid() )
    set_image_pixels( fft_amp_flagged, pix * flag_mask )
  del pix
  del mask
  
  # generate flag table with flags
  if print_info:
    print '... generating flag commands'
  integration_time = restore_parameter( uv, 'integration_time' ) / 86400.
  flag_list = []
  for group in wizardry( flag_uv ):
    uvw = group.uvw
    uv2 = uvw[ 0 ]**2 + uvw[ 1 ]**2
    if ( ( uv2 < min_uv2 ) or ( uv2 > max_uv2 ) ):
      continue
    dus = array( ( uvw[ 0 ] / ( uv_size[ 0 ] * factors ) ).round(), dtype = int64 )
    dvs = array( ( uvw[ 1 ] / ( uv_size[ 1 ] * factors ) ).round(), dtype = int64 )
    csel = awhere( group.visibility[ 0, : , 0, 2 ] > 0. )
    count = len( csel )
    if ( count == 0 ):
      continue
    dus = aget( dus, csel )
    dvs = aget( dvs, csel )
    sel = zeros( ( count, 2 ), dtype = sel.dtype )
    sel[ : , 0 ] = uv_ref[ 0 ] + dus
    sel[ : , 1 ] = uv_ref[ 1 ] + dvs
    flags = aget( flag_mask, sel )
    csel = csel.ravel()
    if ( not alltrue( flags ) ):
      time1 = time_to_dhms( group.time - flag_factor * integration_time )
      time2 = time_to_dhms( group.time + flag_factor * integration_time )
      bl = group.baseline
      bc = 0
      in_range = False
      for c in range( count ):
        if ( ( not in_range ) and ( not flags[ c ] ) ):
          bc = csel[ c ] + 1
          in_range = True
        elif ( in_range and flags[ c ] ):
          flag_list.append( { 'antennas' : bl[ 0 : 1 ], 'baseline' : bl[ 1 : 2 ],
             'bchan' : bc, 'echan' : csel[ c - 1 ] + 1,
             'timerang' : time1 + time2 } )
          in_range = False
      if in_range:
        flag_list.append( { 'antennas' : bl[ 0 : 1 ], 'baseline' : bl[ 1 : 2 ],
            'bchan' : bc, 'echan' : csel[ count - 1 ] + 1,
            'timerang' : time1 + time2 } )
        in_range = False
  del flag_mask
  if ( flag_uv != uv ):
    flag_uv.zap()
  add_flags( uv, flags = flag_list, flag_version = flagver,
      print_info = print_info )
  if print_info:
    fraction = 1. - float( new_count ) / float( old_count )
    print '... flagged %s percent of data' % ( repr( 100. * fraction ) )
  
  return ( fraction > 0. )

###############################################################################

def flag_image_baselines( uv, cutoff = 3., print_info = True, imagr_params = {},
    uv_range = [ 0., 1.e6 ], min_data = -1, apply_flags = True,
    per_baseline = True, per_uv_bin = True, check_baselines = False,
    baseline_cutoff = 4., drop_factor = 0.5, keep_images = False,
    scaling_factor = 1., apply_solutions = True, solution_version = 0,
    flag_version = 0, fft_trunc = 0.6, max_per_bin = 0.05, do_phase = False,
    count_power = 0.25, max_gap = 2.5, flag_factor = 0.6, bin_power = 1.,
    max_image_size = 4096 ):
  
  # check UV range limits in imagr_params
  [ uv_min, uv_max ] = uv_range
  try:
    uvrang = imagr_params[ 'uvrang' ]
    uv_min = max( uv_min, uvrang[ 0 ] * 1.e3 )
    uv_max = min( uv_max, uvrang[ 1 ] * 1.e3 )
  except KeyError:
    pass
  try:
    uvrange = imagr_params[ 'uvrange' ]
    uv_min = max( uv_min, uvrange[ 0 ] * 1.e3 )
    uv_max = min( uv_max, uvrange[ 1 ] * 1.e3 )
  except KeyError:
    pass
  try:
    uvtaper = imagr_params[ 'uvtaper' ]
    uv_max = min( uv_max, min( uvtaper ) * 1.5e3 )
  except KeyError:
    pass
  uv2_min = uv_min**2
  uv2_max = uv_max**2
  if print_info:
    print '... flagging UV range %d to %d lambda' % ( int( uv_min ), int( uv_max ) )
  
  # flag and collapse UV data in frequency
  if print_info:
    print '... making image'
  if ( apply_flags and table_exists( uv, 'FG', 0 ) ):
    flag_uv = apply_flag_table( uv, keep_solutions = True )
  else:
    flag_uv = uv
  if ( apply_solutions and table_exists( uv, 'SN', solution_version ) ):
    docalib = 100
  else:
    docalib = -1
  channel_count = get_channel_count( uv )
  favg_uv = get_aips_file( uv.disk, uv.name, 'FAVG', -1, 'UV' )
  call_aips_task( 'SPLAT', indata = flag_uv, flagver = -1, docalib = docalib,
      gainuse = solution_version, douvcomp = 0, aparm = [ 3, 0, 1, 0, 0, 0, 0 ],
      outdata = favg_uv, stokes = 'I', bchan = 1, echan = 0,
      channel = channel_count, chinc = channel_count )
  if ( flag_uv != uv ):
    flag_uv.zap()
  
  # make large image
  # scale image to lower nearest power of two
  cell_size = scaling_factor * restore_parameter( uv, 'cell_size' )
  image_size = restore_parameter( uv, 'pb_image_size' )
  uv_size = image_size
  image_size = int( 2.**floor( log( image_size ) / log( 2. ) ) )
  image_size = min( image_size, max_image_size )
  image = get_aips_file( uv.disk, 'FOURIER', 'ICL001', -1, 'MA' )
  try:
    call_aips_task( 'IMAGR', indata = favg_uv, nchav = 1, in2disk = uv.disk,
        outdisk = image.disk, overlap = 2, outver = 0, docalib = -1,
        gainuse = 0, outname = image.name, outseq = image.seq,
        cellsize = [ cell_size, cell_size ], do3dimag = 0, niter = 100000,
        flux = 100000., imsize = [ image_size, image_size ], dotv = 0,
        boxfile = '', cmethod = '', minpatch = 0, 
        gain = 0.1, nfield = 1, maxpixel = 0, factor = 0,
        imagrprm = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 ],
        flagver = -1, bcomp = [ 0 for j in range( 64 ) ],
        uvsize = [ uv_size, uv_size ], **imagr_params )
  except RuntimeError:
    return 0.
  beam = get_facet_beam( image )
  beam.zap()
  image.zap_table( 'CC', 0 )
  
  # FFT image
  if print_info:
    print '... Fourier-transforming image'
  fft_real = get_aips_file( image.disk, image.name, 'UVREAL', -1, 'MA' )
  fft_imag = get_aips_file( image.disk, image.name, 'UVIMAG', -1, 'MA' )
  outseq = max( fft_real.seq, fft_imag.seq )
  call_aips_task( 'FFT', indata = image, opcode = 'MARE',
      outdisk = fft_real.disk, outname = fft_real.name, outseq = outseq )
  fft_real = get_aips_file( image.disk, image.name, 'UVREAL', outseq, 'MA' )
  fft_imag = get_aips_file( image.disk, image.name, 'UVIMAG', outseq, 'MA' )
  
  # combine FFT results into amplitude and phase images
  fft_amp = get_aips_file( image.disk, image.name, 'FFTA', -1, 'MA' )
  if do_phase:
    call_aips_task( 'COMB', indata = fft_real, in2data = fft_imag,
        outdata = fft_amp, opcode = 'POLA', aparm = [ 1., 0., 0., 0. ] )
    pix = get_image_pixels( fft_amp, flip = False )
    pix = abs( amodulo( pix + pi, 2. * pi ) - pi )
  else:
    call_aips_task( 'COMB', indata = fft_real, in2data = fft_imag,
        outdata = fft_amp, opcode = 'POLI', aparm = [ 1., 1., 0., 0. ] )
    pix = get_image_pixels( fft_amp, flip = False ) - 1.
  fft_real.zap()
  fft_imag.zap()
  
  # get FFT image info
  for ctype in fft_amp.header.ctype:
    if ( ctype.find( 'UU' ) != - 1 ):
      u_index = fft_amp.header.ctype.index( ctype )
    if ( ctype.find( 'VV' ) != - 1 ):
      v_index = fft_amp.header.ctype.index( ctype )
  uv_ref = [ int( round( fft_amp.header.crpix[ u_index ] ) ) - 1,
      int( round( fft_amp.header.crpix[ v_index ] ) ) - 1 ]
  uv_size = [ fft_amp.header.cdelt[ u_index ], fft_amp.header.cdelt[ v_index ] ]
  factor = get_frequency( favg_uv ) / get_frequency( image )
  du = factor * uv_size[ 0 ]
  dv = factor * uv_size[ 1 ]
  u_min = int( round( 0.5 * ( 1. - fft_trunc ) * pix.shape[ 0 ] ) )
  u_max = pix.shape[ 0 ] - 1 - u_min
  v_min = int( round( 0.5 * ( 1. - fft_trunc ) * pix.shape[ 1 ] ) )
  v_max = pix.shape[ 1 ] - 1 - v_min
  if ( not keep_images ):
    image.zap()
  
  # create UV track mask
  if print_info:
    print '... generating UV baseline masks'
  antenna_count = get_antenna_count( uv )
  group_count = 0
  bl_data = []
  for group in wizardry( favg_uv ):
    group_count = group_count + 1
    uvw = group.uvw
    uv2 = uvw[ 0 ]**2 + uvw[ 1 ]**2
    if ( ( uv2 < uv2_min ) or ( uv2 > uv2_max ) ):
      continue
    time = group.time
    baseline = group.baseline
    index = baseline_to_index( uv, baseline, antenna_count = antenna_count )
    du = int( round( uvw[ 0 ] / ( uv_size[ 0 ] * factor ) ) )
    dv = int( round( uvw[ 1 ] / ( uv_size[ 1 ] * factor ) ) )
    u1 = uv_ref[ 0 ] + du
    u2 = uv_ref[ 0 ] - du
    v1 = uv_ref[ 1 ] + dv
    v2 = uv_ref[ 1 ] - dv
    if ( ( u1 < u_min ) or ( u1 > u_max ) ):
      continue
    if ( ( v1 < v_min ) or ( v1 > v_max ) ):
      continue
    if ( ( u2 < u_min ) or ( u2 > u_max ) ):
      continue
    if ( ( v2 < v_min ) or ( v2 > v_max ) ):
      continue
    bl_data.append( [ time, index, u1,u2,v1,v2 ] )
  favg_uv.zap()
  bl_data.sort( cmp = lambda a, b: cmp( a[ 0 ], b[ 0 ] ) )
  baseline_pixels = []
  baseline_count = get_baseline_count( uv )
  counts = zeros( pix.shape, dtype = int64 )
  for i in range( baseline_count ):
    mask = zeros( pix.shape, dtype = bool )
    for [ time, index, u1,u2,v1,v2 ] in bl_data:
      if ( index == i ):
        mask[ u1, v1 ] = True
        mask[ u2, v2 ] = True
        counts[ u1, v1 ] = counts[ u1, v1 ] + 1
        counts[ u2, v2 ] = counts[ u2, v2 ] + 1
    baseline_pixels.append( awhere( mask * ( pix > 0. ) ) )
  mask = ( counts > 0 )
  if ( len( awhere( mask ) ) == 0 ):
    return 0.
  counts = ( float32( counts )**count_power ) * float32( mask )
  counts = counts + ( float32( mask ) - 1. )
  pix = pix / counts
  del counts
  pix = pix / pix.max()
  if keep_images:
    set_image_pixels( fft_amp, pix, flip = True )
  else:
    fft_amp.zap()
  mask = mask * ( pix > 0. )
  pix = pix * float32( mask )
  orig_count = len( awhere( mask ) )
  flag_mask = mask.copy()
  if ( orig_count == 0 ):
    return 0.
  if ( min_data < 0 ):
    min_data = int( around( sqrt( float( len( awhere( mask ) ) ) ) ) )
  
  # identify amplitude outliers per baseline
  old_count = len( awhere( flag_mask ) )
  if ( per_baseline and ( old_count > 0 ) ):
    if print_info:
      print '... flagging per baseline'
    for i in range( baseline_count ):
      # robust flagging
      while ( True ):
        sel = awhere( aget( flag_mask, baseline_pixels[ i ] ) )
        if ( len( sel ) < 10 ):
          break
        sel = aget( baseline_pixels[ i ], sel )
        sel_pix = aget( pix, sel )
        rms = 1.48 * median( abs( sel_pix ) )
        sel2 = awhere( sel_pix > cutoff * rms )
        if ( len( sel2 ) == 0 ):
          break
        sel2 = aget( sel, sel2 )
        flag_mask = aput( flag_mask, sel2, False )
    if print_info:
      new_count = len( awhere( flag_mask ) )
      fraction = 1. - float( new_count ) / float( old_count )
      print '...... flagged %s percent of data' % ( repr( 100. * fraction ) )
  
  # identify amplitude outliers per UV bin
  old_count = len( awhere( flag_mask ) )
  if ( per_uv_bin and ( old_count > 0 ) ):
    if print_info:
      print '... flagging per UV distance bin'
    sel = awhere( flag_mask )
    data_list = []
    for s in sel:
      u = factor * ( s[ 0 ] - uv_ref[ 0 ] ) * uv_size[ 0 ]
      v = factor * ( s[ 1 ] - uv_ref[ 1 ] ) * uv_size[ 1 ]
      uv2 = u**2 + v**2
      data = pix[ s[ 0 ], s[ 1 ] ]
      data_list.append( [ uv2, s, data ] )
    data_list.sort( cmp = lambda a, b: cmp( a[ 0 ], b[ 0 ] ) )
    sel = array(  [ d[ 1 ] for d in data_list ], dtype = int64 )
    data = array( [ d[ 2 ] for d in data_list ], dtype = float32 )
    ndata = len( sel )
    if ( ndata <= 2. * min_data ):
      nbins = 1
      dbin = ndata
    else:
#      nbins = int( floor( sqrt( ndata / min_data ) ) )
      nbins = int( ceil( sqrt( float( ndata ) ) ) )
      dbin = int( ceil( float( ndata ) / ( float( nbins )**bin_power ) ) )
      dbin = max( min_data, dbin )
      nbins = int( ceil( ( float( ndata ) / float( dbin ) )**( 1. / bin_power ) ) )
    i = -1
    while True:
      i = i + 1
      lim1 = int( round( ( float( i )**bin_power ) * float( dbin ) ) )
      lim2 = int( round( ( float( i + 1 )**bin_power ) * float( dbin ) ) )
      if ( lim1 >= ndata ):
        break
      if ( lim2 > ndata ):
        lim2 = ndata
      bin_sel = sel[ lim1 : lim2 ]
      min_data2 = int( round( ( 1. - max_per_bin ) * float( len( bin_sel ) ) ) )
      while ( True ):
        bin_data = aget( pix, bin_sel )
#        if ( len( bin_data ) <= min_data ):
#          break
        if ( len( bin_data ) <= min_data2 ):
          break
        rms = 1.48 * median( abs( bin_data ) )
        sel2 = awhere( bin_data <= cutoff * rms )
        while ( len( sel2 ) < min_data2 ):
          rms = 1.01 * rms
          sel2 = awhere( bin_data <= cutoff * rms )
        sel2 = awhere( bin_data > cutoff * rms )
        if ( len( sel2 ) == 0 ):
          break
        sel2 = aget( bin_sel, sel2 )
        flag_mask = aput( flag_mask, sel2, False )
        sel2 = awhere( bin_data <= cutoff * rms )
        bin_sel = aget( bin_sel, sel2 )
    if print_info:
      new_count = len( awhere( flag_mask ) )
      fraction = 1. - float( new_count ) / float( old_count )
      print '...... flagged %s percent of data' % ( repr( 100. * fraction ) )
  sel = awhere( flag_mask != mask )
  pix = aput( pix, sel, -1. )
  
  # identify baselines with few points left or persistent high rms
  bl_flag_list = []
  if check_baselines:
    if print_info:
      print '... flagging excessive baselines'
      old_count = len( awhere( flag_mask ) )
    bl_rms_list = []
    for i in range( baseline_count ):
      sel = awhere( aget( pix * float32( flag_mask ), baseline_pixels[ i ] ) > 0. )
      sel = aget( baseline_pixels[ i ], sel )
      sel_pix = aget( pix, sel )
      if ( len( sel_pix ) > 0 ):
        flag_factor = float( len( sel_pix ) ) / float( len( baseline_pixels[ i ] ) )
        if ( flag_factor <= drop_factor ):
          bl_flag_list.append( i )
        continue
      rms = 1.48 * median( abs( sel_pix ) )
      bl_rms_list.append( [ i, rms ] )
    indices = array( [ b[ 0 ] for b in bl_rms_list ], dtype = int64 )
    rmss = array( [ b[ 1 ] for b in bl_rms_list ], dtype = float32 )
    sel = awhere( rmss > baseline_cutoff * median( rmss ) )
    while ( len( sel ) > 0 ):
      sel_indices = aget( indices, sel )
      bl_flag_list = bl_flag_list + sel_indices.tolist()
      sel = awhere( rmss <= baseline_cutoff * median( rmss ) )
      indices = aget( indices, sel )
      rmss = aget( rmss, sel )
      sel = awhere( rmss > baseline_cutoff * median( rmss ) )
    if ( len( bl_flag_list ) > 0 ):
      new_mask = zeros( pix.shape, dtype = bool )
      for i in range( baseline_count ):
        if ( not i in bl_flag_list ):
          new_mask = aput( new_mask, baseline_pixels[ i ], True )
      sel = awhere( new_mask != mask )
      pix = aput( pix, sel, -1. )
      del mask
      mask = new_mask
      flag_mask = flag_mask * mask
    if print_info:
      new_count = len( awhere( flag_mask ) )
      fraction = 1. - float( new_count ) / float( old_count )
      print '...... flagged %s percent of data' % ( repr( 100. * fraction ) )
  if print_info:
    new_count = len( awhere( flag_mask ) )
  
  if keep_images:
    fft_amp_flagged = get_aips_file( image.disk, image.name, 'FFTAF', -1, 'MA' )
    call_aips_task( 'MOVE', indata = fft_amp, outdata = fft_amp_flagged,
        userid = get_aips_userid() )
    set_image_pixels( fft_amp_flagged, pix, flip = True )
  del pix
  
  # generate flag table
  if print_info:
    print '... generating flag table'
  if ( flag_version == -1 ):
    flagver = uv.table_highver( 'FG' ) + 1
  elif ( flag_version == 0 ):
    flagver = max( 1, uv.table_highver( 'FG' ) )
  else:
    flagver = flag_version
  flag_list = []
  # first flag baselines
  if ( len( bl_flag_list ) > 0 ):
    baseline_flags = []
    for index in bl_flag_list:
      [ i, j ] = index_to_baseline( uv, index )
      flag_list.append( { 'antennas' : [ i ], 'baseline' : [ j ] } )
      flag_mask = aput( flag_mask, baseline_pixels[ index ], True )
      mask = aput( mask, baseline_pixels[ index ], True )
  # then flag outlier UV pixels
  flag_mask = ( flag_mask != mask )
  del mask
  
  if any( flag_mask ):
    integration_time = restore_parameter( uv, 'integration_time' ) / 86400.
    for i in range( baseline_count ):
      start_time = -1.
      last_time = -1.
      bl = index_to_baseline( uv, i, antenna_count = antenna_count )
      for [ time, index, u1,u2,v1,v2 ] in bl_data:
        if ( index == i ):
          if ( ( not flag_mask[ u1, v1 ] ) and ( not flag_mask[ u2, v2 ] ) ):
            if ( start_time > -1. ):
              # finish previous flagging
              time1 = time_to_dhms( start_time - flag_factor * integration_time )
              time2 = time_to_dhms( last_time + flag_factor * integration_time )
              flag_list.append( { 'antennas' : bl[ 0 : 1 ], 'baseline' : bl[ 1 : 2 ],
                  'timerang' : time1 + time2 } )
              start_time = -1.
          elif ( abs( time - last_time ) > max_gap * integration_time ):
            if ( start_time > -1. ):
              # finish previous flagging
              time1 = time_to_dhms( start_time - flag_factor * integration_time )
              time2 = time_to_dhms( last_time + flag_factor * integration_time )
              flag_list.append( { 'antennas' : bl[ 0 : 1 ], 'baseline' : bl[ 1 : 2 ],
                  'timerang' : time1 + time2 } )
            start_time = time
          else:
            # start new flagging
            if ( start_time == -1. ):
              start_time = time
          last_time = time
      if ( start_time > -1. ): # finish last flagging
        time1 = time_to_dhms( start_time - flag_factor * integration_time )
        time2 = time_to_dhms( last_time + flag_factor * integration_time )
        flag_list.append( { 'antennas' : bl[ 0 : 1 ], 'baseline' : bl[ 1 : 2 ],
            'timerang' : time1 + time2 } )
  del flag_mask
  add_flags( uv, flags = flag_list, flag_version = flagver, 
      print_info = print_info )
  if print_info:
    fraction = 1. - float( new_count ) / float( orig_count )
    print '... flagged total %s percent of data' % ( repr( 100. * fraction ) )
  
  return ( fraction > 0. )

###############################################################################

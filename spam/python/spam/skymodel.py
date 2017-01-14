###############################################################################

# import Python modules
from os import *
from math import *

# import user modules
from files import *
from aips import *
from sphere import *
from parameter import *
from mpfit import *
from error import *

###############################################################################

def generate_source_list( radec, radius, freq, epoch = 2000.0, use_wenss = True,
    use_nvss = True, assoc_radius = 40. / 3600., spectral_index = None,
    flux_min = None, si_limits = [ -1.5, 0. ], use_vlss = False,
    catalog_path = '${AIPS_ROOT}/TEXT/STARS/' ):
# radius <= 90 degrees
# assoc_radius = 40. / 3600. # e.g. 2004MNRAS.352..909C
  
  catalog_path_e = path.expandvars( catalog_path )
  if ( epoch == 2000.0 ):
    if ( radec[ 1 ] < -35. ):
#      nvss_file_name = catalog_path_e + 'SU00.0006'
      nvss_file_name = catalog_path_e + 'SM00.0006'
    else:
      nvss_file_name = catalog_path_e + 'NV00.0003'
    if ( not file_exists( nvss_file_name ) ):
      if ( radec[ 1 ] < -35. ):
#        nvss_file_name = catalog_path_e + 'SU00.0050'
        nvss_file_name = catalog_path_e + 'SM00.0050'
      else:
        nvss_file_name = catalog_path_e + 'NV00.0030'
    wenss_file_name = catalog_path_e + 'WE00.0000'
    if ( not file_exists( wenss_file_name ) ):
      wenss_file_name = catalog_path_e + 'WE00.0100'
    vlss_file_name = catalog_path_e + 'VL00.0400'
    if ( not file_exists( vlss_file_name ) ):
      vlss_file_name = None
  elif ( epoch == 1950.0 ):
    nvss_file_name = catalog_path_e + 'NV50.0003'
    if ( not file_exists( nvss_file_name ) ):
      nvss_file_name = catalog_path_e + 'NV50.0030'
    wenss_file_name = catalog_path_e + 'WE50.0000'
    if ( not file_exists( wenss_file_name ) ):
      wenss_file_name = catalog_path_e + 'WE50.0100'
    vlss_file_name = None
  if ( ( not file_exists( nvss_file_name ) ) or ( not file_exists( wenss_file_name ) ) ):
    raise error( 'no catalog files available' )
  
  if ( not spectral_index is None ):
    if ( ( spectral_index < si_limits[ 0 ] ) or ( spectral_index > si_limits[ 1 ] ) ):
      raise error( 'default spectral index outside allowed spectral index range' )

  nvss_freq = 1.4e9
  if ( radec[ 1 ] < -35. ):
    nvss_freq = 843.e6
  wenss_freq = 326.e6
  if ( radec[ 1 ] < 0. ):
    wenss_freq = 352.e6
  vlss_freq = 73.8e6

  if ( not flux_min is None ):
    if ( nvss_freq > freq ):
      nvss_flux_min = flux_min * ( nvss_freq / freq )**si_limits[ 0 ]
    else:
      nvss_flux_min = flux_min * ( nvss_freq / freq )**si_limits[ 1 ]
    if ( wenss_freq > freq ):
      wenss_flux_min = flux_min * ( wenss_freq / freq )**si_limits[ 0 ]
    else:
      wenss_flux_min = flux_min * ( wenss_freq / freq )**si_limits[ 1 ]
  else:
    nvss_flux_min = 0.
    wenss_flux_min = 0.

  # determine crude selection criteria
  ra = radec[ 0 ]
  dec = radec[ 1 ]
  ra_range_contains_zero = False
  if ( dec + radius >= 90. ) or ( dec - radius <= -90. ):
    ra_min = 0.
    ra_max = 360.
    if ( dec + radius >= 90. ):
      dec_min = min( [ dec - radius, 180. - ( dec + radius ) ] )
      dec_max = 90.
    else:
      dec_min = - 90.
      dec_max = max( [ dec + radius, - 180. - ( dec - radius ) ] )
  else:
    cos_dec = cos( radians( max( [ abs( dec + radius ), abs( dec - radius ) ] ) ) )
    if ( ( radius / cos_dec ) >= 180. ):
      ra_min = 0.
      ra_max = 360.
    else:
      ra_min = amodulo( ra - ( radius / cos_dec ), 360. )
      ra_max = amodulo( ra + ( radius / cos_dec ), 360. )
      if ( ra_min > ra_max ):
        ra_range_contains_zero = True
    dec_min = dec - radius
    dec_max = dec + radius

  # get NVSS sources
  nvss_source_list = []
  if use_nvss:
    nvss_file = file( nvss_file_name, mode = 'r' )
    for line in nvss_file:
      first_index = 0
      while ( line[ first_index ] == ' ' ):
        first_index = first_index + 1
      line = line[ first_index : ]
      if ( len( line[ first_index : ] ) == 0 ):
        continue
      if ( line[ first_index ] == ';' ):
        continue
      if ( line[ first_index : first_index + 5 ] == 'Found' ):
        continue
      try:
        words = line.split()
        nvss_source = [ [ float( words[ 0 ] ), float( words[ 1 ] ) ],
            float( int( words[ 2 ] ) ) / 1.e3 ]
      except:
        try:
          # F9.5,1X,F9.5,I7,F10.4
          nvss_source = [ [ float( line[ 0 : 9 ] ), float( line[ 10 : 19 ] ) ],
              float( int( line[ 19 : 26 ] ) ) / 1.e3 ]
        except:
          continue
      if ( ( nvss_source[ 0 ][ 1 ] >= dec_min ) and 
          ( nvss_source[ 0 ][ 1 ] <= dec_max ) and
          ( nvss_source[ 1 ] >= nvss_flux_min ) ):
        if ra_range_contains_zero:
          if ( ( nvss_source[ 0 ][ 0 ] >= ra_max ) or
              ( nvss_source[ 0 ][ 0 ] <= ra_min ) ):
            [ nvss_radius, nvss_angle ] = calculate_angular_separation( 
                nvss_source[ 0 ], [ ra, dec ] )
            if ( nvss_radius <= radius ):
              nvss_source_list.append( nvss_source )
        else:
          if ( ( nvss_source[ 0 ][ 0 ] >= ra_min ) and
              ( nvss_source[ 0 ][ 0 ] <= ra_max ) ):
            [ nvss_radius, nvss_angle ] = calculate_angular_separation(
                nvss_source[ 0 ], [ ra, dec ] )
            if ( nvss_radius <= radius ):
              nvss_source_list.append( nvss_source )
    nvss_file.close()
    nvss_source_list.sort( cmp = lambda a, b: cmp( b[ 1 ], a[ 1 ] ) )

  # get WENSS sources
  wenss_source_list = []
  if use_wenss:
    wenss_file = file( wenss_file_name, mode = 'r' )
    for line in wenss_file:
#      columns = [ column.strip() for column in line.split() ]
#      if ( len( columns ) == 3 ):
#        if ( columns[ 0 ][ 0 ] != ';' ):
#          wenss_source_data = [ float( column ) for column in columns ]
#          wenss_source = [ [ wenss_source_data[ 0 ], wenss_source_data[ 1 ] ],
#              wenss_source_data[ 2 ] / 1.e3 ]
      first_index = 0
      while ( line[ first_index ] == ' ' ):
        first_index = first_index + 1
      line = line[ first_index : ]
      if ( len( line[ first_index : ] ) == 0 ):
        continue
      if ( line[ first_index ] == ';' ):
        continue
      if ( line[ first_index : first_index + 5 ] == 'Found' ):
        continue
      try:
        words = line.split()
        wenss_source = [ [ float( words[ 0 ] ), float( words[ 1 ] ) ],
            float( int( words[ 2 ] ) ) / 1.e3 ]
      except:
        try:
          # F9.5,1X,F9.5,I7,F10.4
          wenss_source = [ [ float( line[ 0 : 9 ] ), float( line[ 10 : 19 ] ) ],
              float( int( line[ 19 : 26 ] ) ) / 1.e3 ]
        except:
          continue
      if ( ( wenss_source[ 0 ][ 1 ] >= dec_min ) and
          ( wenss_source[ 0 ][ 1 ] <= dec_max ) and
          ( wenss_source[ 1 ] >= wenss_flux_min ) ):
        if ra_range_contains_zero:
          if ( ( wenss_source[ 0 ][ 0 ] >= ra_max ) or
              ( wenss_source[ 0 ][ 0 ] <= ra_min ) ):
            [ wenss_radius, wenss_angle ] = calculate_angular_separation(
                wenss_source[ 0 ], [ ra, dec ] )
            if ( wenss_radius <= radius ):
              wenss_source_list.append( wenss_source )
        else:
          if ( ( wenss_source[ 0 ][ 0 ] >= ra_min ) and
              ( wenss_source[ 0 ][ 0 ] <= ra_max ) ):
            [ wenss_radius, wenss_angle ] = calculate_angular_separation(
                wenss_source[ 0 ], [ ra, dec ] )
            if ( wenss_radius <= radius ):
              wenss_source_list.append( wenss_source )
    wenss_file.close()
    wenss_source_list.sort( cmp = lambda a, b: cmp( b[ 1 ], a[ 1 ] ) )

  # get VLSS sources
  vlss_source_list = []
  if ( use_vlss and ( not vlss_file_name is None ) ):
    vlss_file = file( vlss_file_name, mode = 'r' )
    for line in vlss_file:
#      columns = [ column.strip() for column in line.split() ]
#      if ( len( columns ) == 3 ):
#        if ( columns[ 0 ][ 0 ] != ';' ):
#          vlss_source_data = [ float( column ) for column in columns ]
#          vlss_source = [ [ vlss_source_data[ 0 ], vlss_source_data[ 1 ] ],
#              vlss_source_data[ 2 ] / 1.e3 ]
      first_index = 0
      while ( line[ first_index ] == ' ' ):
        first_index = first_index + 1
      line = line[ first_index : ]
      if ( len( line[ first_index : ] ) == 0 ):
        continue
      if ( line[ first_index ] == ';' ):
        continue
      if ( line[ first_index : first_index + 5 ] == 'Found' ):
        continue
      try:
        words = line.split()
        vlss_source = [ [ float( words[ 0 ] ), float( words[ 1 ] ) ],
            float( int( words[ 2 ] ) ) / 1.e3 ]
      except:
        try:
          # F9.5,1X,F9.5,I7,F10.4
          vlss_source = [ [ float( line[ 0 : 9 ] ), float( line[ 10 : 19 ] ) ],
              float( int( line[ 19 : 26 ] ) ) / 1.e3 ]
        except:
          continue
      if ( ( vlss_source[ 0 ][ 1 ] >= dec_min ) and
          ( vlss_source[ 0 ][ 1 ] <= dec_max ) ):
        if ra_range_contains_zero:
          if ( ( vlss_source[ 0 ][ 0 ] >= ra_max ) or
              ( vlss_source[ 0 ][ 0 ] <= ra_min ) ):
            [ vlss_radius, vlss_angle ] = calculate_angular_separation(
                vlss_source[ 0 ], [ ra, dec ] )
            if ( vlss_radius <= radius ):
              vlss_source_list.append( vlss_source )
        else:
          if ( ( vlss_source[ 0 ][ 0 ] >= ra_min ) and
              ( vlss_source[ 0 ][ 0 ] <= ra_max ) ):
            [ vlss_radius, vlss_angle ] = calculate_angular_separation(
                vlss_source[ 0 ], [ ra, dec ] )
            if ( vlss_radius <= radius ):
              vlss_source_list.append( vlss_source )
    vlss_file.close()
    vlss_source_list.sort( cmp = lambda a, b: cmp( b[ 1 ], a[ 1 ] ) )

  # look for WENSS sources in NVSS catalog
  # extrapolate source fluxes
  # check for additional VLSS source
  # build source list
  source_list = []
  i = 0
  for wenss_source in wenss_source_list:
    i = i + 1
    cross_source_list = []
    ra_min = wenss_source[ 0 ][ 0 ] - ( assoc_radius / cos( radians( wenss_source[ 1 ] ) ) )
    ra_max = wenss_source[ 0 ][ 0 ] + ( assoc_radius / cos( radians( wenss_source[ 1 ] ) ) )
    dec_min = wenss_source[ 0 ][ 1 ] - assoc_radius
    dec_max = wenss_source[ 0 ][ 1 ] + assoc_radius
    for nvss_source in nvss_source_list:
      if ( ( nvss_source[ 0 ][ 0 ] >= ra_min ) and ( nvss_source[ 0 ][ 0 ] <= ra_max ) and
           ( nvss_source[ 0 ][ 1 ] >= dec_min ) and ( nvss_source[ 0 ][ 1 ] <= dec_max ) ):
        [ cross_radius, cross_angle ] = calculate_angular_separation(
            nvss_source[ 0 ], wenss_source[ 0 ] )
        if ( cross_radius <= assoc_radius ):
          cross_source_list.append( nvss_source )
          nvss_source_list.remove( nvss_source )
    if ( len( cross_source_list ) > 0 ):
      nvss_flux = float( sum( [ cross_source[ 1 ] for cross_source in cross_source_list ] ) )
      source_sp_index = log( nvss_flux / wenss_source[ 1 ] ) / log( nvss_freq / wenss_freq )
      if ( source_sp_index < si_limits[ 0 ] ):
        source_sp_index = si_limits[ 0 ]
      elif ( source_sp_index > si_limits[ 1 ] ):
        source_sp_index = si_limits[ 1 ]
      # check for VLSS counterpart
      if ( len( vlss_source_list ) > 0 ):
        for vlss_source in vlss_source_list:
          if ( ( vlss_source[ 0 ][ 0 ] >= ra_min ) and ( vlss_source[ 0 ][ 0 ] <= ra_max ) and
             ( vlss_source[ 0 ][ 1 ] >= dec_min ) and ( vlss_source[ 0 ][ 1 ] <= dec_max ) ):
            [ cross_radius, cross_angle ] = calculate_angular_separation(
                vlss_source[ 0 ], wenss_source[ 0 ] )
            if ( cross_radius <= assoc_radius ):
              cross_source_list.append( vlss_source )
              vlss_source_list.remove( vlss_source )
        if ( len( cross_source_list ) > 0 ):
          vlss_flux = float( sum( [ cross_source[ 1 ] for cross_source in cross_source_list ] ) )
          vlss_model_flux = nvss_flux * pow( vlss_freq / nvss_freq, source_sp_index )
          if ( vlss_flux < vlss_model_flux ):
            coef = 1. / ( ( 1. / vlss_flux ) - ( 1. / vlss_model_flux ) )
            source_flux_1 = nvss_flux * pow( freq / nvss_freq, source_sp_index )
            source_flux_2 = coef * pow( freq / vlss_freq, 2.5 )
            source_flux = 1. / ( ( 1. / source_flux_1 ) + ( 1. / source_flux_2 ) )
          else:
            source_flux = nvss_flux * pow( freq / nvss_freq, source_sp_index )
          source_list.append( [ wenss_source[ 0 ], source_flux ] )
        else: # no VLSS counterfound found
          source_flux = nvss_flux * pow( freq / nvss_freq, source_sp_index )
          source_list.append( [ wenss_source[ 0 ], source_flux ] )
      else: # no VLSS counterfound found
        source_flux = nvss_flux * pow( freq / nvss_freq, source_sp_index )
        source_list.append( [ wenss_source[ 0 ], source_flux ] )
    elif ( not spectral_index is None ): # no NVSS counterpart found
      if ( len( vlss_source_list ) > 0 ):
        for vlss_source in vlss_source_list:
          if ( ( vlss_source[ 0 ][ 0 ] >= ra_min ) and ( vlss_source[ 0 ][ 0 ] <= ra_max ) and
             ( vlss_source[ 0 ][ 1 ] >= dec_min ) and ( vlss_source[ 0 ][ 1 ] <= dec_max ) ):
            [ cross_radius, cross_angle ] = calculate_angular_separation(
                vlss_source[ 0 ], wenss_source[ 0 ] )
            if ( cross_radius <= assoc_radius ):
              cross_source_list.append( vlss_source )
              vlss_source_list.remove( vlss_source )
        if ( len( cross_source_list ) > 0 ):
          vlss_flux = float( sum( [ cross_source[ 1 ] for cross_source in cross_source_list ] ) )
          source_sp_index = log( vlss_flux / wenss_source[ 1 ] ) / log( vlss_freq / wenss_freq )
          if ( source_sp_index < si_limits[ 0 ] ):
            source_sp_index = si_limits[ 0 ]
          elif ( spectral_index > si_limits[ 1 ] ):
            source_sp_index = si_limits[ 1 ]
          source_flux = vlss_flux * pow( freq / vlss_freq, source_sp_index )
          source_list.append( [ wenss_source[ 0 ], source_flux ] )
        else: # no NVSS & VLSS counterpart found
          source_flux = wenss_source[ 1 ] * pow( freq / wenss_freq, spectral_index )
          source_list.append( [ wenss_source[ 0 ], source_flux ] )
      else: # no NVSS & VLSS counterpart found
        source_flux = wenss_source[ 1 ] * pow( freq / wenss_freq, spectral_index )
        source_list.append( [ wenss_source[ 0 ], source_flux ] )
  for vlss_source in vlss_source_list: # no WENSS counterpart found
    i = i + 1
    cross_source_list = []
    ra_min = vlss_source[ 0 ][ 0 ] - ( assoc_radius / cos( radians( vlss_source[ 1 ] ) ) )
    ra_max = vlss_source[ 0 ][ 0 ] + ( assoc_radius / cos( radians( vlss_source[ 1 ] ) ) )
    dec_min = vlss_source[ 0 ][ 1 ] - assoc_radius
    dec_max = vlss_source[ 0 ][ 1 ] + assoc_radius
    for nvss_source in nvss_source_list:
      if ( ( nvss_source[ 0 ][ 0 ] >= ra_min ) and ( nvss_source[ 0 ][ 0 ] <= ra_max ) and
           ( nvss_source[ 0 ][ 1 ] >= dec_min ) and ( nvss_source[ 0 ][ 1 ] <= dec_max ) ):
        [ cross_radius, cross_angle ] = calculate_angular_separation(
            nvss_source[ 0 ], vlss_source[ 0 ] )
        if ( cross_radius <= assoc_radius ):
          cross_source_list.append( nvss_source )
          nvss_source_list.remove( nvss_source )
    if ( len( cross_source_list ) > 0 ):
      nvss_flux = float( sum( [ cross_source[ 1 ] for cross_source in cross_source_list ] ) )
      source_sp_index = log( nvss_flux / vlss_source[ 1 ] ) / log( nvss_freq / vlss_freq )
      if ( source_sp_index < si_limits[ 0 ] ):
        source_sp_index = si_limits[ 0 ]
      elif ( source_sp_index > si_limits[ 1 ] ):
        source_sp_index = si_limits[ 1 ]
      source_flux = nvss_flux * pow( freq / nvss_freq, source_sp_index )
      source_list.append( [ vlss_source[ 0 ], source_flux ] )
    else: # no WENSS & NVSS counterpart
      pass
  for nvss_source in nvss_source_list: # no VLSS & WENSS counterpart found
    source_flux = nvss_source[ 1 ]
    if ( not spectral_index is None ):
      source_flux = source_flux * pow( freq / nvss_freq, spectral_index )
    source_list.append( [ nvss_source[ 0 ], source_flux ] )
  if ( len( source_list ) == 0 ):
    return source_list
  
  # apply minimum flux criteria
  source_list.sort( cmp = lambda a, b: cmp( b[ 1 ], a[ 1 ] ) )
  if ( not flux_min is None ):
    for source in source_list:
      if ( source[ 1 ] < flux_min ):
        break
    index = source_list.index( source )
    source_list = source_list[ 0 : index ]
  
  return source_list

###############################################################################

def generate_source_catalog( radec, radius, freq, source_file_name, epoch = 2000.0,
    flux_min = None, append = False, use_nvss = True, use_wenss = True,
    use_vlss = False ):
# radius <= 90 degrees
  
  # generate source list from WENSS & NVSS
  source_list = generate_source_list( radec, radius, freq, epoch = epoch,
      flux_min = flux_min, spectral_index = -0.8, use_nvss = use_nvss,
      use_wenss = use_wenss, use_vlss = use_vlss )
  
  # check for old source list
  old_source_list = []
  source_file_name_e = path.expandvars( source_file_name )
  if file_exists( source_file_name_e ):
    if append:
      source_file = file( source_file_name_e, mode = 'r' )
      for line in source_file:
        words = line.split()
        if ( len( words ) >= 3 ):
          if ( words[ 0 ][ 0 ] != ';' ):
            source_radec = [ float( words[ 0 ] ), float( words[ 1 ] ) ]
            source_flux = float( int( words[ 2 ] ) ) / 1.e3
            old_source_list.append( [ source_radec, source_flux ] )
      source_file.close()
  
  # remove double entries
  copy_source_list = [ x for x in source_list ]
  if ( append and ( len( old_source_list ) > 0 ) ):
    for source1 in old_source_list:
      [ [ ra1, dec1 ], flux1 ] = source1
      for source2 in copy_source_list:
        [ [ ra2, dec2 ], flux2 ] = source2
        if ( ( abs( ra2 - ra1 ) < 2.e-5 ) and ( abs( dec2 - dec1 ) < 2.e-5 ) and
            ( abs( flux2 - flux1 ) < 2.e-3 ) ):
          source_list.remove( source2 )
          break
    if ( len( source_list ) > 0 ):
      source_list = source_list + old_source_list
      source_list.sort( cmp = lambda a, b: cmp( b[ 1 ], a[ 1 ] ) )
  
  # write source list to disk
  if ( len( source_list ) > 0 ):
    if file_exists( source_file_name_e ):
      remove_file( source_file_name_e )
    source_file = file( source_file_name_e, mode = 'w' )
    source_file.write( ';RA(%s) Dec(%s)   Flux\n' % ( repr( int( epoch ) ), 
        repr( int( epoch ) ) ) )
    source_file.write( '; deg       deg        mJy\n' )
    for source in source_list:
      [ source_radec, source_flux ] = source
      # F9.5,1X,F9.5,I7,F10.4
      source_file.write( '%9.5f %9.5f%7d\n' % ( source_radec[ 0 ], source_radec[ 1 ],
          min( int( source_flux * 1.e3 ), 9999999 ) ) )
    source_file.close()
  if ( len( source_list ) > 0 ):
    return ( source_list[ 0 ][ 1 ] )
  else:
    return ( 0. )

###############################################################################

def make_pb_model_from_source_list( uv, source_list, pbm_facet_size = 32.,
    imagr_params = {}, beam_correction = True, flux_min = 0., cutoff = 0.3 ):
  
  if ( source_list == [] ):
    return None
  
  # loop over source list
  radec_list = []
  flux_list = []
  for source in source_list:
    [ source_radec, source_flux ] = source
    radec_list.append( source_radec )
    flux_list.append( source_flux )
  
  # loop over candidate list
  if beam_correction:
    pb_source_list = []
    A_list = get_primary_beam_attenuations( uv, radec_list, cutoff = cutoff )
    for i in range( len( radec_list ) ):
      if ( A_list[ i ] >= cutoff ):
        source_flux = A_list[ i ] * flux_list[ i ]
        if ( source_flux >= flux_min ):
          pb_source_list.append( [ radec_list[ i ], source_flux ] )
  else:
    pb_source_list = [ [ x for x in y ] for y in source_list ]
  pb_source_list.sort( cmp = lambda a, b: cmp( b[ 1 ], a[ 1 ] ) )
  
  # create dummy facet 
  pbm_facets = get_aips_file( uv.disk, 'PBM', 'IMO001', -1, 'MA' )
  empty_facet = get_aips_file( uv.disk, 'EMPTY', 'ICL001', -1, 'MA' )
  channel_count = get_channel_count( uv )
  cell_size = restore_parameter( uv, 'cell_size' )
  call_aips_task( 'IMAGR', indata = uv, nchav = channel_count, flagver = -1,
      outdisk = empty_facet.disk, outname = empty_facet.name, outseq = empty_facet.seq,
      outver = 0, cellsize = [ cell_size, cell_size ], flux = 100000., gain = 0.1,
      imsize = [ pbm_facet_size, pbm_facet_size ], do3dimag = 1, niter = 100000,
      fldsize = [ [ pbm_facet_size, pbm_facet_size ] ], nfield = 1, overlap = 2,
      cmethod = '', minpatch = pbm_facet_size - 1, dotv = 0, **imagr_params )
  fill_image( empty_facet, do_all = True, value = 0. )
  empty_facet.zap_table( 'CC', 0 )
  empty_beam = get_facet_beam( empty_facet )
  empty_beam.zap()
  
  # loop over source list
  if True:
#  if ( source_list == [] ):
#    # generate one facet with zero flux source model (should we??)
#    source_count = 1
#    pbm_facet_i = get_facet( pbm_facets, source_count )
#    call_aips_task( 'MOVE', indata = empty_facet, outdata = pbm_facet_i, userid = get_aips_userid(), opcode = '' )
#    facet_pixel_ref = get_pixel_reference( pbm_facet_i )
#    call_aips_task( 'CCMOD', indata = pbm_facet_i, invers = -1, opcode = 'POIN', flux = 0., pixxy = facet_pixel_ref )
#    pbm_facet_count = 1
#  else:
    i = 0
    for source in pb_source_list:
      i = i + 1
      [ source_radec, source_flux ] = source
      pbm_facet_i = get_facet( pbm_facets, i )
      call_aips_task( 'MOVE', indata = empty_facet, outdata = pbm_facet_i,
          userid = get_aips_userid() )
      set_radec( pbm_facet_i, source_radec )
      facet_pixel_ref = get_pixel_reference( pbm_facet_i )
      call_aips_task( 'CCMOD', indata = pbm_facet_i, invers = -1,
          opcode = 'POIN', flux = source_flux, pixxy = facet_pixel_ref )
    pbm_facet_count = i
  
  # delete dummy facet
  empty_facet.zap()
  
  # store parameters
  store_parameter( uv, 'pbm_facet_count', pbm_facet_count )
  store_parameter( pbm_facets, 'facet_count', pbm_facet_count )
  
  return pbm_facets

###############################################################################

def make_pb_model( uv, flux_min, pbm_facet_size = 32, pb_factor = 1.2,
    use_nvss = True, use_wenss = True, epoch_correction = False, cutoff = 0.3,
    imagr_params = {}, use_vlss = False ):
  
  # generate source catalog
  radec = get_radec( uv )
  field_size = restore_parameter( uv, 'field_size' )
  search_radius = pb_factor * 0.5 * field_size
  frequency = get_central_frequency( uv )
  epoch = get_epoch( uv )
  source_list = generate_source_list( radec, search_radius, frequency, epoch = epoch,
      assoc_radius = 40. / 3600., spectral_index = -0.8, use_nvss = use_nvss,
      use_wenss = use_wenss, use_vlss = use_vlss )
  
  # adjust source_list to observing epoch
  if epoch_correction:
    obs_epoch = get_observing_epoch( uv )
    obs_radec = convert_radec_from_j2000( radec, obs_epoch )
    [ radius, angle ] = calculate_angular_separation( obs_radec, radec )
    for i in range( len( source_list ) ):
      [ source_radec, source_flux ] = source_list[ i ]
      obs_source_radec = convert_radec_from_j2000( source_radec, obs_epoch )
      new_source_radec = calculate_offset_position( obs_source_radec, radius, angle )
      source_list[ i ] = [ new_source_radec, source_flux ]
  
  # loop over source list
  pb_source_list = []
  for source in source_list:
    [ source_radec, source_flux ] = source
    if ( source_flux >= flux_min ):
      [ radius, angle ] = calculate_angular_separation( radec, source_radec )
      if ( radius <= search_radius ):
        pb_source_list.append( source )
  
  # generate facets containing source models multiplied with primary beam shape
  pbm_facets = make_pb_model_from_source_list( uv, pb_source_list, 
      pbm_facet_size = pbm_facet_size, imagr_params = imagr_params,
      flux_min = flux_min, cutoff = cutoff )
  
  return pbm_facets

###############################################################################

def define_pb_facets( uv, max_smearing = 0.1, max_facet_size = 5.,
    central_field_ratio = 0.25, overlap = 10, file_prefix = '${DAT}/' ):
# max_smearing is max. bandwidth/time averaging smearing in units of beamwidth
# max_facet_size is max. facet size in degrees

  # calculate size of FOV
  cell_size = restore_parameter( uv, 'cell_size' )
  field_size = restore_parameter( uv, 'field_size' )
  pb_image_size = int( ceil( 3600. * field_size / cell_size ) ) # pixels
  cpb_image_size = int( ceil( central_field_ratio * 3600. * field_size / cell_size ) ) # pixels

  # determine facet size in arcsec
  # limiting factors for facet size:
  # ( 1 ) bandwidth smearing
  # ( 2 ) time averaging smearing
  # ( 3 ) flat sky approximation (through max_facet_size parameter)
  # ( 4 ) ionospheric patch size (through max_facet_size parameter)
  # ( 5 ) ionospheric changes time scale (not incorporated)
#  beam_size = restore_parameter( uv, 'beam_size' )
#  integration_time = restore_parameter( uv, 'integration_time' )
#  earth_rotation_rate = ( ( 2. * pi ) / 86144. ) # rad / s
#  bw_limit = ( max_smearing * beam_size ) / ( ( get_channel_width( uv ) / frequency ) * ( cell_size / 2. ) )
#  ta_limit = ( max_smearing * beam_size ) / ( ( earth_rotation_rate * integration_time ) * ( cell_size / 2. ) )
#  fs_ip_limit = max_facet_size * 3600. / cell_size
#  print 'facet limits (bandw.,time,flat sky/ion. patch) = ', repr( [ bw_limit, ta_limit, fs_ip_limit ] )
#  limit = min( [ bw_limit, ta_limit, fs_ip_limit ] )
#  pb_facet_size = pow( 2., floor( log10( limit ) / log10( 2. ) ) ) # pixels
#  w_limit = ( 3600. / ( cell_size / 2. ) ) * ( 360. / ( 2. * pi ) ) * sqrt( 2. * max_smearing * w_max )
  fs_ip_limit = max_facet_size * 3600. / cell_size
  w_max = restore_parameter( uv, 'w_max' )
  w_limit = ( 3600. / cell_size ) * degrees( sqrt( 8. * max_smearing / w_max ) )
  limit = min( [ fs_ip_limit, w_limit ] )
  pb_facet_size = pow( 2., floor( log10( limit ) / log10( 2. ) ) ) # pixels

  # define facets for imaging primary beam
  pb_radius = 0.5 * field_size
  pb_facet_file_name = file_prefix + aips_file_name_to_string( uv ) + '.PB.BOX'
  pb_facet_file_name_e = path.expandvars( pb_facet_file_name )
  if file_exists( pb_facet_file_name_e ):
    remove_file( pb_facet_file_name_e )
  outputs = call_aips_task( 'SETFC', indata = uv, cellsize = [ cell_size, cell_size ],
      imsize = [ pb_facet_size, pb_facet_size ], bparm = [ pb_radius, overlap, 0, 0, 0, 0, 0, 0, 0, 0 ],
      boxfile = pb_facet_file_name_e, outputs = [ 'nfield' ] )
  pb_facet_count = int( around( outputs[ 0 ] ) )

  # replace clean box definitions from current boxfile with dummy clean boxes
  temp_facet_file_name = pb_facet_file_name_e + '.TEMP'
  if file_exists( temp_facet_file_name ):
    remove_file( temp_facet_file_name )
  facet_list = range( 1, pb_facet_count + 1 )
  extract_facet_definitions( pb_facet_file_name_e, facet_list, temp_facet_file_name,
      include_clean_boxes = False )
  remove_file( pb_facet_file_name_e )
  move_file( temp_facet_file_name, pb_facet_file_name_e )
  for i in range( 1, pb_facet_count + 1 ):
    add_rectangular_clean_box( pb_facet_file_name_e, i, [ 0, 0 ], [ 0, 0 ] )

  # count facets in central part of primary beam (used for noise measurements)
  cpb_radius = central_field_ratio * pb_radius
  cpb_facet_file_name = pb_facet_file_name_e + '.CPB'
  if file_exists( cpb_facet_file_name ):
    remove_file( cpb_facet_file_name )
  outputs = call_aips_task( 'SETFC', indata = uv, cellsize = [ cell_size, cell_size ],
      imsize = [ pb_facet_size, pb_facet_size ], bparm = [ cpb_radius, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
      boxfile = cpb_facet_file_name, outputs = [ 'nfield' ] )
  cpb_facet_count = int( around( outputs[ 0 ] ) )
  remove_file( cpb_facet_file_name )

  # store parameters with UV data
#  store_parameter( uv, 'field_size', field_size )
  store_parameter( uv, 'pb_image_size', pb_image_size )
  store_parameter( uv, 'pb_facet_size', pb_facet_size )
  store_parameter( uv, 'pb_facet_file_name', pb_facet_file_name )
  store_parameter( uv, 'pb_facet_count', pb_facet_count )
  store_parameter( uv, 'cpb_image_size', cpb_image_size )
  store_parameter( uv, 'cpb_facet_count', cpb_facet_count )

  return

###############################################################################

# outliers just outside primary beam
def define_o_facets( uv, flux_min, field_ratio = 2., o_facet_size = None,
    append = False, inner_facet_overlap = 0.25, use_nvss = True,
    use_wenss = True, use_vlss = False, print_info = True,
    file_prefix = '${DAT}/' ): #overlap = 10, 
  
  # generate source catalog
  radec = get_radec( uv )
  file_prefix_e = path.expandvars( file_prefix )
  catalog_name = file_prefix_e + aips_file_name_to_string( uv ) + '.OM.CAT'
  if file_exists( catalog_name ):
    remove_file( catalog_name )
  o_facet_file_name = file_prefix + aips_file_name_to_string( uv ) + '.O.BOX'
  o_facet_file_name_e = path.expandvars( o_facet_file_name )
  field_size = restore_parameter( uv, 'field_size' )
  outer_radius = 0.5 * field_size * field_ratio
  frequency = get_central_frequency( uv )
  max_flux = generate_source_catalog( radec, outer_radius, frequency, catalog_name,
      flux_min = flux_min, use_nvss = use_nvss, use_wenss = use_wenss,
      use_vlss = use_vlss )
  if ( o_facet_size is None ):
    o_facet_size = restore_parameter( uv, 'pb_facet_size' )
  if ( ( not file_exists( catalog_name ) ) and ( not append ) ):
    make_file( o_facet_file_name_e )
    store_parameter( uv, 'o_facet_file_name', o_facet_file_name )
    store_parameter( uv, 'o_facet_count', 0 )
    store_parameter( uv, 'o_facet_size', o_facet_size )
    return
  
  # define facets for imaging bright outlier sources just outside the primary beam
  # (actually allow for a small overlap with primary beam for catching sources on the edge) 
  cell_size = restore_parameter( uv, 'cell_size' )
  pb_facet_size = restore_parameter( uv, 'pb_facet_size' )
  inner_radius = ( 0.5 * field_size ) - ( 
      inner_facet_overlap * pb_facet_size * cell_size / 3600. )
  if file_exists( o_facet_file_name_e ):
    if append:
      bcount = get_facet_count( o_facet_file_name_e, count_gaps = True ) + 1
    else:
      remove_file( o_facet_file_name_e )
      bcount = 1
  else:
    remove_file( o_facet_file_name_e )
    bcount = 1
  cell_size = 0.
  overlap = 0
  flux = 1.e-9
#  outputs = call_aips_task( 'SETFC', indata = uv, cellsize = [ cell_size, cell_size ],
#      bparm = [ - inner_radius, overlap, 0, outer_radius, flux, o_facet_size,
#      0, 1, 0, 0 ], boxfile = o_facet_file_name, outputs = [ 'nfield' ], 
#      inlist = catalog_name, bcount = bcount, imsize = [ o_facet_size, o_facet_size ] )
  # check for old source list
  catalog_file = file( catalog_name, mode = 'r' )
  for line in catalog_file:
    words = line.split()
    if ( len( words ) >= 3 ):
      if ( words[ 0 ][ 0 ] != ';' ):
        source_radec = [ float( words[ 0 ] ), float( words[ 1 ] ) ]
        source_flux = float( int( words[ 2 ] ) ) / 1.e3
        [ r, p ] = calculate_angular_separation( radec, source_radec )
        if ( ( r > inner_radius ) and ( r < outer_radius ) and
            ( source_flux > flux ) ):
          add_facet( o_facet_file_name_e, source_radec, [ o_facet_size, o_facet_size ],
              facet_id = bcount, add_clean_box = True )
          bcount = bcount + 1
  catalog_file.close()
  if print_info:
    print '... added %d nearby outlier facets' % ( bcount )
  
  # remove double entries
  if append:
    removed = []
    facet_list = get_facet_list( o_facet_file_name_e )
    copy_facet_list = [ x for x in facet_list ]
    for i in range( len( copy_facet_list ) ):
      [ ra1, dec1 ] = copy_facet_list[ i ][ 1 ]
      for j in range( i + 1, len( copy_facet_list ) ):
        [ ra2, dec2 ] = copy_facet_list[ j ][ 1 ]
        if ( ( abs( ra2 - ra1 ) < 2.e-5 ) and ( abs( dec2 - dec1 ) < 2.e-5 ) ):
          if ( not ( j + 1 in removed ) ):
            facet_list.remove( copy_facet_list[ j ] )
            remove_facet( o_facet_file_name_e, copy_facet_list[ j ][ 0 ] )
            removed.append( j + 1 )
    if ( len( removed ) > 0 ):
      facet_list = range( 1, 1 + len( copy_facet_list ) )
      for j in removed:
        facet_list.remove( j )
      new_o_facet_file_name = o_facet_file_name_e + '.TEMP'
      extract_facet_definitions( o_facet_file_name_e, facet_list, 
          new_o_facet_file_name )
      remove_file( o_facet_file_name_e )
      move_file( new_o_facet_file_name, o_facet_file_name_e )
    if print_info:
      print '... removed %d double entry facets' % ( len( removed ) )
  o_facet_count = get_facet_count( o_facet_file_name_e )
  
  # store parameters with UV data
  store_parameter( uv, 'o_facet_file_name', o_facet_file_name )
  store_parameter( uv, 'o_facet_count', o_facet_count )
  store_parameter( uv, 'o_facet_size', o_facet_size )
  
  return

###############################################################################

# very bright, far away outliers (A-team) and Sun
def define_a_s_facets( uv, a_facet_size = 128., s_ang_facet_size = 1.,
    file_prefix = '${DAT}/' ):
  
  # generate A-team catalog file
  a_list = [
      [ 23, 23, 24,      58, 48, 54    ],  # 3C461  CasA      (SNR)
      [ 19, 59, 28.357,  40, 44,  2.10 ],  # 3C405  CygA      (RG)
      [  5, 34, 31.95,   22,  0, 52.1  ],  # 3C144  TauA,M1   (SNR)
      [ 12, 30, 49.423,  12, 23, 28.04 ],  # 3C274  VirA,M87  (RG)
      [ 16, 51,  9,       4, 59, 34    ],  # 3C348  HerA      (RG)
      [  9, 18,  5.651, -12,  5, 43.99 ],  # 3C218  HyaA      (RG)
      [  4, 37,  4.375,  29, 40, 13.82 ],  # 3C123  PerB      (RG)
      [ 17, 20, 28.164, -0., 58, 46.60 ],  # 3C353            (RG)
      [ 13, 25, 27.6,   -43,  1,  9    ],  #        CenA      (RG)
      [  5, 19, 49.7,   -45, 46, 44    ] ] #        PicA      (RG)
  a_catalog_name = path.expandvars( file_prefix ) + 'A_TEAM.CAT'
  if file_exists( a_catalog_name ):
    remove_file( a_catalog_name )
  cat_file = file( a_catalog_name, mode = 'w' )
  cat_file.write( ";RA(2000) Dec(2000)   Flux\n" )
  cat_file.write( "; deg       deg        mJy\n" )
  for a in a_list:
    x = hmsdms_to_degdeg( a )
    cat_file.write( "%9.5f  %9.5f  1000\n" % ( x[ 0 ], x[ 1 ] ) )
  cat_file.close()
  
  # define facets for A-team sources and the Sun
  cell_size = restore_parameter( uv, 'cell_size' )
  s_facet_size = pow( 2., ceil( log10( 3600. * s_ang_facet_size / cell_size ) / log10( 2. ) ) ) # pixels
  as_facet_file_name = file_prefix + aips_file_name_to_string( uv ) + '.AS.BOX'
  as_facet_file_name = path.expandvars( as_facet_file_name )
  if file_exists( as_facet_file_name ):
    remove_file( as_facet_file_name )
  outputs = call_aips_task( 'SETFC', indata = uv, cellsize = [ cell_size, cell_size ],
      bparm = [ 0, 0, 0, 200., 0., a_facet_size, s_facet_size, 1, 0, 0 ],
      boxfile = as_facet_file_name, inlist = a_catalog_name, outputs = [ 'nfield' ],
       imsize = [ a_facet_size, a_facet_size ] )
  as_facet_count = int( around( outputs[ 0 ] ) )
  a_facet_count = as_facet_count - 1
  s_facet_count = 1
  
  # split as facet file in a and s facet files and change the clean box sizes
  a_facet_list = range( 1, a_facet_count + 1 )
  a_facet_file_name = file_prefix + aips_file_name_to_string( uv ) + '.A.BOX'
  a_facet_file_name_e = path.expandvars( a_facet_file_name )
  a_clean_radius = int( ( a_facet_size - 10. ) / 2. )
  extract_facet_definitions( as_facet_file_name, a_facet_list, a_facet_file_name_e, 
      new_clean_box_radius = a_clean_radius )
  s_facet_list = [ a_facet_count + 1 ]
  s_facet_file_name = file_prefix + aips_file_name_to_string( uv ) + '.S.BOX'
  s_facet_file_name_e = path.expandvars( s_facet_file_name )
  s_clean_radius = int( ( s_facet_size - 10. ) / 2. )
  extract_facet_definitions( as_facet_file_name, s_facet_list, s_facet_file_name_e, 
      new_clean_box_radius = s_clean_radius )
  remove_file( as_facet_file_name )
  
  # store parameters with UV data
  store_parameter( uv, 'a_facet_file_name', a_facet_file_name )
  store_parameter( uv, 'a_facet_count', a_facet_count )
  store_parameter( uv, 'a_facet_size', a_facet_size )
  store_parameter( uv, 's_facet_file_name', s_facet_file_name )
  store_parameter( uv, 's_facet_count', s_facet_count )
  store_parameter( uv, 's_facet_size', s_facet_size )
  
  return

###############################################################################

def find_source_facets( facets, source_radec, facet_list = [], primary_facet_only = False ):
  facet_source_list = []
  if ( len( facet_list ) > 0 ):
    used_facet_list = facet_list
  else:
    used_facet_list = []
    facet_count = restore_parameter( facets, 'facet_count' )
    for i in range( 1, 1 + facet_count ):
      facet_i = get_facet( facets, i )
      facet_size = get_image_size( facet_i )
      [ x, y ] = calculate_source_position( facet_i, source_radec )
      if ( ( x >= 0.5 ) and ( x <= facet_size[ 0 ] + 0.5 ) and 
           ( y >= 0.5 ) and ( y <= facet_size[ 1 ] + 0.5 ) ):
        used_facet_list = [ i ] + get_facet_overlap( facet_i )
        break
  for i in used_facet_list:
    facet_i = get_facet( facets, i )
    facet_size = get_image_size( facet_i )
    [ x, y ] = calculate_source_position( facet_i, source_radec )
    if ( ( x >= 0.5 ) and ( x <= facet_size[ 0 ] + 0.5 ) and 
         ( y >= 0.5 ) and ( y <= facet_size[ 1 ] + 0.5 ) ):
      facet_source_list.append( [ i, [ x, y ] ] )
  if ( ( primary_facet_only ) and ( len( facet_source_list ) > 1 ) ):
    min_distance_2 = 100000000.
    for facet_source in facet_source_list:
      [ i, [ x, y ] ] = facet_source
      facet_i = get_facet( facets, i )
      pixel_ref = get_pixel_reference( facet_i )
      distance_2 = ( x - pixel_ref[ 0 ] )**2 + ( y - pixel_ref[ 1 ] )**2
      if ( distance_2 < min_distance_2 ):
        min_distance_2 = distance_2
        min_facet_source = facet_source
    facet_source_list = [ min_facet_source ]
  return facet_source_list

###############################################################################

def determine_facet_overlap_old( facets ):
# this function requires facets to be rectangular, non-rotated in RA and DEC
# and < 90 degrees in size
  
  # get necessary info from facets
  facet_count = restore_parameter( facets, 'facet_count' )
  facet_list = range( 1, facet_count + 1 )
  facet_info_list = []
  for i in facet_list:
    facet_i = get_facet( facets, i )
    radec_i = get_radec( facet_i )
    pixel_ref_i = get_pixel_reference( facet_i )
    pixel_size_i = get_pixel_size( facet_i, make_absolute = True )
    facet_size_i = get_image_size( facet_i )
    dx_max_i = max( [ pixel_ref_i[ 0 ] - 1, facet_size_i[ 0 ] - pixel_ref_i[ 0 ] ] )
    dy_max_i = max( [ pixel_ref_i[ 1 ] - 1, facet_size_i[ 1 ] - pixel_ref_i[ 1 ] ] )
    max_radius_i = sqrt( ( dx_max_i * pixel_size_i[ 0 ] )**2 + ( dy_max_i * pixel_size_i[ 1 ] )**2 ) / 3600.
    corners_i = [ calculate_source_radec( facet_i, [ 0., 0. ] ),
                  calculate_source_radec( facet_i, [ facet_size_i[ 0 ] + 1, 0. ] ),
                  calculate_source_radec( facet_i, [ 0., facet_size_i[ 1 ] + 1 ] ),
                  calculate_source_radec( facet_i, [ facet_size_i[ 0 ] + 1, facet_size_i[ 1 ] + 1 ] ) ]
    facet_info_list.append( [ i, radec_i, max_radius_i, corners_i, facet_size_i ] )
  
  overlap_list = []
  for facet_info_i in facet_info_list:
    [ i, radec_i, max_radius_i, corners_i, facet_size_i ] = facet_info_i
    facet_i = get_facet( facets, i )
    new_overlap = [ i ]
    
    # search for already found overlaps with previous facets
    for overlap in overlap_list:
      if ( i in overlap ):
        new_overlap.append( overlap[ 0 ] )
    
    # find overlaps with next facets
    for facet_info_j in facet_info_list:
      facets_overlap = False
      [ j, radec_j, max_radius_j, corners_j, facet_size_j ] = facet_info_j
      if ( j > i ):
        facet_j = get_facet( facets, j )
        [ radius_ij, angle_ij ] = calculate_angular_separation( radec_i, radec_j )
        if ( radius_ij <= ( max_radius_i + max_radius_j ) ):
          for corner_i in corners_i:
            corner_ij = calculate_source_position( facet_j, corner_i )
            if ( ( corner_ij[ 0 ] > 0. ) and ( corner_ij[ 0 ] < facet_size_j[ 0 ] ) and
                 ( corner_ij[ 1 ] > 0. ) and ( corner_ij[ 1 ] < facet_size_j[ 1 ] ) ):
              facets_overlap = True
              break
          if ( not facets_overlap ):
            for corner_j in corners_j:
              corner_ji = calculate_source_position( facet_i, corner_j )
              if ( ( corner_ji[ 0 ] > 0. ) and ( corner_ji[ 0 ] < facet_size_i[ 0 ] ) and
                   ( corner_ji[ 1 ] > 0. ) and ( corner_ji[ 1 ] < facet_size_i[ 1 ] ) ):
                facets_overlap = True
                break
      if facets_overlap:
        new_overlap.append( j )
    overlap_list.append( new_overlap )
  
  # store overlap info to facets
  for i in facet_list:
    overlap_found = False
    facet_i = get_facet( facets, i )
    for overlap in overlap_list:
      if ( overlap[ 0 ] == i ):
        overlap_count = len( overlap[ 1 : ] )
        store_parameter( facet_i, 'overlap_count', overlap_count )
        for j in range( 1, 1 + overlap_count ):
          store_parameter( facet_i, 'overlap%d' % j, overlap[ j ] )
      overlap_found = True
      # don't break yet, as there might be newer overlaps
    if ( not overlap_found ):
      store_parameter( facet_i, 'overlap_count', 0 )
  
  return

###############################################################################

def determine_facet_overlap( facets ):
# this function requires facets to be rectangular, non-rotated in RA and DEC
# and < 90 degrees in size
  
  # get necessary info from facets
  facet_count = restore_parameter( facets, 'facet_count' )
  facet_list = range( 1, facet_count + 1 )
  facet_info_list = []
  for i in facet_list:
    facet_i = get_facet( facets, i )
    radec_i = get_radec( facet_i )
    pixel_ref_i = get_pixel_reference( facet_i )
    pixel_size_i = get_pixel_size( facet_i, make_absolute = False )
    facet_size_i = get_image_size( facet_i )
    limits_i = [ [ 0.5, facet_size_i[ 0 ] + 0.5 ], [ 0.5, facet_size_i[ 1 ] + 0.5 ] ]
    dx_max_i = max( [ pixel_ref_i[ 0 ] - limits_i[ 0 ][ 0 ], limits_i[ 0 ][ 1 ] - pixel_ref_i[ 0 ] ] )
    dy_max_i = max( [ pixel_ref_i[ 1 ] - limits_i[ 1 ][ 0 ], limits_i[ 1 ][ 1 ] - pixel_ref_i[ 1 ] ] )
    max_radius_i = ( sqrt( ( dx_max_i * pixel_size_i[ 0 ] )**2 + 
        ( dy_max_i * pixel_size_i[ 1 ] )**2 ) / 3600. )
    corners_i = [ [ limits_i[ 0 ][ 0 ], limits_i[ 1 ][ 1 ] ], [ limits_i[ 0 ][ 0 ], limits_i[ 1 ][ 0 ] ],
        [ limits_i[ 0 ][ 1 ], limits_i[ 1 ][ 0 ] ], [ limits_i[ 0 ][ 1 ], limits_i[ 1 ][ 1 ] ] ]
    corner_angles_i = [ amodulo( degrees( atan2( pixel_ref_i[ 0 ] - corner[ 0 ],
        corner[ 1 ] - pixel_ref_i[ 1 ] ) ), 360. ) for corner in corners_i ]
    facet_info_list.append( [ i, radec_i, pixel_ref_i, limits_i, max_radius_i, corner_angles_i ] )
  
  overlap_list = []
  for facet_info_i in facet_info_list:
    [ i, radec_i, pixel_ref_i, limits_i, max_radius_i, corner_angles_i ] = facet_info_i
    facet_i = get_facet( facets, i )
    new_overlap = [ i ]
    
    # search for already found overlaps with previous facets
    for overlap in overlap_list:
      if ( i in overlap ):
        new_overlap.append( overlap[ 0 ] )
    
    # find overlaps with next facets
    for facet_info_j in facet_info_list:
      [ j, radec_j, pixel_ref_j, limits_j, max_radius_j, corner_angles_j ] = facet_info_j
      if ( j > i ):
        facet_j = get_facet( facets, j )
        [ radius_ij, angle_ij ] = calculate_angular_separation( radec_i, radec_j )
        if ( radius_ij <= ( max_radius_i + max_radius_j ) ):
          # calculate RADEC of facet edge towards other facet
          angle_ij = amodulo( angle_ij, 360. )
          if ( ( angle_ij > corner_angles_i[ 0 ] ) and ( angle_ij <= corner_angles_i[ 1 ] ) ):
            ex_ij = limits_i[ 0 ][ 0 ]
            ey_ij = ( pixel_ref_i[ 1 ] + 
                ( pixel_ref_i[ 0 ] - limits_i[ 0 ][ 0 ] ) * cos( radians( angle_ij ) ) )
          elif ( ( angle_ij > corner_angles_i[ 1 ] ) and ( angle_ij <= corner_angles_i[ 2 ] ) ):
            ey_ij = limits_i[ 1 ][ 0 ]
            ex_ij = ( pixel_ref_i[ 0 ] - 
                ( pixel_ref_i[ 1 ] - limits_i[ 1 ][ 0 ] ) * sin( radians( angle_ij ) ) )
          elif ( ( angle_ij > corner_angles_i[ 2 ] ) and ( angle_ij <= corner_angles_i[ 3 ] ) ):
            ex_ij = limits_i[ 0 ][ 1 ]
            ey_ij = ( pixel_ref_i[ 1 ] + 
                ( limits_i[ 0 ][ 1 ] - pixel_ref_i[ 0 ] ) * cos( radians( angle_ij ) ) )
          else:
            ey_ij = limits_i[ 1 ][ 1 ]
            ex_ij = ( pixel_ref_i[ 0 ] - 
                ( limits_i[ 1 ][ 1 ] - pixel_ref_i[ 1 ] ) * sin( radians( angle_ij ) ) )
          eradec_ij = calculate_source_radec( facet_i, [ ex_ij, ey_ij ] )
          # check if RADEC lies within other facet
          epos_ji = calculate_source_position( facet_j, eradec_ij )
          if ( ( epos_ji[ 0 ] >= limits_j[ 0 ][ 0 ] ) and ( epos_ji[ 0 ] <= limits_j[ 0 ][ 1 ] ) and
              ( epos_ji[ 1 ] >= limits_j[ 1 ][ 0 ] ) and ( epos_ji[ 1 ] <= limits_j[ 1 ][ 1 ] )  ):
            new_overlap.append( j )
    overlap_list.append( [ o for o in new_overlap ] )
  
  # store overlap info to facets
  for i in facet_list:
    overlap_found = False
    facet_i = get_facet( facets, i )
    for overlap in overlap_list:
      if ( overlap[ 0 ] == i ):
        overlap_count = len( overlap[ 1 : ] )
        store_parameter( facet_i, 'overlap_count', overlap_count )
        for j in range( 1, 1 + overlap_count ):
          store_parameter( facet_i, 'overlap%d' % j, overlap[ j ] )
      overlap_found = True
      # don't break yet, as there might be newer overlaps
    if ( not overlap_found ):
      store_parameter( facet_i, 'overlap_count', 0 )
  
  return

###############################################################################

def get_facet_overlap( facet ):
  overlap = []
  overlap_count = restore_parameter( facet, 'overlap_count' )
  for j in range( 1, 1 + overlap_count ):
    overlap.append( restore_parameter( facet, 'overlap%d' % j ) )
  return overlap

###############################################################################

def update_facet_overlap_old( facets, remove_facet_list = [], add_facet_list = [] ):
# this function assumes that facets already contains the ones specified in added_facet_list 

  if ( ( len( remove_facet_list ) == 0 ) and ( len( add_facet_list ) == 0 ) ):
    determine_facet_overlap( facets )
    return

  # restore current overlap info
  # remove old facets from info
  facet_count = restore_parameter( facets, 'facet_count' )
#  facet_list = range( 1, facet_count + 1 )
  overlap_list = []
#  for i in facet_list:
  for i in range( 1, 1 + facet_count ):
    if ( ( i in remove_facet_list ) or ( i in add_facet_list ) ):
      continue
    facet_i = get_facet( facets, i )
    try:
      overlap_count = restore_parameter( facet_i, 'overlap_count' )
    except:
      continue
    overlap = [ i ]
    for j in range( 1, 1 + overlap_count ):
      overlap_i = restore_parameter( facet_i, 'overlap%d' % ( j ) )
      if ( not ( ( overlap_i in remove_facet_list ) or ( overlap_i in remove_facet_list ) ) ):
        overlap.append( overlap_i )
    overlap_list.append( overlap )

  # get necessary info from all facets
  facet_info_list = []
#  for i in facet_list:
  for i in range( 1, 1 + facet_count ):
    facet_i = get_facet( facets, i )
    radec_i = get_radec( facet_i )
    pixel_ref_i = get_pixel_reference( facet_i )
    pixel_size_i = get_pixel_size( facet_i, make_absolute = True )
    facet_size_i = get_image_size( facet_i )
    dx_max_i = max( [ pixel_ref_i[ 0 ] - 1, facet_size_i[ 0 ] - pixel_ref_i[ 0 ] ] )
    dy_max_i = max( [ pixel_ref_i[ 1 ] - 1, facet_size_i[ 1 ] - pixel_ref_i[ 1 ] ] )
    max_radius_i = sqrt( ( dx_max_i * pixel_size_i[ 0 ] )**2 + ( dy_max_i * pixel_size_i[ 1 ] )**2 ) / 3600.
    corners_i = [ calculate_source_radec( facet_i, [ 0., 0. ] ),
                  calculate_source_radec( facet_i, [ facet_size_i[ 0 ], 0. ] ),
                  calculate_source_radec( facet_i, [ 0., facet_size_i[ 1 ] ] ),
                  calculate_source_radec( facet_i, facet_size_i ) ]
    facet_info_list.append( [ i, radec_i, max_radius_i, corners_i, facet_size_i ] )

  # find overlaps for new facets
  for facet_info_i in facet_info_list:
    [ i, radec_i, max_radius_i, corners_i, facet_size_i ] = facet_info_i
    if ( not ( i in add_facet_list ) ):
      continue
    facet_i = get_facet( facets, i )
    new_overlap = [ i ]

    # find overlaps with all facets
    for facet_info_j in facet_info_list:
      facets_overlap = False
      [ j, radec_j, max_radius_j, corners_j, facet_size_j ] = facet_info_j
      if ( j != i ):
        facet_j = get_facet( facets, j )
        [ radius_ij, angle_ij ] = calculate_angular_separation( radec_i, radec_j )
        if ( radius_ij <= ( max_radius_i + max_radius_j ) ):
          for corner_i in corners_i:
            corner_ij = calculate_source_position( facet_j, corner_i )
            if ( ( corner_ij[ 0 ] > 0. ) and ( corner_ij[ 0 ] < facet_size_j[ 0 ] ) and
                 ( corner_ij[ 1 ] > 0. ) and ( corner_ij[ 1 ] < facet_size_j[ 1 ] ) ):
              facets_overlap = True
              break
          if ( not facets_overlap ):
            for corner_j in corners_j:
              corner_ji = calculate_source_position( facet_i, corner_j )
              if ( ( corner_ji[ 0 ] > 0. ) and ( corner_ji[ 0 ] < facet_size_i[ 0 ] ) and
                   ( corner_ji[ 1 ] > 0. ) and ( corner_ji[ 1 ] < facet_size_i[ 1 ] ) ):
                facets_overlap = True
                break
      if facets_overlap:
        new_overlap.append( j )

    for k in range( len( overlap_list ) ):
      if ( overlap_list[ k ][ 0 ] in new_overlap[ 1 : ] ):
        overlap_list[ k ].append( i )
    overlap_list.append( new_overlap )

  # store overlap info to facets
  for i in range( 1, 1 + facet_count ):
    overlap_found = False
    facet_i = get_facet( facets, i )
    for overlap in overlap_list:
      if ( overlap[ 0 ] == i ):
        overlap_count = len( overlap[ 1 : ] )
        store_parameter( facet_i, 'overlap_count', overlap_count )
        for j in range( 1, 1 + overlap_count ):
          store_parameter( facet_i, 'overlap%d' % j, overlap[ j ] )
      overlap_found = True
      # don't break yet, as there might be newer overlaps
    if ( not overlap_found ):
      store_parameter( facet_i, 'overlap_count', 0 )

  return

###############################################################################

def update_facet_overlap( facets, remove_facet_list = [], add_facet_list = [] ):
# this function assumes that facets already contains the ones specified in added_facet_list 
  
  if ( ( len( remove_facet_list ) == 0 ) and ( len( add_facet_list ) == 0 ) ):
    determine_facet_overlap( facets )
    return
  
  # restore current overlap info
  # remove old facets from info
  facet_count = restore_parameter( facets, 'facet_count' )
  overlap_list = []
  for i in range( 1, 1 + facet_count ):
    if ( ( i in remove_facet_list ) or ( i in add_facet_list ) ):
      continue
    facet_i = get_facet( facets, i )
    try:
      overlap_count = restore_parameter( facet_i, 'overlap_count' )
    except:
      continue
    overlap = [ i ]
    for j in range( 1, 1 + overlap_count ):
      overlap_i = restore_parameter( facet_i, 'overlap%d' % ( j ) )
      if ( not ( ( overlap_i in remove_facet_list ) or ( overlap_i in remove_facet_list ) ) ):
        overlap.append( overlap_i )
    overlap_list.append( overlap )
  
  # get necessary info from all facets
  facet_info_list = []
  for i in range( 1, 1 + facet_count ):
    facet_i = get_facet( facets, i )
    radec_i = get_radec( facet_i )
    pixel_ref_i = get_pixel_reference( facet_i )
    pixel_size_i = get_pixel_size( facet_i, make_absolute = False )
    facet_size_i = get_image_size( facet_i )
    limits_i = [ [ 0.5, facet_size_i[ 0 ] + 0.5 ], [ 0.5, facet_size_i[ 1 ] + 0.5 ] ]
    dx_max_i = max( [ pixel_ref_i[ 0 ] - limits_i[ 0 ][ 0 ], limits_i[ 0 ][ 1 ] - pixel_ref_i[ 0 ] ] )
    dy_max_i = max( [ pixel_ref_i[ 1 ] - limits_i[ 1 ][ 0 ], limits_i[ 1 ][ 1 ] - pixel_ref_i[ 1 ] ] )
    max_radius_i = ( sqrt( ( dx_max_i * pixel_size_i[ 0 ] )**2 + 
        ( dy_max_i * pixel_size_i[ 1 ] )**2 ) / 3600. )
    corners_i = [ [ limits_i[ 0 ][ 0 ], limits_i[ 1 ][ 1 ] ], [ limits_i[ 0 ][ 0 ], limits_i[ 1 ][ 0 ] ],
        [ limits_i[ 0 ][ 1 ], limits_i[ 1 ][ 0 ] ], [ limits_i[ 0 ][ 1 ], limits_i[ 1 ][ 1 ] ] ]
    corner_angles_i = [ amodulo( degrees( atan2( pixel_ref_i[ 0 ] - corner[ 0 ],
        corner[ 1 ] - pixel_ref_i[ 1 ] ) ), 360. ) for corner in corners_i ]
    facet_info_list.append( [ i, radec_i, pixel_ref_i, limits_i, max_radius_i, corner_angles_i ] )
  
  # find overlaps for new facets
  for facet_info_i in facet_info_list:
    [ i, radec_i, pixel_ref_i, limits_i, max_radius_i, corner_angles_i ] = facet_info_i
    if ( not ( i in add_facet_list ) ):
      continue
    facet_i = get_facet( facets, i )
    new_overlap = [ i ]
    
    # find overlaps with all facets
    for facet_info_j in facet_info_list:
      facets_overlap = False
      [ j, radec_j, pixel_ref_j, limits_j, max_radius_j, corner_angles_j ] = facet_info_j
      if ( j != i ):
        facet_j = get_facet( facets, j )
        [ radius_ij, angle_ij ] = calculate_angular_separation( radec_i, radec_j )
        if ( radius_ij <= ( max_radius_i + max_radius_j ) ):
          # calculate RADEC of facet edge towards other facet
          angle_ij = amodulo( angle_ij, 360. )
          if ( ( angle_ij > corner_angles_i[ 0 ] ) and ( angle_ij <= corner_angles_i[ 1 ] ) ):
            ex_ij = limits_i[ 0 ][ 0 ]
            ey_ij = ( pixel_ref_i[ 1 ] + 
                ( pixel_ref_i[ 0 ] - limits_i[ 0 ][ 0 ] ) * cos( radians( angle_ij ) ) )
          elif ( ( angle_ij > corner_angles_i[ 1 ] ) and ( angle_ij <= corner_angles_i[ 2 ] ) ):
            ey_ij = limits_i[ 1 ][ 0 ]
            ex_ij = ( pixel_ref_i[ 0 ] - 
                ( pixel_ref_i[ 1 ] - limits_i[ 1 ][ 0 ] ) * sin( radians( angle_ij ) ) )
          elif ( ( angle_ij > corner_angles_i[ 2 ] ) and ( angle_ij <= corner_angles_i[ 3 ] ) ):
            ex_ij = limits_i[ 0 ][ 1 ]
            ey_ij = ( pixel_ref_i[ 1 ] + 
                ( limits_i[ 0 ][ 1 ] - pixel_ref_i[ 0 ] ) * cos( radians( angle_ij ) ) )
          else:
            ey_ij = limits_i[ 1 ][ 1 ]
            ex_ij = ( pixel_ref_i[ 0 ] - 
                ( limits_i[ 1 ][ 1 ] - pixel_ref_i[ 1 ] ) * sin( radians( angle_ij ) ) )
          eradec_ij = calculate_source_radec( facet_i, [ ex_ij, ey_ij ] )
          # check if RADEC lies within other facet
          epos_ji = calculate_source_position( facet_j, eradec_ij )
          if ( ( epos_ji[ 0 ] >= limits_j[ 0 ][ 0 ] ) and ( epos_ji[ 0 ] <= limits_j[ 0 ][ 1 ] ) and
              ( epos_ji[ 1 ] >= limits_j[ 1 ][ 0 ] ) and ( epos_ji[ 1 ] <= limits_j[ 1 ][ 1 ] )  ):
            new_overlap.append( j )
    for k in range( len( overlap_list ) ):
      if ( overlap_list[ k ][ 0 ] in new_overlap[ 1 : ] ):
        overlap_list[ k ].append( i )
    overlap_list.append( new_overlap )
  
  # store overlap info to facets
  for i in range( 1, 1 + facet_count ):
    overlap_found = False
    facet_i = get_facet( facets, i )
    for overlap in overlap_list:
      if ( overlap[ 0 ] == i ):
        overlap_count = len( overlap[ 1 : ] )
        store_parameter( facet_i, 'overlap_count', overlap_count )
        for j in range( 1, 1 + overlap_count ):
          store_parameter( facet_i, 'overlap%d' % j, overlap[ j ] )
      overlap_found = True
      # don't break yet, as there might be newer overlaps
    if ( not overlap_found ):
      store_parameter( facet_i, 'overlap_count', 0 )
  
  return

###############################################################################

def get_pbparms( uv ):
  pbparm3 = restore_parameter( uv, 'pbparm3' )
  pbparm4 = restore_parameter( uv, 'pbparm4' )
  pbparm5 = restore_parameter( uv, 'pbparm5' )
  pbparm6 = restore_parameter( uv, 'pbparm6' )
  pbparm7 = restore_parameter( uv, 'pbparm7' )
  return [ pbparm3, pbparm4, pbparm5, pbparm6, pbparm7 ]

###############################################################################

def get_primary_beam_attenuations( uv, radec_list, cutoff = 0.3 ):
  A_list = []
  pbparms = [ cutoff, 1. ] + get_pbparms( uv )
  radec = get_radec( uv )
  freq = get_central_frequency( uv )
  for source_radec in radec_list:
    [ radius, angle ] = calculate_angular_separation( radec, source_radec )
    A_list.append( calculate_pbparm_attenuation( freq, radius, pbparms ) )
  return A_list

###############################################################################

def get_3C48_PT_model():
  # order = 4 P L C X U K Q
  freqs = [ 73.8e6, 320.e6, 1.51e9, 4.75e9, 8.4e9, 14.9e9, 23.e9, 45.e9 ]
  fluxes = [ 71.0, 42.0, 16.5, 5.48, 3.25, 1.78, 1.13, 0.64 ]
  dfluxes = [ 0.05 * f for f in fluxes ]
  parms = [ 1.31752, -0.74090, -0.16708, 0.01525 ]
  [ new_parms, new_dparms ] = fit_PT_calibrator_spectrum( freqs, fluxes,
      dfluxes, parameters = parms )
  return new_parms

def get_3C138_PT_model():
  # order = 4 P L C X U Q
  freqs = [ 73.8e6, 320.e6, 1.51e9, 4.75e9, 8.4e9, 14.9e9, 45.e9 ]
  fluxes = [ 14.6, 17.5, 8.47, 3.78, 2.52, 1.56, 0.40 ]
  dfluxes = [ 0.05 * f for f in fluxes ]
  parms = [ 1.00761, -0.55629, -0.11134, -0.01460 ]
  [ new_parms, new_dparms ] = fit_PT_calibrator_spectrum( freqs, fluxes,
      dfluxes, parameters = parms )
  return new_parms

def get_3C147_PT_model():
  # order = 4 P L C X U Q
  freqs = [ 73.8e6, 320.e6, 1.51e9, 4.75e9, 8.4e9, 14.9e9, 45.e9 ]
  fluxes = [ 55.0, 52.3, 22.5, 7.94, 4.84, 2.78, 0.91 ]
  dfluxes = [ 0.05 * f for f in fluxes ]
  parms = [ 1.44856, -0.67252, -0.21124, 0.04077 ]
  [ new_parms, new_dparms ] = fit_PT_calibrator_spectrum( freqs, fluxes,
      dfluxes, parameters = parms )
  return new_parms

def get_3C286_PT_model():
  # order = 4 P L C X U K Q
  freqs = [ 73.8e6, 320.e6, 1.51e9, 4.75e9, 8.4e9, 14.9e9, 23.e9, 45.e9 ]
  fluxes = [ 30.3, 26.0, 15.0, 7.47, 5.23, 3.40, 2.59, 1.45 ]
  dfluxes = [ 0.05 * f for f in fluxes ]
  parms = [ 1.23734, -0.43276, -0.14223, 0.00345 ]
  [ new_parms, new_dparms ] = fit_PT_calibrator_spectrum( freqs, fluxes,
      dfluxes, parameters = parms )
  return new_parms

###############################################################################

def plot_PT_model_fit( freqs, fluxs, parms, new_parms ):
  loglog( array( freqs, dtype = float64 ) / 1.e9,
      array( fluxs, dtype = float64 ), 'ko' )
  freq_array = 10.**( ( arange( 101. ) / 100. ) * ( 11. - 7. ) + 7. )
  logf = log10( freq_array / 1.e9 )
  logs = ( parms[ 0 ] + ( parms[ 1 ] * logf ) + ( parms[ 2 ] * logf**2 ) +
      ( parms[ 3 ] * logf**3 ) )
  model = 10.**logs
  logs = ( new_parms[ 0 ] + ( new_parms[ 1 ] * logf ) + 
      ( new_parms[ 2 ] * logf**2 ) + ( new_parms[ 3 ] * logf**3 ) )
  new_model = 10.**logs
  loglog( freq_array / 1.e9, model, 'b' )
  loglog( freq_array / 1.e9, new_model, 'r' )
  show()
  return

###############################################################################

def fit_perley_taylor( p, dojac = None, x = None, y = None, err = None ):
  if ( ( not dojac is None ) or ( x is None ) or ( y is None ) or
      ( err is None ) ):
    return -1, None, None
  logf = log10( x / 1.e9 )
  logs = p[ 0 ] + ( p[ 1 ] * logf ) + ( p[ 2 ] * logf**2 ) + ( p[ 3 ] * logf**3 )
  f = 10.**logs
  return 0, ( ( y - f ) / err ), None

def fit_PT_calibrator_spectrum( frequencies, fluxes, flux_errors, 
    parameters = [] ):
  freqs = array( frequencies, dtype = float64 )
  fluxs = array( fluxes, dtype = float64 )
  flux_errs = array( flux_errors, dtype = float64 )
  if ( len( parameters ) == 4 ):
    parms = array( parameters, dtype = float64 )
  else:
    parms = zeros( ( 4 ), dtype = float64 )
  # calculate power-law fits to source fluxes from catalogs
  function_keywords = { 'x' : freqs, 'y' : fluxs, 'err' : flux_errs }
  parameter_info = [ 
      { 'parname' : 'A', 'value' : parms[ 0 ] },
      { 'parname' : 'B', 'value' : parms[ 1 ] },
      { 'parname' : 'C', 'value' : parms[ 2 ] },
      { 'parname' : 'D', 'value' : parms[ 3 ] } ]
  fit = mpfit( fit_perley_taylor, functkw = function_keywords,
      parinfo = parameter_info, quiet = True, autoderivative = True,
      debug = False, fastnorm = False, nocovar = False, dblprec = True )
  if ( not ( fit.status in [ 1, 2, 3 ] ) ):
    raise error( 'ERROR: fit failed, returned %s: %s' %( repr( fit_status ),
        fit.errmsg ) )
  new_parms = fit.params.copy()
  new_dparms = fit.perror
  return [ new_parms.tolist(), new_dparms.tolist() ]

###############################################################################

def get_SH_model_flux_spix( cal_name, frequency, model = None ):
  if ( model is None ):
    calibrators = [ '3C48', '3C147', '3C196', '3C286', '3C295', '3C380', '3C468.1' ]
    models = [
        [ log10( 64.768 ), -0.387, -0.420,  0.181,      0 ],  # 3C48
        [ log10( 66.738 ), -0.022, -1.012,  0.549,      0 ],  # 3C147
        [ log10( 83.084 ), -0.699, -0.110,      0,      0 ],  # 3C196
        [ log10( 27.477 ), -0.158,  0.032, -0.180,      0 ],  # 3C286
        [ log10( 97.763 ), -0.582, -0.298,  0.583, -0.363 ],  # 3C295
        [ log10( 77.352 ), -0.767,      0,      0,      0 ],  # 3C380
        [ log10( 39.688 ), -0.420, -0.807,  0.365, -0.077 ] ] # 3C468.1
    if ( not cal_name in calibrators ):
      raise error( 'calibrator %s not in Scaife & Heald (2012)' % ( cal_name ) )
    index = calibrators.index( cal_name )
    selected_model = array( models[ index ] )
  else:
    if ( len( model ) != 5 ):
      raise error( 'model should have 5 constants' )
    try:
      selected_model = array( model, dtype = float64 )
    except:
      raise error( 'model should have 5 floating point constants' )
  log_v = log10( frequency / 150.e6 )
  flux_vector = array( [ 1., log_v, log_v**2, log_v**3, log_v**4 ] )
  spix_vector = array( [ 0., 1., 2. * log_v, 3. * log_v**2, 4. * log_v**3 ] )
  flux = 10.**( ( flux_vector * selected_model ).sum() )
  spix = ( spix_vector * selected_model ).sum()
  return [ flux, spix ]

def get_SH_model_flux_spix_array( cal_name, frequency, times = 10000 ):
  calibrators = [ '3C48', '3C147', '3C196', '3C286', '3C295', '3C380', '3C468.1' ]
  models = [
      [ 64.768, -0.387, -0.420,  0.181,  0     ],  # 3C48
      [ 66.738, -0.022, -1.012,  0.549,  0     ],  # 3C147
      [ 83.084, -0.699, -0.110,  0,      0     ],  # 3C196
      [ 27.477, -0.158,  0.032, -0.180,  0     ],  # 3C286
      [ 97.763, -0.582, -0.298,  0.583, -0.363 ],  # 3C295
      [ 77.352, -0.767,  0,      0,      0     ],  # 3C380
      [ 39.688, -0.420, -0.807,  0.365, -0.077 ] ] # 3C468.1
  dmodels = [
      [ 1.761, 0.039, 0.031, 0.060, 0     ],  # 3C48
      [ 2.490, 0.030, 0.167, 0.170, 0     ],  # 3C147
      [ 1.862, 0.014, 0.024, 0,     0     ],  # 3C196
      [ 0.746, 0.033, 0.043, 0.052, 0     ],  # 3C286
      [ 2.787, 0.045, 0.085, 0.116, 0.137 ],  # 3C295
      [ 1.164, 0.013, 0,     0,     0     ],  # 3C380
      [ 1.577, 0.031, 0.107, 0.149, 0.055 ] ] # 3C468.1
  if ( not cal_name in calibrators ):
    raise error( 'calibrator %s not in Scaife & Heald (2012)' % ( cal_name ) )
  index = calibrators.index( cal_name )
  A0 = log10( draw_from_gaussian_distribution( models[ index ][ 0 ],
      dmodels[ index ][ 0 ] ) )
  selected_model = [ A0 ]
  for i in range( 1, 5 ):
    Ai = draw_from_gaussian_distribution( models[ index ][ i ],
      dmodels[ index ][ i ] )
    selected_model.append( Ai )
  selected_model = array( selected_model )
  log_v = log10( frequency / 150.e6 )
  flux_vector = array( [ 1., log_v, log_v**2, log_v**3, log_v**4 ] )
  spix_vector = array( [ 0., 1., 2. * log_v, 3. * log_v**2, 4. * log_v**3 ] )
  flux = 10.**( dot( flux_vector, selected_model ) )
  spix = dot( spix_vector, selected_model )
  return [ flux, spix ]

###############################################################################

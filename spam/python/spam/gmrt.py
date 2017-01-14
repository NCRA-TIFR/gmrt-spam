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

def fix_fq_table( uv ):
  old_fq = uv.table( 'FQ', 1 )
  for old_row in old_fq:
    break
  try:
    dummy = old_row.side_band
  except:
    pass
  else:
    new_fq = new_table( uv, 'FQ', 2 )
    for key in old_fq.keywords:
      new_fq.keywords[ key ] = old_fq.keywords[ key ]
    new_row = new_table_row( new_fq )
    for old_row in old_fq:
      new_row.ch_width = old_row.ch_width
      new_row.if_freq = old_row.if_freq
      new_row.total_bandwidth = old_row.total_bandwidth
      new_row.frqsel = old_row.frqsel
      new_row.sideband = old_row.side_band
      new_fq.append( new_row )
    new_fq.close()
    uv.zap_table( 'FQ', 1 )
    call_aips_task( 'TACOP', indata = uv, inext = 'FQ', invers = 2, ncount = 1,
        outdata = uv, outvers = 1 )
    uv.zap_table( 'FQ', 2 )
  return

###############################################################################

def fix_time_offset( uv, time_offset ):
  # time offset in seconds

  # initial checks
#  if uv.header.object.strip() == 'MULTI':
#    raise error( 'only single source data is currently supported' )
  if uv.header.sortord.strip() != 'TB':
    raise error( 'only Time-Baseline sorted UV data is currently supported' )

  # axis definitions and sizes
  naxis = uv.header.naxis[ 0 : uv.header.ndim ]
  ctype = [ ctype.strip() for ctype in uv.header.ctype[ 0 : uv.header.ndim ] ]
  crval = uv.header.crval[ 0 : uv.header.ndim ]
  cdelt = uv.header.cdelt[ 0 : uv.header.ndim ]
  crpix = uv.header.crpix[ 0 : uv.header.ndim ]
  crota = uv.header.crota[ 0 : uv.header.ndim ]
  ra_index = ctype.index( 'RA' )
  dec_index = ctype.index( 'DEC' )
  if_index = ctype.index( 'IF' )
  freq_index = ctype.index( 'FREQ' )
  stokes_index = ctype.index( 'STOKES' )
  complex_index = ctype.index( 'COMPLEX' )

  # group parameters definition (is not used)
  ptype = uv.header.ptype[ 0 : uv.header.nrparm ]
  u_index = [ ( pt.find( 'UU' ) != -1  ) for pt in ptype ].index( True )
  v_index = [ ( pt.find( 'VV' ) != -1  ) for pt in ptype ].index( True )
  w_index = [ ( pt.find( 'WW' ) != -1  ) for pt in ptype ].index( True )
  baseline_index = [ ( pt.find( 'BASELINE' ) != -1  ) for pt in ptype ].index( True )
  time_index = [ ( pt.find( 'TIME' ) != -1  ) for pt in ptype ].index( True )

  # do some data checks
  if ( crval[ ra_index ] == 0 ) or ( crval[ dec_index ] == 0 ):
    raise error( 'RA and/or DEC not present in header' )
  if crval[ if_index ] != 1:
    raise error( 'only single IF data is currently supported' )
  if crval[ freq_index ] == 0:
    raise error( 'FREQ not present in header' )
  if ( complex_index != 0 ) or ( stokes_index != 1 ) or ( freq_index != 2 ) or ( if_index != 3 ):
    raise error( 'UV data groups not organised as COMPLEX-STOKES-FREQ-IF' )

  # create channel frequency table and scale table
#  freq = array( type = Float32, shape = [ naxis[ freq_index ] ] )
#  scale = array( type = Float32, shape = [ naxis[ freq_index ] ] )
#  for freq_i in range( naxis[ freq_index ] ):
#    freq[ freq_i ] = crval[ freq_index ] + ( float( freq_i ) + 1. - crpix[ freq_index ] ) * cdelt[ freq_index ]
#    scale[ freq_i ] = freq[ freq_i ] / crval[ freq_index ]

###########################################################################
# The following code is based on Dave Green's AIPS task UVFXT
###########################################################################

  # Calculate UVW to XYZ and ROTATION matrices (using reference
  # epoch RA and DEC, converted deg to rad; this is not strictly correct
  # due to precession)
  alpha = crval[ ra_index ] * ( pi / 180. )
  delta = crval[ dec_index ] * ( pi / 180. )
  theta = ( time_offset / 86400. ) * ( 2. * pi )
  uvw2xyz = array ( [ [ - sin( alpha ), - cos( alpha ) * sin( delta ), cos( alpha ) * cos( delta ) ],
                      [   cos( alpha ), - sin( alpha ) * sin( delta ), sin( alpha ) * cos( delta ) ],
                      [             0.,                  cos( delta ),                sin( delta ) ] ], 
                    dtype = float64 )
  rot_z = array( [ [ cos( theta ), - sin( theta ), 0. ],
                   [ sin( theta ),   cos( theta ), 0. ],
                   [           0.,             0., 1. ] ], dtype = float64 )

  # convert UVW and time in database
  mod = dot( transpose( uvw2xyz ), dot( rot_z, uvw2xyz ) )
  wiz_uv = wizardry( uv )
  for group in wiz_uv:
    old_uvw = array( [ [ group.uvw[ i ] ] for i in range( 3 ) ], dtype = float64 )
    new_uvw = dot( mod, old_uvw )
    group.uvw = [ new_uvw[ i ][ 0 ] for i in range( 3 ) ]
    new_time = group.time + ( time_offset / 86400. )
    group.time = new_time
    group.update()

  # make sure last update was also written (was bug in ParselTongue)
  for group in wiz_uv:
    break

  return

###############################################################################

def calculate_haslam_tsys_correction_old( uvim, Tsys, cal_name, beta = [ -2.5, 0.1 ], 
    print_info = True, apply_correction = False ):
  
  if print_info:
    print 'checking Tsky'
  
  # set receiver and ground temperatures
  if ( uvim.header.telescop == 'GMRT' ):
    freq = get_frequency( uvim )
    if ( ( freq > 100.e6 ) and ( freq < 200.e6 ) ):
      Trg = 295. + 12.
    elif ( ( freq > 200.e6 ) and ( freq < 300.e6 ) ):
      Trg = 106. + 32.
    elif ( ( freq > 300.e6 ) and ( freq < 400.e6 ) ):
      Trg = 53. + 13.
    elif ( ( freq > 550.e6 ) and ( freq < 700.e6 ) ):
      Trg = 60. + 32.
    elif ( ( freq > 900.e6 ) and ( freq < 1700.e6 ) ):
      Trg = 45. + 24.
    else:
      raise error( 'frequency not recognized' )
  else:
    raise error( 'telescope not recognized' )
  
  # calculate calibrator system temperature
  if ( cal_name == '3C48' ):
    Tsky_cal = draw_from_gaussian_distribution( 24.1, 0.3 )
  elif ( cal_name == '3C147' ):
    Tsky_cal = draw_from_gaussian_distribution( 35.3, 0.4 )
  elif ( cal_name == '3C196' ):
    Tsky_cal = draw_from_gaussian_distribution( 21.9, 0.4 )
  elif ( cal_name == '3C286' ):
    Tsky_cal = draw_from_gaussian_distribution( 21.1, 0.2 )
  elif ( cal_name == '3C295' ):
    Tsky_cal = draw_from_gaussian_distribution( 21.3, 0.3 )
  else:
    raise error( 'calibrator name not recognized' )
  beta_cal = draw_from_gaussian_distribution( -2.5, 0.1 )
  Tsys_cal = Trg + Tsky_cal * ( freq / 408.e6 )**beta_cal
  
  # calculate target system temperature
  Tsky_target = draw_from_gaussian_distribution( Tsys[ 0 ], Tsys[ 1 ] )
  beta_target = draw_from_gaussian_distribution( beta[ 0 ], beta[ 1 ] )
  Tsys_target = Trg + Tsky_target * ( freq / 408.e6 )**beta_target
  
  # calculate flux correction
  corrections = get_robust_mean_deviations( Tsys_target / Tsys_cal )
  if print_info:
    print '... Tsys-based correction = %f +/- %f' % ( corrections[ 0 ],
        max( corrections[ 1 ], -corrections[ 2 ] ) )
  
  # correct visibilities or image (or not)
  if ( not apply_correction ):
    return corrections
  if is_image( uvim ):
    cor_uvim = get_aips_file( uvim.disk, uvim.name, 'TCOR', -1, 'MA' )
    units = im.header.bunit
    call_aips_task( 'MATHS', indata = im, outdata = cor_im,
        opcode = 'POLY', cparm = [ 0, corrections[ 0 ], 0,0 ] )
    set_header_keyword( cor_im, 'bunit', units )
  else:
    sol_version = uvim.table_highver( 'SN' )
    if ( sol_version == 0 ):
      reference_antenna = 1
    else:
      reference_antenna = get_reference_antenna( uvim, sol_version )
    sn_version = generate_fluxscale_solutions( uvim, reference_antenna,
        corrections[ 0 ] )
    cal_uv = apply_solution_table( uvim, version = sn_version )
    uv.zap_table( 'SN', sn_version )
    cor_uvim = get_aips_file( uvim.disk, uvim.name, 'TCOR', -1, 'UV' )
    cal_uv.rename( name = cor_uvim.name, klass = cor_uvim.klass,
        seq = cor_uvim.seq )
    for version in range( 1, 1 + sol_version ):
      if table_exists( uvim, 'SN', version ):
        call_aips_task( 'TACOP', indata = uvim, outdata = cor_uvim, inext = 'SN',
            invers = version, outvers = version, ncount = 1 )
  return cor_uvim

###############################################################################

def check_astrometry( im, radec_list, ref_radec_list = None, beam_factor = 2.,
    min_points = 10, plot_offsets = False, rejection = 3., print_info = True,
    apply_correction = False, use_nvss = True ):
  
  if print_info:
    print 'checking astrometry'
  
  # gather image info
  frequency = get_frequency( im )
  beam_size = array( get_beam_size( im ) )
  assoc_radius = max( beam_size ) * beam_factor
  epoch = get_epoch( im )
  image_size = array( get_image_size( im ) )
  pixel_size = array( get_pixel_size( im ) )
  radius = max( image_size * pixel_size ) / ( 2. * 3600. )
  radec = calculate_source_radec( im, ( image_size / 2. ).tolist() )
  
  # build catalogs
  cat = [ [ 'RA', 'DEC', 'FLUX' ] ]
  ref_cat = [ [ 'RA', 'DEC', 'FLUX' ] ]
  for [ ra, dec ] in radec_list:
    cat.append( [ ra, dec, 1. ] )
  if ( ref_radec_list is None ):
    source_list = generate_source_list( radec, radius, frequency, epoch = epoch,
        use_nvss = use_nvss, use_wenss = ( not use_nvss ) )
    for [ [ ra, dec ], flux ] in source_list:
      ref_cat.append( [ ra, dec, flux ] )
  else:
    for [ ra, dec ] in ref_radec_list:
      ref_cat.append( [ ra, dec, 1. ] )
  
  # check astrometry
  [ match_cat, match_ref_cat ] = match_source_catalogs( cat, ref_cat, 
      assoc_radius = assoc_radius )
  if ( len( match_cat ) <= min_points ):
    raise error( 'there are too few matches to determine RA DEC offset' )
  dradec_list = []
  for i in range( 1, len( match_cat ) ):
    source = match_cat[ i ] 
    ref_source = match_ref_cat[ i ]
    source_radec = source[ 0 : 2 ]
    ref_radec = ref_source[ 0 : 2 ]
    [ r, p ] = calculate_angular_separation( ref_radec, source_radec )
    dra = 3600. * r * sin( radians( p ) )
    ddec = 3600. * r * cos( radians( p ) )
    dradec_list.append( [ dra, ddec ] )
  dradec_array = array( dradec_list, dtype = float64 )
  while ( len( dradec_array ) >= min_points ):
    dra_median = median( dradec_array[ : , 0 ] )
    ddec_median = median( dradec_array[ : , 1 ] )
    dradec_median = array( [ dra_median, ddec_median ], dtype = float64 )
    dr_array = sqrt( add.reduce( ( dradec_array - dradec_median )**2, 1 ) )
    dr_std = 1.4826 * median( dr_array )
    sel = awhere( dr_array < rejection * dr_std )
    if ( len( sel ) == len( dr_array ) ):
      break
    dradec_array = aget( dradec_array, sel )
  dradec_array = array( dradec_list, dtype = float64 )
  
  if print_info:
    print '... median dRA,dDEC = %f, %f arcsec' % ( dra_median, ddec_median )
    print '... dRADEC std = %f arcsec' % ( dr_std )
    print '... based on %d source matches' % len( sel )
  
  # plot result
  if plot_offsets:
    figure( figsize = [ 6.,5.7 ] )
#    max_dra = max( abs( dradec_array[ : , 0 ] ) )
#    max_ddec = max( abs( dradec_array[ : , 1 ] ) )
#    max_dradec = max( max_dra, max_ddec )
    max_dradec = assoc_radius
    cp = arange( 0., 360. )
    cra = dr_std * cos( radians( cp ) )
    cdec = dr_std * sin( radians( cp ) )
    plot( dradec_array[ : , 0 ], dradec_array[ : , 1 ], 'k.' )
    plot( [ -max_dradec, max_dradec ], [ 0.,0. ], 'k' )
    plot( [ 0.,0. ], [ -max_dradec, max_dradec ], 'k' )
    plot( [ -max_dradec, max_dradec ], [ ddec_median, ddec_median ], 'k:' )
    plot( [ dra_median, dra_median ], [ -max_dradec, max_dradec ], 'k:' )
    plot( dra_median + cra, ddec_median + cdec, 'k:' )
    axis( [ -max_dradec, max_dradec, -max_dradec, max_dradec ] )
    plt.gca().invert_xaxis()
    xlabel( r'delta RA ["]' )
    ylabel( r'delta DEC ["]' )
    title( r'median offset = (%4.2f",%4.2f"), std = %5.2f"' % ( dra_median, ddec_median,
        dr_std ) )
    show()
  
  # correct image (or not)
  if ( not apply_correction ):
    return [ dra_median, ddec_median, dr_std ]
  cor_im = get_aips_file( im.disk, im.name, 'POSCOR', -1, 'MA' )
  call_aips_task( 'MOVE', indata = im, outdata = cor_im, userid = get_aips_userid() )
  offset = [ -dra_median / 3600., -ddec_median / 3600. ]
  cor_radec = calculate_offset_position_dradec( get_radec( cor_im ), offset  )
  set_radec( cor_im, cor_radec )
  return cor_im

###############################################################################

def fix_su_table( uv, tolerance = 1. ):
  # read source info, delete qualifiers and calcodes
#  source_list = []
  idno_list = []
  radec_list = []
  change_list = []
  su_table = uv.table( 'SU', 1 )
  for row in su_table:
#    source = row.source.strip()
#    if ( source in source_list ):
#      index = source_list.index( source )
#      idno = idno_list[ index ]
#      change_list.append( [ row.id__no, idno ] )
#    else:
#      source_list.append( source )
#      idno_list.append( row.id__no )
    source_radec = [ row.raobs, row.decobs ]
    source_found = False
    for radec in radec_list:
      [ r, p ] = calculate_angular_separation( source_radec, radec )
      if ( r < tolerance / 3600. ):
        source_found = True
        break
    if source_found:
      index = radec_list.index( radec )
      idno = idno_list[ index ]
      change_list.append( [ row.id__no, idno ] )
    else:
      radec_list.append( source_radec )
      idno_list.append( row.id__no )
  # change source ids in visibilities
  changes = len( change_list )
  for i in range( 0, len( change_list ), 10 ):
    aparm = [ x[ 0 ] for x in change_list[ i : min( changes, i + 10 ) ] ]
    bparm = [ x[ 1 ] for x in change_list[ i : min( changes, i + 10 ) ] ]
    call_aips_task( 'DSORC', indata = uv, outdata = uv, aparm = aparm, bparm = bparm )
  # clobber unused source names
  idno_list = [ x[ 0 ] for x in change_list ]
  wiz_uv = wizardry( uv )
  su_table = wiz_uv.table( 'SU', 1 )
  for row in su_table:
    row.calcode = '    '
    row.qual = 0
    if ( row.id__no in idno_list ):
      new_source = ''
      new_source = new_source.join( choice( ascii_uppercase + digits )
          for i in range( len( row.source ) ) )
      row.source = new_source.strip().ljust( len( row.source ) )
    row.update()
  su_table.close()
  su_table = wiz_uv.table( 'SU', 1 )
  del wiz_uv
  return

###############################################################################

def read_raw_archive_text( raw_file_name ):
  # read in raw text
  lines = []
  raw_file = file( raw_file_name, mode = 'r' )
  header_found = False
  for raw_line in raw_file:
    raw_words = raw_line.split()
    if ( not header_found ):
      if ( len( raw_words ) == 0 ):
        continue
      if ( not ( ( raw_words[ 0 ] == 'Project' ) and ( raw_words[ 1 ] == 'code' ) ) ):
        continue
      header_found = True
      blank_line = True
      line = []
    # read in lines
    if ( len( raw_words ) == 0 ):
      if ( blank_line == True ):
        line.append( '' )
      blank_line = True
      continue
    word = raw_words[ 0 ]
    for raw_word in raw_words[ 1 : ]:
      word = word + '_' + raw_word
    if blank_line:
      blank_line = False
    else:
      lines.append( line )
      line = []
      if ( raw_words[ 0 ][ 0 : 2 ] == '<<' ):
        break
    line.append( word )
  raw_file.close()
  
  # clean up table
  for i in range( len( lines ) ):
    for j in range( len( lines[ i ] ) ):
      word = lines[ i ][ j ]
      if ( word == '' ):
        lines[ i ][ j ] = '-'
        continue
      if ( ( 'm' in word ) and ( 's' in word ) ):
        if ( 'h' in word ):
          h_index = word.index( 'h' )
          m_index = word.index( 'm' )
          s_index = word.index( 's' )
          if ( not ( ( h_index < m_index ) and ( m_index < s_index ) ) ):
            continue
          try:
            h = int( word[ 0 : h_index ] )
            m = int( word[ h_index + 1 : m_index ] )
            s = float( word[ m_index + 1 : s_index ] )
          except ValueError:
            break
            continue
          ra = hmsdms_to_degdeg( [ h,m,s, 0,0,0 ] )[ 0 ]
          lines[ i ][ j ] = ra
          continue
        elif ( 'd' in word ):
          d_index = word.index( 'd' )
          m_index = word.index( 'm' )
          s_index = word.index( 's' )
          if ( not ( ( d_index < m_index ) and ( m_index < s_index ) ) ):
            continue
          try:
            d = int( word[ 0 : d_index ] )
            m = int( word[ d_index + 1 : m_index ] )
            s = float( word[ m_index + 1 : s_index ] )
          except ValueError:
            break
            continue
          dec = hmsdms_to_degdeg( [ 0,0,0, d,m,s ] )[ 1 ]
          lines[ i ][ j ] = dec
          continue
        else:
          continue
      if ( '.' in word ):
        try:
          word = float( word )
        except ValueError:
          continue
        else:
          lines[ i ][ j ] = word
      else:
        try:
          word = int( word )
        except ValueError:
          pass
        else:
          lines[ i ][ j ] = word
  return lines

###############################################################################

def condense_archive_catalog( cat ):
  observation_index = cat[ 0 ].index( 'Observation_No' )
  time_index = cat[ 0 ].index( 'Time_on_src(Mins)' )
  new_cat = []
  time_list = []
  for source in cat[ 1 : ]:
    new_source = [ x for x in source ]
    time = source[ time_index ]
    new_source[ observation_index ] = -1
    new_source[ time_index ] = -1
    if ( not new_source in new_cat ):
      new_cat.append( new_source )
      time_list.append( time )
    else:
      source_index = new_cat.index( new_source )
      time_list[ source_index ] = time_list[ source_index ] + time
  for i in range( len( new_cat ) ):
    new_cat[ i ][ time_index ] = time_list[ i ]
  new_cat = cat[ 0 : 1 ] + new_cat
  return new_cat

###############################################################################

def calculate_haslam_tsys_correction_old( uvim, cal_name, beta = [ -2.5, 0.1 ], 
    print_info = True, apply_correction = False, radius_factor = 1.5,
    haslam_fits_file = '${AIPS_ROOT}/FITS/HASLAM_408MHZ.FITS' ):
  
  if print_info:
    print 'calculating Tsky correction'
  
  # calculate calibrator system temperature
  if ( cal_name == '3C48' ):
    cal_radec = hmsdms_to_degdeg( [ 1,37,41.299431, 33,9,35.132990 ] )
  elif ( cal_name == '3C147' ):
    cal_radec = hmsdms_to_degdeg( [ 5,42,36.137916, 49,51,7.233560 ] )
  elif ( cal_name == '3C196' ):
    cal_radec = hmsdms_to_degdeg( [ 8,13,36.0518, 48,13,2.262 ] )
  elif ( cal_name == '3C286' ):
    cal_radec = hmsdms_to_degdeg( [ 13,31,8.287984, 30,30,32.958850 ] )
  elif ( cal_name == '3C295' ):
    cal_radec = hmsdms_to_degdeg( [ 14,11,20.6477, 52,12,9.141 ] )
  else:
    raise error( 'calibrator name not recognized' )
  cal_gal = equatorial_to_galactic( cal_radec )
  
  # set antenna, receiver and ground temperatures
  if ( uvim.header.telescop == 'GMRT' ):
    frequency = get_frequency( uvim )
    if ( ( frequency > 100.e6 ) and ( frequency < 200.e6 ) ):
      Trg = 295. + 12.
      Tcal = 0.33 * get_SH_model_flux_spix( cal_name, frequency )[ 0 ]
      beam_radius = ( 186. / 120. )
    elif ( ( frequency > 200.e6 ) and ( frequency < 300.e6 ) ):
      Trg = 106. + 32.
      Tcal = 0.33 * get_SH_model_flux_spix( cal_name, frequency )[ 0 ]
      beam_radius = ( 114. / 120. )
    elif ( ( frequency > 300.e6 ) and ( frequency < 400.e6 ) ):
      Trg = 53. + 13.
      Tcal = 0.32 * get_SH_model_flux_spix( cal_name, frequency )[ 0 ]
      beam_radius = ( 81. / 120. )
    elif ( ( frequency > 550.e6 ) and ( frequency < 700.e6 ) ):
      Trg = 60. + 32.
      Tcal = 0.32 * get_SH_model_flux_spix( cal_name, frequency )[ 0 ]
      beam_radius = ( 43. / 120. )
    elif ( ( frequency > 900.e6 ) and ( frequency < 1700.e6 ) ):
      Trg = 45. + 24.
      Tcal = 0.22 * get_SH_model_flux_spix( cal_name, frequency )[ 0 ]
      beam_radius = ( 24. / 120. ) * ( 1.4e9 / frequency )
    else:
      raise error( 'frequency not recognized' )
  else:
    raise error( 'telescope not recognized' )
  
  # read in Haslam et al. map
  haslam_fits_file_e = path.expandvars( haslam_fits_file )
  if ( not file_exists( haslam_fits_file_e ) ):
    raise error( 'Haslam et al. 408 MHz all-sky map is missing' )
  haslam_image = get_aips_file( uvim.disk, 'HASLAM', 'IMAGE', -1, 'MA' )
  read_fits_image( haslam_fits_file_e, haslam_image )
  haslam_frequency = get_frequency( haslam_image )
  haslam_pixels = get_image_pixels( haslam_image )
  ll_index = haslam_image.header.ctype.index( 'LL' )
  mm_index = haslam_image.header.ctype.index( 'MM' )
  image_size = [ haslam_image.header.naxis[ ll_index ], haslam_image.header.naxis[ mm_index ] ]
  pixel_ref = [ haslam_image.header.crpix[ ll_index ], haslam_image.header.crpix[ mm_index ] ]
  pixel_size = [ haslam_image.header.cdelt[ ll_index ], haslam_image.header.cdelt[ mm_index ] ]
  haslam_image.zap()
  
  # calculate coordinate grid
  coord_ref = [ haslam_image.header.crval[ ll_index ], haslam_image.header.crval[ mm_index ] ]
  coord_grid = array( [ [ [ coord_ref[ 0 ] + ( float( x ) + 1.0 - pixel_ref[ 0 ] ) * pixel_size[ 0 ],
      coord_ref[ 1 ] + ( float( y ) + 1.0 - pixel_ref[ 1 ] ) * pixel_size[ 1 ] ] 
      for y in range( image_size[ 1 ] ) ] for x in range( image_size[ 0 ] ) ] )
  
  # determine Tsys towards calibrator
  radial_grid = array( [ [ calculate_angular_separation( cal_gal, coord_grid[ x, y ] )[ 0 ] 
      for y in range( image_size[ 1 ] ) ] for x in range( image_size[ 0 ] ) ] )
#  sel = awhere( array( radial_grid ) < radius )
  sel = awhere( ( radial_grid > 1.2 ) & ( radial_grid < 2.5 ) )
  cal_pixels = aget( haslam_pixels, sel )
  median_pixels = median( cal_pixels )
  sigma_pixels = 1.4826 * median( abs( cal_pixels - median_pixels ) )
  Tsky_cal = draw_from_gaussian_distribution( median_pixels, sigma_pixels )
  beta_cal = draw_from_gaussian_distribution( beta[ 0 ], beta[ 1 ] )
  Tsys_cal = Tcal + Trg + Tsky_cal * ( frequency / haslam_frequency )**beta_cal
  if print_info:
    Tsys = get_robust_mean_deviations( Tsys_cal )
    print '... Tsys for calibrator %s = %f +/- %f' % ( cal_name, Tsys[ 0 ],
        max( Tsys[ 1 ], -Tsys[ 2 ] ) )
  
  # determine Tsys towards target
  target_radec = get_radec( uvim )
  target_gal = equatorial_to_galactic( target_radec )
  radial_grid = array( [ [ calculate_angular_separation( target_gal, coord_grid[ x, y ] )[ 0 ] 
      for y in range( image_size[ 1 ] ) ] for x in range( image_size[ 0 ] ) ] )
  sel = awhere( radial_grid < beam_radius * radius_factor )
  target_pixels = aget( haslam_pixels, sel )
  median_pixels = median( target_pixels )
  sigma_pixels = 1.4826 * median( abs( target_pixels - median_pixels ) )
  Tsky_target = draw_from_gaussian_distribution( median_pixels, sigma_pixels )
  beta_target = draw_from_gaussian_distribution( beta[ 0 ], beta[ 1 ] )
  Tsys_target = Trg + Tsky_target * ( frequency / haslam_frequency )**beta_target
  if print_info:
    Tsys = get_robust_mean_deviations( Tsys_target )
    print '... Tsys for target = %f +/- %f' % ( Tsys[ 0 ], max( Tsys[ 1 ], -Tsys[ 2 ] ) )
  
  # calculate flux correction
  corrections = get_robust_mean_deviations( Tsys_target / Tsys_cal )
  if print_info:
    print '... Tsys-based correction = %f +/- %f' % ( corrections[ 0 ],
        max( corrections[ 1 ], -corrections[ 2 ] ) )
  
  # correct visibilities or image (or not)
  if ( not apply_correction ):
    return corrections
  
  if print_info:
    print '... applying correction'
  if is_image( uvim ):
    cor_uvim = get_aips_file( uvim.disk, uvim.name, 'TCOR', -1, 'MA' )
    units = uvim.header.bunit
    call_aips_task( 'MATHS', indata = uvim, outdata = cor_uvim,
        opcode = 'POLY', cparm = [ 0, corrections[ 0 ], 0,0 ] )
    set_header_keyword( cor_uvim, 'bunit', units )
  else:
    sol_version = uvim.table_highver( 'SN' )
    if ( sol_version == 0 ):
      reference_antenna = 1
    else:
      reference_antenna = get_reference_antenna( uvim, sol_version )
    sn_version = generate_fluxscale_solutions( uvim, reference_antenna,
        corrections[ 0 ] )
    cal_uv = apply_solution_table( uvim, version = sn_version )
    uvim.zap_table( 'SN', sn_version )
    cor_uvim = get_aips_file( uvim.disk, uvim.name, 'TCOR', -1, 'UV' )
    cal_uv.rename( name = cor_uvim.name, klass = cor_uvim.klass,
        seq = cor_uvim.seq )
    for version in range( 1, 1 + sol_version ):
      if table_exists( uvim, 'SN', version ):
        call_aips_task( 'TACOP', indata = uvim, outdata = cor_uvim, inext = 'SN',
            invers = version, outvers = version, ncount = 1 )
  store_parameter( cor_uvim, 'tsys_correction', corrections[ 0 ] )
  store_parameter( cor_uvim, 'tsys_dcorrection', max( corrections[ 1 ], -corrections[ 2 ] ) )
  return cor_uvim

###############################################################################

def calculate_haslam_tsys_correction_old2( uvim, cal_name, beta = [ -2.5, 0.1 ], 
    print_info = True, apply_correction = False, haslam_beam = 0.5,
    haslam_fits_file = '${AIPS_ROOT}/FITS/HASLAM_408MHZ.FITS' ):
  
  if print_info:
    print 'calculating Tsky correction'
  
  # calculate calibrator system temperature
  if ( cal_name == '3C48' ):
    cal_radec = hmsdms_to_degdeg( [ 1,37,41.299431, 33,9,35.132990 ] )
  elif ( cal_name == '3C147' ):
    cal_radec = hmsdms_to_degdeg( [ 5,42,36.137916, 49,51,7.233560 ] )
  elif ( cal_name == '3C196' ):
    cal_radec = hmsdms_to_degdeg( [ 8,13,36.0518, 48,13,2.262 ] )
  elif ( cal_name == '3C286' ):
    cal_radec = hmsdms_to_degdeg( [ 13,31,8.287984, 30,30,32.958850 ] )
  elif ( cal_name == '3C295' ):
    cal_radec = hmsdms_to_degdeg( [ 14,11,20.6477, 52,12,9.141 ] )
  else:
    raise error( 'calibrator name not recognized' )
  cal_gal = equatorial_to_galactic( cal_radec )
  
  # set antenna, receiver and ground temperatures
  if ( uvim.header.telescop == 'GMRT' ):
    frequency = get_frequency( uvim )
    if ( ( frequency > 100.e6 ) and ( frequency < 200.e6 ) ):
      Trg = 295. + 12.
      Tcal = 0.33 * get_SH_model_flux_spix( cal_name, frequency )[ 0 ]
      Tcal_factor = 0.4
      beam_diameter = ( 186. / 60. )
    elif ( ( frequency > 200.e6 ) and ( frequency < 300.e6 ) ):
      Trg = 106. + 32.
      Tcal = 0.33 * get_SH_model_flux_spix( cal_name, frequency )[ 0 ]
      Tcal_factor = 0.95
      beam_diameter = ( 114. / 60. )
    elif ( ( frequency > 300.e6 ) and ( frequency < 400.e6 ) ):
      Trg = 53. + 13.
      Tcal = 0.32 * get_SH_model_flux_spix( cal_name, frequency )[ 0 ]
      Tcal_factor = 1.15
      beam_diameter = ( 81. / 60. )
    elif ( ( frequency > 550.e6 ) and ( frequency < 700.e6 ) ):
      Trg = 60. + 32.
      Tcal = 0.32 * get_SH_model_flux_spix( cal_name, frequency )[ 0 ]
      Tcal_factor = 0.7
      beam_diameter = ( 43. / 60. )
    elif ( ( frequency > 900.e6 ) and ( frequency < 1700.e6 ) ):
      Trg = 45. + 24.
      Tcal = 0.22 * get_SH_model_flux_spix( cal_name, frequency )[ 0 ]
      Tcal_factor = 0.25
      beam_diameter = ( 24. / 60. ) * ( 1.4e9 / frequency )
    else:
      raise error( 'frequency not recognized' )
  else:
    raise error( 'telescope not recognized' )
  Tcal = Tcal_factor * Tcal
  beam_radius = sqrt( haslam_beam**2 + beam_diameter**2 ) / 2.
  
  # read in Haslam et al. map
  haslam_fits_file_e = path.expandvars( haslam_fits_file )
  if ( not file_exists( haslam_fits_file_e ) ):
    raise error( 'Haslam et al. 408 MHz all-sky map is missing' )
  haslam_image = get_aips_file( uvim.disk, 'HASLAM', 'IMAGE', -1, 'MA' )
  read_fits_image( haslam_fits_file_e, haslam_image )
  haslam_frequency = get_frequency( haslam_image )
  haslam_pixels = get_image_pixels( haslam_image )
  ll_index = haslam_image.header.ctype.index( 'LL' )
  mm_index = haslam_image.header.ctype.index( 'MM' )
  image_size = [ haslam_image.header.naxis[ ll_index ], haslam_image.header.naxis[ mm_index ] ]
  pixel_ref = [ haslam_image.header.crpix[ ll_index ], haslam_image.header.crpix[ mm_index ] ]
  pixel_size = [ haslam_image.header.cdelt[ ll_index ], haslam_image.header.cdelt[ mm_index ] ]
  haslam_image.zap()
  
  # calculate coordinate grid
  coord_ref = [ haslam_image.header.crval[ ll_index ], haslam_image.header.crval[ mm_index ] ]
  coord_grid = array( [ [ [ coord_ref[ 0 ] + ( float( x ) + 1.0 - pixel_ref[ 0 ] ) * pixel_size[ 0 ],
      coord_ref[ 1 ] + ( float( y ) + 1.0 - pixel_ref[ 1 ] ) * pixel_size[ 1 ] ] 
      for y in range( image_size[ 1 ] ) ] for x in range( image_size[ 0 ] ) ] )
  
  # determine Tsys towards calibrator
  radial_grid = array( [ [ calculate_angular_separation( cal_gal, coord_grid[ x, y ] )[ 0 ] 
      for y in range( image_size[ 1 ] ) ] for x in range( image_size[ 0 ] ) ] )
  sel = awhere( ( radial_grid > 1.5 ) & ( radial_grid < 3. ) )
  cal_pixels = aget( haslam_pixels, sel )
  median_pixels = median( cal_pixels )
  sigma_pixels = 1.4826 * median( abs( cal_pixels - median_pixels ) )
  Tsky_cal = draw_from_gaussian_distribution( median_pixels, sigma_pixels )
  beta_cal = draw_from_gaussian_distribution( beta[ 0 ], beta[ 1 ] )
  Tsys_cal = Tcal + Trg + Tsky_cal * ( frequency / haslam_frequency )**beta_cal
  if print_info:
    Tsys = get_robust_mean_deviations( Tsys_cal )
    print '... Tsys for calibrator %s = %f +/- %f' % ( cal_name, Tsys[ 0 ],
        max( Tsys[ 1 ], -Tsys[ 2 ] ) )
  
  # determine Tsys towards target
  target_radec = get_radec( uvim )
  target_gal = equatorial_to_galactic( target_radec )
  radial_grid = array( [ [ calculate_angular_separation( target_gal, coord_grid[ x, y ] )[ 0 ] 
      for y in range( image_size[ 1 ] ) ] for x in range( image_size[ 0 ] ) ] )
  sel = awhere( radial_grid < beam_radius )
  target_pixels = aget( haslam_pixels, sel )
  median_pixels = median( target_pixels )
  sigma_pixels = 1.4826 * median( abs( target_pixels - median_pixels ) )
  Tsky_target = draw_from_gaussian_distribution( median_pixels, sigma_pixels )
  beta_target = draw_from_gaussian_distribution( beta[ 0 ], beta[ 1 ] )
  Tsys_target = Trg + Tsky_target * ( frequency / haslam_frequency )**beta_target
  if print_info:
    Tsys = get_robust_mean_deviations( Tsys_target )
    print '... Tsys for target = %f +/- %f' % ( Tsys[ 0 ], max( Tsys[ 1 ], -Tsys[ 2 ] ) )
  
  # calculate flux correction
  corrections = get_robust_mean_deviations( Tsys_target / Tsys_cal )
  if print_info:
    print '... Tsys-based correction = %f +/- %f' % ( corrections[ 0 ],
        max( corrections[ 1 ], -corrections[ 2 ] ) )
  
  # correct visibilities or image (or not)
  if ( not apply_correction ):
    return corrections
  
  if print_info:
    print '... applying correction'
  if is_image( uvim ):
    cor_uvim = get_aips_file( uvim.disk, uvim.name, 'TCOR', -1, 'MA' )
    units = uvim.header.bunit
    call_aips_task( 'MATHS', indata = uvim, outdata = cor_uvim,
        opcode = 'POLY', cparm = [ 0, corrections[ 0 ], 0,0 ] )
    set_header_keyword( cor_uvim, 'bunit', units )
  else:
    sol_version = uvim.table_highver( 'SN' )
    if ( sol_version == 0 ):
      reference_antenna = 1
    else:
      reference_antenna = get_reference_antenna( uvim, sol_version )
    sn_version = generate_fluxscale_solutions( uvim, reference_antenna,
        corrections[ 0 ] )
    cal_uv = apply_solution_table( uvim, version = sn_version )
    uvim.zap_table( 'SN', sn_version )
    cor_uvim = get_aips_file( uvim.disk, uvim.name, 'TCOR', -1, 'UV' )
    cal_uv.rename( name = cor_uvim.name, klass = cor_uvim.klass,
        seq = cor_uvim.seq )
    for version in range( 1, 1 + sol_version ):
      if table_exists( uvim, 'SN', version ):
        call_aips_task( 'TACOP', indata = uvim, outdata = cor_uvim, inext = 'SN',
            invers = version, outvers = version, ncount = 1 )
  store_parameter( cor_uvim, 'tsys_correction', corrections[ 0 ] )
  store_parameter( cor_uvim, 'tsys_dcorrection', max( corrections[ 1 ], -corrections[ 2 ] ) )
  return cor_uvim

###############################################################################

def calculate_haslam_tsys_correction( uvim, cal_name, beta = [ -2.5, 0.1 ], 
    print_info = True, apply_correction = False, haslam_beam = 0.85,
    haslam_fits_file = '${AIPS_ROOT}/FITS/HASLAM_408MHZ.FITS',
    use_cal_model = False, beam_factor = 3. ):
  
  if print_info:
    print 'calculating Tsky correction'
  
  # read in Haslam et al. map
  haslam_fits_file_e = path.expandvars( haslam_fits_file )
  if ( not file_exists( haslam_fits_file_e ) ):
    raise error( 'Haslam et al. 408 MHz all-sky map is missing' )
  haslam_image = get_aips_file( uvim.disk, 'HASLAM', 'IMAGE', -1, 'MA' )
  read_fits_image( haslam_fits_file_e, haslam_image )
  haslam_frequency = get_frequency( haslam_image )
  haslam_pixels = get_image_pixels( haslam_image )
  ll_index = haslam_image.header.ctype.index( 'LL' )
  mm_index = haslam_image.header.ctype.index( 'MM' )
  image_size = [ haslam_image.header.naxis[ ll_index ], haslam_image.header.naxis[ mm_index ] ]
  pixel_ref = [ haslam_image.header.crpix[ ll_index ], haslam_image.header.crpix[ mm_index ] ]
  pixel_size = [ haslam_image.header.cdelt[ ll_index ], haslam_image.header.cdelt[ mm_index ] ]
  haslam_image.zap()
  
  # calculate coordinate grid
  coord_ref = [ haslam_image.header.crval[ ll_index ], haslam_image.header.crval[ mm_index ] ]
  coord_grid = array( [ [ [ coord_ref[ 0 ] + ( float( x ) + 1.0 - pixel_ref[ 0 ] ) * pixel_size[ 0 ],
      coord_ref[ 1 ] + ( float( y ) + 1.0 - pixel_ref[ 1 ] ) * pixel_size[ 1 ] ] 
      for y in range( image_size[ 1 ] ) ] for x in range( image_size[ 0 ] ) ] )
  
  # calculate calibrator system temperature
  if ( cal_name == '3C48' ):
    cal_radec = hmsdms_to_degdeg( [ 1,37,41.299431, 33,9,35.132990 ] )
  elif ( cal_name == '3C147' ):
    cal_radec = hmsdms_to_degdeg( [ 5,42,36.137916, 49,51,7.233560 ] )
  elif ( cal_name == '3C196' ):
    cal_radec = hmsdms_to_degdeg( [ 8,13,36.0518, 48,13,2.262 ] )
  elif ( cal_name == '3C286' ):
    cal_radec = hmsdms_to_degdeg( [ 13,31,8.287984, 30,30,32.958850 ] )
  elif ( cal_name == '3C295' ):
    cal_radec = hmsdms_to_degdeg( [ 14,11,20.6477, 52,12,9.141 ] )
  elif ( cal_name == '3C468.1' ):
    cal_radec = hmsdms_to_degdeg( [ 23,50,54.849, 64,40,19.54 ] )
  else:
    raise error( 'calibrator name not recognized' )
  cal_gal = equatorial_to_galactic( cal_radec )
  
  # set telescope specific parameters
  if ( uvim.header.telescop == 'GMRT' ):
    frequency = get_frequency( uvim )
    if ( ( frequency > 100.e6 ) and ( frequency < 200.e6 ) ):
      Trg = 295. + 12.
      gain = 0.33
      beam_diameter = ( 186. / 60. )
    elif ( ( frequency > 200.e6 ) and ( frequency < 300.e6 ) ):
      Trg = 106. + 32.
      gain = 0.33
      beam_diameter = ( 114. / 60. )
    elif ( ( frequency > 300.e6 ) and ( frequency < 400.e6 ) ):
      Trg = 53. + 13.
      gain = 0.32
      beam_diameter = ( 81. / 60. )
    elif ( ( frequency > 550.e6 ) and ( frequency < 700.e6 ) ):
      Trg = 60. + 32.
      gain = 0.32
      beam_diameter = ( 43. / 60. )
    elif ( ( frequency > 900.e6 ) and ( frequency < 1700.e6 ) ):
      Trg = 45. + 24.
      gain = 0.22
      beam_diameter = ( 24. / 60. ) * ( 1.4e9 / frequency )
    else:
      raise error( 'frequency not recognized' )
  else:
    raise error( 'telescope not recognized' )
  
  # determine Tsys towards calibrator
  beta_cal = draw_from_gaussian_distribution( beta[ 0 ], beta[ 1 ] )
  radial_grid = array( [ [ calculate_angular_separation( cal_gal, coord_grid[ x, y ] )[ 0 ] 
      for y in range( image_size[ 1 ] ) ] for x in range( image_size[ 0 ] ) ] )
  if use_cal_model:
    sel = awhere( ( radial_grid > 1.5 ) & ( radial_grid < 3. ) )
    weights = cos( aradians( aget( coord_grid, sel )[ : , 1 ] ) )
    weight_grid = azeros( radial_grid )
    weight_grid = aput( weight_grid, sel, weights )
    Tsky_cal = ( weight_grid * haslam_pixels ).sum() / weight_grid.sum()
    Tcal = gain * get_SH_model_flux_spix( cal_name, frequency )[ 0 ]
    Tsys_cal = Tcal + Trg + Tsky_cal * ( frequency / haslam_frequency )**beta_cal
  else:
    if ( beam_diameter <= haslam_beam ):
      sigma = max( abs( array( pixel_size ) ) )
    else:
      sigma = sqrt( ( beam_diameter / sqrt( 2. * log( 2. ) ) )**2 - haslam_beam**2 ) / 2.
    sel = awhere( radial_grid < beam_factor * sigma )
    weights = ( cos( aradians( aget( coord_grid, sel )[ : , 1 ] ) ) *
        exp( - aget( radial_grid, sel )**2 / ( 2. * sigma**2 ) ) )
    weight_grid = azeros( radial_grid )
    weight_grid = aput( weight_grid, sel, weights )
    Tsky_cal = ( weight_grid * haslam_pixels ).sum() / weight_grid.sum()
    Tsys_cal = Trg + Tsky_cal * ( frequency / haslam_frequency )**beta_cal
  if print_info:
    Tsys = get_robust_mean_deviations( Tsys_cal )
    print '... Tsys for calibrator %s = %f +/- %f' % ( cal_name, Tsys[ 0 ],
        max( Tsys[ 1 ], -Tsys[ 2 ] ) )
  
  # determine Tsys towards target
  target_radec = get_radec( uvim )
  target_gal = equatorial_to_galactic( target_radec )
  beta_target = draw_from_gaussian_distribution( beta[ 0 ], beta[ 1 ] )
  radial_grid = array( [ [ calculate_angular_separation( target_gal, coord_grid[ x, y ] )[ 0 ] 
      for y in range( image_size[ 1 ] ) ] for x in range( image_size[ 0 ] ) ] )
  if ( beam_diameter <= haslam_beam ):
    sigma = max( abs( array( pixel_size ) ) )
  else:
    sigma = sqrt( ( beam_diameter / sqrt( 2. * log( 2. ) ) )**2 - haslam_beam**2 ) / 2.
  sel = awhere( radial_grid < beam_factor * sigma )
  weights = ( cos( aradians( aget( coord_grid, sel )[ : , 1 ] ) ) *
      exp( - aget( radial_grid, sel )**2 / ( 2. * sigma**2 ) ) )
  weight_grid = azeros( radial_grid )
  weight_grid = aput( weight_grid, sel, weights )
  Tsky_target = ( weight_grid * haslam_pixels ).sum() / weight_grid.sum()
  Tsys_target = Trg + Tsky_target * ( frequency / haslam_frequency )**beta_target
  if print_info:
    Tsys = get_robust_mean_deviations( Tsys_target )
    print '... Tsys for target = %f +/- %f' % ( Tsys[ 0 ], max( Tsys[ 1 ], -Tsys[ 2 ] ) )
  
  # calculate flux correction
  corrections = get_robust_mean_deviations( Tsys_target / Tsys_cal )
  if print_info:
    print '... Tsys-based correction = %f +/- %f' % ( corrections[ 0 ],
        max( corrections[ 1 ], -corrections[ 2 ] ) )
  
  # correct visibilities or image (or not)
  if ( not apply_correction ):
    return corrections
  
  if print_info:
    print '... applying correction'
  if is_image( uvim ):
    cor_uvim = get_aips_file( uvim.disk, uvim.name, 'TCOR', -1, 'MA' )
    units = uvim.header.bunit
    call_aips_task( 'MATHS', indata = uvim, outdata = cor_uvim,
        opcode = 'POLY', cparm = [ 0, corrections[ 0 ], 0,0 ] )
    set_header_keyword( cor_uvim, 'bunit', units )
  else:
    sol_version = uvim.table_highver( 'SN' )
    if ( sol_version == 0 ):
      reference_antenna = 1
    else:
      reference_antenna = get_reference_antenna( uvim, sol_version )
    sn_version = generate_fluxscale_solutions( uvim, reference_antenna,
        corrections[ 0 ] )
    cal_uv = apply_solution_table( uvim, version = sn_version )
    uvim.zap_table( 'SN', sn_version )
    cor_uvim = get_aips_file( uvim.disk, uvim.name, 'TCOR', -1, 'UV' )
    cal_uv.rename( name = cor_uvim.name, klass = cor_uvim.klass,
        seq = cor_uvim.seq )
    for version in range( 1, 1 + sol_version ):
      if table_exists( uvim, 'SN', version ):
        call_aips_task( 'TACOP', indata = uvim, outdata = cor_uvim, inext = 'SN',
            invers = version, outvers = version, ncount = 1 )
  store_parameter( cor_uvim, 'flux_calibrator', cal_name )
  store_parameter( cor_uvim, 'tsys_correction', corrections[ 0 ] )
  store_parameter( cor_uvim, 'tsys_dcorrection', max( corrections[ 1 ], -corrections[ 2 ] ) )
  return cor_uvim

###############################################################################

def convert_lta_to_uvfits( lta_file_name, uvfits_file_name = None, target_list = [],
    stokes_list = [ 'RR', 'LL' ], keep_all_cals = True, max_scans = 128,
    scan_offset = 0, correction_list = [] ):
  # correction_list: [ [ scan#, old_source, new_source ] ]
  
  # hard-coded
  pad_char_offset = 127
  record_length = 336
  name_length = 12
  name_offset_start = -5 * 16
  name_offset_end = 14 * 16
  
  # check presence of files
  if ( not file_exists( lta_file_name ) ):
    raise error( 'LTA file %s does not exist' % ( lta_file_name ) )
  if ( uvfits_file_name is None ):
    fits_path = path.expandvars( '${FIT}/' )
    uvfits_file_name = fits_path + lta_file_name.split( '/' )[ -1 ] + '.UVFITS'
  if file_exists( uvfits_file_name ):
    raise error( 'UVFITS file %s already exist' % ( uvfits_file_name ) )
  
  # list contents of LTA file
  file_name = lta_file_name.split( '/' )[ -1 ]
  extension = file_name.split( '.' )[ -1 ]
  if ( extension in [ 'lta', 'ltb' ] ):
    extension = ''
  else:
    extension = '.' + extension
    file_name = file_name[ 0 : - len( extension ) ]
  ltab = file_name.split( '.' )[ -1 ]
  if ( not ltab in [ 'lta', 'ltb' ] ):
    raise error( 'LTA file name %s extension not recognized' % ( lta_file_name ) )
  file_name = file_name[ : file_name.index( '.' + ltab ) ]
  lsin_file_name = file_name + '.log'
  if file_exists( lsin_file_name ):
    remove_file( lsin_file_name )
  plan_file_name = file_name + '.plan'
  if file_exists( plan_file_name ):
    remove_file( plan_file_name )
  if ( sys.stdout.name == '<stdout>' ):
    result = system( 'listscan-1.17 %s' % ( lta_file_name ) )
  else:
    sys.stdout.flush()
    sys.stderr.flush()
    result = system( 'listscan-1.17 %s >> %s 2>&1' % ( lta_file_name, sys.stdout.name ) )
    sys.stdout.flush()
    sys.stderr.flush()
  if ( result != 0 ):
    raise error( 'LTA to UVFITS conversion failed for %s' % ( lta_file_name ) )
#  pdb.set_trace()
  
  # get overview of scans
  scan_list = []
  source_list = []
  lsin_file = file( lsin_file_name, mode = 'r' )
  lsin_list = []
  for line in lsin_file:
    lsin_list.append( line )
    words = line.split()
    if ( len( words ) == 0 ):
      continue
    if ( words[ 0 ] != 'Scan' ):
      continue
    scan_list.append( int( words[ 1 ] ) )
    source_list.append( words[ 2 ] )
  lsin_file.close()
  plan_file = file( plan_file_name, mode = 'rb' )
  plan_data = plan_file.read()
  plan_file.close()
  index_list = []
  index = -name_length
  pad_char = None
  for source in source_list:
    dindex = plan_data[ index + name_length : ].find( source )
    if ( index == -1 ):
      raise error( 'mismatch between listscan log and plan files' )
    index = index + name_length + dindex
    index_list.append( index )
    if ( pad_char == None ):
      if ( plan_data[ index + name_length - 1 ] in [ '\x00', ' ' ] ):
        pad_char = plan_data[ index + name_length - 1 ]
  # fix listscan files if needed
  if ( len( correction_list ) > 0 ):
    for [ scan_id, old_name, new_name ] in correction_list:
      if ( ( scan_id >= 0 ) and ( not scan_id in scan_list ) ):
        raise error( 'scan id %d not found in observation' % ( scan_id ) )
#        print( 'scan id %d not found in observation' % ( scan_id ) )
#        correction_list = []
#        break
      if ( scan_id >= 0 ):
        if ( source_list[ scan_list.index( scan_id ) ] != old_name ):
          raise error( 'source name %s does not correspond to scan id %d' %
              ( old_name, scan_id ) )
#          print( 'source name %s does not correspond to scan id %d' %
#              ( old_name, scan_id ) )
#          correction_list = []
#          break
      else:
        if ( not old_name in source_list ):
          raise error( 'source name %s not found in observation' % ( old_name ) )
#          print( 'source name %s not found in observation' % ( old_name ) )
#          correction_list = []
#          break
  if ( len( correction_list ) > 0 ):
    for [ scan_id, old_name, new_name ] in correction_list:
      if ( scan_id >= 0 ):
        for j, line in enumerate( lsin_list ):
          words = line.split()
          if ( len( words ) == 0 ):
            continue
          if ( words[ 0 ] != 'Scan' ):
            continue
          if ( int( words[ 1 ] ) != scan_id ):
            continue
          index = line.index( words[ 2 ] )
          lsin_list[ j ] = ( line[ : index ] + new_name.ljust( name_length ) + 
              line[ index + name_length : ] )
          break
        old_index = index_list[ scan_list.index( scan_id ) ]
        new_index = index_list[ source_list.index( new_name ) ]
        plan_data = ( plan_data[ : old_index + name_offset_start ] +
            plan_data[ new_index + name_offset_start : new_index + name_offset_end ] +
            plan_data[ old_index + name_offset_end : ] )
        source_list[ scan_id ] = new_name
      else:
        i = -1
        while ( old_name in source_list[ i + 1 : ] ):
          i = i + 1 + source_list[ i + 1 : ].index( old_name )
          for j, line in enumerate( lsin_list ):
            words = line.split()
            if ( len( words ) == 0 ):
              continue
            if ( words[ 0 ] != 'Scan' ):
              continue
            if ( int( words[ 1 ] ) != scan_list[ i ] ):
              continue
            index = line.index( words[ 2 ] )
            lsin_list[ j ] = ( line[ : index ] + new_name.ljust( name_length ) + 
                line[ index + name_length : ] )
            break
          old_index = index_list[ i ]
          new_index = index_list[ source_list.index( new_name ) ]
          plan_data = ( plan_data[ : old_index + name_offset_start ] +
              plan_data[ new_index + name_offset_start : new_index + name_offset_end ] +
              plan_data[ old_index + name_offset_end : ] )
          source_list[ scan_list[ i ] ] = new_name
    plan_file = file( plan_file_name, mode = 'wb' )
    plan_data = plan_file.write( plan_data )
    plan_file.close()
  
  # filter contents of LTA file
  lsout_file_name = lsin_file_name + '.sel'
  if file_exists( lsout_file_name ):
    remove_file( lsout_file_name )
  file_name = ( file_name + ltab[ -1 : ] ).upper() + '.UVFITS' + extension
  if file_exists( file_name ):
    remove_file( file_name )
  lsout_file = file( lsout_file_name, mode = 'w' )
  source_list = []
  scan_in_range = False
  for line in lsin_list:
    words = line.split()
    if ( len( words ) == 0 ):
      continue
    if ( words[ 0 ] == 'FITS' ):
      line = line[ : line.index( words[ 1 ] ) ]
      line = line + file_name + '\n'
    elif ( words[ 0 ] == 'STOKES' ):
      line = line[ : line.index( words[ 1 ] ) ]
      for stoke in stokes_list:
        if ( stoke in words[ 1 : ] ):
          line = line + ' ' + stoke
      line = line + '\n'
    elif ( words[ 0 ] == 'Scan' ):
      if ( not int( words[ 1 ] ) in range( scan_offset, scan_offset + max_scans ) ):
        continue
      scan_in_range = True
      keep_scan = False
      if ( ( not keep_scan ) and ( words[ 2 ] in [ '3C48', '3C147', '3C286' ] ) ):
        keep_scan = True
      if ( ( not keep_scan ) and ( words[ 2 ] in target_list ) ):
        keep_scan = True
      if ( ( not keep_scan ) and ( len( target_list ) == 0 ) ):
        if ( len( words[ 2 ] ) == 8 ):
          try:
            dummy = int( words[ 2 ][ 0 : 4 ] )
            dummy = int( words[ 2 ][ 5 : 8 ] )
          except:
            keep_scan = True
          else:
            if ( words[ 2 ][ 4 ] in [ '+', '-' ] ):
              if keep_all_cals:
                keep_scan = True
            else:
              keep_scan = True
        else:
          keep_scan = True
      if ( not keep_scan ):
        continue
      if ( not words[ 2 ] in source_list ):
        source_list.append( words[ 2 ] )
    lsout_file.write( line )
  lsout_file.close()
  if ( not scan_in_range ):
    remove_file( lsin_file_name )
    remove_file( lsout_file_name )
    return None
  
  # convert to UVFITS and clean up
  if ( len( source_list ) > 0 ):
    if ( sys.stdout.name == '<stdout>' ):
      result = system( 'gvfits-1.17 %s' % ( lsout_file_name ) )
    else:
      sys.stdout.flush()
      sys.stderr.flush()
      result = system( 'gvfits-1.17 %s >> %s 2>&1' % ( lsout_file_name, sys.stdout.name ) )
      sys.stdout.flush()
      sys.stderr.flush()
    if ( result != 0 ):
      raise error( 'LTA to UVFITS conversion failed for %s' % ( lta_file_name ) )
    move_file( file_name, uvfits_file_name )
    remove( 'gvfits.log' )
  remove_file( lsin_file_name )
  remove_file( lsout_file_name )
  remove_file( plan_file_name )
  
  return source_list

###############################################################################

def pre_calibrate_targets( uvfits_file_name, flags_file_name = None,
    target_names = [], bw_smearing = 0.7, ta_smearing = 0.35, 
    keep_channel_one = False, do_tb_sort = False, force_date = False,
    force_stokes = True, force_frequency = True, do_plots = True,
    bpcal_names = [ '3C48','3C147','3C286' ], bpcal_snr_power = 1.,
    bpcal_per_scan = True, initial_reference_antenna = 0, amplitude_limit = 20.,
    channel_range = None, bpcal_scan = None, bpcal_min_scan_length = 5. ):
  
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
  
  # get list of FLAGS files
  flags_file_names = []
  if ( not flags_file_name is None ):
    if ( type( flags_file_name ) == type( [] ) ):
      for file_name in flags_file_name:
        flags_file_names = flags_file_names + glob.glob( file_name )
    else:
      flags_file_names = glob.glob( flags_file_name )
    if ( len( flags_file_names ) == 0 ):
      raise error( 'FLAGS file(s) %s does not exist' % ( repr( flags_file_name ) ) )
  
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
  project_name = 'GMRT'
  for uvfits_file_name in uvfits_file_names:
    
    # read UV data
    multi_uv = get_aips_file( aips_disk, project_name, 'UV', -1, 'UV' )
    read_fits_uv( uvfits_file_name, multi_uv )
    if ( multi_uv.header.telescop != 'GMRT' ):
      raise error( 'UVFITS file %s is not a GMRT observation' % ( uvfits_file_name ) )
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
  
  # determine processing parameters
  integration_time = around( find_integration_time( multi_uv ), decimals = 2 )
  dfrequency = get_channel_width( multi_uv )
  frequency = mean( freqs_list[ 0 ] )
  channel_count = len( freqs_list[ 0 ] )
  if keep_channel_one:
    channel_increment = int( round( ( freqs_list[ 0 ][ 0 ] / dfrequency ) * 
        bw_smearing * ( ( 2. * 45. ) / 25.e3 ) ) )
  else:
    channel_increment = int( round( ( frequency / dfrequency ) * 
        bw_smearing * ( ( 2. * 45. ) / 25.e3 ) ) )
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
  central_channels = [ central_channel - 4, central_channel + 4 ]
  time_increment = max( 1, int( round( ( 1. / integration_time ) *
      ta_smearing * ( ( 2. * 45. ) / 25.e3 ) * ( 86400. / ( 2. * pi ) ) ) ) )
  if ( frequency > 140.e6 ) and ( frequency < 170.e6 ):
    uv_range = [ 0.5,1.e6 ]
    bp_dfluxes = [ 50., 30., 20. ]
    clip_level = 1000.
    extension = '_GMRT150' + extension
  elif ( frequency > 200.e6 ) and ( frequency < 280.e6 ):
    uv_range = [ 0.75,1.e6 ]
    bp_dfluxes = [ 25.,15.,10. ]
    clip_level = 200.
    extension = '_GMRT235' + extension
  elif ( frequency > 300.e6 ) and ( frequency < 350.e6 ):
    uv_range = [ 1.,1.e6 ]
    bp_dfluxes = [ 12.,7.,5. ]
    clip_level = 40.
    extension = '_GMRT325' + extension
  elif ( frequency > 580.e6 ) and ( frequency < 650.e6 ):
    uv_range = [ 2.,1.e6 ]
    bp_dfluxes = [ 5.,3.,2. ]
    clip_level = 8.
    extension = '_GMRT610' + extension
  elif ( frequency > 900.e6 ) and ( frequency < 1600.e6 ):
    uv_range = [ 4.,1.e6 ]
    bp_dfluxes = [ 5.,3.,2. ]
    clip_level = 4.
    extension = '_GMRT1200' + extension
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
  fix_uv = get_aips_file( aips_disk, project_name, 'UVCOR', -1, 'UV' )
  call_aips_task( 'UVFIX', indata = multi_uv, outdata = fix_uv )
  multi_uv.zap()
  fix_uv.rename( name = multi_uv.name, klass = multi_uv.klass, seq = multi_uv.seq )
  
  # get overview of raw UV content
  call_aips_task( 'INDXR', indata = multi_uv,
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
  
  # apply online flags
  for flags_file_name in flags_file_names:
    online_flags = read_gmrt_flag_file( multi_uv, flags_file_name,
        day_offset = -5.5 / 24., skip_flags = [], ignore = True )
    if ( len( online_flags ) > 0 ):
      for i in range( len( online_flags ) )[ : : -1 ]:
        if ( online_flags[ i ][ 'timerang' ][ 4 : 8 ] != [ 999,0,0,0 ] ):
          add_flags( multi_uv, flags = online_flags[ : i + 1 ], flag_version = 1 )
          break
  
  # create separate UV files for different sources
  if ( len( target_names ) == 0 ):
    target_names = get_source_names( multi_uv )
    for bpcal_name in bpcal_names:
      if ( bpcal_name in target_names ):
        target_names.remove( bpcal_name )
  source_names = bpcal_names + target_names
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
      
      # initial flagging
      flags_bpcal = [
          { 'bchan' : 1, 'echan' : 1 },
          { 'antennas' : [ 31,32 ] }]
      add_flags( main_bp_uv, flags = flags_bpcal, flag_version = 1 )
      
      # process per scan
      scan_list = find_scans( main_bp_uv )
      if ( ( not bpcal_per_scan ) and ( len( scan_list ) > 1 ) ):
        scan_list = [ [ scan_list[ 0 ][ 0 ], scan_list[ -1 ][ 1 ] ] ]
    except:
      continue
    
    for scan_id, scan_times in enumerate( scan_list ):
      
      if ( scan_times[ 1 ] - scan_times[ 0 ] < bpcal_min_scan_length / 1440. ):
        continue
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
        call_aips_task( 'UVCOP', indata = bp_cal_uv, outdata = temp_uv, flagver = 1 )
        
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
          call_aips_task( 'UVCOP', indata = temp_uv, outdata = temp_uv2, flagver = 1 )
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
  for i, bpcal_snr in enumerate( bpcal_snrs ):
    [ bpcal_name, scan_id, total_snr, bp_uv, bp_cal_uv, bp_ical_uv ] = bpcal_snr
    print '%02d: calibrator %s scan %d has an SNR of %d' % (
        i + 1, bpcal_name, scan_id + 1, int( total_snr ) )
  total_snrs = [ x[ 2 ] for x in bpcal_snrs ]
  if ( bpcal_scan is None ):
    bpcal_index = total_snrs.index( max( total_snrs ) )
  else:
    bpcal_index = bpcal_scan - 1
  [ bpcal_name, scan_id, total_snr, bp_uv, bp_cal_uv, bp_ical_uv ] = bpcal_snrs[ bpcal_index ]
  print '... selecting calibrator %s scan %d' % ( bpcal_name, scan_id + 1 )
  
  # make gain amplitude plots
  if do_plots:
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
  if do_plots:
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
  if do_plots:
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

def combine_usb_lsb( uvfits_file_name_1, uvfits_file_name_2,
    uvfits_file_name, overwrite = True, time_offset = 0 ):
  
  # basic checking
  if ( not file_exists( uvfits_file_name_1 ) ):
    raise error( 'UVFITS file %s does not exist' % ( uvfits_file_name_1 ) )
  if ( not file_exists( uvfits_file_name_2 ) ):
    raise error( 'UVFITS file %s does not exist' % ( uvfits_file_name_2 ) )
  if ( ( not overwrite ) and file_exists( uvfits_file_name ) ):
    raise error( 'UVFITS file %s already exists' % ( uvfits_file_name ) )
  
  # find free AIPS disk
  aips_disk = 0
  for i in range( 1, len( AIPS.disks ) ):
    if ( len( AIPSCat( i )[ i ] ) == 0 ):
      aips_disk = i
      break
  if ( aips_disk == 0 ):
    raise error( 'no free AIPS disk available' )
  
  # read and check UV data
  uv1 = get_aips_file( aips_disk, 'UV1', 'UV', -1, 'UV' )
  read_fits_uv( uvfits_file_name_1, uv1 )
  if ( uv1.header.telescop != 'GMRT' ):
    raise error( 'UVFITS file %s is not a GMRT observation' % ( uvfits_file_name_1 ) )
  uv2 = get_aips_file( aips_disk, 'UV2', 'UV', -1, 'UV' )
  read_fits_uv( uvfits_file_name_2, uv2 )
  if ( uv2.header.telescop != 'GMRT' ):
    raise error( 'UVFITS file %s is not a GMRT observation' % ( uvfits_file_name_2 ) )
  date1 = get_reference_date( uv1 ).replace( '-', '' )
  date2 = get_reference_date( uv2 ).replace( '-', '' )
  if ( date1 != date2 ):
    raise error( 'UVFITS files are from different observating days' )
  stokes1 = get_stokes( uv1 )
  stokes2 = get_stokes( uv2 )
  if ( stokes1 != stokes2 ):
    raise error( 'UVFITS files have different stokes' )
  freqs1 = array( get_frequency_list( uv1 ) )
  freqs2 = array( get_frequency_list( uv2 ) )
  dfreq1 = median( freqs1[ 1 : ] - freqs1[ : -1 ] )
  dfreq2 = median( freqs2[ 1 : ] - freqs2[ : -1 ] )
  if ( dfreq1 * dfreq2 > 0. ):
    raise error( 'UVFITS files are not USB and LSB' )
  if ( dfreq1 + dfreq2 > 1.e-3 * abs( dfreq1 ) ):
    raise error( 'USB and LSB have different channel widths' )
  if ( dfreq1 > 0. ):
    usb_uv = uv1
    lsb_uv = uv2
    usb_freqs = freqs1
    lsb_freqs = freqs2
    dfreq = dfreq1
  else:
    usb_uv = uv2
    lsb_uv = uv1
    usb_freqs = freqs2
    lsb_freqs = freqs1
    dfreq = dfreq2
  channel_count_usb = len( usb_freqs )
  channel_count_lsb = len( lsb_freqs )
  channel_gap = ( usb_freqs[ 0 ] - lsb_freqs[ 0 ] ) / dfreq
  if ( abs( channel_gap - round( channel_gap ) ) > 1.e-2 ):
    raise error( 'USB and LSB are on different frequency grids' )
  channel_gap = int( round( channel_gap ) ) - 1
  if ( channel_gap < 0 ):
    raise error( 'USB and LSB overlap' )
  if ( channel_gap > 10 ):
    raise error( 'USB and LSB are more than 10 channels apart' )
  
  # combine USB and LSB
  flip_lsb_uv = get_aips_file( aips_disk, 'LSB', 'FLIP', -1, 'UV' )
  call_aips_task( 'FLOPM', indata = lsb_uv, outdata = flip_lsb_uv )
  lsb_uv.zap()
  bloat_lsb_uv = get_aips_file( aips_disk, 'LSB', 'BLOAT', -1, 'UV' )
  call_aips_task( 'BLOAT', indata = flip_lsb_uv, outdata = bloat_lsb_uv,
      aparm = [ channel_count_lsb + channel_count_usb + channel_gap, 
      1, 1, channel_count_lsb, 0, 0 ] )
  flip_lsb_uv.zap()
  bloat_usb_uv = get_aips_file( aips_disk, 'USB', 'BLOAT', -1, 'UV' )
  call_aips_task( 'BLOAT', indata = usb_uv, outdata = bloat_usb_uv,
      aparm = [ channel_count_lsb + channel_count_usb + channel_gap,
      channel_count_lsb + channel_gap + 1, 1, channel_count_usb, 0, 0 ] )
  usb_uv.zap()
  if ( time_offset != 0 ):
    change_uv_reference_date( bloat_usb_uv, bloat_usb_uv, offset = 0 )
    change_uv_reference_date( bloat_lsb_uv, bloat_usb_uv, offset = time_offset )
  dbcon_uv = get_aips_file( aips_disk, 'USBLSB', 'DBCON', -1, 'UV' )
  call_aips_task( 'DBCON', indata = bloat_lsb_uv, in2data = bloat_usb_uv,
      outdata = dbcon_uv, doarray = 1 )
  bloat_lsb_uv.zap()
  bloat_usb_uv.zap()
  tavg_uv = time_average_uv_data( dbcon_uv, 1 )
  dbcon_uv.zap()
  
  # save result to UVFITS, clean-up and exit
  write_fits_uv( tavg_uv, uvfits_file_name )
  clear_aips_disk( aips_disk )
  return

###############################################################################

def combine_stokes( uvfits_file_name_1, uvfits_file_name_2, uvfits_file_name,
    overwrite = True, time_offset = 0 ):
  
  # basic checking
  if ( not file_exists( uvfits_file_name_1 ) ):
    raise error( 'UVFITS file %s does not exist' % ( uvfits_file_name_1 ) )
  if ( not file_exists( uvfits_file_name_2 ) ):
    raise error( 'UVFITS file %s does not exist' % ( uvfits_file_name_2 ) )
  if ( ( not overwrite ) and file_exists( uvfits_file_name ) ):
    raise error( 'UVFITS file %s already exists' % ( uvfits_file_name ) )
  
  # find free AIPS disk
  aips_disk = 0
  for i in range( 1, len( AIPS.disks ) ):
    if ( len( AIPSCat( i )[ i ] ) == 0 ):
      aips_disk = i
      break
  if ( aips_disk == 0 ):
    raise error( 'no free AIPS disk available' )
  
  # read and check UV data
  uv1 = get_aips_file( aips_disk, 'UV1', 'UV', -1, 'UV' )
  read_fits_uv( uvfits_file_name_1, uv1 )
  if ( uv1.header.telescop != 'GMRT' ):
    raise error( 'UVFITS file %s is not a GMRT observation' % ( uvfits_file_name_1 ) )
  uv2 = get_aips_file( aips_disk, 'UV2', 'UV', -1, 'UV' )
  read_fits_uv( uvfits_file_name_2, uv2 )
  if ( uv2.header.telescop != 'GMRT' ):
    raise error( 'UVFITS file %s is not a GMRT observation' % ( uvfits_file_name_2 ) )
  date1 = get_reference_date( uv1 ).replace( '-', '' )
  date2 = get_reference_date( uv2 ).replace( '-', '' )
  if ( date1 != date2 ):
    raise error( 'UVFITS files are from different observating days' )
  freqs1 = array( get_frequency_list( uv1 ) )
  freqs2 = array( get_frequency_list( uv2 ) )
  dfreq1 = median( freqs1[ 1 : ] - freqs1[ : -1 ] )
  dfreq2 = median( freqs2[ 1 : ] - freqs2[ : -1 ] )
  if ( abs( dfreq1 - dfreq2 ) > 1.e-3 * abs( dfreq1 ) ):
    raise error( 'USB and LSB have different channel widths' )
  if ( sometrue( abs( freqs1 - freqs2 ) > 1.e-3 * abs( dfreq1 ) ) ):
    raise error( 'UVFITS files have different frequency setup' )
  stokes1 = get_stokes( uv1 )
  stokes2 = get_stokes( uv2 )
  if ( stokes1 != stokes2 ):
    raise error( 'UVFITS files have different stokes' )
  
  # combine USB and LSB
  if ( time_offset != 0 ):
    change_uv_reference_date( uv1, uv1, offset = 0 )
    change_uv_reference_date( uv2, uv1, offset = time_offset )
  dbcon_uv = get_aips_file( aips_disk, 'UV12', 'DBCON', -1, 'UV' )
  call_aips_task( 'DBCON', indata = uv1, in2data = uv2,
      outdata = dbcon_uv, doarray = 1 )
  uv1.zap()
  uv2.zap()
  tavg_uv = time_average_uv_data( dbcon_uv, 1 )
  dbcon_uv.zap()
  
  # save result to UVFITS, clean-up and exit
  write_fits_uv( tavg_uv, uvfits_file_name )
  clear_aips_disk( aips_disk )
  return

###############################################################################

def split_scans( uvfits_file_name, overwrite = True, max_gap = 60. ):
  
  # basic checking
  if ( not file_exists( uvfits_file_name ) ):
    raise error( 'UVFITS file %s does not exist' % ( uvfits_file_name_1 ) )
  
  # find free AIPS disk
  aips_disk = 0
  for i in range( 1, len( AIPS.disks ) ):
    if ( len( AIPSCat( i )[ i ] ) == 0 ):
      aips_disk = i
      break
  if ( aips_disk == 0 ):
    raise error( 'no free AIPS disk available' )
  
  # read UV data
  uv = get_aips_file( aips_disk, 'UV', 'UV', -1, 'UV' )
  read_fits_uv( uvfits_file_name, uv )
  scan_list = find_scans( uv, max_gap = max_gap )
  for scan_id, scan_times in enumerate( scan_list ):
    scan_uvfits_file_name = uvfits_file_name + '.%03d' % ( scan_id + 1 )
    if ( ( not overwrite ) and file_exists( scan_uvfits_file_name ) ):
      clear_aips_disk( aips_disk )
      raise error( 'UVFITS file %s already exists' % ( scan_uvfits_file_name ) )
  scan_uvfits_file_names = []
  for scan_id, scan_times in enumerate( scan_list ):
    scan_uv = get_aips_file( aips_disk, 'SCAN%03d' % ( scan_id + 1 ), 'UV', -1, 'UV' )
    time_range = timetime_to_dhmsdhms( scan_times )
    call_aips_task( 'SPLAT', indata = uv, docalib = -1, doband = -1,
        flagver = 0, blver = -1, douvcomp = 0, aparm = [ 3,0,1,0,0,0,0 ],
        stokes = '', outdata = scan_uv, timerang = time_range )
    scan_uvfits_file_name = uvfits_file_name + '.%03d' % ( scan_id + 1 )
    write_fits_uv( scan_uv, scan_uvfits_file_name )
    scan_uvfits_file_names.append( scan_uvfits_file_name )
  
  # clean-up
  clear_aips_disk( aips_disk )  
  return scan_uvfits_file_names

###############################################################################


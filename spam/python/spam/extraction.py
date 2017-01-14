###############################################################################

# import Python modules
from sys import *
from os import *
from math import *
import random

# import 3rd party modules
from numpy import *
import pyfits as pf

# import user modules
from acalc import *
from files import *
from aips import *
from sphere import *
from parameter import *
from skymodel import *
from error import *

###############################################################################

def extract_source_catalog_table( im, catalog_version = 0, sigma_min = 5., 
    iterations = 5, max_gaussians_per_island = 4, rms_boxsize = 256,
    rms_stepsize = 32, make_residual_map = True, keep_rms_map = False ):
  
  # handle inputs
  if ( iterations > 10 ):
    raise error( 'number of iterations may not exceed 10' )
  if ( max_gaussians_per_island > 4 ):
    raise error( 'number of gaussians per island may not exceed 4' )
  if make_residual_map:
    res_image = get_aips_file( im.disk, im.name, 'RESID', - 1, 'MA' )
    doresid = 1
  else:
    doresid = 0
  if ( max_gaussians_per_island == 1 ):
    doall = 0
  else:
    doall = 1
  dowidth = []
  for i in range( max_gaussians_per_island ):
    dowidth.append( [ 1, 1, 1 ] )
  for i in range( max_gaussians_per_island, 4 ):
    dowidth.append( [ 0, 0, 0 ] )
  
  # remove blanked borders of image
#  image = crop_image( im )
#  [ flux_max, pos_max ] = get_image_maximum( image )
#  if ( flux_min >= flux_max ):
#    raise error( 'peak flux in image does not exceed flux threshold' )
  image = im
  
  # calculate limits for iterations
  cparm = [ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. ]
  for i in range( iterations ):
    cparm[ i ] =  2.**float( iterations - i - 1 ) * sigma_min
  
  # calculate rms image
  rms_image = get_aips_file( im.disk, im.name, 'RMS', - 1, 'MA' )
  call_aips_task( 'RMSD', indata = image, outdata = rms_image,
     imsize = [ rms_boxsize, rms_boxsize ], xinc = rms_stepsize,
     yinc = rms_stepsize )
  
  # do source extraction
  call_aips_task( 'SAD', indata = image, invers = 0, outvers = - 1,
      in2data = rms_image, docrt = 0, ngauss = 40000, cparm = cparm,
      icut = 0., bwsmear = 0, sort = 'S', doall = doall, dowidth = dowidth,
      gain = 0.1, dparm = [ sigma_min, 0, 0, 0, 0, 0, 0, 0, 2, 0 ],  
      doresid = doresid, outdata = res_image )
  
  # copy catalog table to original image
#  call_aips_task( 'TACOP', indata = image, inext = 'MF', invers = 0,
#      ncount = 1, outdata = im, outvers = catalog_version )
#  image.zap()
  if ( not keep_rms_map ):
    rms_image.zap()
  
  if make_residual_map:
    return res_image
  else:
    return

###############################################################################

# format example
#
# Title: AIPS MODEL FIT COMPONENTS TABLE
# Created by      SAD   on 07-DEC-2007 14:19:27
# Last written by SAD   on 07-DEC-2007 14:19:27
# Ncol  39  Nrow     109    Sort cols:
#     Table has     7 keyword-value pairs:
#    REVISION =            3
#    DEPTH1   =            1
#    DEPTH2   =            1
#    DEPTH3   =            1
#    DEPTH4   =            1
#    DEPTH5   =            1
#    REALRMS  =  8.7719038128853D-02
#    Table can be written as a FITS ASCII table
# 
# COL. NO.      1          2          3          4           5          6
#      ROW  PLANE      PEAK INT   I FLUX     DELTAX     DELTAY      MAJOR AX
#   NUMBER             JY/BEAM    JY         DEGREES    DEGREES     DEGREES
#        1  1.000E+00  4.439E+01  4.480E+01  2.695E+00  -1.503E+00  2.240E-02
# 
# COL. NO.      7          8          9          10         11         12
#      ROW  MINOR AX   POSANGLE   Q FLUX     U FLUX     V FLUX     ERR PEAK
#   NUMBER  DEGREES    DEGREES    JY         JY         JY         JY/BEAM
#        1  2.228E-02  1.698E+01  0.000E+00  0.000E+00  0.000E+00  8.772E-02
# 
# COL. NO.      13         14         15         16         17         18
#      ROW  ERR FLUX   ERR DLTX   ERR DLTY   ERR MAJA   ERR MINA   ERR PA
#   NUMBER  JY         DEGREES    DEGREES    DEGREES    DEGREES    DEGREES
#        1  1.529E-01  1.868E-05  1.878E-05  4.425E-05  4.397E-05  1.235E+01
# 
# COL. NO.      19         20         21         22         23         24
#      ROW  ERR QFLX   ERR UFLX   ERR VFLX   TYPE MOD   D0 MAJOR   D0 MINOR
#   NUMBER  JY         JY         JY                    DEGREES    DEGREES
#        1  0.000E+00  0.000E+00  0.000E+00  1.000E+00  2.786E-03  1.142E-03
# 
# COL. NO.      25         26         27         28         29         30
#      ROW  D0 POSAN   D- MAJOR   D- MINOR   D- POSAN   D+ MAJOR   D+ MINOR
#   NUMBER  DEGREES    DEGREES    DEGREES    DEGREES    DEGREES    DEGREES
#        1  1.698E+01  2.405E-03  0.000E+00  5.718E+00  3.122E-03  1.809E-03
# 
# COL. NO.      31         32         33         34         35         36
#      ROW  D+ POSAN   RES RMS    RES PEAK   RES FLUX   CENTER X   CENTER Y
#   NUMBER  DEGREES    JY/BEAM    JY/BEAM    JY         PIXELS     PIXELS
#        1  2.933E+01  1.130E-01  3.329E-01  1.960E-01  3.214E+03  3.605E+03
# 
# COL. NO.          37             38             39
#      ROW      MAJ AXIS       MIN AXIS       PIXEL PA
#   NUMBER      PIXELS         PIXELS         DEGREES
#        1      7.330E+00      7.282E+00      1.807E+01

def read_aips_model_fit_table( image, version = 0 ):
  
  if ( not table_exists( image, 'MF', version ) ):
    raise error( 'source catalog table does not exist' )
  wiz_im = wizardry( image )
  cat_table = wiz_im.table( 'MF', version )
  
  cat_list = [ [ 'ID', 'RA', 'DEC', 'PEAK', 'FLUX', 'MAJ', 'MIN', 'PA',
      'PEAK_RES', 'RA_ERROR', 'DEC_ERROR', 'PEAK_ERROR', 'FLUX_ERROR',
      'MAJ_ERROR', 'MIN_ERROR', 'PA_ERROR' ] ]
  id = 0
  for row in cat_table:
    id = id + 1
    [ ra, dec ] = calculate_source_radec( image, [ row.center_x, row.center_y ] )
    peak = row.peak_int
    flux = row.i_flux
    maj = row.major_ax
    min = row.minor_ax
    pa = row.posangle
    peak_res = row.res_peak
    ra_error = row.err_dltx
    dec_error = row.err_dlty
    peak_error = row.err_peak
    flux_error = row.err_flux
    maj_error = row.err_maja
    min_error = row.err_mina
    pa_error = row.err_pa
    cat_list.append( [ id, ra, dec, peak, flux, maj, min, pa, peak_res,
        ra_error, dec_error, peak_error, flux_error, maj_error, min_error, pa_error ] )
  
  return cat_list


###############################################################################

# format example
#             RA           Dec          Peak    Flux    IRMS  Fit Maj Fit min   PA    res. RMS res Peak  PixX    PixY   Field
#    1  12 28 46.5232  -21 26 49.634  1726.27  2047.85 210.435  94.902  79.999 -122.4    0.047    0.099  3233.7  1914.4 1300-208
#                7.20           6.29   217.24   257.71          12.913   9.310   27.5 

def read_obit_source_catalog( cat_file_name ):

  if ( not file_exists( cat_file_name ) ):
    raise error( 'source catalog %s was not found' % cat_file_name )

  cat_list = [ [ 'ID', 'RA', 'DEC', 'PEAK', 'FLUX', 'MAJ', 'MIN', 'PA', 'PEAK_RES',
      'RA_ERROR', 'DEC_ERROR', 'PEAK_ERROR', 'FLUX_ERROR', 'MAJ_ERROR', 'MIN_ERROR', 'PA_ERROR' ] ]
  cat_file = file( cat_file_name, mode = 'r' )

  header_found = False
  error_line = False
  for line in cat_file:
    words = [ word.strip() for word in line.split() ]
    if ( len( words ) == 0 ):
      continue
    if ( not header_found ):
      if ( words[ 0 ] == 'RA' ):
        if ( line != '             RA           Dec          Peak    Flux    IRMS  ' + 
            'Fit Maj Fit min   PA    res. RMS res Peak  PixX    PixY   Field\n' ):
          raise error( 'unexpected row format' )
        header_found = True
    elif ( not error_line ):
      id = int( words[ 0 ] )
      [ ra, dec ] = hmsdms_to_degdeg( [ float( words[ 1 ] ), float( words[ 2 ] ), float( words[ 3 ] ),
          float( words[ 4 ] ), float( words[ 5 ] ), float( words[ 6 ] ) ] )
      peak = float( words[ 7 ] ) / 1.e3
      flux = float( words[ 8 ] ) / 1.e3
      maj = float( words[ 10 ] )
      min = float( words[ 11 ] )
      pa = float( words[ 12 ] )
      peak_res = float( words[ 14 ] ) * peak
      error_line = True
    else:
      ra_error = float( words[ 0 ] ) / 3600.
      dec_error = float( words[ 1 ] ) / 3600.
      peak_error = float( words[ 2 ] ) / 1.e3
      flux_error = float( words[ 3 ] ) / 1.e3
      maj_error = float( words[ 4 ] )
      min_error = float( words[ 5 ] )
      pa_error = float( words[ 6 ] )
      cat_list.append( [ id, ra, dec, peak, flux, maj, min, pa, peak_res, 
          ra_error, dec_error, peak_error, flux_error, maj_error, min_error, pa_error ] )
      error_line = False

  cat_file.close()

  return cat_list

###############################################################################

# format example
# src# flag code tot_Jy err CenPkJy err MaxPkJy err CenRA err CenDec err MaxRA err MaxDec err bmaj_asec_fw err bmin_asec_fw err bpa_deg err deconv_bmaj_bmin_bpa_asec_fw  &errors rms_isl num_gaus
# fmt 76 "(2(i5,1x),a4,1x,6(1Pe11.3,1x),8(0Pf13.9,1x),12(0Pf10.5,1x),1(1Pe11.3,1x),i3)"
#    1     0 S      4.412E+01   1.387E-01   4.401E+01   8.014E-02   4.401E+01   8.014E-02 197.913289400   0.001773634 -22.278580020   0.026014849 197.913289400   0.001773634 -22.278580020   0.026014849   80.55641    0.14753   79.65808    0.14426   10.85062    6.54240    9.45851    0.14753    0.00000    0.14426   10.85062    6.54240   1.823E-01    1

# gaul# island# flag tot_Jy err peak_Jy err   RA err DEC err  xpos_pix err ypos_pix err bmaj_asec_fw err bmin_asec_fw err bpa_deg err deconv_bmaj_bmin_bpa_asec_fw &errors src_rms src_av isl_rms isl_av spin e_spin src#  blc1  blc2  trc1 trc2





def read_bdsm_source_catalog( cat_file_name, use_centroid = True ):
  if ( not file_exists( cat_file_name ) ):
    raise error( 'source catalog %s was not found' % cat_file_name )
  if use_centroid:
    col_list = [ 'src#', 'flag', 'code', 'CenRA', 'CenDec', 'CenPkJy', 'tot_Jy',
        'bmaj_asec_fw', 'bmin_asec_fw', 'bpa_deg' ]
  else:
    col_list = [ 'src#', 'flag', 'code', 'MaxRA', 'MaxDec', 'MaxPkJy', 'tot_Jy',
        'bmaj_asec_fw', 'bmin_asec_fw', 'bpa_deg' ]
  cat_list = [ [ 'ID', 'RA', 'DEC', 'PEAK', 'FLUX', 'MAJ', 'MIN', 'PA',
      'RA_ERROR', 'DEC_ERROR', 'PEAK_ERROR', 'FLUX_ERROR', 'MAJ_ERROR', 'MIN_ERROR', 'PA_ERROR' ] ]
  cat_file = file( cat_file_name, mode = 'r' )
  header_found = False
  for line in cat_file:
    words = [ word.strip() for word in line.split() ]
    if ( len( words ) == 0 ):
      continue
    if ( not header_found ):
      if ( words[ 0 ] == col_list[ 0 ] ):
        try:
          index_list = [ words.index( col ) for col in col_list ]
        except:
          raise error( 'unexpected row format' )
        index_count = len( index_list )
        err_list = [ words[ index + 1 ] == 'err' for index in index_list[ 3 : index_count ] ]
        if ( False in err_list ):
          raise error( 'unexpected row format' )
        for index in index_list[ 3 : index_count ]:
          index_list.append( index + 1 )
      elif ( words[ 0 ] == 'fmt' ):
        header_found = True
    else:
      if ( int( words[ index_list[ 1 ] ] ) > 0 ): # flag
        continue
#      if ( words[ index_list[ 2 ] ] != 'S' ): # code
#        continue
      id = int( words[ index_list[ 0 ] ] )
      cat_list.append( [ id ] + [ float( words[ index ] ) for index in index_list[ 3 : ] ] )
  cat_file.close()

  return cat_list

###############################################################################

# Gaussian list made by bdsm.f
# Image_name BOOTES_GMRT_150_PBCOR
# Image_size_x  2737
# Image_size_y  2737
# Island_list_name BOOTES_GMRT_150_PBCOR.524
# Number_of_islands  637
# Number_of_sources  672
# Number_of_gaussians  742
# Max_gaussians_per_island  6
# RMS_map map
# Sigma   0.00188698388
# Detect_threshold   5.
# src# isl# flag code tot_Jy err CenPkJy err MaxPkJy err CenRA err CenDec err MaxRA err MaxDec err bmaj_asec_fw err bmin_asec_fw err bpa_deg err deconv_bmaj_bmin_bpa_asec_fw  &errors rms_isl num_gaus blc1 blc2 trc1 trc2
# fmt 85 "(3(i5,1x),a4,1x,6(1Pe11.3,1x),8(0Pf13.9,1x),12(0Pf10.5,1x),1(1Pe11.3,1x),i3,4(1x,i5))"
#     1     1     0 C      4.192E+00   9.394E-03   4.172E+00   5.408E-03   4.172E+00   5.408E-03 219.704046843   0.000440722  33.837557532   0.004750058 219.704046843   0.000440722  33.837557532   0.004750058   26.04024    0.03666   22.06899    0.02634   81.86934    0.29090    3.93741    0.03666    0.00000    0.02634  125.97794    0.29090   1.444E-01   1   337  1050   393  1068

def read_raw_bdsm_source_catalog( cat_file_name, drop_flagged = True, 
    print_info = True ):
  if ( not file_exists( cat_file_name ) ):
    raise error( 'source catalog %s was not found' % cat_file_name )
  cat_file = file( cat_file_name, mode = 'r' )
  header_found = False
  cat_list = []
  dropped_count = 0
  for line in cat_file:
    words = [ word.strip() for word in line.split() ]
    if ( len( words ) == 0 ):
      continue
    if ( not header_found ):
      if ( 'src#' in words ):
        flag_index = words.index( 'flag' )
        deconv_index = words.index( 'deconv_bmaj_bmin_bpa_asec_fw' )
        words = words[ 0 : deconv_index ] + [ 'deconv_bmaj_asec_fw', 'err', 
            'deconv_bmin_asec_fw', 'err', 'deconv_bpa_deg', 'err' ] + words[ deconv_index + 2 : ]
        for i in range( len( words ) ):
          if ( words[ i ] == 'err' ):
            words[ i ] = words[ i - 1 ] + '_err'
        cat_list.append( words )
      elif ( words[ 0 ] == 'fmt' ):
        header_found = True
    else:
      if ( len( words ) != len( cat_list[ 0 ] ) ):
         raise error( 'data row length deviates from header row length' )
      if ( drop_flagged and ( int( words[ flag_index ] ) > 0 ) ):
        dropped_count = dropped_count + 1
        continue
      values = []
      for word in words:
        if ( '.' in word ):
          try:
            value = float( word )
          except:
            value = word
        else:
          try:
            value = int( word )
          except:
            value = word
        values.append( value )
      cat_list.append( [ value for value in values ] )
  cat_file.close()
  if ( print_info and ( dropped_count > 0 ) ):
    print '... dropped %s flagged source components' % ( repr( dropped_count ) )
  return cat_list

###############################################################################

# Gaussian list made by bdsm.f
# Image_name A2256.MOSAIC
# Image_size_x  6000
# Image_size_y  6000
# Island_list_name A2256.MOSAIC.5_3_4
# Number_of_islands  1280
# Number_of_gaussians  2074
# Max_gaussians_per_island  105
# RMS_map map
# Sigma   0.00275386184
# Detect_threshold   5.
# gaul# island# flag tot_Jy err peak_Jy err   RA err DEC err  xpos_pix err ypos_pix err bmaj_asec_fw err bmin_asec_fw err bpa_deg err deconv_bmaj_bmin_bpa_asec_fw &errors src_rms src_av isl_rms isl_av spin e_spin src#  blc1  blc2  trc1 trc2  dumr1 dumr2 dumr3 dumr4 dumr5 dumr6
# fmt 107 "(3(i5,1x),4(1Pe11.3,1x),4(0Pf13.9,1x),16(0Pf11.5,1x),4(1Pe11.3,1x),0Pf7.3,1x,0Pf7.3,5(1x,i5),6(1x,1Pe11.3))"

def read_raw_bdsm_gaussian_catalog( cat_file_name, drop_flagged = True,
    fix_position_errors = False, print_info = True ):
  if ( not file_exists( cat_file_name ) ):
    raise error( 'gaussian catalog %s was not found' % cat_file_name )
  cat_file = file( cat_file_name, mode = 'r' )
  header_found = False
  cat_list = []
  dropped_count = 0
  for line in cat_file:
    words = [ word.strip() for word in line.split() ]
    if ( len( words ) == 0 ):
      continue
    if ( not header_found ):
      if ( 'gaul#' in words ):
        flag_index = words.index( 'flag' )
        deconv_index = words.index( 'deconv_bmaj_bmin_bpa_asec_fw' )
        words = ( words[ 0 : deconv_index ] + [ 'deconv_bmaj_asec_fw', 'err', 
            'deconv_bmin_asec_fw', 'err', 'deconv_bpa_deg', 'err' ] + 
            words[ deconv_index + 2 : ] )
        for i in range( len( words ) ):
          if ( words[ i ] == 'err' ):
            words[ i ] = words[ i - 1 ] + '_err'
        cat_list.append( words )
      elif ( words[ 0 ] == 'fmt' ):
        header_found = True
    else:
      if ( len( words ) != len( cat_list[ 0 ] ) ):
         raise error( 'data row length deviates from header row length' )
      if ( drop_flagged and ( int( words[ flag_index ] ) > 0 ) ):
        dropped_count = dropped_count + 1
        continue
      values = []
      for word in words:
        if ( '.' in word ):
          try:
            value = float( word )
          except:
            value = word
        else:
          try:
            value = int( word )
          except:
            value = word
        values.append( value )
      cat_list.append( [ value for value in values ] )
  cat_file.close()
  if ( print_info and ( dropped_count > 0 ) ):
    print '... dropped %s flagged gaussian components' % ( repr( dropped_count ) )
  if fix_position_errors:
    indices = [ cat_list[ 0 ].index( x ) for x in 
        [ 'RA_err', 'DEC_err', 'xpos_pix_err', 'ypos_pix_err', 'bpa_deg' ] ]
    for s in cat_list[ 1 : ]:
      [ dra, ddec, dx, dy, bpa ] = [ s[ x ] for x in indices ]
      dradx = dra / dx
      ddecdy = ddec / dy
      bpar = radians( bpa )
      new_dx = sqrt( dx**2 * sin( bpar )**2 + dy**2 * cos( bpar )**2 )
      new_dy = sqrt( dx**2 * cos( bpar )**2 + dy**2 * sin( bpar )**2 )
      new_dra = dradx * new_dx
      new_ddec = dradx * new_dy
      new_values = [ new_dra, new_ddec, new_dx, new_dy ]
      for i in range( len( new_values ) ):
        s[ indices[ i ] ] = new_values[ i ]
  return cat_list

###############################################################################

# gaul# island# flag tot_Jy err peak_Jy err   RA err DEC err  xpos_pix err ypos_pix err bmaj_asec_fw err bmin_asec_fw err bpa_deg err deconv_bmaj_bmin_bpa_asec_fw &errors src_rms src_av isl_rms isl_av spin e_spin src#  blc1  blc2  trc1 trc2

def read_bdsm_gaussian_catalog( cat_file_name, fix_position_errors = False ):
  if ( not file_exists( cat_file_name ) ):
    raise error( 'source catalog %s was not found' % cat_file_name )
  col_list = [ 'gaul#', 'flag', 'RA', 'DEC', 'peak_Jy', 'tot_Jy',
      'bmaj_asec_fw', 'bmin_asec_fw', 'bpa_deg', 'xpos_pix', 'ypos_pix' ]
  cat_list = [ [ 'ID', 'RA', 'DEC', 'PEAK', 'FLUX', 'MAJ', 'MIN', 'PA', 'X', 'Y',
      'RA_ERROR', 'DEC_ERROR', 'PEAK_ERROR', 'FLUX_ERROR', 'MAJ_ERROR', 'MIN_ERROR', 'PA_ERROR',
      'X_ERROR', 'Y_ERROR' ] ]
  cat_file = file( cat_file_name, mode = 'r' )
  header_found = False
  for line in cat_file:
    words = [ word.strip() for word in line.split() ]
    if ( len( words ) == 0 ):
      continue
    if ( not header_found ):
      if ( words[ 0 ] == col_list[ 0 ] ): # gaul#
        try:
          index_list = [ words.index( col ) for col in col_list ]
        except:
          raise error( 'unexpected row format' )
        index_count = len( index_list )
        err_list = [ words[ index + 1 ] == 'err' for index in index_list[ 2 : index_count ] ]
        if ( False in err_list ):
          raise error( 'unexpected row format' )
        for index in index_list[ 2 : index_count ]:
          index_list.append( index + 1 )
      elif ( words[ 0 ] == 'fmt' ):
        header_found = True
    else:
      if ( int( words[ index_list[ 1 ] ] ) > 0 ): # flag
        continue
      id = int( words[ index_list[ 0 ] ] )
      cat_list.append( [ id ] + [ float( words[ index ] ) for index in index_list[ 2 : ] ] )
  cat_file.close()
  if fix_position_errors:
    indices = [ cat_list[ 0 ].index( x ) for x in 
        [ 'RA_ERROR', 'DEC_ERROR', 'X_ERROR', 'Y_ERROR', 'PA' ] ]
    for s in cat_list[ 1 : ]:
      [ dra, ddec, dx, dy, bpa ] = [ s[ x ] for x in indices ]
      dradx = dra / dx
      ddecdy = ddec / dy
      bpar = radians( bpa )
      new_dx = sqrt( dx**2 * sin( bpar )**2 + dy**2 * cos( bpar )**2 )
      new_dy = sqrt( dx**2 * cos( bpar )**2 + dy**2 * sin( bpar )**2 )
      new_dra = dradx * new_dx
      new_ddec = dradx * new_dy
      new_values = [ new_dra, new_ddec, new_dx, new_dy ]
      for i in range( len( new_values ) ):
        s[ indices[ i ] ] = new_values[ i ]
  return cat_list

###############################################################################

def associate_sources( radec_list_1, radec_list_2, assoc_radius = 100.,
    allow_multiple = False, pre_filter = 5. ):
# assoc_radius in arcsec
  
  assoc_list = []
  if ( ( len( radec_list_1 ) == 0 ) or ( len( radec_list_2 ) == 0 ) ):
    return assoc_list
  
  radius_list = []
  for i in range( len( radec_list_1 ) ):
    radec_i = radec_list_1[ i ]
    if ( not pre_filter is None ):
      [ ra, dec ] = radec_i
      radius = pre_filter * assoc_radius / 3600.
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
    for j in range( len( radec_list_2 ) ):
      radec_j = radec_list_2[ j ]
      if ( pre_filter is None ):
        [ r, phi ] = calculate_angular_separation( radec_i, radec_j )
        if ( 3600. * r < assoc_radius ):
          radius_list.append( [ i, j, 3600. * r ] )
        continue
      [ ra, dec ] = radec_j
      if ( ( dec < dec_min ) or ( dec > dec_max ) ):
        continue
      if ( ra_range_contains_zero and ( ra > ra_min ) and ( ra < ra_max ) ):
        continue
      if ( ( not ra_range_contains_zero ) and ( ( ra < ra_min ) or ( ra > ra_max ) ) ):
        continue
      [ r, phi ] = calculate_angular_separation( radec_i, radec_j )
      if ( 3600. * r < assoc_radius ):
        radius_list.append( [ i, j, 3600. * r ] )
  i_array = array( [ x[ 0 ] for x in radius_list ] )
  j_array = array( [ x[ 1 ] for x in radius_list ] )
  r_array = array( [ x[ 2 ] for x in radius_list ] )
  
  for i in range( len( radec_list_1 ) ):
    j_matches = aget( j_array, awhere( i_array == i ) )
    if ( len( j_matches ) == 0 ):
      assoc_list.append( [ i, -1 ] )
      continue
    if allow_multiple:
      for j in j_matches.ravel():
        assoc_list.append( [ i, j ] )
    else:
      r_matches = aget( r_array, awhere( i_array == i ) )
      j_min = aget( j_matches, awhere( r_matches == r_matches.min() ) )[ 0 ]
      i_matches = aget( i_array, awhere( j_array == j_min ) )
      r_matches = aget( r_array, awhere( j_array == j_min ) )
      i_min = aget( i_matches, awhere( r_matches == r_matches.min() ) )[ 0 ]
      if ( i_min == i ):
        assoc_list.append( [ i, j_min ] )
      else:
        assoc_list.append( [ i, -1 ] )
  assoc_array = array( assoc_list )
  for j in range( len( radec_list_2 ) ):
    if ( len( awhere( assoc_array[ : , 1 ] == j ) ) == 0 ):
      assoc_list.append( [ -1, j ] )
  
  return assoc_list

###############################################################################

def match_source_catalogs( cat_list_1, cat_list_2, assoc_radius = 100.,
    allow_multiple = False, pre_filter = 2., difference = False ):
  ra_index_1 = [ 'RA' in word for word in cat_list_1[ 0 ] ].index( True )
  dec_index_1 = [ 'DE' in word for word in cat_list_1[ 0 ] ].index( True )
  ra_index_2 = [ 'RA' in word for word in cat_list_2[ 0 ] ].index( True )
  dec_index_2 = [ 'DE' in word for word in cat_list_2[ 0 ] ].index( True )
  radec_list_1 = [ [ s[ ra_index_1 ], s[ dec_index_1 ] ] for s in cat_list_1[ 1 : ] ]
  radec_list_2 = [ [ s[ ra_index_2 ], s[ dec_index_2 ] ] for s in cat_list_2[ 1 : ] ]
  assoc_list = associate_sources( radec_list_1, radec_list_2,
      assoc_radius = assoc_radius, allow_multiple = allow_multiple,
      pre_filter = pre_filter )
  sub_cat_list_1 = cat_list_1[ 0 : 1 ]
  sub_cat_list_2 = cat_list_2[ 0 : 1 ]
  if difference:
    for assoc in assoc_list:
      if ( assoc[ 1 ] == -1 ):
        sub_cat_list_1.append( [ x for x in cat_list_1[ assoc[ 0 ] + 1 ] ] )
      elif ( assoc[ 0 ] == -1 ):
        sub_cat_list_2.append( [ x for x in cat_list_2[ assoc[ 1 ] + 1 ] ] )
  else:
    for assoc in assoc_list:
      if ( not -1 in assoc ):
        sub_cat_list_1.append( [ x for x in cat_list_1[ assoc[ 0 ] + 1 ] ] )
        sub_cat_list_2.append( [ x for x in cat_list_2[ assoc[ 1 ] + 1 ] ] )
  return [ sub_cat_list_1, sub_cat_list_2 ]

###############################################################################

def add_point_sources_to_image( image, flux_min, flux_list, radec_list = [], check_radius = 10,
    add_sources = True, print_info = False ):

  if ( len( radec_list ) > 0 ):
    if ( len( radec_list ) != len( flux_list ) ):
      raise error( 'lengths of flux list and RA DEC list do not match' )

  # make list of valid source pixel positions
  new_radec_list = []
  pos_list = []
  pixel_map = get_image_pixels( image )
  image_size = list( pixel_map.shape )
  sub_image_mask = zeros( ( 2 * check_radius + 1, 2 * check_radius + 1 ), dtype = bool )
  [ X, Y ] = indices( ( 2 * check_radius + 1, 2 * check_radius + 1 ), dtype = int32 )
  R2 = ( X - float( check_radius ) )**2 + ( Y - float( check_radius ) )**2
  putmask( sub_image_mask, R2 <= check_radius**2, True )
  if ( len( radec_list ) > 0 ):
    for radec in radec_list:
      pos = calculate_source_position( image, radec )
      if ( ( pos[ 0 ] < 0.5 + float( check_radius ) ) or 
           ( pos[ 0 ] > float( image_size[ 0 ] ) + 0.5 - float( check_radius ) ) or
           ( pos[ 1 ] < 0.5 + float( check_radius ) ) or 
           ( pos[ 1 ] > float( image_size[ 1 ] ) + 0.5 - float( check_radius ) ) ):
        pos = [ - 1, - 1 ]
      else:
        int_pos = [ int( around( pos[ 0 ] ) ) - 1, int( around( pos[ 1 ] ) ) - 1 ]
        sub_image = pixel_map[ int_pos[ 0 ] - check_radius : int_pos[ 0 ] + check_radius + 1, 
                               int_pos[ 1 ] - check_radius : int_pos[ 1 ] + check_radius + 1 ]
        putmask( sub_image, sub_image_mask == False, 0. )
        sel = awhere( ( sub_image == get_aips_magic_value() ) | ( abs( sub_image ) > flux_min ) )
        if ( len( sel ) == 0 ):
          new_radec_list.append( radec )
          # blank selected area in image to prevent re-selection
          sub_image = pixel_map[ int_pos[ 0 ] - check_radius : int_pos[ 0 ] + check_radius + 1, 
                                 int_pos[ 1 ] - check_radius : int_pos[ 1 ] + check_radius + 1 ]
          putmask( sub_image, sub_image_mask == True, get_aips_magic_value() )
          pixel_map[ int_pos[ 0 ] - check_radius : int_pos[ 0 ] + check_radius + 1, 
                     int_pos[ 1 ] - check_radius : int_pos[ 1 ] + check_radius + 1 ] = sub_image
        else:
          pos = [ - 1, - 1 ]
          new_radec_list.append( [ -1., -1. ] )
          if print_info:
            print '... source position %s has been rejected' % ( repr( degdeg_to_hmsdms( radec ) ) )
      pos_list.append( pos )
  else:
    for i in range( len( flux_list ) ):
      pos_found = False
      while ( not pos_found ):
        pos = [ random.random() * float( image_size[ 0 ] ) + 0.5, 
                random.random() * float( image_size[ 1 ] ) + 0.5 ]
        if ( ( pos[ 0 ] < 0.5 + float( check_radius ) ) or 
             ( pos[ 0 ] > float( image_size[ 0 ] ) + 0.5 - float( check_radius ) ) or
             ( pos[ 1 ] < 0.5 + float( check_radius ) ) or 
             ( pos[ 1 ] > float( image_size[ 1 ] ) + 0.5 - float( check_radius ) ) ):
           continue
        int_pos = [ int( around( pos[ 0 ] ) ) - 1, int( around( pos[ 1 ] ) ) - 1 ]
        sub_image = pixel_map[ int_pos[ 0 ] - check_radius : int_pos[ 0 ] + check_radius + 1, 
                               int_pos[ 1 ] - check_radius : int_pos[ 1 ] + check_radius + 1 ]
        putmask( sub_image, sub_image_mask == False, 0. )
        sel = awhere( ( sub_image == get_aips_magic_value() ) | ( abs( sub_image ) > flux_min ) )
        if ( len( sel ) > 0 ):
          continue
        pos_list.append( pos )
        radec = calculate_source_radec( image, pos )
        new_radec_list.append( radec )
        pos_found = True
        # blank selected area in image to prevent re-selection
        sub_image = pixel_map[ int_pos[ 0 ] - check_radius : int_pos[ 0 ] + check_radius + 1, 
                               int_pos[ 1 ] - check_radius : int_pos[ 1 ] + check_radius + 1 ]
        putmask( sub_image, sub_image_mask == True, get_aips_magic_value() )
        pixel_map[ int_pos[ 0 ] - check_radius : int_pos[ 0 ] + check_radius + 1, 
                   int_pos[ 1 ] - check_radius : int_pos[ 1 ] + check_radius + 1 ] = sub_image

  if ( not add_sources ):
    return [ None, new_radec_list ]

  # add model sources to image
  temp_image = get_aips_file( image.disk, image.name, 'TEMP', - 1, 'MA' )
  call_aips_task( 'MOVE', indata = image, outdata = temp_image, userid = get_aips_userid() )
  pixel_ref = get_pixel_reference( temp_image )
  pixel_scale = get_pixel_size( temp_image, make_absolute = False )
  model_table = new_table( temp_image, 'CC', 0 )
  model_row = new_table_row( model_table )
  added_source_count = 0
  for i in range( len( flux_list ) ):
    pos = pos_list[ i ]
    if ( pos != [ - 1, - 1 ] ):
      model_row.deltax = ( pos[ 0 ] - pixel_ref[ 0 ] ) * pixel_scale[ 0 ] / 3600.
      model_row.deltay = ( pos[ 1 ] - pixel_ref[ 1 ] ) * pixel_scale[ 1 ] / 3600.
      model_row.flux = flux_list[ i ]
      model_table.append( model_row )
      added_source_count = added_source_count + 1
  model_table.close()
  if ( added_source_count > 0 ):
    new_image = get_aips_file( image.disk, image.name, 'ADDSRC', - 1, 'MA' )
    model_version = temp_image.table_highver( 'CC' )
    call_aips_task( 'CCRES', indata = temp_image, in2data = temp_image, invers = model_version,
        outdata = new_image, optype = 'ADD' )
    new_image.zap_table( 'CC', 0 )
  else:
    new_image = None
  temp_image.zap()

  return [ new_image, new_radec_list ]

###############################################################################

def add_point_sources_to_noise_image( image, noise, flux_list, radec_list = [],
    check_radius = 10., print_info = False ):

  if ( len( radec_list ) > 0 ):
    if ( len( radec_list ) != len( flux_list ) ):
      raise error( 'lengths of flux list and RA DEC list do not match' )

  # make list of valid source pixel positions
  new_radec_list = []
  pos_list = []
  pixel_map = get_image_pixels( image )
  image_size = list( pixel_map.shape )
  sub_image_mask = zeros( ( 2 * check_radius + 1, 2 * check_radius + 1 ), dtype = bool )
  [ X, Y ] = indices( ( 2 * check_radius + 1, 2 * check_radius + 1 ), dtype = int32 )
  R2 = ( X - float( check_radius ) )**2 + ( Y - float( check_radius ) )**2
  putmask( sub_image_mask, R2 <= check_radius**2, True )
  if ( len( radec_list ) > 0 ):
    for radec in radec_list:
      pos = calculate_source_position( image, radec )
      if ( ( pos[ 0 ] < 0.5 + float( check_radius ) ) or 
           ( pos[ 0 ] > float( image_size[ 0 ] ) + 0.5 - float( check_radius ) ) or
           ( pos[ 1 ] < 0.5 + float( check_radius ) ) or 
           ( pos[ 1 ] > float( image_size[ 1 ] ) + 0.5 - float( check_radius ) ) ):
        pos = [ - 1, - 1 ]
      else:
        int_pos = [ int( around( pos[ 0 ] ) ) - 1, int( around( pos[ 1 ] ) ) - 1 ]
        sub_image = pixel_map[ int_pos[ 0 ] - check_radius : int_pos[ 0 ] + check_radius + 1, 
                               int_pos[ 1 ] - check_radius : int_pos[ 1 ] + check_radius + 1 ]
        putmask( sub_image, sub_image_mask, 0. )
        sel = awhere( sub_image == get_aips_magic_value() )
        if ( len( sel ) == 0 ):
          new_radec_list.append( radec )
        else:
          pos = [ - 1, - 1 ]
          new_radec_list.append( [ -1., -1. ] )
          if print_info:
            print '... source position %s has been rejected' % ( repr( degdeg_to_hmsdms( radec ) ) )
      pos_list.append( pos )
  else:
    for i in range( len( flux_list ) ):
      pos_found = False
      while ( not pos_found ):
        pos = [ random.random() * float( image_size[ 0 ] ) + 0.5, 
                random.random() * float( image_size[ 1 ] ) + 0.5 ]
        if ( ( pos[ 0 ] < 0.5 + float( check_radius ) ) or 
             ( pos[ 0 ] > float( image_size[ 0 ] ) + 0.5 - float( check_radius ) ) or
             ( pos[ 1 ] < 0.5 + float( check_radius ) ) or 
             ( pos[ 1 ] > float( image_size[ 1 ] ) + 0.5 - float( check_radius ) ) ):
           continue
        int_pos = [ int( around( pos[ 0 ] ) ) - 1, int( around( pos[ 1 ] ) ) - 1 ]
        sub_image = pixel_map[ int_pos[ 0 ] - check_radius : int_pos[ 0 ] + check_radius + 1, 
                               int_pos[ 1 ] - check_radius : int_pos[ 1 ] + check_radius + 1 ]
        putmask( sub_image, sub_image_mask, 0. )
        sel = awhere( sub_image == get_aips_magic_value() )
        if ( len( sel ) > 0 ):
          continue
        pos_list.append( pos )
        radec = calculate_source_radec( image, pos )
        new_radec_list.append( radec )
        pos_found = True

  # add model sources to image
  new_image = get_aips_file( image.disk, image.name, 'ADDSRC', - 1, 'MA' )
  call_aips_task( 'MOVE', indata = image, outdata = new_image, userid = get_aips_userid() )
  sub_flux_list = []
  sub_pos_list = []
  sub_beam_list = []
  sub_count = 0
  flux = noise
  factor = 0.
  for i in range( len( flux_list ) ):
    pos = pos_list[ i ]
    if ( pos != [ - 1, - 1 ] ):
      sub_flux_list.append( flux_list[ i ] )
      sub_pos_list.append( pos )
      sub_beam_list.append( convert_beam_size( image, to_pixel = True ) )
      sub_count = sub_count + 1
      if ( sub_count == 4 ):
        temp_image = get_aips_file( image.disk, image.name, 'TEMP', - 1, 'MA' )
        call_aips_task( 'IMMOD', indata = new_image, outdata = temp_image, factor = factor,
            flux = flux, opcode = 'GAUS', ngaus = sub_count, fmax = sub_flux_list,
            fpos = sub_pos_list, fwidth = sub_beam_list )
        new_image.zap()
        temp_image.rename( name = new_image.name, klass = new_image.klass, seq = new_image.seq )
        sub_flux_list = []
        sub_pos_list = []
        sub_beam_list = []
        sub_count = 0
        flux = 0.
        factor = 1.
  if ( sub_count != 0 ):
    temp_image = get_aips_file( image.disk, image.name, 'TEMP', - 1, 'MA' )
    call_aips_task( 'IMMOD', indata = new_image, outdata = temp_image, factor = factor,
        flux = flux, opcode = 'GAUS', ngaus = sub_count, fmax = sub_flux_list,
        fpos = sub_pos_list, fwidth = sub_beam_list )
    new_image.zap()
    temp_image.rename( name = new_image.name, klass = new_image.klass, seq = new_image.seq )

  return [ new_image, new_radec_list ]

###############################################################################

def write_catalog( cat, file_name ):
  row_length = len( cat[ 0 ] )
  for row in cat:
    if ( len( row ) != row_length ):
      raise error( 'catalog should have rows of equal length' )
  if file_exists( file_name ):
    remove_file( file_name )
  cat_file = file( file_name, mode = 'w' )
  line = ''
  for field in cat[ 0 ]:
    line = line + field + '  '
  line = line[ : - 2 ] + '\n'
  cat_file.write( line )
  for row in cat[ 1 : ]:
    line = ''
    for field in row:
      line = line + repr( field ) + '  '
    line = line[ : - 2 ] + '\n'
    cat_file.write( line )
  cat_file.close()
  return

###############################################################################

def read_catalog( file_name ):
  if ( not file_exists( file_name ) ):
    raise error( 'catalog file %s does not exists' % ( file_name ) )
  cat = []
  cat_file = file( file_name, mode = 'r' )
  line_count = 0
  for line in cat_file:
    words = [ word.strip() for word in line.split() ]
    if ( len( words ) == 0 ):
      continue
    if ( words[ 0 ][ 0 ] == '#' ):
      continue
    line_count = line_count + 1
    if ( line_count == 1 ):
      row_length = len( words )
      cat.append( words )
    else:
      if ( len( words ) != row_length ):
        raise error( 'catalog should have rows of equal length' )
      values = []
      for word in words:
        if ( '.' in word ):
          values.append( float( word ) )
        else:
          values.append( int( word ) )
      cat.append( values )
  cat_file.close()
  return cat

###############################################################################

def convert_catalog_epoch_from_b1950_to_j2000( cat_list ):
  new_cat_list = [ cat_list[ 0 ] ]
  ra_index = cat_list[ 0 ].index( 'RA' )
  dec_index = cat_list[ 0 ].index( 'DEC' )
  for cat in cat_list[ 1 : ]:
    new_cat = [ field for field in cat ]
    radec = [ cat[ ra_index ], cat[ dec_index ] ]
    new_radec = convert_b1950_to_j2000( radec )
    new_cat[ ra_index ] = new_radec[ 0 ]
    new_cat[ dec_index ] = new_radec[ 1 ]
    new_cat_list.append( new_cat )
  return new_cat_list

###############################################################################

def convert_catalog_epoch_from_j2000_to_b1950( cat_list ):
  new_cat_list = [ cat_list[ 0 ] ]
  ra_index = cat_list[ 0 ].index( 'RA' )
  dec_index = cat_list[ 0 ].index( 'DEC' )
  for cat in cat_list[ 1 : ]:
    new_cat = [ field for field in cat ]
    radec = [ cat[ ra_index ], cat[ dec_index ] ]
    new_radec = convert_j2000_to_b1950( radec )
    new_cat[ ra_index ] = new_radec[ 0 ]
    new_cat[ dec_index ] = new_radec[ 1 ]
    new_cat_list.append( new_cat )
  return new_cat_list

###############################################################################

def correct_catalog_for_pb_attenuation( uvim, cat_list, cutoff = 0.1,
    invert = False, offset_rp = None ):
  new_cat_list = [ cat_list[ 0 ] ]
  ra_index = cat_list[ 0 ].index( 'RA' )
  dec_index = cat_list[ 0 ].index( 'DEC' )
  try:
    peak_index = cat_list[ 0 ].index( 'PEAK' )
  except:
    peak_index = - 1
  try:
    flux_index = cat_list[ 0 ].index( 'FLUX' )
  except:
    flux_index = - 1
  try:
    peak_err_index = cat_list[ 0 ].index( 'PEAK_ERROR' )
  except:
    peak_err_index = - 1
  try:
    flux_err_index = cat_list[ 0 ].index( 'FLUX_ERROR' )
  except:
    flux_err_index = - 1
  radec_list = []
  for cat in cat_list[ 1 : ]:
    radec = [ cat[ ra_index ], cat[ dec_index ] ]
    if ( not offset_rp is None ):
      [ r, p ] = offset_rp
      radec = calculate_offset_position( radec, r, p )
    radec_list.append( radec )
  A_list = get_primary_beam_attenuations( uvim, radec_list )
  for i in range( len( A_list ) ):
    A = A_list[ i ]
    if ( A >= cutoff ):
      if invert:
        A = 1. / A
      new_cat = [ field for field in cat_list[ i + 1 ] ]
      if ( peak_index > - 1 ):
        new_cat[ peak_index ] = new_cat[ peak_index ] / A
      if ( flux_index > - 1 ):
        new_cat[ flux_index ] = new_cat[ flux_index ] / A
      if ( peak_err_index > - 1 ):
        new_cat[ peak_err_index ] = new_cat[ peak_err_index ] / A
      if ( flux_err_index > - 1 ):
        new_cat[ flux_err_index ] = new_cat[ flux_err_index ] / A
      new_cat_list.append( new_cat )
  return new_cat_list

###############################################################################

def generate_power_law_flux_list( count, alpha = - 2., flux_limits = [ 0.1, 1. ] ):
  C0 = flux_limits[ 0 ]**( alpha + 1. )
  C1 = flux_limits[ 1 ]**( alpha + 1. ) - C0
  flux_list = []
  for i in range( count ):
    r = random.random()
    flux = ( C0 + C1 * r )**( 1. / ( alpha + 1. ) )
    flux_list.append( flux )
  flux_list.sort( cmp = lambda a, b: cmp( b, a ) )
  return flux_list

###############################################################################

def correct_catalog_for_b1950_to_j2000( cat_list ):
  new_cat_list = [ cat_list[ 0 ] ]
  ra_index = cat_list[ 0 ].index( 'RA' )
  dec_index = cat_list[ 0 ].index( 'DEC' )
  for cat in cat_list[ 1 : ]:
    radec = [ cat[ ra_index ], cat[ dec_index ] ]
    new_radec = convert_b1950_to_j2000( radec )
    new_cat = [ field for field in cat ]
    new_cat[ ra_index ] = new_radec[ 0 ]
    new_cat[ dec_index ] = new_radec[ 1 ]
    new_cat_list.append( new_cat )
  return new_cat_list

###############################################################################

def write_starbase_catalog( cat, file_name ):
  row_length = len( cat[ 0 ] )
  for row in cat:
    if ( len( row ) != row_length ):
      raise error( 'catalog should have rows of equal length' )
  if file_exists( file_name ):
    remove_file( file_name )
  cat_file = file( file_name, mode = 'w' )
  line = file_name.split( '/' )[ - 1 ] + '\n'
  cat_file.write( line )
  line = ''
  for field in cat[ 0 ]:
    line = line + field + '\t'
  line = line[ : - 1 ] + '\n'
  cat_file.write( line )
  line = '-----\n'
  cat_file.write( line )
  for row in cat[ 1 : ]:
    line = ''
    for field in row:
      if ( type( field ) == type( '' ) ):
        line = line + field + '\t'
      else:
        line = line + repr( field ) + '\t'
    line = line[ : - 1 ] + '\n'
    cat_file.write( line )
  cat_file.close()
  return

###############################################################################

def read_starbase_catalog( file_name, drop_strings = False, new_header_row = None ):
  if ( not file_exists( file_name ) ):
    raise error( 'catalog file %s does not exists' % ( file_name ) )
  cat = []
  cat_file = file( file_name, mode = 'r' )
  line_count = 0
  table_found = False
  for line in cat_file:
    words = [ word.strip() for word in line.split( '\t' ) ]
    if ( len( words ) == 0 ):
      continue
    if ( len( words[ 0 ] ) == 0 ):
      continue
    if ( words[ 0 ][ 0 ] == '#' ):
      continue
    if ( words[ 0 ][ 0 ] == '[EOD]' ):
      break
    if ( not table_found ):
      if ( words[ 0 ][ 0 ] == '-' ):
        row_length = len( last_words )
        cat.append( last_words )
        table_found = True
      else:
        last_words = [ w for w in words ]
    else:
      if ( len( words ) != row_length ):
        raise error( 'catalog should have rows of equal length' )
      values = []
      for word in words:
        try:
          if ( '.' in word ):
            values.append( float( word ) )
          else:
            values.append( int( word ) )
        except ValueError:
          values.append( word )
      cat.append( values )
  cat_file.close()
  
  if ( drop_strings and ( len( cat ) > 1 ) ):
    indices = range( len( cat[ 0 ] ) )
    for line in cat[ 1 : ]:
      for i in range( len( line ) ):
        if ( i in indices ):
          if ( type( line[ i ] ) == type( '' ) ):
            indices.remove( i )
    cat2 = []
    for line in cat:
      cat2.append( [ line[ i ] for i in indices ] )
    cat = cat2
  
  if ( not new_header_row is None ):
    if ( len( new_header_row ) == len( cat[ 0 ] ) ):
      cat[ 0 ] = [ x for x in new_header_row ]
    elif ( drop_strings and ( len( new_header_row ) == row_length ) ):
      cat[ 0 ] = [ new_header_row[ i ] for i in indices ]
    else:
      raise error( 'new header row does not match number of columns' )
  
  return cat

###############################################################################

def write_lions_catalog( cat, file_name ):
  row_length = len( cat[ 0 ] )
  for row in cat:
    if ( len( row ) != row_length ):
      raise error( 'catalog should have rows of equal length' )
  if file_exists( file_name ):
    remove_file( file_name )
  cat_file = file( file_name, mode = 'w' )
  id = 1
  ra_index = cat[ 0 ].index( 'RA' )
  dec_index = cat[ 0 ].index( 'DEC' )
  flux_index = cat[ 0 ].index( 'FLUX' )
  for row in cat[ 1 : ]:
    line = ''
    for field in [ id, row[ ra_index ], row[ dec_index ], row[ flux_index ] ]:
      line = line + repr( field ).ljust( 20 )
    line = line.strip() + '\n'
    cat_file.write( line )
    id = id + 1
  cat_file.close()
  return

###############################################################################

def merge_catalog_components( cat_list, max_separation = 120. ):
# max_separation in arcsec
  
  ra_index = cat_list[ 0 ].index( 'RA' )
  dec_index = cat_list[ 0 ].index( 'DEC' )
  flux_index = cat_list[ 0 ].index( 'FLUX' )
  flux_err_index = cat_list[ 0 ].index( 'FLUX_ERROR' )
  try:
    peak_index = cat_list[ 0 ].index( 'PEAK' )
  except:
    peak_index = flux_index
  
  # group multiples
  merge_ids = [ 0 for l in cat_list ]
  merge_id = 0
  for i in range( 1, len( cat_list ) ):
    radec1 = [ cat_list[ i ][ ra_index ], cat_list[ i ][ dec_index ] ]
    for j in range( i + 1, len( sp5sc_cat ) ):
      radec2 = [ cat_list[ j ][ ra_index ], cat_list[ j ][ dec_index ] ]
      [ r, p ] = calculate_angular_separation( radec1, radec2 )
      if ( r < max_separation / 3600. ):
        if ( merge_ids[ i ] == 0 ):
          if ( merge_ids[ j ] == 0 ):
            merge_id = merge_id + 1
            merge_ids[ i ] = merge_id
            merge_ids[ j ] = merge_id
          else:
            merge_ids[ i ] = merge_ids[ j ]
        else:
          if ( merge_ids[ j ] == 0 ):
            merge_ids[ j ] = merge_ids[ i ]
          else:
            if ( merge_ids[ i ] != merge_ids[ j ] ):
              old_merge_id = merge_ids[ j ]
              for k in range( 1, len( cat_list ) ):
                if ( merge_ids[ k ] == old_merge_id ):
                  merge_id[ k ] = merge_id[ i ]
  
  # make merged catalog
  if ( merge_id == 0 ):
    merged_cat_list = [ [ x for x in s ] for s in cat_list ]
  else:
    merged_cat_list = [ [ x for x in cat_list[ 0 ] ] ]
    for i in range( 1, len( cat_list ) ):
      if ( merge_ids[ i ] == 0 ):
        merged_cat_list.append( [ x for x in cat_list[ i ] ] )
    for k in range( 1, merge_id + 1 ):
      try:
        i = merge_ids.index( k )
      except:
        continue
      else:
        l = [ x for x in cat_list[ i ] ]
        for j in range( i + 1, len( cat_list ) ):
          if ( merge_ids[ j ] == k ):
            m = [ x for x in cat_list[ j ] ]
            if ( l[ peak_index ] > m[ peak_index ] ):
              l[ flux_index ] = l[ flux_index ] + m[ flux_index ]
              l[ flux_err_index ] = sqrt( l[ flux_err_index ]**2 + m[ flux_err_index ]**2 )
            else:
              m[ flux_index ] = m[ flux_index ] + l[ flux_index ]
              m[ flux_err_index ] = sqrt( m[ flux_err_index ]**2 + l[ flux_err_index ]**2 )
              l = [ x for x in m ]
        merged_cat_list.append( l )
  
  # sort on peak flux
  merged_cat_list.sort( cmp = lambda a, b: cmp( b[ peak_index ], a[ peak_index ] ) )
  
  return merged_cat_list

###############################################################################

def read_ska_catalog( file_name, header_row = False ):
  if ( not file_exists( file_name ) ):
    raise error( 'catalog file %s does not exists' % ( file_name ) )
  cat = []
  row_length = 4
  if header_row:
    cat.append( [ 'ID', 'RA', 'DEC', 'FLUX' ] )
  cat_file = file( file_name, mode = 'r' )
  line_count = 0
  table_found = False
  for line in cat_file:
    words = [ word.strip() for word in line.split( '\t' ) ]
    if ( len( words ) == 0 ):
      continue
    if ( len( words[ 0 ] ) == 0 ):
      continue
    if ( words[ 0 ][ 0 ] == '#' ):
      continue
    if ( len( words ) != row_length ):
      raise error( 'catalog should have rows of equal length' )
    values = []
    for word in words:
      try:
        if ( '.' in word ):
          values.append( float( word ) )
        else:
          values.append( int( word ) )
      except ValueError:
        values.append( word )
    cat.append( values )
  cat_file.close()
  return cat

###############################################################################

# format example
#
#Title: AIPS VL
#Created by      FITLD on 03-DEC-2010 10:38:11
#Last written by FITLD on 03-DEC-2010 10:38:11
#Ncol  18  Nrow 1810672    Sort cols:   1 (ASCEND)
#    Table has    30 keyword-value pairs:
#   REVISION =            1
#   BM_MAJOR =  1.2500000186265D-02
#   BM_MINOR =  1.2500000186265D-02
#   BM_PA    =  0.0000000000000D+00
#   SORTORT  =            1
#   NUM_INDE =      1810672
#   INDEX00  =            1
#   INDEX01  =        74608
#   INDEX02  =       148837
#   INDEX03  =       223936
#   INDEX04  =       299369
#   INDEX05  =       373984
#   INDEX06  =       450586
#   INDEX07  =       525193
#   INDEX08  =       600541
#   INDEX09  =       674935
#   INDEX10  =       749192
#   INDEX11  =       823509
#   INDEX12  =       898254
#   INDEX13  =       972989
#   INDEX14  =      1047472
#   INDEX15  =      1121930
#   INDEX16  =      1195172
#   INDEX17  =      1267782
#   INDEX18  =      1351430
#   INDEX19  =      1436113
#   INDEX20  =      1512483
#   INDEX21  =      1590676
#   INDEX22  =      1663460
#   INDEX23  =      1735545
#   Table can be written as a FITS ASCII table
# 
#COL. NO.        1             2            3          4          5           6           7
#     ROW  RA(2000)      DEC(2000)      PEAK INT   MAJOR AX   MINOR AX   POSANGLE    Q CENTER
#  NUMBER  DEGREE        DEGREE         JY/BEAM    DEGREE     DEGREE     DEGREEE     JY/BEAM
#       1  3.999606D-04  -3.411935D+01  2.483E-03  1.447E-02  1.250E-02  -1.259E+01  -4.436E-04
#       2  5.354161D-04  -3.844127D+01  2.388E-03  2.212E-02  1.250E-02  -2.727E+01  -1.512E-06
#
#COL. NO.       8          9          10         11         12         13         14          15
#     ROW  U CENTER    P FLUX     I RMS      POL RMS    RES RMS    RES PEAK   RES FLUX    CENTER X
#  NUMBER  JY/BEAM     JY         JY/BEAM    JY/BEAM    JY/BEAM    JY/BEAM    JY          PIXEL
#       1  -7.046E-04  7.322E-04  4.803E-04  3.008E-04  3.477E-04  8.829E-04  -5.658E-04  5.119E+02
#       2   1.209E-04  7.179E-04  4.837E-04  3.070E-04  3.629E-04  8.158E-04  -5.483E-04  5.119E+02
#
#COL. NO.          16            17            18
#     ROW      CENTER Y       FIELD         JD PROCE
#  NUMBER      PIXEL                        DAY
#       1      9.643E+02      C0000M36       2450823
#       2      8.870E+02      C0000M40       2450823

def read_aips_catalog_table( image, version = 0 ):
  # WARNING: NOT YET DONE

  
  if ( not table_exists( image, 'VL', catalog_version ) ):
    raise error( 'source catalog table does not exist' )
  wiz_im = wizardry( image )
  cat_table = wiz_im.table( 'VL', catalog_version )

  cat_list = [ [ 'ID', 'RA', 'DEC', 'PEAK', 'FLUX', 'MAJ', 'MIN', 'PA', 'RMS' ] ]
  id = 0
  for row in cat_table:
    id = id + 1
    ra = row[ 'ra(2000)' ]
    dec = row[ 'dec(2000)' ]
    peak = row.peak_int
    flux = row.i_flux
    maj = row.major_ax
    min = row.minor_ax
    pa = row.posangle
    peak_res = row.res_peak
    ra_error = row.err_dltx
    dec_error = row.err_dlty
    peak_error = row.err_peak
    flux_error = row.err_flux
    maj_error = row.err_maja
    min_error = row.err_mina
    pa_error = row.err_pa
    cat_list.append( [ id, ra, dec, peak, flux, maj, min, pa, peak_res,
        ra_error, dec_error, peak_error, flux_error, maj_error, min_error, pa_error ] )

  return cat_list


###############################################################################

def run_bdsm( fits_file_name, work_path = '/users/hintema/data2/bdsm/',
#    bdsm_exe = '/export/data_2/hintema/archive/bdsm/bdsm',
    bdsm_exe = 'bdsm',
    run_id = 'def', verbose = False, do_gaussians = True, do_shapelets = False,
    pixel_sigma = 5., island_sigma = 3., island_pixels = 4, image_plane = 0,
    stack_size_gb = 4., cleanup = True, keep_residual = False,
    box_size = 0, step_size = 0, keep_metadata = True,
    fitfreely = 'true', iniguess = 'default', flagsmallsrc = 'true', 
    mean_map = 'default' ):
  
  if ( not file_exists( fits_file_name ) ):
    raise error( 'FITS file %s does not exists' % ( fits_file_name ) )
  if ( fits_file_name[ 0 ] == '/' ):
    full_fits_name = fits_file_name
  else:
    current_path = os.getenv( 'PWD' )
    full_fits_name = current_path + '/' + fits_file_name
  index = len( full_fits_name ) - full_fits_name[ : : -1 ].find( '/' )
  fits_name = full_fits_name[ index : ]
  fits_file_path = full_fits_name[ 0 : index ]
  if ( not directory_exists( work_path ) ):
    mkdir( work_path )
  fits_path = work_path + 'fits/'
  scratch_path = work_path + 'scratch/'
  srl_path = work_path + 'srl/'
  plot_path = work_path + 'plots/'
  if ( not directory_exists( fits_path ) ):
    mkdir( fits_path )
  if ( not directory_exists( scratch_path ) ):
    mkdir( scratch_path )
  if ( not directory_exists( srl_path ) ):
    mkdir( srl_path )
#  if ( not directory_exists( plot_path ) ):
#    mkdir( plot_path )
  
  # make paradefine file
  par_file_name = 'paradefine'
  if file_exists( par_file_name ):
    remove_file( par_file_name )
  par_file = file( par_file_name, mode = 'w' )
  par_file.write( "fitsdir = '%s'\n" % ( fits_path ) )
  par_file.write( "scratch = '%s'\n" % ( scratch_path ) )
  par_file.write( "srldir = '%s'\n" % ( srl_path ) )
  par_file.write( "plotdir = '%s'\n" % ( scratch_path ) )
  par_file.write( "fitsname = '%s'\n" % ( fits_name ) )
  par_file.write( "solnname = '%s'\n" % ( run_id ) )
  if verbose:
    par_file.write( "runcode = 'av'\n" )
  else:
    par_file.write( "runcode = 'aq'\n" )
  switch = ''
  if do_gaussians:
    switch = switch + '1'
  else:
    switch = switch + '0'
  if do_shapelets:
    switch = switch + '1'
  else:
    switch = switch + '0'
  par_file.write( "gausshap = '%s'\n" % ( switch ) )
  par_file.write( "bmpersrc_th = 0\n" ) # 0\n" )
  par_file.write( "boxsize_th = %s\n" % ( repr( box_size ) ) )
  par_file.write( "stepsize_th = %s\n" % ( repr( step_size ) ) )
  par_file.write( "thresh = 'hard'\n" ) # 'default'\n" )
  par_file.write( "fdr_alpha = 0\n" )
  par_file.write( "thresh_pix = %s\n" % ( repr( pixel_sigma ) ) )
  par_file.write( "thresh_isl = %s\n" % ( repr( island_sigma ) ) )
  par_file.write( "minpix_isl = %s\n" % ( repr( island_pixels ) ) )
  par_file.write( "maxsize_beam = 40\n" ) # 0\n" )
  par_file.write( "rms_map = 'default'\n" )
  par_file.write( "ch0plane = %s\n" % ( repr( image_plane ) ) )
  par_file.write( "takemeanclip = 'default'\n" )
  par_file.write( "BMAJ = 0.0\n" )
  par_file.write( "BMIN = 0.0\n" )
  par_file.write( "BPA = 0.0\n" )
  par_file.write( "editor = 'gedit'\n" )
  # new paradefine entries
  par_file.write( "mean_map = '%s'\n" % ( mean_map ) )
  par_file.write( "fitfreely = '%s'\n" % ( fitfreely ) )
  par_file.write( "iniguess = '%s'\n" % ( iniguess ) )
  par_file.write( "flagsmallsrc = '%s'\n" % ( flagsmallsrc ) )
  par_file.write( "dummy = 'dummy'\n" )
  par_file.close()
  
  # create link to fits file
  if file_exists( fits_path + fits_name ):
    remove_file( fits_path + fits_name )
  system( "ln -s %s %s" % ( full_fits_name, fits_path + fits_name ) )
  
  # run BDSM
  ret = system( "ulimit -s %d; %s" % (
      int( stack_size_gb * 1024**2 ), bdsm_exe ) )
  if ( ret != 0 ):
    print 'BDSM returns exit status = %s' % ( repr( ret ) )
  
  # copy relevant output files
  index = fits_name[ : : -1 ].find( '.' )
  name = fits_name
  if ( index > -1 ):
    index = len( fits_name ) - index
    if ( fits_name[ index : index + 3 ].lower() == 'fit' ):
      name = fits_name[ 0 : index - 1 ]
  if keep_metadata:
    bstat_name = scratch_path + name + '.' + run_id + '.bstat'
    copy_file( bstat_name, name + '.' + run_id + '.bstat' )
    bparms_name = scratch_path + name + '.' + run_id + '.bparms'
    copy_file( bparms_name, name + '.' + run_id + '.bparms' )
  srl_name = srl_path + name + '.' + run_id + '.srl'
  copy_file( srl_name, name + '.' + run_id + '.srl' )
  if do_gaussians:
    gaul_name = srl_path + name + '.' + run_id + '.gaul'
    copy_file( gaul_name, name + '.' + run_id + '.gaul' )
  
  #cleanup
  if cleanup:
    remove_file( fits_path + fits_name )
    prefix = scratch_path + name + '.'
    for pf in [ 'header', 'img', 'qc_cc' ]:
      remove_file( prefix + pf )
    prefix = prefix + run_id + '.'
    for pf in [ 'bparms', 'bstat', 'ch0.img', 'header', 'islandlist', 'islandlist.asc',
        'mean.img', 'rank.img', 'rank.nmg', 'resid.gaus.img', 'rmsd.img', 'scratch',
        'snum', 'srcim.gaus.img', 'norm.img' ]:
      remove_file( prefix + pf )
    prefix = fits_path + name + '.' + run_id + '.'
    for pf in [ 'resid.gaus.FITS' ]:
      if keep_residual:
        res_fits_file_name = fits_file_path + name + '.' + run_id + '.' + pf
        copy_file( prefix + pf, res_fits_file_name )
      remove_file( prefix + pf )
    prefix = srl_path + name + '.' + run_id + '.'
    for pf in [ 'gaul', 'gaul.bin', 'gaul.FITS', 'gaul.reg', 'gaul.star', 'ini.gaul',
        'ini.gaul.bin', 'outputfiles', 'shapcoef.c', 'shapelet.c', 'srl', 'srl.bin',
        'srl.FITS', 'srl.star' ]:
      remove_file( prefix + pf )
  
  return ret

###############################################################################

def read_pybdsm_ascii_catalog( cat_file_name ):
  if ( not file_exists( cat_file_name ) ):
    raise error( 'source catalog %s was not found' % cat_file_name )
  cat_file = file( cat_file_name, mode = 'r' )
  cat_list = []
  column_count = -1
  header_found = False
  for line in cat_file:
    words = [ word.strip() for word in line.split() ]
    if ( len( words ) == 0 ):
      continue
    if ( not header_found ):
      if ( ( words[ 0 ] == '#' ) and ( 
          ( words[ 1 ] == 'Source_id' ) or ( words[ 1 ] == 'Gaus_id' ) ) ):
        column_count = len( words ) - 1
        cat_list.append( words[ 1 : ] )
        header_found = True
    else:
      if ( len( words ) != column_count ):
        continue
      new_words = []
      for word in words:
        try:
          dummy = float( word )
        except ValueError:
          new_words.append( word )
          continue
        try:
          dummy = int( word )
        except ValueError:
          new_words.append( float( word ) )
        else:
          new_words.append( int( word ) )
      cat_list.append( new_words )
  cat_file.close()
  return cat_list

###############################################################################

def combine_pybdsm_srl_gaul_catalogs( srl_cat, gaul_cat, group_islands = None ):
  
  # initialize
  isl_list = []
  src_list = []
  gau_list = []
  
  # get lists of ids from gaussian cat
  isl_gindex = gaul_cat[ 0 ].index( 'Isl_id' )
  src_gindex = gaul_cat[ 0 ].index( 'Source_id' )
  gau_gindex = gaul_cat[ 0 ].index( 'Gaus_id' )
  for i in range( 1, len( gaul_cat ) ):
    g = gaul_cat[ i ]
    isl_id = g[ isl_gindex ]
    if ( not isl_id in isl_list ):
      isl_list.append( isl_id )
      src_list.append( [] )
      gau_list.append( [] )
    isl_index = isl_list.index( isl_id )
    src_id = g[ src_gindex ] + 100000L
    if ( not src_id in src_list[ isl_index ] ):
      src_list[ isl_index ].append( src_id )
      gau_list[ isl_index ].append( [] )
    src_index = src_list[ isl_index ].index( src_id )
    gau_id = g[ gau_gindex ]
    gau_list[ isl_index ][ src_index ].append( i )
  
  # check lists of ids with source cat
  isl_sindex = srl_cat[ 0 ].index( 'Isl_id' )
  src_sindex = srl_cat[ 0 ].index( 'Source_id' )
  for i in range( 1, len( srl_cat ) ):
    s = srl_cat[ i ]
    isl_id = s[ isl_sindex ]
    if ( not isl_id in isl_list ):
      raise error( 'island ID %d in srl_cat not found in gaul_cat' % ( isl_id )  )
    isl_index = isl_list.index( isl_id )
    src_id = s[ src_sindex ] + 100000L
    if ( not src_id in src_list[ isl_index ] ):
      raise error( 'source ID %d in srl_cat not found in gaul_cat island' %
          ( src_id - 100000L )  )
    src_index = src_list[ isl_index ].index( src_id )
    src_list[ isl_index ][ src_index ] = i
  
  # group islands manually
  if ( not group_islands is None ):
    for group in group_islands:
      main_isl_id = group[ 0 ]
      main_isl_index = isl_list.index( main_isl_id )
      for isl_id in group[ 1 : ]:
        isl_index = isl_list.index( isl_id )
        isl_list.remove( isl_id )
        src_list[ main_isl_index ] = src_list[ main_isl_index ] + src_list[ isl_index ]
        src_list.remove( src_list[ isl_index ] )
        gau_list[ main_isl_index ] = gau_list[ main_isl_index ] + gau_list[ isl_index ]
        gau_list.remove( gau_list[ isl_index ] )
  
  # combine gaussians into islands
  new_cat = [ [ 'ID', 'name', 'RA', 'dRA', 'DEC', 'dDEC', 'Si', 'dSi', 'Sp', 'dSp',
      'noise', 'Ngauss' ] ]
  rms_sindex = srl_cat[ 0 ].index( 'Isl_rms' )
  sp_sindex = srl_cat[ 0 ].index( 'Peak_flux' )
  dsp_sindex = srl_cat[ 0 ].index( 'E_Peak_flux' )
  si_gindex = gaul_cat[ 0 ].index( 'Total_flux' )
  dsi_gindex = gaul_cat[ 0 ].index( 'E_Total_flux' )
  ra_gindex = gaul_cat[ 0 ].index( 'RA' )
  dra_gindex = gaul_cat[ 0 ].index( 'E_RA' )
  dec_gindex = gaul_cat[ 0 ].index( 'DEC' )
  ddec_gindex = gaul_cat[ 0 ].index( 'E_DEC' )
  for i in range( len( isl_list ) ):
    ref_radec = None
    ra_array = None
    dec_array = None
    si_array = None
    sp = None
    dsp = None
    noise = None
    ngauss = 0
    for j in range( len( src_list[ i ] ) ):
      s = srl_cat[ src_list[ i ][ j ] ]
      if ( noise is None ):
        noise = s[ rms_sindex ]
      elif ( s[ rms_sindex ] > noise ):
        noise = s[ rms_sindex ]
      if ( sp is None ):
        sp = s[ sp_sindex ]
        dsp = s[ dsp_sindex ]
      elif ( s[ sp_sindex ] > sp ):
        sp = s[ sp_sindex ]
        dsp = s[ dsp_sindex ]
      for k in range( len( gau_list[ i ][ j ] ) ):
        g = gaul_cat[ gau_list[ i ][ j ][ k ] ]
        ngauss = ngauss + 1
        si = draw_from_gaussian_distribution( g[ si_gindex ], g[ dsi_gindex ] )
        if ( si_array is None ):
          si_array = si.copy()
        else:
          si_array = si_array + si
        radec = [ g[ ra_gindex ], g[ dec_gindex ] ]
        dradec = [ g[ dra_gindex ], g[ ddec_gindex ] ]
        if ( ref_radec is None ):    
          ref_radec = radec
          radec = [ 0., 0. ]
        else:
          [ r, p ] = calculate_angular_separation( ref_radec, radec )
          radec = calculate_offset_position( [ 0., 0. ], r, p )
          radec[ 0 ] = amodulo( radec[ 0 ] + 180., 360. ) - 180.
        ra = draw_from_gaussian_distribution( radec[ 0 ], dradec[ 0 ] )
        dec = draw_from_gaussian_distribution( radec[ 1 ], dradec[ 1 ] )
        if ( ra_array is None ):
          ra_array = si * ra
          dec_array = si * dec
        else:
          ra_array = ra_array + si * ra
          dec_array = dec_array + si * dec
    si = get_robust_mean_deviations( si_array )
    dsi = max( si[ 1 ], -si[ 2 ] )
    si = si[ 0 ]
    ra = get_robust_mean_deviations( ra_array / si_array )
    dra = max( ra[ 1 ], -ra[ 2 ] )
    ra = ra[ 0 ]
    dec = get_robust_mean_deviations( dec_array / si_array )
    ddec = max( dec[ 1 ], -dec[ 2 ] )
    dec = dec[ 0 ]
    [ r, p ] = calculate_angular_separation( [ 0., 0. ], [ ra, dec ] )
    radec = radec = calculate_offset_position( ref_radec, r, p )
    name = 'J' + radec_to_string( radec, decimals = [ 0, 0 ], 
        separators = [ '', '', '', '', '', '' ] )
    new_cat.append( [ isl_list[ i ], name, radec[ 0 ], dra, radec[ 1 ], ddec,
        si, dsi, sp, dsp, noise, ngauss ] )
  
  return new_cat

###############################################################################

def create_source_list_from_catalog( cat, use_peak_flux = False ):
  if ( len( cat ) <= 1 ):
    raise error( 'catalog has no content' )
  ra_index = -1
  for ra_label in [ '_RAJ2000', 'RA' ]: # TODO: to be expanded
    try:
      ra_index = cat[ 0 ].index( ra_label )
    except:
      continue
    else:
      break
  if ( ra_index < 0 ):
    raise error( 'RA column not found in catalog' )
  dec_index = -1
  for dec_label in [ '_DEJ2000', 'DEC', 'DE', 'Dec' ]: # TODO: to be expanded
    try:
      dec_index = cat[ 0 ].index( dec_label )
    except:
      continue
    else:
      break
  if ( dec_index < 0 ):
    raise error( 'DEC column not found in catalog' )
  if use_peak_flux:
    flux_labels = [ 'Peak_flux', 'Sp','Peak', 'peak_Jy' 'PEAK' ] # TODO: to be expanded
  else:
    flux_labels = [ 'Total_flux', 'St', 'Total', 'total_Jy', 'TOTAL' ] # TODO: to be expanded
  flux_index = -1
  for flux_label in flux_labels: 
    try:
      flux_index = cat[ 0 ].index( flux_label )
    except:
      continue
    else:
      break
  if ( flux_index < 0 ):
    raise error( 'FLUX column not found in catalog' )
  source_list = []
  for s in cat[ 1 : ]:
    source_list.append( [ [ float( s[ ra_index ] ), float( s[ dec_index ] ) ],
        s[ flux_index ] ] )
  return source_list

###############################################################################

def read_atlas_dr2_catalog( file_name, drop_strings = False, new_header_row = None ):
  if ( not file_exists( file_name ) ):
    raise error( 'catalog file %s does not exists' % ( file_name ) )
  cat = []
  cat_file = file( file_name, mode = 'r' )
  line_count = 0
  header_found = False
  table_found = False
  for line in cat_file:
    line_count = line_count + 1
    words = [ word.strip() for word in line.split() ]
    if ( not header_found ):
      if ( len( words ) <= 2 ):
        continue
      if ( len( words[ 0 ] ) == 0 ):
        continue
      if ( words[ 0 ][ 0 ] != '#' ):
        raise error( 'no header found' )
      if ( len( words[ 0 ] ) > 1 ):
        words[ 0 ] = words[ 0 ][ 1 : ]
      else:
        words = words[ 1 : ]
      if ( words[ 0 ] != 'ID' ):
        continue
      cat.append( words )
      column_count = len( words )
      header_found = True
      continue
    if ( not table_found ):
      if ( len( words ) == 0 ):
        continue
      if ( len( words[ 0 ] ) == 0 ):
        continue
      if ( words[ 0 ][ 0 ] == '#' ):
        continue
      table_found = True
    values = []
    if ( len( words ) != column_count ):
      raise error( 'number of columns does not match header in line %d' % 
          ( line_count ) )
    for word in words:
      try:
        if ( '.' in word ):
          values.append( float( word ) )
        else:
          values.append( int( word ) )
      except ValueError:
        values.append( word )
    cat.append( values )
  cat_file.close()
  
  if ( drop_strings and ( len( cat ) > 1 ) ):
    indices = range( len( cat[ 0 ] ) )
    for line in cat[ 1 : ]:
      for i in range( len( line ) ):
        if ( i in indices ):
          if ( type( line[ i ] ) == type( '' ) ):
            indices.remove( i )
    cat2 = []
    for line in cat:
      cat2.append( [ line[ i ] for i in indices ] )
    cat = cat2
  
  if ( not new_header_row is None ):
    if ( len( new_header_row ) == len( cat[ 0 ] ) ):
      cat[ 0 ] = [ x for x in new_header_row ]
    elif ( drop_strings and ( len( new_header_row ) == row_length ) ):
      cat[ 0 ] = [ new_header_row[ i ] for i in indices ]
    else:
      raise error( 'new header row does not match number of columns' )
  
  return cat

###############################################################################

def read_atlas_dr3_catalog( file_name, drop_strings = False, new_header_row = None ):
  if ( not file_exists( file_name ) ):
    raise error( 'catalog file %s does not exists' % ( file_name ) )
  cat = []
  cat_file = file( file_name, mode = 'r' )
  line_count = 0
  header_found = False
  table_found = False
  for line in cat_file:
    line_count = line_count + 1
    if ( not header_found ):
      words = [ word.strip() for word in line.split( '|' ) ]
      if ( len( words ) <= 1 ):
        continue
      if ( len( words[ 0 ] ) == 0 ):
        continue
      if ( words[ 0 ][ 0 ] != '#' ):
        raise error( 'no header found' )
      words[ 0 ] = words[ 0 ][ 1 : ].strip()
      cat.append( words )
      column_count = len( words )
      header_found = True
      continue
    words = [ word.strip() for word in line.split() ]
    if ( not table_found ):
      if ( len( words ) == 0 ):
        continue
      if ( len( words[ 0 ] ) == 0 ):
        continue
      if ( words[ 0 ][ 0 ] == '#' ):
        continue
      table_found = True
    values = []
    if ( len( words ) != column_count ):
      raise error( 'number of columns does not match header in line %d' % 
          ( line_count ) )
    for word in words:
      try:
        if ( '.' in word ):
          values.append( float( word ) )
        else:
          values.append( int( word ) )
      except ValueError:
        values.append( word )
    cat.append( values )
  cat_file.close()
  
  if ( drop_strings and ( len( cat ) > 1 ) ):
    indices = range( len( cat[ 0 ] ) )
    for line in cat[ 1 : ]:
      for i in range( len( line ) ):
        if ( i in indices ):
          if ( type( line[ i ] ) == type( '' ) ):
            indices.remove( i )
    cat2 = []
    for line in cat:
      cat2.append( [ line[ i ] for i in indices ] )
    cat = cat2
  
  if ( not new_header_row is None ):
    if ( len( new_header_row ) == len( cat[ 0 ] ) ):
      cat[ 0 ] = [ x for x in new_header_row ]
    elif ( drop_strings and ( len( new_header_row ) == row_length ) ):
      cat[ 0 ] = [ new_header_row[ i ] for i in indices ]
    else:
      raise error( 'new header row does not match number of columns' )
  
  return cat

###############################################################################

def compress_atlas_dr3_catalog( cat, group_sources = None ):
  
  # build source component list
  src_list = []
  cmp_list = []
  index_list = []
  id_index = cat[ 0 ].index( 'ID' )
  for i in range( 1, len( cat ) ):
    c = cat[ i ]
    cname = c[ id_index ]
    cat_id = cname[ 0 : 2 ]
    if ( not cat_id in [ 'EI' ] ):
      raise error( 'unknown component ID format: %s' % ( cname ) )
    try:
      src_id = int( cname[ 2 : 6 ] )
    except:
      raise error( 'unknown component ID format: %s' % ( cname ) )
    cmp_id = 1
    if ( len( cname ) > 6 ):
      if ( cname[ 6 ] != 'C' ):
        raise error( 'unknown component ID format: %s' % ( cname ) )
      try:
        cmp_id = int( cname[ 7 : ] )
      except:
        raise error( 'unknown component ID format: %s' % ( cname ) )
    if ( not src_id in src_list ):
      src_list.append( src_id )
      cmp_list.append( [] )
      index_list.append( [] )
    src_index = src_list.index( src_id )
    if ( cmp_id in cmp_list ):
      raise error( 'double entry for source component %s' % ( cname ) )
    cmp_list[ src_index ].append( cmp_id )
    index_list[ src_index ].append( i )
  
  # group sources manually
  if ( not group_sources is None ):
    for group in group_sources:
      main_src_id = group[ 0 ]
      main_src_index = src_list.index( main_src_id )
      for src_id in group[ 1 : ]:
        src_index = src_list.index( src_id )
        src_list.remove( src_id )
        index_list[ main_src_index ] = index_list[ main_src_index ] + index_list[ src_index ]
        index_list.remove( index_list[ src_index ] )
  
  # merge components and build new cat
  new_cat = [ [ 'ID', 'name', 'RA', 'dRA', 'DEC', 'dDEC', 'Si', 'dSi', 'Sp', 'dSp',
      'noise', 'Ncomp', 'freq', 'SPIX', 'dSPIX' ] ]
  rms_index = cat[ 0 ].index( 'RMS' ) # mJy/beam
  sp_index = cat[ 0 ].index( 'Sp' ) # mJy
  dsp_index = cat[ 0 ].index( 'Sp_ERR' ) # mJy
  si_index = cat[ 0 ].index( 'S' ) # mJy
  dsi_index = cat[ 0 ].index( 'S_ERR' ) # mJy
  ra_index = cat[ 0 ].index( 'RA' )
  dra_index = cat[ 0 ].index( 'RA_ERR' ) # arcsec
  dec_index = cat[ 0 ].index( 'DEC' )
  ddec_index = cat[ 0 ].index( 'DEC_ERR' ) # arcsec
  freq_index = cat[ 0 ].index( 'OBS_FREQ' )
  spix_index = cat[ 0 ].index( 'SINDEX' )
  dspix_index = cat[ 0 ].index( 'INDEX_ERR' ) # arcsec
  for i in range( len( src_list ) ):
    ref_radec = None
    ra_array = None
    dec_array = None
    si_array = None
    sp = None
    dsp = None
    noise = None
    ncomp = 0
    freq = None
    spix = None
    dspix = None
    for j in range( len( index_list[ i ] ) ):
      ncomp = ncomp + 1
      c = cat[ index_list[ i ][ j ] ]
      if ( noise is None ):
        noise = c[ rms_index ] / 1.e3
      elif ( c[ rms_index ] / 1.e3 > noise ):
        noise = c[ rms_index ] / 1.e3
      if ( sp is None ):
        sp = c[ sp_index ] / 1.e3
        dsp = c[ dsp_index ] / 1.e3
      elif ( c[ sp_index ] / 1.e3 > sp ):
        sp = c[ sp_index ] / 1.e3
        dsp = c[ dsp_index ] / 1.e3
      if ( freq is None ):
        freq = c[ freq_index ] * 1.e6
      if ( spix is None ):
        spix = c[ spix_index ]
        dspix = c[ dspix_index ]
      si = draw_from_gaussian_distribution( c[ si_index ] / 1.e3, c[ dsi_index ] / 1.e3 )
      if ( si_array is None ):
        si_array = si.copy()
      else:
        si_array = si_array + si
      radec = [ c[ ra_index ], c[ dec_index ] ]
      dradec = [ c[ dra_index ] / 3600., c[ ddec_index ] / 3600. ]
      if ( ref_radec is None ):    
        ref_radec = radec
        radec = [ 0., 0. ]
      else:
        [ r, p ] = calculate_angular_separation( ref_radec, radec )
        radec = calculate_offset_position( [ 0., 0. ], r, p )
        radec[ 0 ] = amodulo( radec[ 0 ] + 180., 360. ) - 180.
      ra = draw_from_gaussian_distribution( radec[ 0 ], dradec[ 0 ] )
      dec = draw_from_gaussian_distribution( radec[ 1 ], dradec[ 1 ] )
      if ( ra_array is None ):
        ra_array = si * ra
        dec_array = si * dec
      else:
        ra_array = ra_array + si * ra
        dec_array = dec_array + si * dec
    si = get_robust_mean_deviations( si_array )
    dsi = max( si[ 1 ], -si[ 2 ] )
    si = si[ 0 ]
    ra = get_robust_mean_deviations( ra_array / si_array )
    dra = max( ra[ 1 ], -ra[ 2 ] )
    ra = ra[ 0 ]
    dec = get_robust_mean_deviations( dec_array / si_array )
    ddec = max( dec[ 1 ], -dec[ 2 ] )
    dec = dec[ 0 ]
    [ r, p ] = calculate_angular_separation( [ 0., 0. ], [ ra, dec ] )
    radec = radec = calculate_offset_position( ref_radec, r, p )
    src_id = 'EI%04d' % ( src_list[ i ] )
    name = 'ATLASDR3_J' + radec_to_string( radec, decimals = [ 3, 3 ], 
        separators = [ '', '', '', '', '', '' ] ) + '_I'
    new_cat.append( [ src_id, name, radec[ 0 ], dra, radec[ 1 ], ddec, si, dsi,
        sp, dsp, noise, ncomp, freq, spix, dspix ] )
  
  return new_cat

###############################################################################

def read_fits_table( fits_table_file_name, table_id = 1 ):
  fits_table = pf.open( fits_table_file_name )
  data = fits_table[ table_id ].data
  fits_table.close()
  column_names = list( data.dtype.names )
  column_types = [ data.dtype[ i ].name for i in range( len( column_names ) ) ]
  column_list = []
  for column_name in column_names:
    column_list.append( data[ column_name ].tolist() )
  fits_table = [ column_names ] + list( map( lambda *a: list( a ), *column_list ) )
  return fits_table

###############################################################################

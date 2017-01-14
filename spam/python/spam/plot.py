###############################################################################

# import Python modules
from sys import *
from os import *
from datetime import *
from math import *
from thread import *

# import 3rd party modules
from numpy import *
from pylab import *
#from matplotlib import *
#from RO.DS9 import *

# import user modules
from error import *
from files import *
from aips import *
from skymodel import *
from acalc import *
from ionosphere import *

###############################################################################

def plot_image( im, xy_range = None, z_range = None, file_name = None, dpi = 100,
    xyz_labels = None, plot_title = None, marker_list = None ):

  pixels = transpose( get_image_pixels( im ) )
  if ( xy_range is None ):
    im_size = get_image_size( im )
    im_pixel_ref = get_pixel_reference( im )
    im_xy_range = [ - im_pixel_ref[ 0 ] + 0.5, im_size[ 0 ] - im_pixel_ref[ 0 ] + 0.5,
                    - im_pixel_ref[ 1 ] + 0.5, im_size[ 1 ] - im_pixel_ref[ 1 ] + 0.5 ]
  else:
    im_xy_range = xy_range
  if ( z_range is None ):
    im_min = get_image_minimum( im )
    im_max = get_image_maximum( im )
    im_z_range = [ im_min[ 0 ], im_max[ 0 ] ]
  else:
    im_z_range = z_range
  if ( xyz_labels is None ):
    xyz_labels = [ 'pixels', 'pixels', 'Jy' ]

  clf()
  imshow( pixels, vmin = im_z_range[ 0 ], vmax = im_z_range[ 1 ], extent = im_xy_range,
      origin = 'lower', interpolation = 'nearest' )
  axis( im_xy_range, 'scaled' )
  xlabel( xyz_labels[ 0 ] )
  ylabel( xyz_labels[ 1 ] )
  if ( not plot_title is None ):
    title( plot_title )
  hot()
  colorbar()
  hold( True )

  if ( not marker_list is None ):
    for marker in marker_list:
      [ x, y, r, s ] = marker
      plot( [ x - im_pixel_ref[ 0 ] ], [ y - im_pixel_ref[ 1 ] ], s,
          markeredgewidth = 2., markersize = 2. * r, markeredgecolor = s[ 0 : 1 ], 
#          markerfacecolor = None ) # can be uncommented with new version of matplotlib/gtk
          markerfacecolor = s[ 0 : 1 ] )
      axis( im_xy_range, 'scaled' )

  if ( file_name is None ):
#    hold( True )
#    plot_process_id = start_new_thread( show, () )
    show()
  else:
    savefig( file_name, dpi = dpi )
  return

###############################################################################

def make_movie_from_images( image_file_name_list, movie_file_name, format = 'avi',
    frames_per_second = 10 ):

# This routine requires mencoder and ImageMagick's convert to be installed

  # currently only 10 fps supported
  if ( frames_per_second != 10 ):
    raise error( 'currently only 10 fps supported' )

  # currently only png input images supported
  extension = image_file_name_list[ 0 ][ - 4 : ]
  if ( extension.lower() != '.png' ):
    raise error( 'currently only PNG format images supported' )

  # check output format
  # check output filename extension
  if ( ( format == 'avi' ) or ( format == 'mpg' ) or ( format == 'gif' ) ):
    extension = movie_file_name[ - 4 : ]
    if ( extension.lower() == '.' + format ):
      output_file_name = movie_file_name
    else:
      output_file_name = movie_file_name + '.' + format
  else:
    raise error( 'unknown output format: ' + repr( format ) )

  # make consecutive series of images
  j = 0
  for image_file_name in image_file_name_list:
    if ( format == 'mpg' ):
      # copy every image 3 times
      copy_file( image_file_name, output_file_name + '_temp_%05d.png' % ( j ) )
      j = j + 1
      copy_file( image_file_name, output_file_name + '_temp_%05d.png' % ( j ) )
      j = j + 1
      copy_file( image_file_name, output_file_name + '_temp_%05d.png' % ( j ) )
      j = j + 1
    else:
      copy_file( image_file_name, output_file_name + '_temp_%05d.png' % ( j ) )
      j = j + 1
  image_count = j

  # combine separate frame files into one animation
  if file_exists( output_file_name ):
    remove_file( output_file_name )
  if ( format == 'avi' ):
    system( 'mencoder "mf://' + output_file_name + '_temp_*.png" -mf type=png:fps=' + repr( int( frames_per_second ) ) +
        ' -really-quiet -ovc lavc -ffourcc DX50 -noskip -oac copy -o ' + output_file_name )
#  elif ( format == 'mpg' ):
#    system( 'mencoder "mf://' + output_file_name + '_temp_*.png" -mf type=png:fps=' + 
#        repr( int( 3. * frames_per_second ) ) + ' -really-quiet -of mpeg -mpegopts ' +
#        'format=mpeg1:tsaf -ovc lavc -lavcopts vcodec=mpeg1video:vmax_b_frames=0:vmax_p_frames=8 ' +
#        '-oac lavc -lavcopts acodec=mp2 -o ' + output_file_name +
#        ' -ofps ' + repr( int( 3. * frames_per_second ) ) )
  elif ( format == 'mpg' ):
#    system( 'mencoder "mf://' + output_file_name + '_temp_*.png" -mf type=png:fps=' + 
#        repr( int( 3. * frames_per_second ) ) + ' -really-quiet -of mpeg -mpegopts ' +
#        'format=mpeg1:tsaf -ovc lavc -lavcopts vcodec=mpeg1video:vmax_b_frames=0 ' +
#        '-oac copy -o ' + output_file_name +
#        ' -ofps ' + repr( int( 3. * frames_per_second ) ) )
    system( 'mencoder "mf://' + output_file_name + '_temp_*.png" -mf type=png:fps=' + 
        repr( int( 3. * frames_per_second ) ) + ' -really-quiet -of mpeg -mpegopts ' +
        'format=mpeg2:tsaf -ovc lavc -lavcopts vcodec=mpeg2video:vmax_b_frames=0 ' +
        '-oac copy -o ' + output_file_name +
        ' -ofps ' + repr( int( 3. * frames_per_second ) ) )
  else: # ( format == 'gif' ):
    system( '/usr/bin/convert -loop 0 -delay ' + repr( int( 100. / float( frames_per_second ) ) ) + ' ' + 
        output_file_name + '_temp_*.png ' + output_file_name )

  # delete consecutive image files
  for j in range( image_count ):
    image_file_name = output_file_name + '_temp_%05d.png' % ( j )
    remove_file( image_file_name )

  return

###############################################################################

def plot_source( facets, radec, z_range = None ):
  [ i, pos ] = find_source_facets( facets, radec, primary_facet_only = True )[ 0 ]
  facet = get_facet( facets, i )
  plot_image( facet, z_range = z_range, marker_list = [ pos + [ 15., 'bo' ] ] )
  return

###############################################################################

def plot_fit_errors( uv, version = 0 ):
  ni_table = wizardry( uv ).table( 'NI', version )
  time_list = []
  weight_list = []
  for row in ni_table:
    time_list.append( row.time )
    weight_list.append( row.weight )
  time_array = array( time_list )
  weight_array = array( weight_list )
  sel = awhere( weight_array <= 0. )
  weight_array = aput( weight_array, sel, 1. )
  error_array = 1. / weight_array
  error_array = aput( error_array, sel, 0. )
  plot( time_array, error_array, 'k.' )
  show()
  return

###############################################################################

def make_mkl_phase_screen_images( uv, antenna, facets = None, facet_list = [],
    time_steps = [], print_info = False, field_factor = 1., grid_size = 101,
    fit_version = 0, image_prefix = '${PRT}/' ):
  
  # define grid
  radec = get_radec( uv )
  field_size = field_factor * restore_parameter( uv, 'field_size' )
  cell_size = field_size / float( grid_size )
  radec_array = zeros( ( grid_size, grid_size, 2 ), dtype = float64 )
  for y in range( grid_size ):
    for x in range( grid_size ):
      xx = float( x ) - ( float( grid_size - 1 ) / 2. )
      yy = float( y ) - ( float( grid_size - 1 ) / 2. )
      r = sqrt( xx**2 + yy**2 ) * cell_size
      if ( r < field_size / 2. ):
        phi = degrees( atan2( -xx, yy ) )
        radec_array[ y, x, 0 : 2 ] = array( 
            calculate_offset_position( radec, r, phi ), dtype = float64 )
  grid_sel = awhere( alltrue( radec_array[ : , : ] != array( [ 0., 0. ],
      dtype = float64 ), axis = -1 ) )
  radec_table = aget( radec_array, grid_sel ).reshape( len( grid_sel ), 2 ).tolist()
  grid_max = ( float( grid_size ) / 2. ) * cell_size
  grid_min = - grid_max

  # read model fit information
  wiz_uv = wizardry( uv )
  ni_table = wiz_uv.table( 'NI', fit_version )
  model = ni_table.keywords[ 'MODEL' ]
  if ( type( model ) == type( [] ) ):
    model = model[ 0 ]
  model = model.strip()
  if ( model != 'mkl' ):
    raise error( 'unknown model: %s' % ( model ) )
  order = ni_table.keywords[ 'NUM_COEF' ]
  reference_frequency = float32( ni_table.keywords[ 'REF_FREQ' ] )
  beta = float32( ni_table.keywords[ 'BETA' ] )
  r_0 = float32( ni_table.keywords[ 'R_0' ] )
  layer_count = int( ni_table.keywords[ 'LAYERS' ] )
  layer_heights = []
  layer_weights = []
  for l in range( layer_count ):
    layer_heights.append( float32( ni_table.keywords[ 'HEIGHT%d' % ( l + 1 ) ] ) )
    layer_weights.append( float32( ni_table.keywords[ 'WEIGHT%d' % ( l + 1 ) ] ) )
  if ( not facets is None ):
    if ( len( facet_list ) == 0 ):
      used_facet_list = range( 1, 1 + restore_parameter( facets, 'facet_count' ) )
    else:
      used_facet_list = facet_list
  else:
    used_facet_list = []
  ob_table = wiz_uv.table( 'OB', ni_table.version )
  try:
    iterations = ob_table.keywords[ 'ITER' ]
  except:
    iterations = 4
  
  # read ionospheric fit table
  fit_time_list = []
  fit_count_list = []
  fit_coef_table = []
  fit_weight_list = []
  for row in ni_table:
    fit_time_list.append( float32( row.time ) )
    fit_count_list.append( row.antenna_no )
    fit_coef_table.append( [ float32( coef ) for coef in row.coef ] )
    fit_weight_list.append( float32( row.weight ) )
  fit_time_count = len( fit_time_list )
  fit_coef_table = array( fit_coef_table, dtype = float64 )
  
  # remove bad fits
  fit_weight_array = array( fit_weight_list, dtype = float64 )
  sel = awhere( fit_weight_array > 0. )
  sel = sel.ravel().tolist()
  
  # read ionospheric pierce point table
  pierce_time_list = []
  pierce_X_table = []
  pierce_za_list = []
  pierce_index_table = []
  for row in ob_table:
    pierce_time_list.append( float32( row.time ) )
    pierce_X_table.append( [ float32( x ) for x in row.orbxyz[ 0 : 2 ] ] )
    pierce_za_list.append( float32( row.orbxyz[ 2 ] ) )
    pierce_index_table.append( [ int( round( row.orientation ) ) - 1, 
        row.subarray - 1, row.antenna_no - 1 ] )
  pierce_time_count = len( pierce_time_list )
  pierce_time_array = array( pierce_time_list, dtype = float64 )
  pierce_index_array = array( pierce_index_table, dtype = int64 )
  pierce_X_array = array( pierce_X_table, dtype = float64 )
  pierce_za_array = array( pierce_za_list, dtype = float64 )
  
  # get other relevant data
  calibration_data = get_phase_calibration_data( uv, facets,
      time_info = True, source_info = True, antenna_info = True,
      calibration_info = False, facet_list = used_facet_list,
      print_info = print_info )
  time_table = calibration_data[ 0 ]
  center_table = calibration_data[ 1 ]
  source_table = calibration_data[ 2 ]
  array_table = calibration_data[ 3 ]
  antenna_table = calibration_data[ 4 ]
  source_count = len( source_table )
  antenna_count = len( antenna_table )
  time_count = len( time_table )
  
  # generate time table
  fit_time_table = []
  fit_gst_list = get_gst_list( uv, fit_time_list )
  for n in range( fit_time_count ):
    fit_time_table.append( [ fit_time_list[ n ], fit_gst_list[ n ] ] )

  # calculate source positions in plot
  dsource_table = []
  for source in source_table:
    [ r, phi ] = calculate_angular_separation( radec, source )
    dsource_table.append( [ -r * sin( radians( phi ) ), r * cos( radians( phi ) ) ] )
  dsource_array = array( dsource_table, dtype = float64 )
#???  
#  # generate reference antenna table
#  reference_list = []
#  for k in range( source_count ):
#    reference_list.append( r )

  # loop over time stamps
  for nn in range( time_count ):

    if ( len( time_steps ) > 0 ):
      if ( not nn in time_steps ):
        continue

    if print_info:
      print '... time step n = ', nn
    
    # save fit weight
    fit_available = True
    try:
      n = fit_time_list.index( time_table[ nn ][ 0 ] )
    except ValueError:
      fit_available = False
    
    # only process non-rejected fits
    if ( not n in sel ):
      fit_available = False

    weight = 1. / 360.
    if fit_available:
    
      # get model fit parameters
      P = fit_coef_table[ n ]
      weight = fit_weight_list[ n ]
      
      if print_info:
        print '...... calculating base vectors'
      
      active_antennas = []
      Xpl_table = []
      pzal_table = []
      Bl_table = []
      for l in range( layer_count ):
        Xp_table = []
        pza_table = []
      
        # get pierce points 
        sel2 = awhere( ( pierce_time_array == fit_time_list[ n ] ) & 
            ( pierce_index_array[ : , 0 ] == l ) )
        Xp_table = Xp_table + aget( pierce_X_array, sel2 ).tolist()
        pza_table = pza_table + aget( pierce_za_array, sel2 ).tolist()
        for i in aget( pierce_index_array[ : , 2 ], sel2 ).tolist():
          if ( not i in active_antennas ):
             active_antennas.append( i )
  
        Xp_table = array( Xp_table, dtype = float64 )
        pza_table = array( pza_table, dtype = float64 )
        
        # calculate structure matrix
        p_count = len( Xp_table )
        if ( p_count != fit_count_list[ n ] ):
          raise error( 'pierce count does not match' )
        D_table = resize( Xp_table, ( p_count, p_count, 2 ) )
        D_table = transpose( D_table, ( 1, 0, 2 ) ) - D_table
        D_table = add.reduce( D_table**2, 2 )
        D_table = ( D_table / ( r_0**2 ) )**( beta / 2. )
        
        # calculate covariance matrix C
        # calculate partial product for interpolation B
        # reforce symmetry
        C_table = - D_table / 2.
        C_table = transpose( C_table - ( add.reduce( C_table, 0 ) / 
            float( p_count ) ) )
        B_table = add.reduce( C_table, 0 ) / float( p_count )
        C_table = C_table - B_table
        C_table = ( C_table + transpose( C_table ) ) / 2.
        
        # incorporate airmass functions and layer weights
        # add layer C to total C
        A_table = resize( 1. / cos( aradians( pza_table ) ), ( p_count, p_count ) )
        A_table = A_table * transpose( A_table )
        if ( l == 0 ):
          Cl_table = C_table * A_table * ( layer_weights[ l ]**2 )
        else:
          Cl_table = Cl_table + C_table * A_table * ( layer_weights[ l ]**2 )
        
        # save tables per height layer
        Xpl_table.append( Xp_table )
        pzal_table.append( pza_table )
        Bl_table.append( B_table )
      
      # convert to arrays
      Xpl_table = array( Xpl_table, dtype = float64 )
      pzal_table = array( pzal_table, dtype = float64 )
      Bl_table = array( Bl_table, dtype = float64 )
      
      # eigenvalue decomposition
      # reforce symmetry
      # select subset of base vectors
      [ U_table, S, Ut_table ] = linalg.svd( Cl_table )
      U_table = ( U_table + transpose( Ut_table ) ) / 2.
      U_table = U_table[ : , 0 : order ]
      S = S[ 0 : order ]
      
      # calculate interpolation matrix
      F_table = dot( U_table, P / S )
      
      if print_info:
        print '...... calculating pierce point coordinates'
    
      Xl_table = []
      zal_table = []
      ref_list = []
      ref_table = []
      for l in range( layer_count ):
        
        X_table = []
        za_table = []
        
        # get pierce point coordinates
        pierce_table = calculate_pierce_coordinates( fit_time_table[ n ], 
            center_table, radec_table, array_table,
            antenna_table[ antenna - 1 : antenna ], 
            height = layer_heights[ l ], iterations = iterations )
        
        # put all new pierce points into one array
        j = 0
        for pierce_info in pierce_table:
          [ X, za, [ k, i ] ] = pierce_info
          X_table.append( X )
          za_table.append( za )
        Xl_table.append( X_table )
        zal_table.append( za_table )
        
      Xl_table = array( Xl_table, dtype = float64 )
      zal_table = array( zal_table, dtype = float64 )
      
      if print_info:
        print '...... generating solutions'
      
      # calculate pierce point model solutions
      phi_table = phi_mkl_model( layer_weights, Xl_table, zal_table, Xpl_table, 
          pzal_table, Bl_table, F_table, beta = beta, r_0 = r_0 )

    if print_info:
      print '...... generating plot'
    phi_grid = zeros( ( grid_size, grid_size ), dtype = float64 )
    if fit_available:
      phi_grid = aput( phi_grid, grid_sel, phi_table )
      phi_grid = amodulo( phi_grid + 180., 360. ) - 180.

    # remove gradient
    if False:
      l = layer_weights.index( max( layer_weights ) )
      Xp_table = Xpl_table[ l ] / cos( aradians( pzal_table[ l ] ) ).reshape( 
        ( p_count, 1 ) )
      F_table = dot( U_table, P / S )
      phi_offset = phi_mkl_model( layer_weights,
          azeros( Xpl_table[ : , 0 : 1, : ] ), aones( pzal_table[ : , 0 : 1 ] ),
          Xpl_table, pzal_table, Bl_table, F_table, beta = beta, r_0 = r_0 )
      T1 = linalg.inv( dot( transpose( Xp_table ), Xp_table ) )
      T2 = dot( transpose( Xp_table ), dot( U_table, P ) - phi_offset )
      G = dot( T1, T2 )

    # plot phase grid below measured phases
    clf()
    imshow( phi_grid, extent = ( grid_min, grid_max, grid_min, grid_max ),
        interpolation = 'nearest', vmin = -180., vmax = 180., origin = 'lower' )
    axis( [ grid_min, grid_max, grid_min, grid_max ], 'scaled' )
    hsv()
    hold( True )

    plot( dsource_array[ : , 0 ], dsource_array[ : , 1 ], 'k+' )
    axis( [ grid_min, grid_max, grid_min, grid_max ], 'scaled' )
    hsv()
    xlabel( r'$\Delta$RA [deg]' )
    ylabel( r'$\Delta$DEC [deg]' )
    title( r'ant %2d, n=%3d (%s), $\sigma_{\phi}$=%4.1f deg' % (
        antenna, nn, time_to_string( time_table[ nn ][ 0 ] ), 1. / weight ) )
    cb = colorbar()
    cb.ax.set_ylabel( r'phase [deg]' )

    # save plot to image file
#    show()
    image_prefix_e = path.expandvars( image_prefix )
    image_file_name = image_prefix + '%05d.png' % ( nn )
    if file_exists( image_file_name ):
      remove_file( image_file_name )
    savefig( image_file_name, dpi = 75 )
  
  return

###############################################################################

def make_pmkl_phase_screen_images( uv, antenna, facets = None, facet_list = [],
    time_steps = [], print_info = False, field_factor = 1., grid_size = 101,
    fit_version = 0, image_prefix = '${PRT}/', zero_average = False,
    plot_structure = [ True, True, True ] ):
# plot_structure = [ gradient, 2nd order, KL ]
  
  # define grid
  radec = get_radec( uv )
  field_size = field_factor * restore_parameter( uv, 'field_size' )
  cell_size = field_size / float( grid_size )
  radec_array = zeros( ( grid_size, grid_size, 2 ), dtype = float64 )
  for y in range( grid_size ):
    for x in range( grid_size ):
      xx = float( x ) - ( float( grid_size - 1 ) / 2. )
      yy = float( y ) - ( float( grid_size - 1 ) / 2. )
      r = sqrt( xx**2 + yy**2 ) * cell_size
      if ( r < field_size / 2. ):
        phi = degrees( atan2( -xx, yy ) )
        radec_array[ y, x, 0 : 2 ] = array( 
            calculate_offset_position( radec, r, phi ), dtype = float64 )
  grid_sel = awhere( alltrue( radec_array[ : , : ] != array( [ 0., 0. ],
      dtype = float64 ), axis = -1 ) )
  radec_table = aget( radec_array, grid_sel ).reshape( len( grid_sel ), 2 ).tolist()
  grid_max = ( float( grid_size ) / 2. ) * cell_size
  grid_min = - grid_max
  
  # read model fit information
  wiz_uv = wizardry( uv )
  ni_table = wiz_uv.table( 'NI', fit_version )
  model = ni_table.keywords[ 'MODEL' ]
  if ( type( model ) == type( [] ) ):
    model = model[ 0 ]
  model = model.strip()
  if ( model != 'pmkl' ):
    raise error( 'unknown model: %s' % ( model ) )
  order = ni_table.keywords[ 'NUM_COEF' ] - 5
  reference_frequency = float32( ni_table.keywords[ 'REF_FREQ' ] )
  beta = float32( ni_table.keywords[ 'BETA' ] )
  r_0 = float32( ni_table.keywords[ 'R_0' ] )
  layer_count = int( ni_table.keywords[ 'LAYERS' ] )
  layer_heights = []
  layer_weights = []
  for l in range( layer_count ):
    layer_heights.append( float32( ni_table.keywords[ 'HEIGHT%d' % ( l + 1 ) ] ) )
    layer_weights.append( float32( ni_table.keywords[ 'WEIGHT%d' % ( l + 1 ) ] ) )
  l_max = layer_weights.index( max( layer_weights ) )
  if ( not facets is None ):
    if ( len( facet_list ) == 0 ):
      used_facet_list = range( 1, 1 + restore_parameter( facets, 'facet_count' ) )
    else:
      used_facet_list = facet_list
  else:
    used_facet_list = []
  ob_table = wiz_uv.table( 'OB', ni_table.version )
  try:
    iterations = ob_table.keywords[ 'ITER' ]
  except:
    iterations = 4
  
  # read ionospheric fit table
  fit_time_list = []
  fit_count_list = []
  fit_poly_table = []
  fit_coef_table = []
  fit_weight_list = []
  for row in ni_table:
    fit_time_list.append( float32( row.time ) )
    fit_count_list.append( row.antenna_no )
    fit_poly_table.append( [ float32( pol ) for pol in row.coef[ 0 : 5 ] ] )
    fit_coef_table.append( [ float32( coef ) for coef in row.coef[ 5 : ] ] )
    fit_weight_list.append( float32( row.weight ) )
  fit_time_count = len( fit_time_list )
  fit_poly_table = array( fit_poly_table, dtype = float64 )
  fit_coef_table = array( fit_coef_table, dtype = float64 )
  
  # remove bad fits
  fit_weight_array = array( fit_weight_list, dtype = float64 )
  sel = awhere( fit_weight_array > 0. )
  sel = sel.ravel().tolist()
  
  # read ionospheric pierce point table
  pierce_time_list = []
  pierce_X_table = []
  pierce_za_list = []
  pierce_index_table = []
  for row in ob_table:
    pierce_time_list.append( float32( row.time ) )
    pierce_X_table.append( [ float32( x ) for x in row.orbxyz[ 0 : 2 ] ] )
    pierce_za_list.append( float32( row.orbxyz[ 2 ] ) )
    pierce_index_table.append( [ int( round( row.orientation ) ) - 1, 
        row.subarray - 1, row.antenna_no - 1 ] )
  pierce_time_count = len( pierce_time_list )
  pierce_time_array = array( pierce_time_list, dtype = float64 )
  pierce_index_array = array( pierce_index_table, dtype = int64 )
  pierce_X_array = array( pierce_X_table, dtype = float64 )
  pierce_za_array = array( pierce_za_list, dtype = float64 )
  
  # get other relevant data
  calibration_data = get_phase_calibration_data( uv, facets,
      time_info = True, source_info = True, antenna_info = True,
      calibration_info = False, facet_list = used_facet_list,
      print_info = print_info )
  time_table = calibration_data[ 0 ]
  center_table = calibration_data[ 1 ]
  source_table = calibration_data[ 2 ]
  array_table = calibration_data[ 3 ]
  antenna_table = calibration_data[ 4 ]
  source_count = len( source_table )
  antenna_count = len( antenna_table )
  time_count = len( time_table )
  
  # generate time table
  fit_time_table = []
  fit_gst_list = get_gst_list( uv, fit_time_list )
  for n in range( fit_time_count ):
    fit_time_table.append( [ fit_time_list[ n ], fit_gst_list[ n ] ] )
  
  # calculate source positions in plot
  dsource_table = []
  for source in source_table:
    [ r, phi ] = calculate_angular_separation( radec, source )
    dsource_table.append( [ -r * sin( radians( phi ) ), r * cos( radians( phi ) ) ] )
  dsource_array = array( dsource_table, dtype = float64 )
#???  
#  # generate reference antenna table
#  reference_list = []
#  for k in range( source_count ):
#    reference_list.append( r )
  
  # loop over time stamps
  for nn in range( time_count ):

    if ( len( time_steps ) > 0 ):
      if ( not nn in time_steps ):
        continue
    
    if print_info:
      print '... time step n = ', nn
    
    # save fit weight
    fit_available = True
    try:
      n = fit_time_list.index( time_table[ nn ][ 0 ] )
    except ValueError:
      fit_available = False
    
    # only process non-rejected fits
    if ( not n in sel ):
      fit_available = False
    
    weight = 1. / 360.
    if fit_available:
      
      # get model fit parameters
      poly = fit_poly_table[ n ]
      P = fit_coef_table[ n ]
      weight = fit_weight_list[ n ]
      
      if print_info:
        print '...... calculating base vectors'
      
      active_antennas = []
      Xpl_table = []
      pzal_table = []
      Bl_table = []
      for l in range( layer_count ):
        Xp_table = []
        pza_table = []
        
        # get pierce points 
        sel2 = awhere( ( pierce_time_array == fit_time_list[ n ] ) & 
            ( pierce_index_array[ : , 0 ] == l ) )
        Xp_table = Xp_table + aget( pierce_X_array, sel2 ).tolist()
        pza_table = pza_table + aget( pierce_za_array, sel2 ).tolist()
        for i in aget( pierce_index_array[ : , 2 ], sel2 ).tolist():
          if ( not i in active_antennas ):
             active_antennas.append( i )
        
        Xp_table = array( Xp_table, dtype = float64 )
        pza_table = array( pza_table, dtype = float64 )
        
        # calculate structure matrix
        p_count = len( Xp_table )
        if ( p_count != fit_count_list[ n ] ):
          raise error( 'pierce count does not match' )
        if ( plot_structure[ 2 ] ):
          D_table = resize( Xp_table, ( p_count, p_count, 2 ) )
          D_table = transpose( D_table, ( 1, 0, 2 ) ) - D_table
          D_table = add.reduce( D_table**2, 2 )
          D_table = ( D_table / ( r_0**2 ) )**( beta / 2. )
          
          # calculate covariance matrix C
          # calculate partial product for interpolation B
          # reforce symmetry
          C_table = - D_table / 2.
          C_table = transpose( C_table - ( add.reduce( C_table, 0 ) / 
              float( p_count ) ) )
          B_table = add.reduce( C_table, 0 ) / float( p_count )
          C_table = C_table - B_table
          C_table = ( C_table + transpose( C_table ) ) / 2.
          
          # incorporate airmass functions and layer weights
          # add layer C to total C
          A_table = resize( 1. / cos( aradians( pza_table ) ), ( p_count, p_count ) )
          A_table = A_table * transpose( A_table )
          if ( l == 0 ):
            Cl_table = C_table * A_table * ( layer_weights[ l ]**2 )
          else:
            Cl_table = Cl_table + C_table * A_table * ( layer_weights[ l ]**2 )
        
        # save tables per height layer
        Xpl_table.append( Xp_table )
        pzal_table.append( pza_table )
        if ( plot_structure[ 2 ] ):
          Bl_table.append( B_table )
      
      # convert to arrays
      Xpl_table = array( Xpl_table, dtype = float64 )
      pzal_table = array( pzal_table, dtype = float64 )
      if ( plot_structure[ 2 ] ):
        Bl_table = array( Bl_table, dtype = float64 )
        
        # eigenvalue decomposition
        # reforce symmetry
        # select subset of base vectors
        [ U_table, S, Ut_table ] = linalg.svd( Cl_table )
        U_table = ( U_table + transpose( Ut_table ) ) / 2.
        U_table = U_table[ : , 0 : order ]
        S = S[ 0 : order ]
        
        # calculate interpolation matrix
        F_table = dot( U_table, P / S )
      
      if print_info:
        print '...... calculating pierce point coordinates'
      
      Xl_table = []
      zal_table = []
      ref_list = []
      ref_table = []
      for l in range( layer_count ):
        
        X_table = []
        za_table = []
        
        # get pierce point coordinates
        pierce_table = calculate_pierce_coordinates( fit_time_table[ n ], 
            center_table, radec_table, array_table,
            antenna_table[ antenna - 1 : antenna ], 
            height = layer_heights[ l ], iterations = iterations )
        
        # put all new pierce points into one array
        j = 0
        for pierce_info in pierce_table:
          [ X, za, [ k, i ] ] = pierce_info
          X_table.append( X )
          za_table.append( za )
        Xl_table.append( X_table )
        zal_table.append( za_table )
      
      Xl_table = array( Xl_table, dtype = float64 )
      zal_table = array( zal_table, dtype = float64 )
      
      if print_info:
        print '...... generating solutions'
      
      # calculate pierce point model solutions
      if ( not plot_structure[ 0 ] ):
        poly[ 0 : 2 ] = [ 0., 0. ]
      if ( not plot_structure[ 1 ] ):
        poly[ 2 : 5 ] = [ 0., 0., 0. ]
      phi_table = phi_poly_model( Xl_table[ l_max ], poly ) / cos( 
          aradians( zal_table[ l_max ] ) )
      if ( plot_structure[ 2 ] ):
        phi_mkl_table = phi_mkl_model( layer_weights, Xl_table, zal_table,
            Xpl_table, pzal_table, Bl_table, F_table, beta = beta, r_0 = r_0 )
        phi_table = phi_table + phi_mkl_table
    if zero_average:
      phi_table = phi_table - phi_table.mean()
    
    if print_info:
      print '...... generating plot'
    phi_grid = zeros( ( grid_size, grid_size ), dtype = float64 )
    if fit_available:
      phi_grid = aput( phi_grid, grid_sel, phi_table )
      phi_grid = amodulo( phi_grid + 180., 360. ) - 180.
    
    # plot phase grid below measured phases
    clf()
    imshow( phi_grid, extent = ( grid_min, grid_max, grid_min, grid_max ),
        interpolation = 'nearest', vmin = -180., vmax = 180., origin = 'lower' )
    axis( [ grid_min, grid_max, grid_min, grid_max ], 'scaled' )
    hsv()
    hold( True )
    
    plot( dsource_array[ : , 0 ], dsource_array[ : , 1 ], 'k+' )
    axis( [ grid_min, grid_max, grid_min, grid_max ], 'scaled' )
    hsv()
    xlabel( r'$\Delta$RA [deg]' )
    ylabel( r'$\Delta$DEC [deg]' )
    title( r'ant %2d, n=%3d (%s), $\sigma_{\phi}$=%4.1f deg' % (
        antenna, nn, time_to_string( time_table[ nn ][ 0 ] ), 1. / weight ) )
    cb = colorbar()
    cb.ax.set_ylabel( r'phase [deg]' )
    
    # save plot to image file
#    show()
    image_prefix_e = path.expandvars( image_prefix )
    image_file_name = image_prefix_e + '%05d.png' % ( nn )
    if file_exists( image_file_name ):
      remove_file( image_file_name )
    savefig( image_file_name, dpi = 75 )
  
  return

###############################################################################

def plot_overlap_region( frequency = 153.e6, radius = 1.6, offset1 = [ 0, 2 ],
    separation = 3.2 / sqrt( 3. ), amplitudes = [ 1., 1. ], offset2 = [ 0, 0 ],
    pb_parameters = [ 0., 1., -4.04, 76.2, -68.8, 22.03, 0. ] ):

  cell_size = radius / 100.
  x_limit = int( ceil( radius / cell_size ) )
  y_limit = int( ceil( ( radius + ( separation / 2. ) ) / cell_size ) )
  x1 = x2 = x_limit
  x1a = x_limit + offset1[ 0 ]
  x2a = x_limit + offset2[ 0 ]
  y_offset = int( ceil( ( separation / 2. ) / cell_size ) )
  y1 = y_limit - y_offset
  y2 = y_limit + y_offset
  y1a = y_limit - y_offset + offset1[ 1 ]
  y2a = y_limit + y_offset + offset2[ 1 ]
  grid1 = zeros( [ 2 * x_limit + 1, 2 * y_limit + 1 ], dtype = float32 )
  grid2 = zeros( [ 2 * x_limit + 1, 2 * y_limit + 1 ], dtype = float32 )
  for i in range( 2 * x_limit + 1 ):
    for j in range( 2 * y_limit + 1 ):
      r1 = sqrt( float( ( i - x1 )**2 + ( j - y1 )**2 ) ) * cell_size
      if ( r1 < radius ):
        r1a = sqrt( float( ( i - x1a )**2 + ( j - y1a )**2 ) ) * cell_size
        A1 = calculate_pbparm_attenuation( frequency, r1, pb_parameters )
        A1a = calculate_pbparm_attenuation( frequency, r1a, pb_parameters )
        grid1[ i, j ] = A1a / A1
      r2 = sqrt( float( ( i - x2 )**2 + ( j - y2 )**2 ) ) * cell_size
      if ( r2 < radius ):
        r2a = sqrt( float( ( i - x2a )**2 + ( j - y2a )**2 ) ) * cell_size
        A2 = calculate_pbparm_attenuation( frequency, r2, pb_parameters )
        A2a = calculate_pbparm_attenuation( frequency, r2a, pb_parameters )
        grid2[ i, j ] = A2a / A2
  grid12 = grid1 + grid2
  sel = awhere( ( grid1 > 0. ) & (grid2 > 0. ) )
  grid12 = aput( grid12, sel, aget( grid1, sel ) / aget( grid2, sel ) )
  sel = awhere( grid12 == 0. )
  grid12 = aput( grid12, sel, 1. )
  imshow( transpose( grid12 ), origin = 'lower' ); colorbar(); show()

###############################################################################

def convert_facets_to_ds9_region_file( uv, blank_edge = 6, color = 'blue' ):
  
  # get relevant info
  facet_file_name = restore_parameter( uv, 'pb_facet_file_name' )
  facet_file_name = path.expandvars( facet_file_name )
  facet_list = get_facet_list( facet_file_name )
  pixel_size = restore_parameter( uv, 'cell_size' )
  facet_size = restore_parameter( uv, 'pb_facet_size' )
  facet_radius = ( ( facet_size / 2 ) - blank_edge - 1 ) * pixel_size # arcsec
  
  # create regions file
  regions_file_name = facet_file_name + '.facet.reg'
  if file_exists( regions_file_name ):
    remove_file( regions_file_name )
  regions_file = file( regions_file_name, mode = 'w' )
  regions_file.write( '# Region file format: DS9 version 4.1\n' )
  regions_file.write( 'global color=%s dashlist=8 3 width=1 ' % ( color ) +
      'font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 ' + 
      'edit=1 move=1 delete=1 include=1 source=1\nfk5\n' )
  for [ i, radec ] in facet_list:
    regions_file.write( 'circle(%s,%7.3f") # text={%s}\n' % (
        radec_to_string( box_radec, decimals = [ 3,2 ],
        separators = [ ':',':',',',':',':','' ] ), facet_radius, i ) )
  regions_file.close()
  
  return

###############################################################################

def convert_boxes_to_ds9_region_file( uv, color = 'green' ):
  
  # get relevant info
  facet_file_name = restore_parameter( uv, 'pb_facet_file_name' )
  facet_file_name = path.expandvars( facet_file_name )
  facet_list = get_facet_list( facet_file_name )
  clean_box_list = get_clean_boxes( facet_file_name )
  pixel_size = restore_parameter( uv, 'cell_size' )
  facet_size = restore_parameter( uv, 'pb_facet_size' )
  
  # create regions file
  regions_file_name = facet_file_name + '.box.reg'
  if file_exists( regions_file_name ):
    remove_file( regions_file_name )
  regions_file = file( regions_file_name, mode = 'w' )
  regions_file.write( '# Region file format: DS9 version 4.1\n' )
  regions_file.write( 'global color=%s dashlist=8 3 width=1 ' % ( color ) +
      'font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 ' + 
      'edit=1 move=1 delete=1 include=1 source=1\nfk5\n' )
  for [ i, radec ] in facet_list:
    for [ j, dx, dy, x, y ] in clean_box_list:
      if ( i != j ):
        continue
      x = - ( x - ( facet_size / 2 ) ) * pixel_size / 3600.
      y = ( y - ( 1 + facet_size / 2 ) ) * pixel_size / 3600.
      r = sqrt( x**2 + y**2 )
      p = adegrees( atan2( x, y ) )
      box_radec = calculate_offset_position( radec, r, p )
      if ( dx < 0 ):
        box_radius = dy * pixel_size
        regions_file.write( 'circle(%s,%7.3f")\n' % ( radec_to_string( box_radec,
            decimals = [ 3,2 ], separators = [ ':',':',',',':',':','' ] ), box_radius ) )
      else:
        dx = dx * pixel_size
        dy = dy * pixel_size
        regions_file.write( 'box(%s,%7.3f",%7.3f",0.0)\n' % ( radec_to_string( radec,
            decimals = [ 3,2 ], separators = [ ':',':',',',':',':','' ] ), dx, dy ) )
  regions_file.close()
  
  return

###############################################################################


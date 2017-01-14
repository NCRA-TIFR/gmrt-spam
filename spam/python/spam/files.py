###############################################################################

# import Python modules
from os import *
from shutil import *

# import user modules
from error import *
#from aips import *

###############################################################################

def file_exists( file_path ):
  return path.isfile( path.expandvars( file_path ) )

###############################################################################

def directory_exists( dir_path ):
  return path.isdir( path.expandvars( dir_path ) )

###############################################################################

def remove_file( file_path ):
  status = True
  if file_exists( file_path ):
    try:
      remove( path.expandvars( file_path ) )
    except:
      status = False
  return status

###############################################################################

def move_file( src_path, dest_path, overwrite = True ):
  status = True
  if ( ( not file_exists( src_path ) ) or ( ( overwrite == False ) and ( file_exists( dest_path ) ) ) ):
    status = False
  else:
    if ( ( overwrite == True ) and ( file_exists( dest_path ) ) ):
      status = remove_file( dest_path )
  if ( status == True ):
    try:
      move( path.expandvars( src_path ), path.expandvars( dest_path ) )
    except:
      status = False
  return status

###############################################################################

def copy_file( src_path, dest_path, overwrite = True ):
  status = True
  if ( ( not file_exists( src_path ) ) or ( ( overwrite == False ) and ( file_exists( dest_path ) ) ) ):
    status = False
  else:
    if ( ( overwrite == True ) and ( file_exists( dest_path ) ) ):
      status = remove_file( dest_path )
  if ( status == True ):
    try:
      copy( path.expandvars( src_path ), path.expandvars( dest_path ) )
    except:
      status = False
  return status

###############################################################################

def make_directory( dir_path ):
  try:
    mkdir( path.expandvars( dir_path ) )
  except OSError:
    raise error( 'creation of directory %s failed' % ( dir_path ) )
  return

###############################################################################

def remove_directory( dir_path ):
  try:
    rmdir( path.expandvars( dir_path ) )
  except OSError:
    raise error( 'removal of directory %s failed' % ( dir_path ) )
  return

###############################################################################

def make_file( file_path, overwrite = True ):
  if ( overwrite and file_exists( file_path ) ):
    remove_file( file_path )
  if ( not file_exists( file_path ) ):
    new_file = file( path.expandvars( file_path ), mode = 'w' )
    new_file.close()
  return

###############################################################################

def get_directory( dir_path ):
  if ( not directory_exists( dir_path ) ):
    raise error( 'directory %s does not exist' % ( dir_path ) )
  directory = listdir( path.expandvars( dir_path ) )
  return directory

###############################################################################

def read_fits_header( file_path, include_blanks = True ):
  if ( not file_exists( file_path ) ):
    raise error( 'fits file %s does not exist' % ( file_path ) )
  infile = file( path.expandvars( file_path ), mode = 'rb' )
  header = []
  while ( True ):
    line = infile.read( 80 )
    if ( ( not include_blanks ) and ( line.strip() == '' ) ):
      continue
    if ( line.strip() == 'END' ):
      break
    header.append( line.strip() )
  infile.close()
  return header

###############################################################################

def write_fits_header( infile_path, header, outfile_path = None,
    overwrite = True, block_size = 1000000L ):
  if ( not file_exists( infile_path ) ):
    raise error( 'source fits file %s does not exist' % ( infile_path ) )
  if ( not outfile_path is None ):
    if file_exists( outfile_path ):
      if overwrite:
        remove_file( outfile_path )
      else:
        raise error( 'destination fits file %s does already exist' % ( outfile_path ) )
    outfile_path_e = path.expandvars( outfile_path )
    replace = False
  else:
    outfile_path_e = path.expandvars( infile_path + '.temp' )
    replace = True
  # skip input header
  infile = file( path.expandvars( infile_path ), mode = 'rb' )
  end_found = False
  byte_count = 0L
  while ( True ):
    line = infile.read( 80 )
    byte_count = byte_count + 80L
    if ( not end_found ):
      if ( line.strip() == 'END' ):
        end_found = True
    if end_found:
      if ( byte_count % 2880L == 0L ):
        break
  # write output header
  outfile = file( outfile_path_e, mode = 'wb' )
  byte_count = 0L
  for line in header:
    outfile.write( line.strip().ljust( 80 )[ 0 : 80 ] )
    byte_count = byte_count + 80L
  line = 'END'
  outfile.write( line.strip().ljust( 80 )[ 0 : 80 ] )
  byte_count = byte_count + 80L
  line = ''
  while ( byte_count % 2880L != 0L ):
    outfile.write( line.strip().ljust( 80 )[ 0 : 80 ] )
    byte_count = byte_count + 80L
  outfile.flush()
  # copy binary part
  while ( True ):
    bin_data = infile.read( block_size )
    outfile.write( bin_data )
    outfile.flush()
    byte_count = byte_count + len( bin_data )
    if ( len( bin_data ) < block_size ):
      break
  infile.close()
  outfile.close()
  if replace:
    move_file( outfile_path_e, path.expandvars( infile_path ) )
  if ( byte_count % 2880L != 0L ):
    raise error( 'byte count %d in output fits file is not a multiple of 2880' %
        ( byte_count ) )
  return

###############################################################################

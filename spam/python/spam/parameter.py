###############################################################################

# import Python modules
from sys import *
from datetime import *

# import 3rd party modules
from numpy import *

# import user modules
from aips import *
from error import *

###############################################################################

def read_history( uvim, strings = [], count = 1, word_match = False ):
# strings = search terms
# count = 0: return all occurrences
# count = N: return first N occurrences
# count = -N: return last N occurrences in reverse order

  if not table_exists( uvim, 'HI', 1 ):
    raise error( 'AIPS HI table does not exist' )
  hi_table = uvim.history
  hi_list = []
  if ( count >= 0 ):
    found_count = 0
    for line in hi_table:
      if word_match:
        words = [ word.strip() for word in line.split() ]
      if ( len( strings ) == 0 ):
        hi_list.append( line )
        found_count = found_count + 1
      else:
        found = True
        for string in strings:
          if ( line.find( string ) == - 1 ):
            found = False
            break
          elif word_match:
            try:
              dummy = words.index( string )
            except:
              found = False
              break
        if found:
          hi_list.append( line )
          found_count = found_count + 1
      if ( found_count == count ):
        break
  else:
    history_size = get_history_size( uvim )
    found_count = 0
    for index in range( history_size - 1, - 1, - 1 ):
      line = hi_table[ index ]
      if word_match:
        words = [ word.strip() for word in line.split() ]
      if ( len( strings ) == 0 ):
        hi_list.append( line )
        found_count = found_count + 1
      else:
        found = True
        for string in strings:
          if ( line.find( string ) == - 1 ):
            found = False
            break
          elif word_match:
            try:
              dummy = words.index( string )
            except:
              found = False
              break
        if found:
            hi_list.append( line )
            found_count = found_count + 1
      if ( found_count == - count ):
        break

  return hi_list

###############################################################################

def store_parameter_old( uvim, par_name, par_value ):
  if len( par_name.split() ) > 1:
    raise error( 'parameter name %s contains white spaces' % ( par_name ) )
  if isinstance( par_value, bool ):
    par_type = 'BOOL'
  elif isinstance( par_value, int ):
    par_type = 'INT'
  elif isinstance( par_value, float ):
    par_type = 'FLOAT'
  else:
    par_type = 'STR'
  par_string = 'PARAM %s %s = %s' % ( par_type, par_name, repr( par_value ) )
  if par_type != 'STR':
    if ( len( par_string ) <= 70 ):
      write_history( uvim, [ par_string ] )
    else:
      raise error( 'parameter name+value of %s is too long' % par_name )
  else:
    if ( len( par_string ) <= 70 ):
      write_history( uvim, [ par_string ] )
    else:
      par_strings = [ par_string[ 0 : 70 ] ]
      i = 0
      while ( ( i < 99 ) and ( len( par_string ) > 70 ) ):
        i = i + 1
        par_string
        par_string = 'PAR%02d %s %s = %s' % ( i, par_type, par_name, par_string[ 70 : ] )
        if ( len( par_string ) <= 70 ):
          par_strings.append( par_string )
        else:
          par_strings.append( par_string[ 0 : 70 ] )
      if ( len( par_string ) <= 70 ):
        write_history( uvim, par_strings )
      else:
        raise error( 'parameter name+value of %s is too long' % par_name )
  return

###############################################################################

def restore_parameter_old( uvim, par_name ):
  if ( len( par_name.split() ) > 1 ):
    raise error( 'parameter name %s contains white spaces' % ( par_name ) )
  lines = read_history( uvim, strings = [ 'PARAM', par_name ], count = -1, word_match = True )
  if ( len( lines ) == 0 ):
    raise error( 'parameter name %s was not found' % ( par_name ) )
  line = lines[ 0 ]
  columns = [ column.strip() for column in line.split() ]
  par_type = columns[ 1 ]
  par_value_string = columns[ 4 ]
  if ( par_type == 'BOOL' ):
    if ( par_value_string == repr( True ) ):
      par_value = True
    elif ( par_value_string == repr( False ) ):
      par_value = False
    else:
      raise error( 'invalid value %s for parameter %s' % ( par_value_string, par_name ) )
  elif ( par_type == 'INT' ):
    par_value = int( par_value_string )
  elif ( par_type == 'FLOAT' ):
    par_value = float( par_value_string )
  elif ( par_type == 'STR' ):
    if ( par_value_string[ 0 ] != "'" ):
      raise error( 'invalid value %s for parameter %s' % ( par_value_string, par_name ) )
    start_index = line.find( "'" ) + 1
    stop_index = line[ start_index : ].find( "'" )
    if ( stop_index > - 1 ):
      par_value = line[ start_index : start_index + stop_index ]
    else:
      par_value = line[ start_index : ]
      i = 0
      while ( i < 99 ):
        i = i + 1
        lines = read_history( uvim, strings = [ 'PAR%02d' % ( i ), 'STR', par_name ], count = -1, word_match = True )
        if ( len( lines ) == 0 ):
          raise error( 'parameter name %s is an unterminated string' % ( par_name ) )
        line = lines[ 0 ]
        start_index = line.find( "=" ) + 2
        stop_index = line[ start_index : ].find( "'" )
        if ( stop_index > - 1 ):
          par_value = par_value + line[ start_index : start_index + stop_index ]
          break
        else:
          par_value = par_value + line[ start_index : ]
      if ( stop_index == - 1 ):
        raise error( 'parameter name %s is an unterminated string' % ( par_name ) )
  else:
    raise error( 'unknown parameter type %s' % ( par_type ) )

  return par_value

###############################################################################

def write_ps_row( uvim, field_string, comment_string ):
  if ( len( field_string ) > 8 ):
    raise error( 'length of field_string is too large' )
  if ( len( comment_string ) > 80 ):
    raise error( 'length of comment_string is too large' )
  if ( not table_exists( uvim, 'PS', 0 ) ):
    ps_table = new_table( uvim, 'PS', 1 )
  else:
    wiz_uvim = wizardry( uvim )
    ps_table = wiz_uvim.table( 'PS', 0 )
  [ date, time ] = datetime.now().isoformat( ' ' ).split()
  date = date.replace( '-', '' )
  time = time.replace( ':', '' )
  time = time.replace( '.', '' )
  ps_row = new_table_row( ps_table )
  ps_row.field_name = field_string
  ps_row.comments = comment_string
  ps_row.obsdate = date[ 0 : 8 ]
  ps_row.status = time[ 0 : 8 ]
  ps_table.append( ps_row )
  ps_table.close()
  return

###############################################################################

def read_ps_rows( uvim, field_string ):
  lines = []
  if ( not table_exists( uvim, 'PS', 0 ) ):
    raise error( 'PS table does not exist' )
  wiz_uvim = wizardry( uvim )
  ps_table = wiz_uvim.table( 'PS', 0 )
  for ps_row in ps_table:
    if ( ps_row.field_name.strip() == field_string ):
      lines.append( ps_row.comments )
  ps_table.close()
  return ( lines )

###############################################################################

def store_parameter( uvim, par_name, par_value ):
  if ( len( par_name.split() ) > 1 ):
    raise error( 'parameter name %s contains white spaces' % ( par_name ) )
  if ( type( par_value ) in [ bool, bool_, bool8 ] ):
    par_type = 'BOOL'
  elif ( type( par_value ) in [ int, intc, intp, int_, int0, int8, int16, int32, int64,
      uintc, uintp, uint, uint0, uint8, uint16, uint32, uint64 ] ):
    par_type = 'INT'
  elif ( type( par_value ) in [ float, float_, float32, float64 ] ): # , float96 ] ):
    par_type = 'FLOAT'
  else:
    par_type = 'STR'
  par_string = '%s %s = %s' % ( par_type, par_name, repr( par_value ) )
  if ( par_type != 'STR' ):
    if ( len( par_string ) <= 80 ):
      write_ps_row( uvim, 'PARAM', par_string )
    else:
      raise error( 'parameter name+value of %s is too long' % par_name )
  else:
    if ( len( par_string ) <= 80 ):
      write_ps_row( uvim, 'PARAM', par_string )
    else:
      write_ps_row( uvim, 'PARAM', par_string[ 0 : 80 ] )
      par_string = par_string[ 80 : ]
      while ( len( par_string ) > 80 ):
        write_ps_row( uvim, 'PARAM', par_string[ 0 : 80 ] )
        par_string = par_string[ 80 : ]
      if ( len( par_string ) > 0 ):
        write_ps_row( uvim, 'PARAM', par_string )
  return

###############################################################################

#def restore_parameter( uvim, par_name ):
#  try:
#    par_value = restore_parameter_new( uvim, par_name )
#  except:
#    par_value = restore_parameter_old( uvim, par_name )
#  return par_value

###############################################################################

#def restore_parameter_new( uvim, par_name ):
def restore_parameter( uvim, par_name ):
  if len( par_name.split() ) > 1:
    raise error( 'parameter name %s contains white spaces' % ( par_name ) )
  lines = read_ps_rows( uvim, 'PARAM' )
  if ( len( lines ) == 0 ):
    raise error( 'parameter name %s was not found' % ( par_name ) )
  index_list = []
  for i in range( len( lines ) ):
    line = lines[ i ]
    columns = [ column.strip() for column in line.split() ]
    if ( len( columns ) > 1 ):
      if ( columns[ 1 ] == par_name ):
        index_list.append( i )
  if ( len( index_list ) == 0 ):
    raise error( 'parameter name %s was not found' % ( par_name ) )
  line = lines[ index_list[ - 1 ] ]
  columns = [ column.strip() for column in line.split() ]
  par_type = columns[ 0 ]
  par_value_string = columns[ 3 ]
  if ( par_type == 'BOOL' ):
    if ( par_value_string == repr( True ) ):
      par_value = True
    elif ( par_value_string == repr( False ) ):
      par_value = False
    else:
      raise error( 'invalid value %s for parameter %s' % ( par_value_string, par_name ) )
  elif ( par_type == 'INT' ):
    par_value = int( par_value_string )
  elif ( par_type == 'FLOAT' ):
    par_value = float( par_value_string )
  elif ( par_type == 'STR' ):
    if ( par_value_string[ 0 ] != "'" ):
      raise error( 'invalid value %s for parameter %s' % ( par_value_string, par_name ) )
    start_index = line.find( "'" ) + 1
    stop_index = line[ start_index : ].find( "'" )
    if ( stop_index > -1 ):
      par_value = line[ start_index : start_index + stop_index ]
    else:
      par_value = line[ start_index : ]
      index = index_list[ - 1 ] + 1
      while ( index < len( lines ) ):
        line = lines[ index ]
        stop_index = line.find( "'" )
        if ( stop_index > - 1 ):
          par_value = par_value + line[ 0 : stop_index ]
          break
        else:
          par_value = par_value + line
        index = index + 1
      if ( stop_index == - 1 ):
        raise error( 'parameter name %s is an unterminated string' % ( par_name ) )
  else:
    raise error( 'unknown parameter type %s' % ( par_type ) )

  return par_value

###############################################################################


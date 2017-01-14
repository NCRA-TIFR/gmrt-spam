###############################################################################

# import Python modules
from os import *
from errno import *
from sys import *
from math import *
from datetime import *

# import 3rd party modules
from numpy import *

# import ParselTongue modules
import Obit, History, OErr
import AIPS
from AIPSTask import *
from AIPSData import *
import Wizardry.AIPSData

# import user modules
from files import *
from sphere import *
from error import *

# HISTORY:
# 20090925 HTI: Replaced numarray by numpy.

###############################################################################

this_aips_id = -1

###############################################################################

class _aips_plane_iter( object ):

  def __init__( self, data, err ):
    self._err = err
    self._data = data
    self._count = -1
    self._desc = self._data.Desc.Dict
    self._size = 1
    for i in range( 0, 2 ):
      self._size = self._size * self._desc[ 'inaxes' ][ i ]
    self._len = 1
    for i in range( 2, self._desc['naxis'] ):
      self._len = self._len * self._desc[ 'inaxes' ][ i ]
    self._blc = [ 1 for i in range( self._desc[ 'naxis' ] ) ]
    self._trc = [ 1 for i in range( self._desc[ 'naxis' ] ) ]
    self._trc[ 0 : 2 ] = self._desc[ 'inaxes' ][ 0 : 2 ]
    self._dirty = False
    return

  def __len__( self ):
    return self._size

  def next( self ):
    self._count = self._count + 1
    if ( self._count >= self._len ):
      if self._dirty:
        self.update()
      raise StopIteration
    self._fill()
    return self

  def _fill( self ):
    if self._dirty:
      self.update()
    dummy = self._count
    for i in range( 2, self._desc[ 'naxis' ] ):
      self._blc[ i ] = ( dummy % self._desc[ 'inaxes' ][ i ] ) + 1
      self._trc[ i ] = self._blc[ i ]
      dummy = int( floor( dummy / self._desc[ 'inaxes' ][ i ] ) )
    self._data.ReadPlane( self._err, blc = self._blc, trc = self._trc )
    if self._err.isErr:
      raise RuntimeError, "Reading image pixels"
    shape = []
    for i in range( 0, 2 ):
      shape.insert( 0, self._desc[ 'inaxes' ][ i ] )
    shape = tuple( shape )
    self._buffer = fromstring( self._data.PixBuf, 
        dtype = float32 ).reshape( shape )
  def update( self ):
    if self._dirty:
      self._data.PixBuf[ : ] = self._buffer.tostring()
      Obit.ImageWrite( self._data.me, self._err.me )
      if self._err.isErr:
        raise RuntimeError, "Writing image pixels"
      self._dirty = False
    return

  def _get_blc( self ):
    return self._blc
  blc = property( _get_blc )

  def _get_trc( self ):
      return self._trc
  trc = property( _get_trc )

  def _get_pixels( self ):
    pixels = array( self._buffer, dtype = float32 )
    return pixels
  def _set_pixels( self, pixels ):
    self._buffer = array( pixels, dtype = float32 )
    self._dirty = True
    return
  pixels = property( _get_pixels, _set_pixels )

###############################################################################

class aips_image( Wizardry.AIPSData.AIPSImage ):
  
  _initialised = False
  _iter = None
  _len = None
  
  def _initialise( self ):
    if not self._initialised:
      self._len = 1
      for i in range( 2, self._data.Desc.Dict[ 'naxis' ] ):
        self._len = self._len * self._data.Desc.Dict[ 'inaxes' ][ i ]
      self._initialised = True
    return
  
  def __len__( self ):
    if not self._initialised:
      self._initialise()
    return self._len
  
  def __getitem__( self, index ):
    if not self._initialised:
      self._initialise()
    if index >= self._len:
      raise IndexError
    if self._iter is None:
      self.__iter__()
    self._iter._count = index
    self._iter._fill()
    return self._iter
  
  def __iter__( self ):
    if not self._initialised:
      self._initialise()
    self._data.Open( 3, self._err )
    if self._err.isErr:
      raise RuntimeError
    self._iter = _aips_plane_iter( self._data, self._err )
    return self._iter
  
#  def update( self ):
#    _AIPSData.update( self )
  
  def add_table( self, name, version, **kwds ):
    if ( not name.startswith('AIPS ') ):
      name = 'AIPS ' + name
    try:
      return self.attach_table( name, version, **kwds )
    except NotImplementedError:
      if ( name == 'AIPS HI' ):
        data = Obit.ImageCastData( self._data.me )
        version = 1
        Obit.TableHistory( data, [ version ], 3, name, self._err.me )
        if self._err.isErr:
          raise RuntimeError
        return Wizardry.AIPSData._AIPSHistory( self._data )
      else:
        msg = 'Attaching %s tables is not implemented yet' % name
        raise NotImplementedError, msg
    else:
      msg = 'Attaching %s tables is not implemented yet' % name
      raise NotImplementedError, msg
#    if ( version == 0 ):
#      version = Obit.ImageGetHighVer( self._data.me, name ) + 1
#    no_parms = 0
#    if ( 'no_parms' in kwds ):
#      no_parms = kwds[ 'no_parms' ]
#    if ( 'no_if' in kwds ):
#      no_if = kwds[ 'no_if' ]
#    if 'no_pol' in kwds:
#      no_pol = kwds[ 'no_pol' ]
#    data = Obit.ImageCastData( self._data.me )
#    if ( name == 'AIPS Hi' ):
#      Obit.TableHistory( data, [ version ], 3, name, self._err.me )
#    elif ( name == 'AIPS PS' ):
#      Obit.TablePS( data, [ version ], 3, name, self._err.me )
#    elif ( name == 'AIPS SN' ):
#      Obit.TableSN( data, [ version ], 3, name, no_pol, no_if, self._err.me )
#    else:
#      msg = 'Attaching %s tables is not implemented yet' % name
#      raise NotImplementedError, msg
#    if self._err.isErr:
#      raise RuntimeError
#    return Wizardry.AIPSData._AIPSTable( self._data, name, version )

###############################################################################

class aips_uv_data( Wizardry.AIPSData.AIPSUVData ):
  
  def add_table( self, name, version, **kwds ):
    if not name.startswith('AIPS '):
      name = 'AIPS ' + name
    try:
      return self.attach_table( name, version, **kwds )
    except NotImplementedError:
      if ( name == 'AIPS OB' ):
        data = Obit.UVCastData( self._data.me )
        if ( version == 0 ):
          version = Obit.UVGetHighVer( self._data.me, name ) + 1
        Obit.TableOB( data, [ version ], 3, name, self._err.me )
        if self._err.isErr:
          raise RuntimeError
        return Wizardry.AIPSData._AIPSTable( self._data, name, version )
      elif ( name == 'AIPS FG' ):
        data = Obit.UVCastData( self._data.me )
        if ( version == 0 ):
          version = Obit.UVGetHighVer( self._data.me, name ) + 1
        Obit.TableFG( data, [ version ], 3, name, self._err.me )
        if self._err.isErr:
          raise RuntimeError
        return Wizardry.AIPSData._AIPSTable( self._data, name, version )
      elif ( name == 'AIPS HI' ):
        data = Obit.UVCastData( self._data.me )
        version = 1
        Obit.TableHistory( data, [ version ], 3, name, self._err.me )
        if self._err.isErr:
          raise RuntimeError
        return Wizardry.AIPSData._AIPSHistory( self._data )
      else:
        msg = 'Attaching %s tables is not implemented yet' % name
        raise NotImplementedError, msg
    else:
      msg = 'Attaching %s tables is not implemented yet' % name
      raise NotImplementedError, msg
#    if ( name != 'AIPS OB' ):
#    if True:
#      return self.attach_table( name, version, **kwds )
#    if ( version == 0 ):
#      version = Obit.UVGetHighVer( self._data.me, name ) + 1
#    no_parms = 0
#    if ( 'no_parms' in kwds ):
#      no_parms = kwds[ 'no_parms' ]
#    if ( 'no_if' in kwds ):
#      no_if = kwds[ 'no_if' ]
#    if 'no_pol' in kwds:
#      no_pol = kwds[ 'no_pol' ]
#    data = Obit.UVCastData( self._data.me )
#    if ( name == 'AIPS OB' ):
#      Obit.TableOB( data, [ version ], 3, name, self._err.me )
#    elif ( name == 'AIPS NI' ):
#      Obit.TableNI( data, [ version ], 3, name, kwds[ 'num_coef' ], self._err.me )
#    elif ( name == 'AIPS PS' ):
#      Obit.TablePS( data, [ version ], 3, name, self._err.me )
#    else:
#      msg = 'Attaching %s tables is not implemented yet' % name
#      raise NotImplementedError, msg
#    if self._err.isErr:
#      raise RuntimeError
#    return Wizardry.AIPSData._AIPSTable( self._data, name, version )

###############################################################################

def get_aips_userid():
  return AIPS.userno

###############################################################################

def get_aips_magic_value():
  mv_array = fromstring( 'INDEFINE', dtype = float32 )
  return mv_array[ 0 ]

###############################################################################

def get_aips_file( aips_disk, aips_name, aips_class, aips_seq, aips_type,
    use_sequence_convention = True ):

  data = None
  if ( ( aips_type == 'UV' ) or ( aips_type == 'MA' ) ):
    # fully specified AIPS file
    if ( ( aips_disk > 0 ) and ( aips_seq > 0 ) ):
      if ( aips_type == 'UV' ):
        data = AIPSUVData( aips_name, aips_class, aips_disk, aips_seq )
      else:
        data = AIPSImage( aips_name, aips_class, aips_disk, aips_seq )
    # find highest sequence on specified disk
    elif ( ( aips_disk > 0 ) and ( aips_seq == 0 ) ):
      seq = 0
      cat_disk = AIPSCat( aips_disk )[ aips_disk ]
      for cat in cat_disk:
        if ( [ cat.type, cat.klass, cat.name ] == [ aips_type, aips_class, aips_name ] ):
          if ( cat.seq > seq ):
            seq = cat.seq
            if ( aips_type == 'UV' ):
              data = AIPSUVData( aips_name, aips_class, aips_disk, seq )
            else:
              data = AIPSImage( aips_name, aips_class, aips_disk, seq )
    # find highest available sequence on specified disk
    elif ( ( aips_disk > 0 ) and ( aips_seq < 0 ) ):
      if use_sequence_convention:
        # use convention seq = ddss with dd = disk and ss = seq
        seq = aips_disk * 100
        cat_disk = AIPSCat( aips_disk )[ aips_disk ]
        for cat in cat_disk:
          if ( [ cat.type, cat.klass, cat.name ] == [ aips_type, aips_class, aips_name ] ):
            if ( ( cat.seq >= aips_disk * 100 ) and ( cat.seq < ( aips_disk + 1 ) * 100 ) and ( cat.seq >= seq ) ):
              seq = cat.seq + 1
        if ( seq < ( aips_disk + 1 ) * 100 ):
          if ( aips_type == 'UV' ):
            data = AIPSUVData( aips_name, aips_class, aips_disk, seq )
          else:
            data = AIPSImage( aips_name, aips_class, aips_disk, seq )
      else:
        # search all disks for lowest available sequence number
        seq = 1
        cat_table = AIPSCat()
        for disk in range( 1, len( AIPS.disks ) ):
          cat_disk = cat_table[ disk ]
          for cat in cat_disk:
            if ( [ cat.seq, cat.type, cat.klass, cat.name,  ] == [ aips_seq, aips_type, aips_class, aips_name ] ):
              if ( cat.seq >= seq ):
                seq = cat.seq + 1
        if ( aips_type == 'UV' ):
          data = AIPSUVData( aips_name, aips_class, aips_disk, seq )
        else:
          data = AIPSImage( aips_name, aips_class, aips_disk, seq )
    # find AIPS file on unspecified disk
    elif ( ( aips_disk <= 0 ) and ( aips_seq > 0 ) ):
      if ( use_sequence_convention and ( aips_seq > 100 ) ):
        # use convention seq = ddss with dd = disk and ss = seq
        disk = aips_seq / 100
        if ( disk >= len( AIPS.disks ) ):
          disk = 1
        if ( aips_type == 'UV' ):
          data = AIPSUVData( aips_name, aips_class, disk, aips_seq )
        else:
          data = AIPSImage( aips_name, aips_class, disk, aips_seq )
      else:
        cat_table = AIPSCat()
        for disk in range( len( AIPS.disks ) - 1, 0, - 1 ): # search in reverse order
          cat_disk = cat_table[ disk ]
          for cat in cat_disk:
            if ( [ cat.seq, cat.type, cat.klass, cat.name,  ] == [ aips_seq, aips_type, aips_class, aips_name ] ):
              if ( aips_type == 'UV' ):
                data = AIPSUVData( aips_name, aips_class, disk, aips_seq )
              else:
                data = AIPSImage( aips_name, aips_class, disk, aips_seq )
              break
          if ( not data is None ):
            break
    # find highest sequence on unspecified disk
    elif ( ( aips_disk <= 0 ) and ( aips_seq == 0 ) ):
      seq = 0
      cat_table = AIPSCat()
      for disk in range( 1, len( AIPS.disks ) ):
        cat_disk = cat_table[ disk ]
        for cat in cat_disk:
          if ( [ cat.type, cat.klass, cat.name ] == [ aips_type, aips_class, aips_name ] ):
            if ( cat.seq > seq ):
              seq = cat.seq
              if ( aips_type == 'UV' ):
                data = AIPSUVData( aips_name, aips_class, disk, seq )
              else:
                data = AIPSImage( aips_name, aips_class, disk, seq )
    # find highest available sequence on unspecified disk -> disk = 1
    elif ( ( aips_disk <= 0 ) and ( aips_seq < 0 ) ):
      if use_sequence_convention:
        # use convention seq = ddss with dd = disk and ss = seq
        # allow continuation above seq = 0199, map to 36ss etc
        disk = 1
        seq = disk * 100
        cat_disk = AIPSCat( disk )[ disk ]
        for cat in cat_disk:
          if ( [ cat.type, cat.klass, cat.name ] == [ aips_type, aips_class, aips_name ] ):
            if ( cat.seq >= seq ):
              seq = cat.seq + 1
          if ( ( seq >= ( disk + 1 ) * 100 ) and ( seq < len( AIPS.disks ) * 100 ) ):
            seq = len( AIPS.disks ) * 100
        if ( aips_type == 'UV' ):
          data = AIPSUVData( aips_name, aips_class, disk, seq )
        else:
          data = AIPSImage( aips_name, aips_class, disk, seq )
      else:
        # search all disks for lowest available sequence number
        seq = 1
        cat_table = AIPSCat()
        for disk in range( 1, len( AIPS.disks ) ):
          cat_disk = cat_table[ disk ]
          for cat in cat_disk:
            if ( [ cat.seq, cat.type, cat.klass, cat.name,  ] == [ aips_seq, aips_type, aips_class, aips_name ] ):
              if ( cat.seq >= seq ):
                seq = cat.seq + 1
        if ( aips_type == 'UV' ):
          data = AIPSUVData( aips_name, aips_class, 1, seq )
        else:
          data = AIPSImage( aips_name, aips_class, 1, seq )
  else:
    raise error( 'unknown or unsupported file type:%s' % ( aips_type ) )
  if ( data is None ):
    raise error( 'AIPS file %d.%s.%s.%d.%s could not be found or allocated' % 
        ( aips_disk, aips_name, aips_class, aips_seq, aips_type ) )
  return data

###############################################################################

def allocate_aips_id():
  
  global this_aips_id
  this_pid = os.getpid()
  if ( ( this_aips_id < 1 ) or ( this_aips_id > 35 ) or
      ( not file_exists( '/tmp/AIPS%d.%d' % ( this_aips_id, this_pid ) ) ) ):
    # find/allocate AIPS ID
    done = False
    while ( not done ):
      aips_list = []
      tmp_file_names = get_directory( '/tmp' )
      for file_name in tmp_file_names:
        if ( file_name[ 0 : 4 ] == 'AIPS' ):
          aips_id = ehex_to_decimal( file_name[ 4 ] )
          pid = int( file_name[ 6 : ] )
          # ping processes to see active state
          if ( pid == this_pid ):
            this_aips_id = aips_id
            done = True
            break
          try:
            os.kill( pid, 0 )
          except OSError, err:
            if ( err.errno != EPERM ):
              # process must have died, so try to free AIPS ID
              if ( not remove_file( '/tmp/' + file_name ) ):
                aips_list.append( [ aips_id, pid ] )
          else:
            aips_list.append( [ aips_id, pid ] )
      # check for multiple occurences of same AIPS ID
      if done:
        for [ aips_id, pid ] in aips_list:
          if ( aips_id == this_aips_id ):
            if ( this_pid > pid ):
              remove_file( '/tmp/AIPS%s.%d' % ( decimal_to_ehex( this_aips_id ), this_pid ) )
              this_aips_id = -1
              done = False
              break
      if ( not done ):
        # find first free AIPS ID
        used_aips_ids = list( set( [ x[ 0 ] for x in aips_list ] ) )
        for aips_id in range( 1, 1 + 35 ):
          if ( not aips_id in used_aips_ids ):
            this_aips_id = aips_id
            make_file( '/tmp/AIPS%s.%d' % ( decimal_to_ehex( this_aips_id ), this_pid ) )
            break
        if ( this_aips_id < 0 ):
          raise error( 'no more free AIPS IDs' )
  
  return this_aips_id

###############################################################################

def get_aips_id():
  global this_aips_id
  return this_aips_id

###############################################################################

def free_aips_id():
  global this_aips_id
  this_pid = os.getpid()
  remove_file( '/tmp/AIPS%d.%d' % ( this_aips_id, this_pid ) )
  this_aips_id = -1
  return

###############################################################################

def call_aips_task( task_name, print_info = True, log_to_file = False, 
    retries = 3, keep_messages = False, keep_history = False, **keywords ):
  
#  global this_aips_id
#  this_pid = os.getpid()
#  lock_file_name = '/tmp/AIPS%s.%d' % ( decimal_to_ehex( this_aips_id ), this_pid )
#  if ( not file_exists( lock_file_name ) ):
#    allocate_aips_id()
  task = AIPSTask( task_name )
  if log_to_file:
    task.msgkill = -100
    if file_exists( 'task.log' ):
      remove_file( 'task.log' )
    task.log = file( 'task.log', 'a' )
    aips_log = file( 'aips.log', 'a' )
    aips_log.write( '%s\n' % task_name )
  output_list = None
  for key in keywords.keys():
    if ( key == 'outputs' ):
      output_list = keywords[ key ]
    else:
      keyword = keywords[ key ]
      if isinstance( keyword, ndarray ):
        keyword = keyword.tolist()
      if isinstance( keyword, list ):
        element_list = [ None ]
        for element in keyword:
          if isinstance( element, ndarray ):
            element = element.tolist()
          if isinstance( element, list ):
            sub_element_list = [ None ]
            for sub_element in element:
              if ( type( sub_element ) in [ float_, float32, float64 ] ): # , float96 ] ):
                sub_element = float( sub_element )
              elif ( type( sub_element ) in [ bool, bool_, bool8, intc, intp, int_, int0, int8,
                  int16, int32, int64, uintc, uintp, uint, uint0, uint8, uint16, uint32, uint64 ] ):
                sub_element = int( sub_element )
              sub_element_list.append( sub_element )
            element_list.append( sub_element_list )
          else:
            if ( type( element ) in [ float_, float32, float64 ] ): # , float96 ] ):
              element = float( element )
            elif ( type( element ) in [ bool, bool_, bool8, intc, intp, int_, int0, int8,
                int16, int32, int64, uintc, uintp, uint, uint0, uint8, uint16, uint32, uint64 ] ):
              element = int( element )
            element_list.append( element )
        setattr( task, key, element_list )
        if log_to_file:
          aips_log.write( '  %s = %s\n' % ( key, repr( element_list ) ) )
      else:
        if ( type( keyword ) in [ float_, float32, float64 ] ): # , float96 ] ):
          keyword = float( keyword )
        elif ( type( keyword ) in [ bool, bool_, bool8, intc, intp, int_, int0, int8,
            int16, int32, int64, uintc, uintp, uint, uint0, uint8, uint16, uint32, uint64 ] ):
          keyword = int( keyword )
        setattr( task, key, keyword )
        if log_to_file:
          aips_log.write( '  %s = %s\n' % ( key, repr( keyword ) ) )
  if log_to_file:
    aips_log.flush()
#  remove_file( lock_file_name )
  done = False
  for i in range( retries ):
    if print_info:
      if ( i == 0 ):
        print '*** calling AIPSTask %s' % ( task_name )
      else:
        print '*** calling AIPSTask %s (RETRY %d)' % ( task_name, i )
    try:
      task.go()
    except RuntimeError:
      pass
    else:
      done = True
      break
  if ( not done ):
    if print_info:
      if ( retries == 0 ):
        print '*** calling AIPSTask %s' % ( task_name )
      else:
        print '*** calling AIPSTask %s (RETRY %d, FINAL)' % ( task_name, retries )
    task.go()
#  make_file( lock_file_name )
  if log_to_file:
    task.log.close()
    task_log = file( 'task.log', 'r' )
    for line in task_log:
      aips_log.write( line )
    task.log.close()
  outputs = None
  if ( not output_list is None ):
    outputs = []
    for output in output_list:
      element = getattr( task, output )
      if isinstance( element, list ):
        element_list = []
        for sub_element in element[ 1 : ]:
          if isinstance( sub_element, list ):
            element_list.append( sub_element[ 1 : ] )
          else:
            element_list.append( sub_element )
        outputs.append( element_list )
        if log_to_file:
          aips_log.write( '  %s = %s\n' % ( output, repr( element_list ) ) )
      else:
        outputs.append( element )
        if log_to_file:
          aips_log.write( '  %s = %s\n' % ( output, repr( element ) ) )
  if log_to_file:
    aips_log.close()
  # clear AIPS message file
  if ( not keep_messages ):
    copy_file( get_aips_message_file_name( userid = 1 ), get_aips_message_file_name() )
  # clear AIPS HI table
  if ( not keep_history ):
    indata_found = False
    if ( not indata_found ):
      try:
        uvim = keywords[ 'indata' ]
      except:
        pass
      else:
        if uvim.exists():
          indata_found = True
    if ( not indata_found ):
      try:
        uvim_disk = keywords[ 'indisk' ]
        uvim_name = keywords[ 'inname' ]
        uvim_class = keywords[ 'inclass' ]
        uvim_seq = keywords[ 'inseq' ]
        uvim_type = keywords[ 'intype' ]
        uvim = get_aips_file( uvim_disk, uvim_name, uvim_class, uvim_seq, uvim.type )
      except:
        pass
      else:
        if uvim.exists():
          indata_found = True
    if ( not indata_found ):
      try:
        uvim = get_aips_file( uvim_disk, uvim_name, uvim_class, uvim_seq, 'UV' )
      except:
        pass
      else:
        if uvim.exists():
          indata_found = True
    if ( not indata_found ):
      try:
        uvim = get_aips_file( uvim_disk, uvim_name, uvim_class, uvim_seq, 'MA' )
      except:
        pass
      else:
        if uvim.exists():
          indata_found = True
    if indata_found:
      clear_history( uvim )
  return outputs

###############################################################################

def wizardry( uvim ):
  if is_uv( uvim ):
    wiz = aips_uv_data( uvim )
  elif is_image( uvim ):
    wiz = aips_image( uvim )
  else:
    raise error( 'specified AIPS file not UV nor image' )
  return wiz

###############################################################################

def table_exists( uvim, name, version = 0 ):
  table_found = False
  for table in uvim.tables:
    if ( table[ 1 ] == 'AIPS ' + name ):
      if ( version == 0 ) or ( table[ 0 ] == version ):
       table_found = True
  return table_found

###############################################################################

def get_radec( uvim ):
  for ctype in uvim.header.ctype:
    if ( ctype.find( 'RA' ) != - 1 ):
      ra_index = uvim.header.ctype.index( ctype )
    if ( ctype.find( 'DEC' ) != - 1 ):
      dec_index = uvim.header.ctype.index( ctype )
  return [ uvim.header.crval[ ra_index ], uvim.header.crval[ dec_index ] ]

###############################################################################

def get_glonlat( uvim ):
  for ctype in uvim.header.ctype:
    if ( ctype.find( 'GLON' ) != - 1 ):
      glon_index = uvim.header.ctype.index( ctype )
    if ( ctype.find( 'GLAT' ) != - 1 ):
      glat_index = uvim.header.ctype.index( ctype )
  return [ uvim.header.crval[ glon_index ], uvim.header.crval[ glat_index ] ]

###############################################################################

def set_radec( uvim, radec, shift_model = False, model_version = 0 ):
  for ctype in uvim.header.ctype:
    if ( ctype.find( 'RA' ) != - 1 ):
      ra_index = uvim.header.ctype.index( ctype )
    if ( ctype.find( 'DEC' ) != - 1 ):
      dec_index = uvim.header.ctype.index( ctype )
  wiz_uvim = wizardry( uvim )
  if ( shift_model and table_exists( uvim, 'CC', model_version ) ):
    model_list = []
    pixel_size = get_pixel_size( uvim, make_absolute = False )
    pixel_ref = get_pixel_reference( uvim )
    cc_table = wiz_uvim.table( 'CC', model_version )
    for cc in cc_table:
      cc_x = pixel_ref[ 0 ] + 3600. * cc.deltax / pixel_size[ 0 ]
      cc_y = pixel_ref[ 1 ] + 3600. * cc.deltay / pixel_size[ 1 ]
      cc_radec = calculate_source_radec( uvim, [ cc_x, cc_y ] )
      model_list.append( cc_radec )
    cc_table.close()
  wiz_uvim.header.crval[ ra_index ] = radec[ 0 ]
  wiz_uvim.header.crval[ dec_index ] = radec[ 1 ]
  wiz_uvim.header.obsra = radec[ 0 ]
  wiz_uvim.header.obsdec = radec[ 1 ]
  wiz_uvim.header.update()
  wiz_uvim = wizardry( uvim )
  if ( shift_model and table_exists( uvim, 'CC', model_version ) ):
    cc_table = wiz_uvim.table( 'CC', model_version )
    i = 0
    for cc_row in cc_table:
      cc_radec = model_list[ i ]
      [ cc_x, cc_y ] = calculate_source_position( uvim, cc_radec )
      cc_row.deltax = pixel_size[ 0 ] * ( cc_x - pixel_ref[ 0 ] ) / 3600.
      cc_row.deltay = pixel_size[ 1 ] * ( cc_y - pixel_ref[ 1 ] ) / 3600.
      cc_row.update()
      i = i + 1
    cc_table.close()
  del wiz_uvim
  return

###############################################################################

def set_glonlat( uvim, glonlat ):
  for ctype in uvim.header.ctype:
    if ( ctype.find( 'GLON' ) != - 1 ):
      glon_index = uvim.header.ctype.index( ctype )
    if ( ctype.find( 'GLAT' ) != - 1 ):
      glat_index = uvim.header.ctype.index( ctype )
  wiz_uvim = wizardry( uvim )
  wiz_uvim.header.crval[ glon_index ] = glonlat[ 0 ]
  wiz_uvim.header.crval[ glat_index ] = glonlat[ 1 ]
  wiz_uvim.header.update()
  return

###############################################################################

def get_frequency( uvim ):
  freq_index = uvim.header.ctype.index( 'FREQ' )
  return uvim.header.crval[ freq_index ]

###############################################################################

def set_frequency( uvim, frequency ):
  freq_index = uvim.header.ctype.index( 'FREQ' )
  wiz_uvim = wizardry( uvim )
  wiz_uvim.header.crval[ freq_index ] = frequency
  wiz_uvim.header.update()
  del wiz_uvim
  return

###############################################################################

def get_channel_count( uvim ):
  freq_index = uvim.header.ctype.index( 'FREQ' )
  return uvim.header.naxis[ freq_index ]

###############################################################################

def get_channel_width( uvim ):
  freq_index = uvim.header.ctype.index( 'FREQ' )
  return abs( uvim.header.cdelt[ freq_index ] )

###############################################################################

def get_bandwidth( uvim ):
  return get_channel_count( uvim ) * get_channel_width( uvim )

###############################################################################

def get_frequency_list( uvim ):
  freq_list = []
  freq_index = uvim.header.ctype.index( 'FREQ' )
#  for freq_i in range( uvim.header.naxis[ freq_index ] ):
  for freq_i in range( 1, 1 + uvim.header.naxis[ freq_index ] ):
#    freq_list.append( uvim.header.crval[ freq_index ] + ( float( freq_i ) + 1.5
#        - uvim.header.crpix[ freq_index ] ) * uvim.header.cdelt[ freq_index ] )
    freq_list.append( uvim.header.crval[ freq_index ] + ( float( freq_i )
        - uvim.header.crpix[ freq_index ] ) * uvim.header.cdelt[ freq_index ] )
  return freq_list

###############################################################################

def get_central_frequency( uvim ):
  return float( array( get_frequency_list( uvim ), dtype = float64 ).mean() )

###############################################################################

def get_central_wavelength( uvim ):
  return ( 299792458. / get_central_frequency( uvim ) )

###############################################################################

def get_image_size( im ):
  for ctype in im.header.ctype:
    if ( ctype.find( 'RA--' ) != - 1 ):
      ra_index = im.header.ctype.index( ctype )
    elif ( ctype.find( 'XA--' ) != - 1 ):
      ra_index = im.header.ctype.index( ctype )
    elif ( ctype.find( 'GLON' ) != - 1 ):
      ra_index = im.header.ctype.index( ctype )
    if ( ctype.find( 'DEC-' ) != - 1 ):
      dec_index = im.header.ctype.index( ctype )
    elif ( ctype.find( 'XEC-' ) != - 1 ):
      dec_index = im.header.ctype.index( ctype )
    elif ( ctype.find( 'GLAT' ) != - 1 ):
      dec_index = im.header.ctype.index( ctype )
  return [ im.header.naxis[ ra_index ], im.header.naxis[ dec_index ] ]

###############################################################################

def get_pixel_size( im, make_absolute = True ):
  for ctype in im.header.ctype:
    if ( ctype.find( 'RA--' ) != - 1 ):
      ra_index = im.header.ctype.index( ctype )
    elif ( ctype.find( 'XA--' ) != - 1 ):
      ra_index = im.header.ctype.index( ctype )
    elif ( ctype.find( 'GLON' ) != - 1 ):
      ra_index = im.header.ctype.index( ctype )
    if ( ctype.find( 'DEC-' ) != - 1 ):
      dec_index = im.header.ctype.index( ctype )
    elif ( ctype.find( 'XEC-' ) != - 1 ):
      dec_index = im.header.ctype.index( ctype )
    elif ( ctype.find( 'GLAT' ) != - 1 ):
      dec_index = im.header.ctype.index( ctype )
  if make_absolute:
    pix_size = [ abs( 3600. * im.header.cdelt[ ra_index ] ), abs( 3600. * im.header.cdelt[ dec_index ] ) ]
  else:
    pix_size = [ 3600. * im.header.cdelt[ ra_index ], 3600. * im.header.cdelt[ dec_index ] ]
  return pix_size

###############################################################################

def get_pixel_reference( im ):
  for ctype in im.header.ctype:
    if ( ctype.find( 'RA--' ) != - 1 ):
      ra_index = im.header.ctype.index( ctype )
    elif ( ctype.find( 'XA--' ) != - 1 ):
      ra_index = im.header.ctype.index( ctype )
    elif ( ctype.find( 'GLON' ) != - 1 ):
      ra_index = im.header.ctype.index( ctype )
    if ( ctype.find( 'DEC-' ) != - 1 ):
      dec_index = im.header.ctype.index( ctype )
    elif ( ctype.find( 'XEC-' ) != - 1 ):
      dec_index = im.header.ctype.index( ctype )
    elif ( ctype.find( 'GLAT' ) != - 1 ):
      dec_index = im.header.ctype.index( ctype )
  pix_reference = [ im.header.crpix[ ra_index ], im.header.crpix[ dec_index ] ]
  return pix_reference

###############################################################################

def get_image_rotation( im, epsilon = 1.e-3 ):
  for ctype in im.header.ctype:
    if ( ctype.find( 'RA--' ) != - 1 ):
      ra_index = im.header.ctype.index( ctype )
    elif ( ctype.find( 'XA--' ) != - 1 ):
      ra_index = im.header.ctype.index( ctype )
    elif ( ctype.find( 'GLON' ) != - 1 ):
      ra_index = im.header.ctype.index( ctype )
    if ( ctype.find( 'DEC-' ) != - 1 ):
      dec_index = im.header.ctype.index( ctype )
    elif ( ctype.find( 'XEC-' ) != - 1 ):
      dec_index = im.header.ctype.index( ctype )
    elif ( ctype.find( 'GLAT' ) != - 1 ):
      dec_index = im.header.ctype.index( ctype )
  if ( ( im.header.crota[ ra_index ] != 0. ) or ( im.header.crota[ dec_index ] != 0. ) ):
    if ( im.header.crota[ ra_index ] < epsilon ):
      rot = im.header.crota[ dec_index ]
    else:
      drot = im.header.crota[ ra_index ] - im.header.crota[ dec_index ] 
      drot = amodulo( drot + 180., 360 ) - 180.
      if ( drot < epsilon ):
        rot = im.header.crota[ dec_index ]
      else:
        raise error( "don't know how to handle different rotations on RA and DEC axes" )
  else:
    rot = 0.
  return rot

###############################################################################

def get_beam_size( im ):
  bmaj = 3600. * im.header.bmaj
  bmin = 3600. * im.header.bmin
  return [ bmaj, bmin, im.header.bpa ]

###############################################################################

def set_beam_size( im, beam ):
  [ bmaj, bmin, bpa ] = beam
  wizim = wizardry( im )
  wizim.header.bmaj = bmaj / 3600.
  wizim.header.bmin = bmin / 3600.
  wizim.header.bpa = bpa
  wizim.header.update()
  del wizim
  return

###############################################################################

def get_model_flux_from_position_area( facet, pos, radius, model_version = 0 ):
# radius in pixels
  model_flux = 0.
  [ x, y ] = pos
  pixel_size = get_pixel_size( facet, make_absolute = False )
  pixel_ref = get_pixel_reference( facet )
  radius2 = float( radius )**2
  if table_exists( facet, 'CC', model_version ):
    wiz_facet = wizardry( facet )
    cc_table = wiz_facet.table( 'CC', model_version )
    for cc in cc_table:
      cc_x = pixel_ref[ 0 ] + 3600. * cc.deltax / pixel_size[ 0 ]
      cc_y = pixel_ref[ 1 ] + 3600. * cc.deltay / pixel_size[ 1 ]
      r2 = ( cc_x - pos[ 0 ] )**2 + ( cc_y - pos[ 1 ] )**2
      if ( r2 <= radius2 ):
        model_flux = model_flux + cc.flux
    cc_table.close()
  del wiz_facet
  return model_flux

###############################################################################

def get_model_flux( facet, model_version = 0 ):
  cc_flux = 0.
  if table_exists( facet, 'CC', model_version ):
    wiz_facet = wizardry( facet )
    cc_table = wiz_facet.table( 'CC', model_version )
    for cc in cc_table:
      cc_flux = cc_flux + cc.flux
    cc_table.close()
    del wiz_facet
  return cc_flux

###############################################################################

def get_model_maximum( facet, model_version = 0 ):
  cc_flux = 0.
  if table_exists( facet, 'CC', model_version ):
    wiz_facet = wizardry( facet )
    cc_table = wiz_facet.table( 'CC', model_version )
    for cc in cc_table:
      if ( cc.flux > cc_flux ):
        cc_flux = cc.flux
    cc_table.close()
  del wiz_facet
  return cc_flux

###############################################################################

def get_model_minimum( facet, model_version = 0 ):
  cc_flux = 0.
  if table_exists( facet, 'CC', model_version ):
    wiz_facet = wizardry( facet )
    cc_table = wiz_facet.table( 'CC', model_version )
    for cc in cc_table:
      if ( cc.flux < cc_flux ):
        cc_flux = cc.flux
    cc_table.close()
  del wiz_facet
  return cc_flux

###############################################################################

def get_model_component_count( facet, model_version = 0 ):
  cc_count = 0
  if table_exists( facet, 'CC', model_version ):
    wiz_facet = wizardry( facet )
    cc_table = wiz_facet.table( 'CC', model_version )
    for cc in cc_table:
      cc_count = cc_count + 1
  cc_table.close()
  del wiz_facet
  return cc_count

###############################################################################

def model_table_empty( facet, model_version = 0 ):
  table_empty = True
  if table_exists( facet, 'CC', model_version ):
    wiz_facet = wizardry( facet )
    cc_table = wiz_facet.table( 'CC', model_version )
    for cc in cc_table:
      if ( cc.flux != 0. ):
        table_empty = False
        break
    cc_table.close()
    del wiz_facet
  return table_empty

###############################################################################

def change_source_name( uvim, old_name, new_name, allow_rename = False ):
  wiz_uvim = wizardry( uvim )
  if ( is_uv( uvim ) and table_exists( uvim, 'SU', 1 ) ):
    su_table = wiz_uvim.table( 'SU', 1 )
    for row in su_table:
      if ( row.source.strip() == old_name.strip() ):
        row.source = new_name.strip().ljust( len( row.source ) )
        row.update()
    su_table.close()
    su_table = wiz_uvim.table( 'SU', 1 )
  else:
    if ( uvim.header[ 'object' ] == old_name ):
      wiz_uvim.header[ 'object' ] = new_name
      wiz_uvim.header.update()
    if ( allow_rename and ( uvim.name == old_name ) ):
      if is_uv( uvim ):
        new_uvim = get_aips_file( uvim.disk, new_name, uvim.klass, -1, 'UV' )
      else:
        new_uvim = get_aips_file( uvim.disk, new_name, uvim.klass, -1, 'MA' )
      uvim.rename( name = new_uvim.name, klass = new_uvim.klass, seq = new_uvim.seq )
  del wiz_uvim
  return

###############################################################################

def get_source_names( uvim ):
  source_names = []
  if ( is_uv( uvim ) and table_exists( uvim, 'SU', 1 ) ):
    su_table = uvim.table( 'SU', 1 )
    for row in su_table:
      source_names.append( row.source.strip() )
  else:
    source_names.append( uvim.header[ 'object' ].strip() )
  return source_names

###############################################################################

def add_circular_clean_box( facet_file_name, facet_id, pos, clean_box_radius ):
# clean_box_radius in pixels
  if ( not file_exists( facet_file_name ) ):
    raise error( 'facet file %s does not exist' % facet_file_name )
  facet_file = file( facet_file_name, mode = 'a' )
  facet_file.write( '%05d  %d  %d  %d %d\n' %
      ( facet_id, - 1, int( around( clean_box_radius ) ),
      int( around( pos[ 0 ] ) ), int( around( pos[ 1 ] ) ) ) )
  facet_file.close()
  return

###############################################################################

def add_rectangular_clean_box( facet_file_name, facet_id, pos1, pos2 ):
# clean_box_radius in pixels
  if ( not file_exists( facet_file_name ) ):
    raise error( 'facet file %s does not exist' % facet_file_name )
  facet_file = file( facet_file_name, mode = 'a' )
  facet_file.write( '%05d  %d %d  %d %d\n' %
      ( facet_id, int( around( pos1[ 0 ] ) ), int( around( pos1[ 1 ] ) ),
                  int( around( pos2[ 0 ] ) ), int( around( pos2[ 1 ] ) ) ) )
  facet_file.close()
  return

###############################################################################

def remove_rectangular_clean_box( facet_file_name, facet_id, pos1, pos2 ):
# clean_box_radius in pixels
  if ( not file_exists( facet_file_name ) ):
    raise error( 'facet file %s does not exist' % facet_file_name )
  # scan box file for line(s), keep all other lines
  keep_lines = []
  facet_file = file( facet_file_name, mode = 'r' )
  for line in facet_file:
    line_found = False
    words = [ word.strip() for word in line.split() ]
    if ( len( words ) == 0 ):
      continue
    if ( words[ 0 ][ 0 ] == '#' ):
      continue
    if ( len( words ) == 5 ):
      try:
        line_id = int( words[ 0 ] )
      except:
        pass
      else:
        if ( ( line_id == facet_id ) and 
            ( int( words[ 1 ] ) == int( around( pos1[ 0 ] ) ) ) and
            ( int( words[ 2 ] ) == int( around( pos1[ 1 ] ) ) ) and
            ( int( words[ 3 ] ) == int( around( pos2[ 0 ] ) ) ) and
            ( int( words[ 4 ] ) == int( around( pos2[ 1 ] ) ) ) ):
          line_found = True
    if ( not line_found ):
      keep_lines.append( line )
  # write kept lines to box file
  if ( len( keep_lines ) > 0 ):
    remove_file( facet_file_name )
    facet_file = file( facet_file_name, mode = 'w' )
    for line in keep_lines:
      facet_file.write( line )
    facet_file.close()

  return

###############################################################################

def remove_circular_clean_box( facet_file_name, facet_id, pos, clean_box_radius ):
# clean_box_radius in pixels
  if ( not file_exists( facet_file_name ) ):
    raise error( 'facet file %s does not exist' % facet_file_name )
  # scan box file for line(s), keep all other lines
  keep_lines = []
  facet_file = file( facet_file_name, mode = 'r' )
  for line in facet_file:
    line_found = False
    words = [ word.strip() for word in line.split() ]
    if ( len( words ) == 0 ):
      continue
    if ( words[ 0 ][ 0 ] == '#' ):
      continue
    if ( len( words ) == 5 ):
      try:
        line_id = int( words[ 0 ] )
      except:
        pass
      else:
        if ( ( line_id == facet_id ) and ( int( words[ 1 ] ) == - 1 ) and
            ( int( words[ 2 ] ) == int( around( clean_box_radius ) ) ) and
            ( int( words[ 3 ] ) == int( around( pos[ 0 ] ) ) ) and
            ( int( words[ 4 ] ) == int( around( pos[ 1 ] ) ) ) ):
          line_found = True
    if ( not line_found ):
      keep_lines.append( line )
  # write kept lines to box file
  if ( len( keep_lines ) > 0 ):
    remove_file( facet_file_name )
    facet_file = file( facet_file_name, mode = 'w' )
    for line in keep_lines:
      facet_file.write( line )
    facet_file.close()

  return

###############################################################################

def add_facet( facet_file_name, radec, facet_size, facet_id = 0,
    add_clean_box = False, edge_size = 8, ):

# TODO: take proper care of roundoff towards upper boundary of hmsdms in boxfile
# Note that AIPS doesn't seem to mind

  if file_exists( facet_file_name ):
    if ( facet_id == 0 ):
      new_id = get_facet_count( facet_file_name, count_gaps = True ) + 1
    else:
      # check if field ID is not in use
      facet_file = file( facet_file_name, mode = 'r' )
      for line in facet_file:
        words = [ word.strip() for word in line.split() ]
        if ( len( words ) == 0 ):
          continue
        if ( words[ 0 ][ 0 ] == '#' ):
          continue
        if ( words[ 0 ] == 'C' ) or ( words[ 0 ] == 'F' ):
          facet_number = int( words[ 1 ] )
          if ( facet_number == facet_id ):
            raise error( 'facet ID is already in use' )
      facet_file.close()
      new_id = facet_id
    facet_file = file( facet_file_name, mode = 'a' )
  else:
    if ( facet_id == 0 ):
      new_id = 1
    else:
      new_id = facet_id
    facet_file = file( facet_file_name, mode = 'w' )
  ra_dec = degdeg_to_hmsdms( radec, precision = [ 3, 2 ] )
  if ( asign( ra_dec[ 3 ] ) > 0. ):
    facet_line = ( 'C  %d  %d %d  %02d %02d %06.3f   %02d %02d %05.2f\n' %
        ( new_id, facet_size[ 0 ] - 12, facet_size[ 1 ] - 12,
        ra_dec[ 0 ], ra_dec[ 1 ], ra_dec[ 2 ], ra_dec[ 3 ], ra_dec[ 4 ], ra_dec[ 5 ] ) )
  else:
    facet_line = ( 'C  %d  %d %d  %02d %02d %06.3f  -%02d %02d %05.2f\n' % 
        ( new_id, facet_size[ 0 ] - 12, facet_size[ 1 ] - 12,
        ra_dec[ 0 ], ra_dec[ 1 ], ra_dec[ 2 ], - ra_dec[ 3 ], ra_dec[ 4 ], ra_dec[ 5 ] ) )
  facet_file.write( facet_line )
  facet_file.close()

  if add_clean_box:
    add_circular_clean_box( facet_file_name, new_id,
        [ 1 + int( floor( ( facet_size[ 0 ] - 1. ) / 2. ) ), 1 + int( ceil( ( facet_size[ 1 ] - 1. ) / 2. ) ) ],
        int( around( ( min( facet_size ) - ( 2 * edge_size ) + 1 ) / 2 ) ) )

  return

###############################################################################

def remove_facet( facet_file_name, facet_id ):

  if ( not file_exists( facet_file_name ) ):
    raise error( 'facet file %s does not exist' % facet_file_name )
  # scan box file for line(s), keep all other lines
  keep_lines = []
  facet_file = file( facet_file_name, mode = 'r' )
  for line in facet_file:
    words = [ word.strip() for word in line.split() ]
    if ( len( words ) == 0 ):
      continue
    if ( words[ 0 ][ 0 ] == '#' ):
      continue
    try:
      facet_number = int( words[ 0 ] )
    except:
      try:
        facet_number = int( words[ 1 ] )
      except:
        keep_lines.append( line )
      else:
        if not ( ( words[ 0 ] in [ 'C', 'F' ] ) and ( facet_number == facet_id ) ):
          keep_lines.append( line )
    else:
      if ( facet_number != facet_id ):
        keep_lines.append( line )

  # write kept lines to box file
  remove_file( facet_file_name )
  facet_file = file( facet_file_name, mode = 'w' )
  for line in keep_lines:
    facet_file.write( line )
  facet_file.close()

  return

###############################################################################

def replace_facet( facet_file_name, facet_id, facet_size, radec, keep_boxes = False ):
  
  if ( not file_exists( facet_file_name ) ):
    raise error( 'facet file %s does not exist' % facet_file_name )
  # scan box file for line(s), keep all other lines
  keep_lines = []
  facet_file = file( facet_file_name, mode = 'r' )
  for line in facet_file:
    words = [ word.strip() for word in line.split() ]
    if ( len( words ) == 0 ):
      continue
    if ( words[ 0 ][ 0 ] == '#' ):
      continue
    try:
      facet_number = int( words[ 0 ] )
    except:
      try:
        facet_number = int( words[ 1 ] )
      except:
        keep_lines.append( line )
      else:
        if ( ( words[ 0 ] != 'C' ) or ( facet_number != facet_id ) ):
          keep_lines.append( line )
        else:
          ra_dec = degdeg_to_hmsdms( radec, precision = [ 3, 2 ] )
          if ( asign( ra_dec[ 3 ] ) > 0. ):
            facet_line = ( 'C  %d  %d %d  %02d %02d %06.3f   %02d %02d %05.2f\n' % 
                ( facet_id, facet_size[ 0 ] - 12, facet_size[ 1 ] - 12,
                ra_dec[ 0 ], ra_dec[ 1 ], ra_dec[ 2 ], ra_dec[ 3 ], ra_dec[ 4 ], ra_dec[ 5 ] ) )
          else:
            facet_line = ( 'C  %d  %d %d  %02d %02d %06.3f  -%02d %02d %05.2f\n' % 
                ( facet_id, facet_size[ 0 ] - 12, facet_size[ 1 ] - 12,
                ra_dec[ 0 ], ra_dec[ 1 ], ra_dec[ 2 ], - ra_dec[ 3 ], ra_dec[ 4 ], ra_dec[ 5 ] ) )
          keep_lines.append( facet_line )
    else:
      if ( keep_boxes or ( facet_number != facet_id ) ):
        keep_lines.append( line )
  
  # write kept lines to box file
  if ( len( keep_lines ) > 0 ):
    remove_file( facet_file_name )
    facet_file = file( facet_file_name, mode = 'w' )
    for line in keep_lines:
      facet_file.write( line )
    facet_file.close()
  
  return

###############################################################################

def get_facet_count( facet_file_name, count_gaps = False ):
  facet_count = 0
  facet_list = get_facet_list( facet_file_name )
  if count_gaps:
    for facet in facet_list:
      if ( facet[ 0 ] > facet_count ):
        facet_count = facet[ 0 ]
  else:
    facet_count = len( facet_list )
  return facet_count

###############################################################################

def get_facet_list( facet_file_name ):
  radec_list = []
  if file_exists( facet_file_name ):
    facet_file = file( facet_file_name, mode = 'r' )
    for line in facet_file:
      words = [ word.strip() for word in line.split() ]
      if ( len( words ) == 0 ):
        continue
      if ( words[ 0 ][ 0 ] == '#' ):
        continue
      if ( words[ 0 ] == 'C' ):
        facet_id = int( words[ 1 ] )
        hmsdms = [ float( word ) for word in words[ 4 : 10 ] ]
        radec_list.append( [ facet_id, hmsdms_to_degdeg( hmsdms ) ] )
      elif ( words[ 0 ] == 'F' ):
        raise error( 'interpretation of RA/DEC shift not implemented' )
    facet_file.close()
  return radec_list

###############################################################################

def get_history_size( uvim ):
  hi_version = uvim.table_highver( 'HI' )
  if ( hi_version == 0 ):
    raise error( 'AIPS HI table does not exist' )
  if ( hi_version > 1 ):
    raise error( 'AIPS HI table version > 1' )
  hi_table = uvim.history
  decade = 0.
  decade_found = False
  while ( not decade_found ):
    decade = decade + 1.
    index = int( pow( 10., decade ) )
    try:
      dummy = hi_table[ index ]
    except IndexError:
      decade_found = True
  index_low = int( pow( 10., decade - 1. ) )
  index_high = int( pow( 10., decade ) )
  size_found = False
  while ( not size_found ):
    index = int( around( float( index_low + index_high ) / 2. ) )
    try:
      dummy = hi_table[ index ]
    except IndexError:
      index_high = index
    else:
      index_low = index
    if ( index_high == index_low + 1 ):
      size_found = True

  return index_high

###############################################################################

def write_history( uvim, strings ):
  hi_version = uvim.table_highver( 'HI' )
  if ( hi_version == 0 ):
    raise error( 'AIPS HI table does not exist' )
  if ( hi_version > 1 ):
    raise error( 'AIPS HI table version > 1' )
  wiz_uvim = wizardry( uvim )
  hi_table = wiz_uvim.history
  for string in strings:
    hi_table.append( string )
  hi_table.close()
  del wiz_uvim
  return

###############################################################################

def clear_history_old( uvim ):
  hi_version = uvim.table_highver( 'HI' )
  if ( hi_version == 0 ):
    pass
#    raise error( 'AIPS HI table does not exist' )
  if ( hi_version > 1 ):
    raise error( 'AIPS HI table version > 1' )
  wiz_uvim = wizardry( uvim )
  hi_table = wiz_uvim.history
  History.PZap( hi_table._table, hi_table._err )
  del wiz_uvim
  new_table( uvim, 'HI' )
  write_history( uvim, [ 'CLEARED HISTORY' ] )
  return

###############################################################################

def clear_history( uvim ):
  wiz_uvim = wizardry( uvim )
  hi_table = wiz_uvim.history
  OErr.PInit( hi_table._err, taskLog = 'dummy.log' ) # suppress output
  History.PEdit( hi_table._table, 1, 0, hi_table._err )
  del wiz_uvim
  write_history( uvim, [ 'CLEARED HISTORY' ] )
  remove_file( 'dummy.log' )
  return

###############################################################################

def read_history( uvim, strings = [], count = 1, word_match = False ):
# strings = search terms
# count = 0: return all occurrences
# count = N: return first N occurrences
# count = -N: return last N occurrences in reverse order

  hi_version = uvim.table_highver( 'HI' )
  if ( hi_version == 0 ):
    raise error( 'AIPS HI table does not exist' )
  if ( hi_version > 1 ):
    raise error( 'AIPS HI table version > 1' )
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
    for index in range( history_size - 1, -1, -1 ):
      line = hi_table[ index ]
      if word_match:
        words = [ word.strip() for word in line.split() ]
      if len( strings ) == 0:
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

def get_antenna_positions( uv, antenna_version = 0 ):
  wiz_uv = wizardry( uv )
  an_table = wiz_uv.table( 'AN', antenna_version )
  frame = an_table.keywords[ 'FRAME' ]
  if ( frame == 'ITRF' ):
    rotate = False
  else: # default (UNKNOWN=WGS84)
    rotate = True
  try:
    hand = an_table.keywords[ 'XYZHAND' ]
    hand = hand.strip()
  except:
    hand = 'RIGHT'
    array_name = an_table.keywords[ 'ARRNAM' ]
    if ( type( array_name ) == type( [] ) ):
      array_name = array_name[ 0 ]
    array_name = array_name.strip()
    if ( ( array_name != 'ATCA' ) and ( array_name != 'ATLBA' ) ):
      hand = 'LEFT'
  antenna_list = []
  array_xyz = array( [ an_table.keywords[ 'ARRAYX' ], an_table.keywords[ 'ARRAYY' ],
      an_table.keywords[ 'ARRAYZ' ] ], dtype = float64 )
  if alltrue( array_xyz == azeros( array_xyz ) ):
    for row in an_table:
      antenna_no = row.nosta
      antenna_xyz = array( row.stabxyz, dtype = float64 )
      antenna_list.append( [ antenna_no, antenna_xyz.tolist() ] )
  else:
    array_geo_llh = array( xyz_to_geo_llh( array_xyz.tolist() ), dtype = float64 )
    [ lon, lat ] = aradians( array_geo_llh[ 0 : 2 ] )
    rotation = array( [ [   cos( lon ), sin( lon ), 0. ],
                        [ - sin( lon ), cos( lon ), 0. ],
                        [           0.,         0., 1. ] ], dtype = float64 )
    for row in an_table:
      antenna_no = row.nosta
      dxyz = array( row.stabxyz, dtype = float64 )
      if ( rotate ):
        dxyz = dot( dxyz, rotation )
      antenna_xyz = array_xyz + dxyz
      antenna_list.append( [ antenna_no, antenna_xyz.tolist() ] )
    an_table.close()
  antenna_list.sort( cmp = lambda a, b: cmp( a[ 0 ], b[ 0 ] ) )
  del wiz_uv
  
  # fix Y-coordinate if needed
  if ( hand == 'LEFT' ):
    print 'WARNING: swapping antenna E-W sign'
    for antenna in antenna_list:
      antenna[ 1 ][ 1 ] = -antenna[ 1 ][ 1 ]
  
  return antenna_list

###############################################################################

def get_mean_antenna_position( uv, antenna_version = 0 ):
  antenna_list = get_antenna_positions( uv, antenna_version )
  antenna_array = transpose( array( [ antenna[ 1 ] for antenna in antenna_list ],
      dtype = float64 ) )
  array_xyz = [ float( ( antenna_array[ 0 ] ).mean() ), float( ( antenna_array[ 1 ] ).mean() ),
      float( ( antenna_array[ 2 ] ).mean() ) ]
  return array_xyz

###############################################################################

def get_local_antenna_positions( uv, antenna_version = 0 ):
  antenna_list = []
  array_xyz = array( get_mean_antenna_position( uv,
      antenna_version = antenna_version ), dtype = float64 )
  array_geo_llh = array( xyz_to_geo_llh( array_xyz.tolist() ), dtype = float64 )
  [ lon, lat ] = aradians( array_geo_llh[ 0 : 2 ] )
  rotation = array( [ [ 0., - sin( lat ), cos( lat ) ], 
                      [ 1.,           0.,         0. ],
                      [ 0.,   cos( lat ), sin( lat ) ] ], dtype = float64 )
  wiz_uv = wizardry( uv )
  an_table = wiz_uv.table( 'AN', antenna_version )
  for row in an_table:
    antenna_no = row.nosta
    antenna_xyz = array( row.stabxyz, dtype = float64 )
    local_xyz = dot( antenna_xyz, rotation )
    antenna_list.append( [ antenna_no, local_xyz.tolist() ] )
  an_table.close()
  antenna_list.sort( cmp = lambda a, b: cmp( a[ 0 ], b[ 0 ] ) )
  del wiz_uv
  return antenna_list

###############################################################################

def get_antenna_names( uv, version = 0 ):
  antenna_list = []
  wiz_uv = wizardry( uv )
  an_table = wiz_uv.table( 'AN', version )
  for row in an_table:
    antenna_no = row.nosta
    ant_name = row.anname
    antenna_list.append( [ antenna_no, ant_name ] )
  an_table.close()
  antenna_list.sort( cmp = lambda a, b: cmp( a[ 0 ], b[ 0 ] ) )
  del wiz_uv
  return antenna_list

###############################################################################

def get_gst_list( uv, time_list = [ 0. ] ):
# time in days since reference date
# TODO: check if implementation below is correct
  wiz_uv = wizardry( uv )
  an_table = wiz_uv.table( 'AN', 0 )
  gst0 = an_table.keywords[ 'GSTIA0' ] # degrees
  deg_per_day = an_table.keywords[ 'DEGPDY' ]
  time_system = an_table.keywords[ 'TIMSYS' ]
  if ( type( time_system ) == type( [] ) ):
    time_system = time_system[ 0 ]
  time_system = time_system.strip()
  try:
    dat_utc = an_table.keywords[ 'DATUTC' ] # sec
  except:
    dat_utc = 0.
  if ( time_system == 'UTC' ):
    iat_utc = an_table.keywords[ 'IATUTC' ] # sec
  elif ( ( time_system == 'IAT' ) or ( time_system == 'TAI' ) ):
    iat_utc = 0.
  else:
    raise error( 'time system %s unknown' % ( time_system ) )
  day_offset_array = array( time_list, dtype = float64 ) + ( iat_utc - dat_utc ) / 86400.
  gst_array = amodulo( day_offset_array * deg_per_day + gst0, 360. )
  an_table.close()
  del wiz_uv
  return gst_array.tolist()

###############################################################################

def calculate_rise_transit_set_times( uv, radec = None, alt_limit = 0. ):
  array_xyz = array( get_mean_antenna_position( uv ), dtype = float64 )
  [ lon, lat, h ] = xyz_to_geo_llh( array_xyz.tolist() )
  if ( radec is None ):
    ra_dec = get_radec( uv )
  else:
    ra_dec = radec
  obs_epoch = get_observing_epoch( uv )
  [ ra, dec ] = convert_radec_from_j2000( ra_dec, obs_epoch )
  gst0 = ( ( get_gst_list( uv ) )[ 0 ] ) % 360. # GST at midnight on reference date
  ha0 = ( gst0 - ra + lon ) % 360. # HA of object at midnight on ref date
  time_transit = ( ( - ha0 ) % 360. ) / 360. # transit occurs at hour angle 0 by default
  ha = calculate_hour_angles_at_elevation_limit( lat, dec, elevation_limit = alt_limit )
  time_rise = ( ( ha[ 0 ] - ha0 ) % 360. ) / 360.
  time_set = ( ( ha[ 1 ] - ha0 ) % 360. ) / 360.
  return [ time_rise, time_transit, time_set ]

###############################################################################

def extract_facet_definitions( facet_file_name, facet_list, new_facet_file_name, 
    include_clean_boxes = True, new_clean_box_radius = -1, include_bcomp = True,
    include_weights = True, renumber_facets = True, invert_selection = False ):
  
  if ( not file_exists( facet_file_name ) ):
    raise error( 'facet file %s does not exist' % ( facet_file_name ) )
  facet_file_lines = []
  facet_file = file( facet_file_name, mode = 'r' )
  for line in facet_file:
    facet_file_lines.append( line )
  facet_file.close()
  
  if file_exists( new_facet_file_name ):
    remove_file( new_facet_file_name )
  new_facet_file = file( new_facet_file_name, mode = 'w' )
  for line in facet_file_lines:
    words = [ word.strip() for word in line.split() ]
    if ( len( words ) == 0 ):
      continue
    if ( words[ 0 ][ 0 ] == '#' ):
      continue
    if ( words[ 0 ] == 'C' ):
      facet_number = int( words[ 1 ] )
      if ( ( ( not invert_selection ) and ( facet_number in facet_list ) ) or 
          ( invert_selection and ( not facet_number in facet_list ) ) ):
        if renumber_facets:
          new_facet_number = facet_list.index( facet_number ) + 1
        else:
          new_facet_number = facet_number
        new_line = '%s  %s  %s %s  %s %s %s  %s %s %s\n' % ( words[ 0 ], 
            repr( new_facet_number ), words[ 2 ], words[ 3 ], words[ 4 ],
            words[ 5 ], words[ 6 ], words[ 7 ], words[ 8 ], words[ 9 ] )
        new_facet_file.write( new_line )
    elif ( words[ 0 ] == 'F' ):
      facet_number = int( words[ 1 ] )
      if ( ( ( not invert_selection ) and ( facet_number in facet_list ) ) or 
          ( invert_selection and ( not facet_number in facet_list ) ) ):
        if renumber_facets:
          new_facet_number = facet_list.index( facet_number ) + 1
        else:
          new_facet_number = facet_number
        new_line = '%s  %s  %s %s  %s %s\n' % ( words[ 0 ], repr( new_facet_number ),
            words[ 2 ], words[ 3 ], words[ 4 ], words[ 5 ] )
        new_facet_file.write( new_line )
    elif ( words[ 0 ] == 'B' ):
      if include_bcomp:
        facet_number = int( words[ 1 ] )
        if ( ( ( not invert_selection ) and ( facet_number in facet_list ) ) or 
            ( invert_selection and ( not facet_number in facet_list ) ) ):
          if renumber_facets:
            new_facet_number = facet_list.index( facet_number ) + 1
          else:
            new_facet_number = facet_number
          new_line = '%s  %s  %s\n' % ( words[ 0 ], repr( new_facet_number ), words[ 2 ] )
          new_facet_file.write( new_line )
    elif ( words[ 0 ] == 'W' ):
      if include_weights:
        new_facet_file.write( line )
    else: # must be a clean box definition
      if include_clean_boxes:
        facet_number = int( words[ 0 ] )
        if ( ( ( not invert_selection ) and ( facet_number in facet_list ) ) or 
            ( invert_selection and ( not facet_number in facet_list ) ) ):
          if renumber_facets:
            new_facet_number = facet_list.index( facet_number ) + 1
          else:
            new_facet_number = facet_number
          if ( int( words[ 1 ] ) == - 1 ):
            if ( new_clean_box_radius > 0 ):
              new_line = '%05d  %s  %d  %s %s\n' % ( int( new_facet_number ), 
                  words[ 1 ], int( new_clean_box_radius ), words[ 3 ], words[ 4 ] )
            else:
              new_line = '%05d  %s  %s  %s %s\n' % ( int( new_facet_number ), 
                  words[ 1 ], words[ 2 ], words[ 3 ], words[ 4 ] )
          else:
            new_line = '%05d  %s %s  %s %s\n' % ( int( new_facet_number ), 
                words[ 1 ], words[ 2 ], words[ 3 ], words[ 4 ] )
          new_facet_file.write( new_line )
  new_facet_file.close()
  return

###############################################################################

def merge_facet_definitions_old( facet_file_name_1, facet_file_name_2,
    new_facet_file_name ):
  
  # open facet files
  facet_file_1 = file( facet_file_name_1, mode = 'r' )
  facet_file_2 = file( facet_file_name_2, mode = 'r' )
  if file_exists( new_facet_file_name ):
    remove_file( new_facet_file_name )
  new_facet_file = file( new_facet_file_name, mode = 'w' )
  
  # copy contents of first facet file
  # look for highest facet number
  facet_count_1 = 0
  for line in facet_file_1:
    words = [ word.strip() for word in line.split() ]
    if ( len( words ) == 0 ):
      continue
    if ( words[ 0 ][ 0 ] == '#' ):
      continue
    if ( words[ 0 ] == 'C' ):
      facet_number = int( words[ 1 ] )
      if ( facet_number > facet_count_1 ):
        facet_count_1 = facet_number
    elif ( words[ 0 ] == 'F' ):
      facet_number = int( words[ 1 ] )
      if ( facet_number > facet_count_1 ):
        facet_count_1 = facet_number
    new_facet_file.write( line )
  facet_file_1.close()
  
  # copy contents of second facet file
  # adjust facet numbers
  facet_count_2 = 0
  for line in facet_file_2:
    words = [ word.strip() for word in line.split() ]
    if ( len( words ) == 0 ):
      continue
    if ( words[ 0 ][ 0 ] == '#' ):
      continue
    if ( words[ 0 ] == 'C' ):
      facet_number = int( words[ 1 ] )
      if ( facet_number > facet_count_2 ):
        facet_count_2 = facet_number
      new_facet_number = facet_count_1 + facet_number
      new_line = '%s  %s  %s %s  %s %s %s  %s %s %s\n' % ( words[ 0 ], 
          repr( new_facet_number ), words[ 2 ], words[ 3 ], words[ 4 ],
          words[ 5 ], words[ 6 ], words[ 7 ], words[ 8 ], words[ 9 ] )
    elif ( words[ 0 ] == 'F' ):
      facet_number = int( words[ 1 ] )
      if ( facet_number > facet_count_2 ):
        facet_count_2 = facet_number
      new_facet_number = facet_count_1 + facet_number
      new_line = '%s  %s  %s %s  %s %s\n' % ( words[ 0 ], repr( new_facet_number ),
          words[ 2 ], words[ 3 ], words[ 4 ], words[ 5 ] )
    elif ( ( words[ 0 ] == 'B' ) or ( words[ 0 ] == 'W' ) ):
      new_line = line
    else: # must be a clean box definition
      facet_number = int( words[ 0 ] )
      new_facet_number = facet_count_1 + facet_number
      if ( int( words[ 1 ] ) == - 1 ):
        new_line = '%05d  %s  %s  %s %s\n' % ( int( new_facet_number ), 
            words[ 1 ], words[ 2 ], words[ 3 ], words[ 4 ] )
      else:
        new_line = '%05d  %s %s  %s %s\n' % ( int( new_facet_number ), 
            words[ 1 ], words[ 2 ], words[ 3 ], words[ 4 ] )
    new_facet_file.write( new_line )
  facet_file_2.close()
  
  new_facet_file.close()
  
  return [ facet_count_1, facet_count_2 ]

###############################################################################

def merge_facet_definitions( facet_file_name_1, facet_file_name_2,
    new_facet_file_name, renumber_facets = True, include_clean_boxes = True,
    include_bcomp = True, include_weights = True ):
  
  # check for presence facet files
  if ( not file_exists( facet_file_name_1 ) ):
    raise error( 'facet file %s does not exist' % ( facet_file_name_1 ) )
  if ( not file_exists( facet_file_name_2 ) ):
    raise error( 'facet file %s does not exist' % ( facet_file_name_2 ) )
  
  # read facet files
  facet_file_1_lines = []
  facet_file_1 = file( facet_file_name_1, mode = 'r' )
  for line in facet_file_1:
    facet_file_1_lines.append( line )
  facet_file_1.close()
  facet_file_2_lines = []
  facet_file_2 = file( facet_file_name_2, mode = 'r' )
  for line in facet_file_2:
    facet_file_2_lines.append( line )
  facet_file_2.close()
  new_facet_file_tail_lines = []
  
  # read definitions from first facet file
  # copy lines from first facet file to output file
  facet_ids_1 = []
  facet_defs_1 = []
  box_ids_1 = []
  box_defs_1 = []
  for line in facet_file_1_lines:
    words = [ word.strip() for word in line.split() ]
    if ( len( words ) == 0 ):
      continue
    elif ( words[ 0 ][ 0 ] == '#' ):
      continue
    elif ( words[ 0 ][ 0 ] == 'B' ):
      if include_bcomp:
        new_facet_file_tail_lines.append( line )
      continue
    elif ( words[ 0 ][ 0 ] == 'W' ):
      if include_weights:
        new_facet_file_tail_lines.append( line )
      continue
    elif ( words[ 0 ] == 'C' ):
      facet_id = int( words[ 1 ] )
      facet_def = [ words[ i ] for i in [ 0,2,3,4,5,6,7,8,9 ] ]
      facet_ids_1.append( facet_id )
      facet_defs_1.append( facet_def )
    elif ( words[ 0 ] == 'F' ):
      facet_id = int( words[ 1 ] )
      facet_def = [ words[ i ] for i in [ 0,2,3,4,5 ] ]
      facet_ids_1.append( facet_id )
      facet_defs_1.append( facet_def )
    else:
      try:
        box_id = int( words[ 0 ] )
        box_def = [ int( x ) for x in words[ 1 : 5 ] ]
      except ValueError:
        continue
      else:
        if ( box_def == [ 0,0,0,0 ] ):
          continue
        if ( box_def == [ -1,1,5,5 ] ):
          continue
        box_ids_1.append( box_id )
        box_defs_1.append( box_def )
  
  # read definitions from second facet file
  facet_ids_2 = []
  facet_defs_2 = []
  box_ids_2 = []
  box_defs_2 = []
  for line in facet_file_2_lines:
    words = [ word.strip() for word in line.split() ]
    if ( len( words ) == 0 ):
      continue
    elif ( words[ 0 ][ 0 ] == '#' ):
      continue
    elif ( words[ 0 ][ 0 ] == 'B' ):
      if include_bcomp:
        new_facet_file_tail_lines.append( line )
      continue
    elif ( words[ 0 ][ 0 ] == 'W' ):
      if include_weights:
        new_facet_file_tail_lines.append( line )
      continue
    elif ( words[ 0 ] == 'C' ):
      facet_id = int( words[ 1 ] )
      facet_def = [ words[ i ] for i in [ 0,2,3,4,5,6,7,8,9 ] ]
      facet_ids_2.append( facet_id )
      facet_defs_2.append( facet_def )
    elif ( words[ 0 ] == 'F' ):
      facet_id = int( words[ 1 ] )
      facet_def = [ words[ i ] for i in [ 0,2,3,4,5 ] ]
      facet_ids_2.append( facet_id )
      facet_defs_2.append( facet_def )
    else:
      try:
        box_id = int( words[ 0 ] )
        box_def = [ int( x ) for x in words[ 1 : 5 ] ]
      except ValueError:
        continue
      else:
        if ( box_def == [ 0,0,0,0 ] ):
          continue
        if ( box_def == [ -1,1,5,5 ] ):
          continue
        box_ids_2.append( box_id )
        box_defs_2.append( box_def )
  
  # renumber facets
  if renumber_facets:
    new_facet_ids_1 = range( 1, 1 + len( facet_ids_1 ) )
    new_facet_ids_2 = []
    new_facet_id_2 = 1 + max( new_facet_ids_1 )
    for i in range( len( facet_ids_2 ) ):
      facet_def_2 = facet_defs_2[ i ]
      if ( facet_def_2 in facet_defs_1 ):
        index = facet_defs_1.index( facet_def_2 )
        new_facet_ids_2.append( new_facet_ids_1[ index ] )
      else:
        new_facet_ids_2.append( new_facet_id_2 )
        new_facet_id_2 = new_facet_id_2 + 1
  else:
    # check for continuity of facet numbers
    facet_ids = list( set( facet_ids_1 + facet_ids_2 ) )
    if ( facet_ids != range( 1, 1 + len( facet_ids ) ) ):
      raise error( 'facet numbers from input facet files do not form a continuous series' )
    # check for duplicate facet numbers
    if ( len( facet_ids ) != len( facet_ids_1 + facet_ids_2 ) ):
      for facet_id_2 in facet_ids_2:
        if ( facet_id_2 in facet_ids_1 ):
          facet_def_1 = facet_defs_1[ facet_ids_1.index( facet_id_2 ) ]
          facet_def_2 = facet_defs_2[ facet_ids_2.index( facet_id_2 ) ]
          if ( facet_def_1 != facet_def_2 ):
            raise error( 'duplicate facet numbers detected in input facet files' )
    new_facet_ids_1 = [ x for x in facet_ids_1 ]
    new_facet_ids_2 = [ x for x in facet_ids_2 ]
  
  # build merged facet file
  if file_exists( new_facet_file_name ):
    remove_file( new_facet_file_name )
  new_facet_file = file( new_facet_file_name, mode = 'w' )
  new_facet_count = max( new_facet_ids_1 + new_facet_ids_2 )
  
  # write facet_definitions
  for i in range( 1, 1 + new_facet_count ):
    if ( i in new_facet_ids_1 ):
      idx = new_facet_ids_1.index( i )
      facet_def = facet_defs_1[ idx ]
    else:
      idx = new_facet_ids_2.index( i )
      facet_def = facet_defs_2[ idx ]
    if ( facet_def[ 0 ] == 'C' ):
      new_line = '%s  %d  %s %s  %s %s %s   %s %s %s\n' % ( facet_def[ 0 ], i,
          facet_def[ 1 ], facet_def[ 2 ], facet_def[ 3 ], facet_def[ 4 ],
          facet_def[ 5 ], facet_def[ 6 ], facet_def[ 7 ], facet_def[ 8 ] )
    else: # 'F'
      new_line = '%s  %d  %s %s  %s %s\n' % ( facet_def[ 0 ], i,
          facet_def[ 1 ], facet_def[ 2 ], facet_def[ 3 ], facet_def[ 4 ] )
    new_facet_file.write( new_line )
  
  # write clean box_definitions
  for i in range( 1, 1 + new_facet_count ):
    facet_has_box = False
    if ( include_clean_boxes and ( i in new_facet_ids_1 ) ):
      idx = new_facet_ids_1.index( i )
      facet_id = facet_ids_1[ idx ]
      if ( facet_id in box_ids_1 ):
        facet_has_box = True
        for j in range( len( box_ids_1 ) ):
          if ( box_ids_1[ j ] == facet_id ):
            box_def = box_defs_1[ j ]
            if ( box_def[ 0 ] == -1 ): # circular
              new_line = '%05d  %s  %s  %s %s\n' % ( i,
                  box_def[ 0 ], box_def[ 1 ], box_def[ 2 ], box_def[ 3 ] )
            else:
              new_line = '%05d  %s %s  %s %s\n' % ( i,
                  box_def[ 0 ], box_def[ 1 ], box_def[ 2 ], box_def[ 3 ] )
            new_facet_file.write( new_line )
    if ( include_clean_boxes and ( i in new_facet_ids_2 ) ):
      idx = new_facet_ids_2.index( i )
      facet_id = facet_ids_2[ idx ]
      if ( facet_id in box_ids_2 ):
        facet_has_box = True
        for j in range( len( box_ids_2 ) ):
          if ( box_ids_2[ j ] == facet_id ):
            box_def = box_defs_2[ j ]
            if ( box_def[ 0 ] == -1 ): # circular
              new_line = '%05d  %s  %s  %s %s\n' % ( i,
                  box_def[ 0 ], box_def[ 1 ], box_def[ 2 ], box_def[ 3 ] )
            else:
              new_line = '%05d  %s %s  %s %s\n' % ( i,
                  box_def[ 0 ], box_def[ 1 ], box_def[ 2 ], box_def[ 3 ] )
            new_facet_file.write( new_line )
    if ( not facet_has_box ):
      new_line = '%05d  %d  %d  %d %d\n' % ( i, 0, 0, 0, 0 )
      new_facet_file.write( new_line )
  
  # write remaining lines
  for new_line in new_facet_file_tail_lines:
    new_facet_file.write( new_line )
  new_facet_file.close()
  
  return

###############################################################################

def fill_facet( facet, facet_file_name = '', invert = False, blank_edge = 6,
    do_edge_circle = False, value = get_aips_magic_value() ):

  # copy facet pixel contents
  pixels = get_image_pixels( facet )
  [ X, Y ] = indices( tuple( facet_size ), dtype = int32 )
  mask = zeros( shape = pixels.shape, dtype = bool )
  if ( facet_file_name != '' ):
    # extract facet definition and clean boxes
    i = get_facet_number( facet )
    facet_i_file_name = facet_file_name + '.%03d' % ( i )
    extract_facet_definitions( facet_file_name, [ i ], facet_i_file_name,
        include_bcomp = False, include_weights = False )

    # read clean box definitions from facet file
    # fill mask with ones in clean box areas
    facet_i_file = file( facet_i_file_name, mode = 'r' )
    for line in facet_i_file:
      words = [ word.strip() for word in line.split() ]
      if ( len( words ) == 0 ):
        continue
      if ( words[ 0 ][ 0 ] == '#' ):
        continue
      try:
        facet_number = int( words[ 0 ] )
      except ValueError:
        pass
      else:
        if ( int( words[ 1 ] ) == - 1 ): # circular area
          rc = int( words[ 2 ] )
          xc = int( words[ 3 ] ) - 1
          yc = int( words[ 4 ] ) - 1
          rc2 = float( rc )**2
          R2 = ( X - float( xc ) )**2 + ( Y - float( yc ) )**2
          mask = aput( mask, awhere( R2 <= rc2 ), True )
        else: # rectangular area
          x1 = int( words[ 1 ] ) - 1
          y1 = int( words[ 2 ] ) - 1
          x2 = int( words[ 3 ] ) - 1
          y2 = int( words[ 4 ] ) - 1
          sel = awhere( ( X >= x1 ) & ( X <= x2 ) & ( Y >= y1 ) & ( Y <= y2 ) )
          mask = aput( mask, sel, True )
    facet_i_file.close()
    remove_file( facet_i_file_name )

  # invert mask when requested
  if invert:
    mask = ( mask == False )

  # exclude edge pixels of facet
  if ( blank_edge > 0 ):
    if do_edge_circle:
      # do circle
      xc = float( facet_size[ 0 ] - 1 ) / 2.
      yc = float( facet_size[ 1 ] - 1 ) / 2.
      rc2 = ( max( 0., min( [ xc, yc ] ) - float( blank_edge ) ) )**2
      if ( rc > 0. ):
        rc2 = rc**2
        R2 = ( X - float( xc ) )**2 + ( Y - float( yc ) )**2
        mask = aput( mask, awhere( R2 > rc2 ), True )
      else:
        raise error( 'size of blanking circle is too large' )
    else:
      # do simple rectangle
      mask[ 0 : blank_edge, : ] = True
      mask[ facet_size[ 0 ] - blank_edge : facet_size[ 0 ], : ] = True
      mask[ : , 0 : blank_edge ] = True
      mask[ : , facet_size[ 1 ] - blank_edge : facet_size[ 1 ] ] = True

  # write blanks at masked values
  pixels = aput( pixels, awhere( mask ), value )
  set_image_pixels( facet, pixels )
  del pixels, mask

  return

###############################################################################

def fill_image( image, value = get_aips_magic_value(), do_all = False, 
    do_edge = False, do_edge_circle = False, edge_size = 8 ):
  
  # get image pixels
  pixels = get_image_pixels( image )
  
  # do all pixels
  if do_all:
    pixels[ : , : ] = value

  # do circle
  if do_edge_circle:
    image_size = list( pixels.shape )
    [ X, Y ] = indices( tuple( image_size ), dtype = float64 )
    xc = float( image_size[ 0 ] - 1 ) / 2.
    yc = float( image_size[ 1 ] - 1 ) / 2.
    rc2 = ( max( min( [ xc, yc ] ) - float( edge_size ), 0. ) )**2
    R2 = ( X - float( xc ) )**2 + ( Y - float( yc ) )**2
    pixels = aput( pixels, awhere( R2 >= rc2 ), value )
    del X, Y, R2
  
  # do simple rectangle
  if do_edge:
    image_size = list( pixels.shape )
    pixels[ : min( edge_size, image_size[ 0 ] ), : ] = value
    pixels[ : , : min( edge_size, image_size[ 1 ] ) ] = value
    pixels[ max( image_size[ 0 ] - edge_size, 0 ) : , : ] = value
    pixels[ : , max( image_size[ 1 ] - edge_size, 0 ) : ] = value

  # write pixels back to image
  set_image_pixels( image, pixels )
  del pixels

  return

###############################################################################

def fill_facets( facets, facet_list, invert = False, do_boxes = False, 
    facet_file_name = '', do_edge_circle = False, edge_size = 8,
    value = get_aips_magic_value() ):
  
  # process input parameters
  if do_boxes:
    if ( facet_file_name == '' ):
      raise error( 'facet file name not specified' )
    clean_box_list = get_clean_boxes( facet_file_name, facet_list = facet_list )
  else:
    clean_box_list = []  
  
  # loop over facets
  for i in facet_list:
    
    # create image mask
    facet = get_facet( facets, i )
    facet_size = get_image_size( facet )
    pixels = get_image_pixels( facet )
    [ X, Y ] = indices( tuple( facet_size ), dtype = float64 )
    mask = zeros( shape = pixels.shape, dtype = bool )
    
    # search through clean box list
    for [ j, a, b, c, d ] in clean_box_list:
      if ( j != i ):
        continue  
      # mask circular area
      if ( a == -1 ):
        [ rc2,xc,yc ] = [ float( b )**2, float( c - 1 ), float( d - 1 ) ]
        R2 = ( X - xc )**2 + ( Y - yc )**2
        sel = awhere( R2 <= rc2 )
        mask = aput( mask, sel, True )
        del R2
      # mask rectangular area
      else:
        [ x1,y1,x2,y2 ] = [ float( a - 1 ), float( b - 1 ), float( c - 1 ), float( d - 1 ) ]
        sel = awhere( ( X >= x1 ) & ( X <= x2 ) & ( Y >= y1 ) & ( Y <= y2 ) )
        mask = aput( mask, sel, True )
      
    # invert mask when requested
    if invert:
      mask = ( mask == False )

    # exclude edge pixels of facet
    if ( edge_size > 0 ):
      if ( float( edge_size ) > float( min( facet_size ) - 1 ) / 2. ):
        raise error( 'edge size is larger than image radius' )
      # do edge circle
      if do_edge_circle:
        xc = float( facet_size[ 0 ] - 1 ) / 2.
        yc = float( facet_size[ 1 ] - 1 ) / 2.
        rc2 = ( max( 0., min( [ xc, yc ] ) - float( edge_size ) ) )**2
        R2 = ( X - float( xc ) )**2 + ( Y - float( yc ) )**2
        sel = awhere( R2 >= rc2 )
        mask = aput( mask, sel, True )
        del R2
      else:
        # do simple rectangle
        mask[ 0 : edge_size, : ] = True
        mask[ facet_size[ 0 ] - edge_size : facet_size[ 0 ], : ] = True
        mask[ : , 0 : edge_size ] = True
        mask[ : , facet_size[ 1 ] - edge_size : facet_size[ 1 ] ] = True
    
    # write blanks at masked values
    pixels = aput( pixels, awhere( mask ), value )
    set_image_pixels( facet, pixels )
    del pixels, X, Y, mask
  
  return

###############################################################################

def get_pixel_value( image, pos, to_float = True ):
  image_size = get_image_size( image )
  [ x, y ] = [ int( around( pos[ 0 ] ) ) - 1, int( around( pos[ 1 ] ) ) - 1 ]
  if ( ( x < 0 ) or ( x > image_size[ 0 ] - 1 ) or ( y < 0 ) or ( y > image_size[ 1 ] - 1 ) ):
    value = get_aips_magic_value()
  else:
    pixels = get_image_pixels( image )
    value = pixels[ x, y ]
  if to_float:
    value = float( value )
  return value

###############################################################################

def set_pixel_value( image, pos, value ):
  image_size = get_image_size( image )
  [ x, y ] = [ int( around( pos[ 0 ] ) ) - 1, int( around( pos[ 1 ] ) ) - 1 ]
  if ( ( x < 0 ) or ( x > image_size[ 0 ] - 1 ) or ( y < 0 ) or ( y > image_size[ 1 ] - 1 ) ):
    raise error( 'specified position outside image domain' )
  pixels = get_image_pixels( image )
  pixels[ x, y ] = value
  set_image_pixels( image, pixels )
  return

###############################################################################

def get_image_minimum( im ):
  avg = get_image_mean( im )
  if ( not avg is None ):
    pixels = get_image_pixels( im )
    mask = ( pixels == get_aips_magic_value() )
    if sometrue( mask ):
      pixels = aput( pixels, awhere( mask ), avg )
    im_min = pixels.min()
    xy_min = where( pixels == im_min )
    x_min = xy_min[ 0 ][ 0 ] + 1
    y_min = xy_min[ 1 ][ 0 ] + 1
    return [ float( im_min ), [ x_min, y_min ] ]
  else:
    return None
    
###############################################################################

def get_image_maximum( im ):
  avg = get_image_mean( im )
  if ( not avg is None ):
    pixels = get_image_pixels( im )
    mask = ( pixels == get_aips_magic_value() )
    if sometrue( mask ):
      pixels = aput( pixels, awhere( mask ), avg )
    im_max = pixels.max()
    xy_max = where( pixels == im_max )
    x_max = xy_max[ 0 ][ 0 ] + 1
    y_max = xy_max[ 1 ][ 0 ] + 1
    return [ float( im_max ), [ x_max, y_max ] ]
  else:
    return None

###############################################################################

def get_image_extremum( im, force_positive = False ):
  im_max = get_image_maximum( im )
  if ( ( im_max is None ) or force_positive ):
    extremum = im_max
  else:
    im_min = get_image_minimum( im )
    if ( abs( im_min[ 0 ] ) > abs( im_max[ 0 ] ) ):
      extremum = im_min
    else:
      extremum = im_max
  return extremum

###############################################################################

def extract_model_components( facets, facet_id, facet_file_name, model_version = 0 ):
  
  # process inputs
  facet = get_facet( facets, facet_id )
  if ( not table_exists( facet, 'CC', model_version ) ):
    raise error( 'input model component table does not exist' )
  if ( model_version == 0 ):
    in_version = facet.table_highver( 'CC' )
  else:
    in_version = model_version
  out_version = facet.table_highver( 'CC' ) + 1
  
  # create clean box mask
  # fill mask with ones in clean box areas
  facet_size = get_image_size( facet )
  clean_boxes = get_clean_boxes( facet_file_name, facet_list = [ facet_id ],
      drop_dummy_boxes = True )
  [ X, Y ] = indices( tuple( facet_size ), dtype = int32 )
  mask = zeros( shape = tuple( facet_size ), dtype = bool )
  for [ i,a,b,c,d ] in clean_boxes:
    if ( a == -1 ):
      xc = float( c - 1 )
      yc = float( d - 1 )
      rc2 = float( b )**2
      R2 = ( X - float( xc ) )**2 + ( Y - float( yc ) )**2
      sel = awhere( R2 <= rc2 )
      mask = aput( mask, sel, True )
    else: # rectangular area
      x1 = a - 1
      y1 = b - 1
      x2 = c - 1
      y2 = d - 1
      sel = awhere( ( X >= x1 ) & ( X <= x2 ) & ( Y >= y1 ) & ( Y <= y2 ) )
      mask = aput( mask, sel, True )
  
  # extract clean components within clean box(es)
  pix_scale = get_pixel_size( facet, make_absolute = False )
  pix_ref = get_pixel_reference( facet )
  wiz_facet = wizardry( facet )
  in_cc_table = wiz_facet.table( 'CC', in_version )
  try:
    no_parms = in_cc_table.keywords[ 'NO_PARMS' ]
  except:
    no_parms = 0
  out_cc_table = new_table( facet, 'CC', out_version, no_parms = no_parms )
  for cc in in_cc_table:
    dx = 3600. * cc.deltax / pix_scale[ 0 ]
    dy = 3600. * cc.deltay / pix_scale[ 1 ]
    x = int( around( pix_ref[ 0 ] + dx ) ) - 1
    y = int( around( pix_ref[ 1 ] + dy ) ) - 1
    if mask[ x, y ]:
      out_cc = new_table_row( out_cc_table )
      out_cc.flux = cc.flux
      out_cc.deltax = cc.deltax
      out_cc.deltay = cc.deltay
      if ( no_parms > 0 ):
        out_cc.parms = cc.parms
      out_cc_table.append( out_cc )
  in_cc_table.close()
  out_cc_table.close()
  del wiz_facet
  
  return out_version

###############################################################################

def combine_model_tables( facet, model_versions = [] ):
  
  # make list of tables to merge
  if ( model_versions == [] ):
    in_versions = []
    for model_version in range( 1, 1 + facet.table_highver( 'CC' ) ):
      if table_exists( facet, 'CC', model_version ):
        in_versions.append( model_version )
  else:
    for model_version in model_versions:
      if ( not table_exists( facet, 'CC', model_version ) ):
        raise error( 'input model component table does not exist' )
    in_versions = model_versions
  out_version = facet.table_highver( 'CC' ) + 1
  
  # determine highest no_parms
  max_no_parms = 0
  wiz_facet = wizardry( facet )
  for version in in_versions:
    in_cc_table = wiz_facet.table( 'CC', version )
    try:
      max_no_parms = max( max_no_parms, in_cc_table.keywords[ 'NO_PARMS' ] )
    except:
      pass
    in_cc_table.close()
  
  # copy contents of selected tables to new table
  out_cc_table = new_table( facet, 'CC', out_version, no_parms = max_no_parms )
  for version in in_versions:
    in_cc_table = wiz_facet.table( 'CC', version )
    no_parms = 0
    if ( max_no_parms > 0 ):
      try:
        no_parms = in_cc_table.keywords[ 'NO_PARMS' ]
      except:
        pass
    for cc in in_cc_table:
      if ( cc.flux <= 0. ):
        continue
      out_cc = new_table_row( out_cc_table )
      out_cc.flux = cc.flux
      out_cc.deltax = cc.deltax
      out_cc.deltay = cc.deltay
      if ( no_parms > 0 ):
        out_cc.parms = cc.parms + [ 0. for i in range( no_parms, max_no_parms ) ]
      out_cc_table.append( out_cc )
    in_cc_table.close()
  out_cc_table.close()
  del wiz_facet
  
  return out_version

###############################################################################

def get_time_list( uv ):
  time_list = []
  wiz_uv = wizardry( uv )
  for group in wiz_uv:
    if ( len( time_list ) == 0 ):
      time_list.append( group.time )
    elif ( group.time != time_list[ - 1 ] ):
      time_list.append( group.time )
  del wiz_uv
  return time_list

###############################################################################

def get_reference_antenna( uvim, solution_version = 0 ):
  wiz_uvim = wizardry( uvim )
  sn_table = wiz_uvim.table( 'SN', solution_version )
  for row in sn_table:
    if ( row.weight_1 > 0. ):
      reference_antenna = row.refant_1
      break
  sn_table.close()
  del wiz_uvim
  return reference_antenna

###############################################################################

def convert_beam_size( im, beam = None, to_pixel = False ):
# TODO: add rotation for non-central beams
# TODO: check for mirroring when pixel_size < 0
  if ( beam is None ):
    [ bmaj, bmin, bpa ] = get_beam_size( im ) # beam in arcsec
    if not to_pixel:
      return [ bmaj, bmin, bpa ]
  else:
    [ bmaj, bmin, bpa ] = beam
    if beam == [ 0., 0., 0. ]:
      return beam
  [ scale_x, scale_y ] = get_pixel_size( im, make_absolute = True )
  if to_pixel: # convert from arcsec to pixels, else from pixels to arcsec
    scale_x = 1. / scale_x
    scale_y = 1. / scale_y
  new_bmaj_x = scale_x * sin( radians( bpa ) ) * bmaj
  new_bmaj_y = scale_y * cos( radians( bpa ) ) * bmaj
  new_bmaj = sqrt( new_bmaj_x**2 + new_bmaj_y**2 )
  new_bmin_x = scale_x * cos( radians( bpa ) ) * bmin
  new_bmin_y = scale_y * sin( radians( bpa ) ) * bmin
  new_bmin = sqrt( new_bmin_x**2 + new_bmin_y**2 )
  new_bpa = degrees( atan2( new_bmaj_x, new_bmaj_y ) ) % 180.
  return [ new_bmaj, new_bmin, new_bpa ]

###############################################################################

def put_sub_image( im, subim, blc ):
  im_size = get_image_size( im )
  im_pixels = get_image_pixels( im )
  subim_size = get_image_size( subim )
  subim_pixels = get_image_pixels( subim )
  im_x_min = max( [ 1, blc[ 0 ] ] ) - 1
  im_x_max = min( [ im_size[ 0 ], blc[ 0 ] + subim_size[ 0 ] - 1 ] ) - 1
  im_y_min = max( [ 1, blc[ 1 ] ] ) - 1
  im_y_max = min( [ im_size[ 1 ], blc[ 1 ] + subim_size[ 1 ] - 1 ] ) - 1
  subim_x_min = max( [ 1, - blc[ 0 ] + 2 ] ) - 1
  subim_x_max = min( [ subim_size[ 0 ], im_size[ 0 ] - blc[ 0 ] + 1 ] ) - 1
  subim_y_min = max( [ 1, - blc[ 1 ] + 2 ] ) - 1
  subim_y_max = min( [ subim_size[ 1 ], im_size[ 1 ] - blc[ 1 ] + 1 ] ) - 1
  if ( ( im_x_min > im_x_max ) or ( im_y_min > im_y_max ) or 
      ( subim_x_min > subim_x_max ) or ( subim_y_min > subim_y_max ) ) :
    raise error( 'subimage and target image do not overlap' )
  im_overlap = im_pixels[ im_x_min : im_x_max + 1, im_y_min : im_y_max + 1 ]
  subim_overlap = subim_pixels[ subim_x_min : subim_x_max + 1, subim_y_min : subim_y_max + 1 ]
  mask = ( subim_overlap == get_aips_magic_value() )
  im_pixels[ im_x_min : im_x_max + 1, im_y_min : im_y_max + 1 ] = (
      im_overlap * array( mask, dtype = im_pixels.dtype ) + 
      subim_overlap * array( mask == False, dtype = im_pixels.dtype ) )
  set_image_pixels( im, im_pixels )
  return

###############################################################################

def fill_source( facet, pos, beam = None, invert = False, value = None ):

  # copy facet pixel contents
  facet_size = get_image_size( facet )
  pixels = get_image_pixels( facet )

  # get source position and shape in pixels
  [ xc, yc ] = pos
  [ bmaj, bmin, bpa ] = convert_beam_size( facet, beam = beam, to_pixel = True )

  [ X, Y ] = indices( tuple( facet_size ), dtype = int32 )
  mask = zeros( shape = pixels.shape, dtype = bool )
  NR2 = ( ( ( ( X - float( xc ) ) * cos( radians( bpa ) ) + 
              ( Y - float( yc ) ) * sin( radians( bpa ) ) ) / ( bmin / 2. ) )**2 + 
          ( ( ( Y - float( yc ) ) * cos( radians( bpa ) ) -
              ( X - float( xc ) ) * sin( radians( bpa ) ) ) / ( bmaj / 2. ) )**2 )
  mask = aput( mask, awhere( NR2 <= 1. ), True )

  # invert when requested
  if invert:
    mask = ( mask == False )

  # write blanks at masked values
  if ( value is None ):
    value = get_aips_magic_value()
  pixels = aput( pixels, awhere( mask ), value )
  set_image_pixels( facet, pixels )
  return

###############################################################################

def replace_pixels( im, value, new_value ):
  pixels = get_image_pixels( im )
  pixels = aput( pixels, awhere( pixels == value ), new_value )
  set_image_pixels( im, pixels )
  del pixels
  return

###############################################################################

def invert_image( im ):
  pixels = get_image_pixels( im )
  sel = awhere( pixels != get_aips_magic_value() )
  pixels = aput( pixels, sel, - aget( pixels, sel ) )
  set_image_pixels( im, pixels )
  del pixels
  return

###############################################################################

def get_facet( facets, i ):
  if ( i < 1000 ):
    facet_i = get_aips_file( facets.disk, facets.name, 
        facets.klass[ 0 : 3 ] + '%03d' % ( i ), facets.seq, 'MA' )
  else:
    facet_i = get_aips_file( facets.disk, facets.name,
        facets.klass[ 0 : 2 ] + '%04d' % ( i ), facets.seq, 'MA' )
  return facet_i

###############################################################################

def get_facet_number( facet ):
  try:
    number = int( facet.klass[ 2 : ] )
  except ValueError:
    number = int( facet.klass[ 3 : ] )
  else:
    pass
  return number

###############################################################################

def get_facet_beam( facet_i ):
  number = get_facet_number( facet_i )
  if ( number < 1000 ):
    beam_i = get_aips_file( facet_i.disk, facet_i.name,
        facet_i.klass[ 0 ] + 'BM' + facet_i.klass[ 3 : 6 ], facet_i.seq, 'MA' )
  else:
    beam_i = get_aips_file( facet_i.disk, facet_i.name,
        facet_i.klass[ 0 ] + 'B' + facet_i.klass[ 2 : 6 ], facet_i.seq, 'MA' )
  return beam_i

###############################################################################

def calculate_source_radec( facet, pos ):
  facet_radec = get_radec( facet )
  pixel_ref = get_pixel_reference( facet )
  pixel_size = get_pixel_size( facet, make_absolute = False )
  dx = ( pos[ 0 ] - pixel_ref[ 0 ] ) * pixel_size[ 0 ] / 3600.
  dy = ( pos[ 1 ] - pixel_ref[ 1 ] ) * pixel_size[ 1 ] / 3600.
  radius = sqrt( dx**2 + dy**2 )
  angle = degrees( atan2( dx, dy ) )
  projection = get_image_projection( facet )
  if ( projection == 'ARC' ):
    pass
  elif ( projection == 'SIN' ):
    radius = degrees( asin( radians( radius ) ) )
  elif ( projection == 'TAN' ):
    radius = degrees( atan( radians( radius ) ) )
  else:
    raise error( 'projection %s not supported (yet)' % projection )
  facet_rot = get_image_rotation( facet )
  if ( facet_rot != 0. ):
    angle = angle - facet_rot
  source_radec = calculate_offset_position( facet_radec, radius, angle )
  return source_radec

###############################################################################

def calculate_source_position( facet, radec, to_grid = False ):
  facet_radec = get_radec( facet )
  pixel_ref = get_pixel_reference( facet )
  pixel_size = get_pixel_size( facet, make_absolute = False )
  [ radius, angle ] = calculate_angular_separation( facet_radec, radec )
  projection = get_image_projection( facet )
  if ( projection == 'ARC' ):
    pass
  elif ( projection == 'SIN' ):
    radius = degrees( sin( radians( radius ) ) )
  elif ( projection == 'TAN' ):
    radius = degrees( tan( radians( radius ) ) )
  else:
    raise error( 'projection %s not supported (yet)' % projection )
  facet_rot = get_image_rotation( facet )
  if ( facet_rot != 0. ):
    angle = angle + facet_rot
  dx = 3600. * radius * sin( radians( angle ) ) / pixel_size[ 0 ]
  dy = 3600. * radius * cos( radians( angle ) ) / pixel_size[ 1 ]
  pos_x = pixel_ref[ 0 ] + dx
  pos_y = pixel_ref[ 1 ] + dy
  if to_grid:
    pos_x = float( around( pos_x ) )
    pos_y = float( around( pos_y ) )
  return [ pos_x, pos_y ]

###############################################################################

def get_image_pixels( im, flip = True ):
  wiz_im = wizardry( im )
  plane = wiz_im[ 0 ]
  pixels = plane.pixels.copy()
  if flip:
    pixels = transpose( pixels )
  del wiz_im
  return pixels

###############################################################################

def set_image_pixels( im, pixels, flip = True, units = None ):
  wiz_im = wizardry( im )
  plane = wiz_im[ 0 ]
  plane.pixels = pixels.copy()
  if flip:
    plane.pixels = transpose( plane.pixels )
  plane.update()
  sel = awhere( plane.pixels != get_aips_magic_value() )
  if ( len( sel ) > 0 ):
    wiz_im.header.datamin = aget( plane.pixels, sel ).min()
    wiz_im.header.datamax = aget( plane.pixels, sel ).max()
  else:
    wiz_im.header.datamin = 0.
    wiz_im.header.datamax = 0.
  if ( not units is None ):
    wiz_im.header.bunit = units
  wiz_im.header.update()
  wiz_im = wizardry( im )
  del wiz_im
  return

###############################################################################

def get_epoch( uvim ):
  return uvim.header.epoch

###############################################################################

def set_epoch( uvim, epoch ):
  wiz_uvim = wizardry( uvim )
  wiz_uvim.header.epoch = epoch
  wiz_uvim.header.update()
  wiz_uvim = wizardry( uvim )
  del wiz_uvim
  return

###############################################################################

def get_observer( uvim ):
  return uvim.header.observer

###############################################################################

def set_observer( uvim, observer ):
  wiz_uvim = wizardry( uvim )
  wiz_uvim.header.observer = observer
  wiz_uvim.header.update()
  wiz_uvim = wizardry( uvim )
  del wiz_uvim
  return

###############################################################################

#def calculate_source_shift( uv, radec ):
#  uv_radec = get_radec( uv )
#  [ radius, angle ] = calculate_angular_separation( uv_radec, radec )
#  dx = 3600. * radius * sin( radians( angle ) )
#  dy = 3600. * radius * cos( radians( angle ) )
#  return [ dx, dy ]

###############################################################################

def get_image_mean( im ):
  pixels = get_image_pixels( im )
  mask = ( pixels == get_aips_magic_value() )
  if alltrue( mask ):
    return None
  pixels = aput( pixels, awhere( mask ), 0. )
  count = ( mask == False ).sum()
  mean = float( pixels.sum() / float( count ) )
  return mean

###############################################################################

def get_image_rms( im ):
  pixels = get_image_pixels( im )
  mask = ( pixels == get_aips_magic_value() )
  if alltrue( mask ):
    return None
  pixels = aput( pixels, awhere( mask ), 0. )
  count = ( mask == False ).sum()
  rms = float( sqrt( ( pixels**2 ).sum() / float( count ) ) )
  return rms

###############################################################################

def clip_image( im, limits, values = None, epsilon = 1.e-7 ):
  pixels = get_image_pixels( im )
  sel_min = awhere( pixels < limits[ 0 ] )
  sel_max = awhere( pixels > limits[ 1 ] )
  if ( values is None ):
    pixels = aput( pixels, sel_min, ( 1. - epsilon ) * limits[ 0 ] )
    pixels = aput( pixels, sel_max, ( 1. - epsilon ) * limits[ 1 ] )
  else:
    pixels = aput( pixels, sel_min, values[ 0 ] )
    pixels = aput( pixels, sel_max, values[ 1 ] )
  set_image_pixels( im, pixels )
  del pixels
  return

###############################################################################

def get_image_stddev( im ):
  avg = get_image_mean( im )
  if ( avg is None ):
    return None
  pixels = get_image_pixels( im )
  mask = ( pixels == get_aips_magic_value() )
  if alltrue( mask ):
    return None
  pixels = aput( pixels, awhere( mask ), avg )
  count = ( mask == False ).sum()
  stddev = float( sqrt( ( ( pixels - avg )**2 ).sum() / float( count ) ) )
  return stddev

###############################################################################

def get_observing_epoch( uv ):
  [ year, month, day ] = [ int( a ) for a in uv.header.date_obs.split( '-' ) ]
  days = ( date( year, month, day ) - date( year, 1, 1 ) ).days
  epoch = year + ( float( days ) / 365.25 )
  return epoch

###############################################################################

def transfer_model_components( in_facet, out_facet, model_version = 0 ):
# Note: this function will put clean components at fractional pixel positions.
# This requires the use of cmethod = 'DFT' for AIPS tasks like CALIB or UVSUB
  
  in_pix_ref = get_pixel_reference( in_facet )
  in_pix_scale = get_pixel_size( in_facet, make_absolute = False )
  out_pix_ref = get_pixel_reference( out_facet )
  out_pix_scale = get_pixel_size( out_facet, make_absolute = False )
  out_facet_size = get_image_size( out_facet )
  wiz_in_facet = wizardry( in_facet )
  wiz_out_facet = wizardry( out_facet )

  if ( not table_exists( in_facet, 'CC', model_version ) ):
    raise error( 'input model component table does not exist' )
  if ( model_version == 0 ):
    in_version = in_facet.table_highver( 'CC' )
  else:
    in_version = model_version
  out_version = out_facet.table_highver( 'CC' ) + 1
  
  # transfer CC components to new output CC table
  in_cc_table = wiz_in_facet.table( 'CC', in_version )
  try:
    no_parms = in_cc_table.keywords[ 'NO_PARMS' ]
  except:
    no_parms = 0
  out_cc_table = new_table( out_facet, 'CC', out_version, no_parms = no_parms )
  for cc in in_cc_table:
    x = in_pix_ref[ 0 ] + 3600. * cc.deltax / in_pix_scale[ 0 ]
    y = in_pix_ref[ 1 ] + 3600. * cc.deltay / in_pix_scale[ 1 ]
    radec = calculate_source_radec( in_facet, [ x, y ] )
    [ x, y ] = calculate_source_position( out_facet, radec )
    if ( ( x >= 0.5 ) and ( x <= out_facet_size[ 0 ] + 0.5 ) and 
         ( y >= 0.5 ) and ( y <= out_facet_size[ 1 ] + 0.5 ) ):
      out_cc = new_table_row( out_cc_table )
      out_cc.deltax = ( x - out_pix_ref[ 0 ] ) * out_pix_scale[ 0 ] / 3600.
      out_cc.deltay = ( y - out_pix_ref[ 1 ] ) * out_pix_scale[ 1 ] / 3600.
      out_cc.flux = cc.flux
      if ( no_parms > 0 ):
        out_cc.parms = cc.parms
      out_cc_table.append( out_cc )
  in_cc_table.close()
  out_cc_table.close()
  del wiz_in_facet
  del wiz_out_facet
  
  return out_version

###############################################################################

def new_table( uvim, table, version = 0, **keywords ):
  if ( version == 0 ):
    new_version = uvim.table_highver( table ) + 1
  else:
    new_version = version
  wiz_uvim = wizardry( uvim )
  ntable = wiz_uvim.add_table( table, new_version, **keywords )
  del wiz_uvim
  return ntable

###############################################################################

def new_table_row( table, **keywords ):
  row = Wizardry.AIPSData.AIPSTableRow( table, **keywords )
  return row

###############################################################################

def copy_keywords( in_uvim, out_uvim ):
  wiz_out_uvim = wizardry( out_uvim )
  for keyword in in_uvim.keywords:
    temp = in_uvim.keywords[ keyword ]
    if ( type( temp ) == type( [] ) ):
      temp = temp[ 0 ]
    wiz_out_uvim.keywords[ keyword ] = temp
    wiz_out_uvim.keywords.update()
  del wiz_out_uvim
  return

###############################################################################

def aips_file_name_to_string( uvim ):
  if is_uv( uvim ):
    string = '%s.%s.%s.%s.UV' % ( repr( uvim.disk ), uvim.name, uvim.klass, repr( uvim.seq ) )
  elif is_image( uvim ):
    string = '%s.%s.%s.%s.MA' % ( repr( uvim.disk ), uvim.name, uvim.klass, repr( uvim.seq ) )
  else:
    raise error( 'unknown AIPS file type' )
  return string

###############################################################################

def get_time_list_from_fit_table( uv, fit_version = 0 ):
  wiz_uv = wizardry( uv )
  ni_table = wiz_uv.table( 'NI', fit_version )
  time_list = []
  for row in ni_table:
    time_list.append( row.time )
  time_count = len( time_list )
  del wiz_uv
  return time_list

###############################################################################

def get_time_list_from_solution_table( uv, solution_version = 0 ):
  wiz_uv = wizardry( uv )
  sn_table = wiz_uv.table( 'SN', solution_version )
  time_list = []
  old_time = - 1000000.
  for row in sn_table:
    if ( row.time > old_time ):
      time_list.append( row.time )
      old_time = row.time
  time_count = len( time_list )
  del wiz_uv
  return time_list

###############################################################################

def find_aips_cno( uvim ):
  cno = -1
  aips_disk = uvim.disk
  aips_name = uvim.name
  aips_class = uvim.klass
  aips_seq = uvim.seq
  if is_uv( uvim ):
    aips_type = 'UV'
  elif is_image( uvim ):
    aips_type = 'MA'
  else:
    raise error( 'unknown AIPS file type' )
  cat_disk = AIPSCat( aips_disk )[ aips_disk ]
  for cat in cat_disk:
    if ( [ cat.name, cat.klass, cat.seq, cat.type ] == [ aips_name, aips_class,
        aips_seq, aips_type ] ):
      cno = cat.cno
      break
  return cno

###############################################################################

def get_aips_cno( disk, cno ):
  uvim = None
  cat_disk = AIPSCat( disk )[ disk ]
  for cat in cat_disk:
    if ( cat.cno == cno ):
      uvim = get_aips_file( disk, cat.name, cat.klass, cat.seq, cat.type )
      break
  return uvim

###############################################################################

def read_fits_uv( file_name, uv, compress = False ):
  if ( not file_exists( file_name ) ):
    raise error( 'UV FITS file %s not found' % ( file_name ) )
  if compress:
    douvcomp = 1
  else:
    douvcomp = 0
  call_aips_task( 'FITLD', datain = file_name, optype = 'UV', outdata = uv,
      douvcomp = douvcomp )
  return

###############################################################################

def write_fits_uv( uv, file_name, overwrite = True, fits_table = False, 
    compress = False, keep_history = False ):
  if ( ( not overwrite ) and ( file_exists( file_name ) ) ):
    raise error( 'UV FITS file %s already exists' % ( file_name ) )
  elif ( overwrite and ( file_exists( file_name ) ) ):
    remove_file( file_name )
  if fits_table:
    if compress:
      douvcomp = 1
    else:
      douvcomp = 0
    call_aips_task( 'FITAB', indata = uv, intype = 'UV', dataout = file_name,
        douvcomp = douvcomp, keep_history = keep_history )
  else:
    call_aips_task( 'FITTP', indata = uv, intype = 'UV', dataout = file_name, 
        keep_history = keep_history )
  return

###############################################################################

def read_fits_image( file_name, im ):
  if ( not file_exists( file_name ) ):
    raise error( 'image FITS file %s not found' % ( file_name ) )
  call_aips_task( 'FITLD', datain = file_name, optype = 'IM', outdata = im )
  return

###############################################################################

def write_fits_image( im, file_name, overwrite = True, compress = False,
    include_tables = False, crop = False, output_pixel_size = None,
    output_image_size = None, fits_table = False, keep_history = False ):
  
  if ( ( not overwrite ) and ( file_exists( file_name ) ) ):
    raise error( 'image FITS file %s already exists' % ( file_name ) )
  elif ( overwrite and ( file_exists( file_name ) ) ):
    remove_file( file_name )
  
  # make work copy
  temp_im = get_aips_file( im.disk, im.name, 'TEMP', - 1, 'MA' )
  call_aips_task( 'MOVE', indata = im, outdata = temp_im, userid = get_aips_userid() )
  
  # scale image
  if ( ( not output_pixel_size is None ) and ( not output_image_size is None ) ):
    scale_im = scale_image( temp_im, output_pixel_size, output_image_size )
    temp_im.zap()
  else:
    scale_im = temp_im
  
  # crop image
  if ( crop or ( ( not output_pixel_size is None ) and ( output_image_size is None ) ) ):
    crop_im = crop_image( scale_im )
    scale_im.zap()
  else:
    crop_im = scale_im
  
  # for image compression, replace blanks with zeros
  if compress:
    format = 1
    replace_pixels( crop_im, get_aips_magic_value(), 0. )
  else:
    format = 0
  
  # remove tables
  if ( not include_tables ):
    for table in crop_im.tables:
      table_version = table[ 0 ]
      table_type = table[ 1 ][ -2 : ]
      if ( table_type != 'HI' ):
        crop_im.zap_table( table_type, table_version )
  
  # write image to FITS file
#  call_aips_task( 'FITTP', indata = crop_im, intype = 'MA', dataout = file_name,
#      format = format )
  if fits_table:
    call_aips_task( 'FITAB', indata = crop_im, intype = 'MA', dataout = file_name,
        douvcomp = 0, keep_history = keep_history )
  else:
    call_aips_task( 'FITTP', indata = crop_im, intype = 'MA', dataout = file_name,
        format = format, keep_history = keep_history )
  crop_im.zap()
  return

###############################################################################

def crop_image( im, blank_value = get_aips_magic_value() ):

  pixels = get_image_pixels( im )
  x_size = pixels.shape[ 0 ] 
  y_size = pixels.shape[ 1 ] 
  x_min = 0
  x_max = x_size - 1
  y_min = 0
  y_max = y_size - 1

  # search for lower x boundary
  for x in range( x_size ):
    if sometrue( pixels[ x, : ] != get_aips_magic_value() ):
      x_min = x
      break

  # try a quick scan on upper x boundary, otherwise a thorough search
  if ( ( x_min > 0 ) and ( x_size - x_min > x_min ) ):
    if ( alltrue( pixels[ x_size - x_min, : ] == blank_value ) and
        sometrue( pixels[ x_size - x_min - 1, : ] != blank_value ) ):
      x_max = x_size - x_min - 1
  if ( x_max == x_size - 1 ):
    for x in range( x_size - 1, x_min - 1, - 1 ):
      if sometrue( pixels[ x, : ] != blank_value ):
        x_max = x
        break

  # try a quick scan on lower y boundary, otherwise a thorough search
  if ( ( x_min == 0 ) and sometrue( pixels[ : , 0 ] != blank_value ) ):
    y_min = 0
  else:
    if ( x_min <= y_max ):
      if ( alltrue( pixels[ : , x_min - 1 ] == blank_value ) and
          sometrue( pixels[ : , x_min ] != blank_value ) ):
        y_min = x_min      
    if ( y_min == 0 ):
      for y in range( y_size ):
        if sometrue( pixels[ : ,y ] != blank_value ):
          y_min = y
          break

  # try a quick scan on upper y boundary, otherwise a thorough search
  if ( ( y_min > 0 ) and ( y_size - y_min > y_min ) ):
    if ( alltrue( pixels[ : , y_size - y_min ] == blank_value ) and
        sometrue( pixels[ : , y_size - y_min - 1 ] != blank_value ) ):
      y_max = y_size - y_min - 1
  if ( y_max == y_size - 1 ):
    for y in range( y_size - 1, y_min - 1, - 1 ):
      if sometrue( pixels[ : , y ] != blank_value ):
        y_max = y
        break

  crop_im = get_aips_file( im.disk, im.name, 'CROP', - 1, 'MA' )
  if ( ( x_min != 0 ) or ( x_max != x_size - 1 ) or ( y_min != 0 ) or ( y_max != y_size - 1 ) ):
    call_aips_task( 'SUBIM', indata = im, outdata = crop_im, blc = [ x_min + 1, y_min + 1 ],
        trc = [ x_max + 1, y_max + 1 ] ) 
  else:
    call_aips_task( 'MOVE', indata = im, outdata = crop_im, userid = get_aips_userid() )

  return crop_im

###############################################################################

def scale_image( im, new_pixel_size = None, new_image_size = None, rotate = False ):
  scale_im = get_aips_file( im.disk, im.name, 'SCALE', - 1, 'MA' )
  if ( ( new_pixel_size is None ) and ( new_image_size is None ) ):
    call_aips_task( 'MOVE', indata = im, outdata = scale_im, userid = get_aips_userid() )
    return scale_im
  image_size = get_image_size( im )
  pixel_size = get_pixel_size( im )
  if rotate:
    rotation = - get_image_rotation( im )
  else:
    rotation = 0.
  if ( not new_pixel_size is None ):
    scale = pixel_size[ 0 ] / new_pixel_size[ 0 ]
    dscale = ( pixel_size[ 1 ] / new_pixel_size[ 1 ] ) / scale
  else:
    scale = 0.
    dscale = 0.
  if ( not new_image_size is None ):
    imsize = new_image_size
    if ( new_pixel_size is None ):
      scale = float( image_size[ 0 ] ) / float( new_image_size[ 0 ] )
      dscale = ( float( image_size[ 1 ] ) / float( new_image_size[ 1 ] ) ) / scale
  else:
    if ( not new_pixel_size is None ):
      imsize = [ int( ceil( float( image_size[ 0 ] ) * scale ) ) + 5,
          int( ceil( float( image_size[ 1 ] ) * scale * dscale ) ) + 5 ]
  replace_pixels( im, get_aips_magic_value(), 0. )
  call_aips_task( 'OGEOM', indata = im, outdata = scale_im, imsize = imsize,
      aparm = [ 0, 0, rotation, scale, dscale, 0, 0 ], reweight = [ 3, 0.5 ] )
  replace_pixels( im, 0., get_aips_magic_value() )
  return scale_im

###############################################################################

def scale_model_flux( facet, scale, model_version = 0, threshold = 0. ):
  # threshold (applied before scale):
  #   if > 0: drop components with abs( flux ) < threshold
  #   if < 0: drop components with flux < 0.
  
  if ( not table_exists( facet, 'CC', model_version ) ):
    raise error( 'input model component table does not exist' )
  if ( model_version == 0 ):
    in_version = facet.table_highver( 'CC' )
  else:
    in_version = model_version
  out_version = facet.table_highver( 'CC' ) + 1
  wiz_facet = wizardry( facet )
  in_cc_table = wiz_facet.table( 'CC', in_version )
  try:
    no_parms = in_cc_table.keywords[ 'NO_PARMS' ]
  except:
    no_parms = 0
  out_cc_table = new_table( facet, 'CC', out_version, no_parms = no_parms )
  for cc in in_cc_table:
    if ( ( ( threshold < 0. ) and ( cc.flux < 0. ) ) or
        ( abs( cc.flux ) < threshold ) ):
      continue
    out_cc = new_table_row( out_cc_table )
    out_cc.flux = scale * cc.flux
    out_cc.deltax = cc.deltax
    out_cc.deltay = cc.deltay
    if ( no_parms > 0 ):
      out_cc.parms = cc.parms
    out_cc_table.append( out_cc )
  in_cc_table.close()
  out_cc_table.close()
  del wiz_facet
  
  return out_version

###############################################################################

def get_image_projection( im ):
  for ctype in im.header.ctype:
    if ( ctype.find( 'RA--' ) != - 1 ):
      ra_index = im.header.ctype.index( ctype )
    elif ( ctype.find( 'XA--' ) != - 1 ):
      ra_index = im.header.ctype.index( ctype )
    elif ( ctype.find( 'GLON' ) != - 1 ):
      ra_index = im.header.ctype.index( ctype )
    if ( ctype.find( 'DEC-' ) != - 1 ):
      dec_index = im.header.ctype.index( ctype )
    elif ( ctype.find( 'XEC-' ) != - 1 ):
      dec_index = im.header.ctype.index( ctype )
    elif ( ctype.find( 'GLAT' ) != - 1 ):
      dec_index = im.header.ctype.index( ctype )
  ra_proj = im.header.ctype[ ra_index ][ - 3 : ]
  dec_proj = im.header.ctype[ dec_index ][ - 3 : ]
  if ( ra_proj != dec_proj ):
    raise error( 'RA and DEC have different projections' )
  return ra_proj

###############################################################################

def get_clean_boxes( facet_file_name, facet_list = [], drop_dummy_boxes = True ):
  clean_box_list = []
  if ( not file_exists( facet_file_name ) ):
    raise error( 'facet file %s does not exist' % facet_file_name )
  facet_file = file( facet_file_name, mode = 'r' )
  for line in facet_file:
    words = [ word.strip() for word in line.split() ]
    if ( len( words ) != 5 ):
      continue
    if ( words[ 0 ][ 0 ] == '#' ):
      continue
    try:
      [ i,a,b,c,d ] = [ int( word ) for word in words ]
    except:
      continue
    if ( ( len( facet_list ) > 0 ) and ( not i in facet_list ) ):
      continue
    if ( drop_dummy_boxes and ( 
        ( [ a,b,c,d, ] == [ 0,0,0,0 ] ) or ( [ a,b,c,d, ] == [ -1,1,5,5 ] ) ) ) :
      continue
    clean_box_list.append( [ i,a,b,c,d ] )
  facet_file.close()
  return clean_box_list

###############################################################################

def remove_facets( facets, include_beams = True, max_facet_count = 9999 ):
  for i in range( 1, 1 + max_facet_count ):
    facet = get_facet( facets, i )
    if ( facet.exists() ):
      facet.clrstat()
      facet.zap()
    if include_beams:
      beam = get_facet_beam( facet )
      if ( beam.exists() ):
        beam.clrstat()
        beam.zap()
  return

###############################################################################

def clear_aips_disk( aips_disk ):
  aips_cat = AIPSCat( aips_disk )[ aips_disk ]
  for aips_file in aips_cat:
    uvim = get_aips_file( aips_disk, aips_file.name, aips_file.klass,
        aips_file.seq, aips_file.type )
    try:
      uvim.clrstat()
    except:
      pass
    try:
      uvim.zap()
    except:
      pass
  return

###############################################################################

def copy_facets( facets, new_facets, include_beams = True, max_facet_count = 9999 ):
  for i in range( 1, 1 + max_facet_count ):
    facet = get_facet( facets, i )
    new_facet = get_facet( new_facets, i )
    if ( facet.exists() ):
      call_aips_task( 'MOVE', indata = facet, outdata = new_facet,
          userid = get_aips_userid() )
    if include_beams:
      beam = get_facet_beam( facet )
      new_beam = get_facet_beam( new_facet )
      if ( beam.exists() ):
        call_aips_task( 'MOVE', indata = beam, outdata = new_beam,
            userid = get_aips_userid() )
  return

###############################################################################

def calculate_pbparm_attenuation( freq, radius, pbparms, check_derivative = True ):
  [ cutoff, apply_pbparms, pbparm3, pbparm4, pbparm5, pbparm6, pbparm7 ] = pbparms
  if ( apply_pbparms > 0. ):
    X = ( ( freq / 1.e9 ) * ( radius * 60. ) )**2
    A = ( 1.0 + ( X * pbparm3 / 1.e3 ) + ( X**2 * pbparm4 / 1.e7 ) + 
        ( X**3 * pbparm5 / 1.e10 ) + ( X**4 * pbparm6 / 1.e13 ) + 
        ( X**5 * pbparm7 / 1.e16 ) )
    if ( A < cutoff ):
      A = 0.
    elif check_derivative:
      dXdR = 2 * ( ( freq / 1.e9 ) * ( radius * 60. ) ) * ( ( freq / 1.e9 ) * ( 60. ) )
      dAdX = ( ( pbparm3 / 1.e3 ) + ( 2 * X * pbparm4 / 1.e7 ) +
          ( 3 * X**2 * pbparm5 / 1.e10 ) + ( 4 * X**3 * pbparm6 / 1.e13 ) + 
          ( 5 * X**4 * pbparm7 / 1.e16 ) )
      if ( dAdX * dXdR > 0. ):
        A = 0.
  else:
    A = 1.
  return A

###############################################################################

def fix_antenna_table( uv, version = 0 ):
# This is to fix a bug in AIPS++ ms2uvfix
  
  if ( version == 0 ):
    an_version = uv.table_highver( 'AN' )
  else:
    an_version = version
  
  # make working copy of antenna table
  call_aips_task( 'TACOP', indata = uv, inext = 'AN', invers = an_version, ncount = 1,
      outdata = uv, outvers = 0 )
  wiz_uv = wizardry( uv )
  an_table = wiz_uv.table( 'AN', 0 )
  
  # get array reference position
  array_xyz = array( [ an_table.keywords[ 'ARRAYX' ], an_table.keywords[ 'ARRAYY' ],
      an_table.keywords[ 'ARRAYZ' ] ], dtype = float64 )
  if alltrue( array_xyz == azeros( array_xyz ) ):
    # nothing to be done, so exit
    uv.zap_table( 'AN', 0 )
    return
  array_geo_llh = array( xyz_to_geo_llh( array_xyz.tolist() ), dtype = float64 )
  [ lon, lat ] = aradians( array_geo_llh[ 0 : 2 ] )
  rotation = array( [ [ cos( lon ), - sin( lon ), 0. ],
                      [ sin( lon ),   cos( lon ), 0. ],
                      [         0.,           0., 1. ] ], dtype = float64 )
  
  # rotate delta xyz of individual antennas
  for row in an_table:
    dxyz = dot( array( row.stabxyz, dtype = float64 ), rotation )
    row.stabxyz = dxyz.tolist()
    row.update()
  an_table.close()
  
  # replace original antenna table with working copy
  uv.zap_table( 'AN', an_version )
  call_aips_task( 'TACOP', indata = uv, inext = 'AN', invers = 0, ncount = 1,
      outdata = uv, outvers = an_version )
  uv.zap_table( 'AN', 0 )
  del wiz_uv
  
  return

###############################################################################

def time_to_string( time ):
  dhms = time_to_dhms( time )
  string = '%03dd%02dh%02dm%05.2fs' % ( dhms[ 0 ], dhms[ 1 ], dhms[ 2 ], dhms[ 3 ] )
  return string

###############################################################################

def radec_to_string( radec, decimals = [ 2, 1 ], do_iau = False,
    separators = [ 'h', 'm', 's ', 'd', "'", '"' ] ):
  hmsdms = degdeg_to_hmsdms( radec )
  [ a, b ] = decimals
  if do_iau:
    x = 4
    a = a + x
    b = b + x
  format = '%02d' + separators[ 0 ] + '%02d' + separators[ 1 ] + '%0'
  if ( a > 0 ):
    format = format + repr( 3 + a ) + '.' +  repr( a ) + 'f' + separators[ 2 ]
  else:
    format = format + '2d' + separators[ 2 ]
  string = format % ( hmsdms[ 0 ], hmsdms[ 1 ], hmsdms[ 2 ] )
  if ( asign( hmsdms[ 3 ] ) < 0. ):
    string = string + '-'
    hmsdms[ 3 ] = -hmsdms[ 3 ]
  else:
    string = string + '+'
  format = '%02d' + separators[ 3 ] + '%02d' + separators[ 4 ] + '%0'
  if ( b > 0 ):
    format = format + repr( 3 + b ) + '.' + repr( b ) + 'f' + separators[ 5 ]
  else:
    format = format + '2d' + separators[ 5 ]
  string = string + format % ( hmsdms[ 3 ], hmsdms[ 4 ], hmsdms[ 5 ] )
  if do_iau:
    try:
      index1 = string.index( '+' ) - len( separators[ 2 ] ) - x
    except ValueError:    
      index1 = string.index( '-' ) - len( separators[ 2 ] ) - x
    index2 = len( string ) - len( separators[ 5 ] ) - 4
    string = string[ : index1 ] + string[ index1 + x : index2 ] + string[ index2 + x : ]
    if ( decimals[ 0 ] == 0 ):
      string = string[ : index1 - 1 ] + string[ index1 : ]
    if ( decimals[ 1 ] == 0 ):
      index2 = len( string ) - len( separators[ 5 ] )
      string = string[ : index2 - 1 ] + string[ index2 : ]
  return string

###############################################################################

def get_evla_antennas( uv ):
  ant_names = uv.antennas
  evla_list = []
  for i in range( len( ant_names ) ):
    if ( 'EVLA' in ant_names[ i ] ):
      evla_list.append( i + 1 )
  return evla_list

###############################################################################

def get_vla_antennas( uv ):
  ant_names = uv.antennas
  vla_list = []
  evla_list = get_evla_antennas( uv )
  for i in range( len( ant_names ) ):
    if ( ( not i + 1 in evla_list ) and ( 'VLA' in ant_names[ i ] ) ):
      vla_list.append( i + 1 )
  return vla_list

###############################################################################

def is_uv( uvim ):
  return isinstance( uvim, AIPSUVData )

###############################################################################

def is_image( uvim ):
  return isinstance( uvim, AIPSImage )

###############################################################################

def set_header_keyword( uvim, key, value, index = None ):
  wiz_uvim = wizardry( uvim )
  if ( not index is None ):
    wiz_uvim.header[ key ][ index ] = value
  else:
    wiz_uvim.header[ key ] = value
  wiz_uvim.header.update()
  del wiz_uvim
  return

###############################################################################

def ehex_to_decimal( ehex ):
  conversion = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  decimal = 0
  multiplier = 1
  for char in ehex.upper()[ : : -1 ]:
    decimal = decimal + conversion.index( char ) * multiplier
    multiplier = multiplier * 36
  return decimal

###############################################################################

def decimal_to_ehex( decimal ):
  conversion = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  ehex = ''
  dec = decimal
  while ( dec > 0 ):
    index = dec % 36
    ehex = conversion[ index ] + ehex
    dec = ( dec - index ) / 36
  return ehex

###############################################################################

def zero_pad_ehex( ehex, length ):
  new_ehex = ehex
  for i in range( length - len( ehex) ):
    new_ehex = '0' + new_ehex
  return new_ehex

###############################################################################

def get_aips_file_name( uvim, table = None, version = None ):
  file_path = os.getenv( 'DA' + zero_pad_ehex( decimal_to_ehex( uvim.disk ), 2 ) )
  if ( ( not table is None ) and ( not version is None ) ):
    tab = table[ -2 : ]
    if ( version > 0 ):
      vers = version
    elif ( version == 0 ):
      vers = uvim.table_highver( tab )
    else: # version < 0
      vers = uvim.table_highver( tab ) + 1
  else:
    vers = 1
    if is_uv( uvim ):
      tab = 'UV'
    else: # is_image( uvim )
      tab = 'MA'
  cno = zero_pad_ehex( decimal_to_ehex( find_aips_cno( uvim ) ), 3 )
  vers = zero_pad_ehex( decimal_to_ehex( vers ), 3 )
  user = zero_pad_ehex( decimal_to_ehex( get_aips_userid() ), 3 )
  file_name = file_path + '/' + tab + 'D' + cno + vers + '.' + user + ';'
  return file_name

###############################################################################

def check_aips_disk( disk, user_id = None, delete = False, print_info = True ):
  error_count = 0
  if ( user_id is None ):
    userid = get_aips_userid()
  else:
    userid = user_id
  if print_info:
    print 'checking AIPS disk %s for user %s' % ( repr( disk ), repr( userid ) )
  file_path = os.getenv( 'DA' + zero_pad_ehex( decimal_to_ehex( disk ), 2 ) )
  file_list = os.listdir( file_path )
  file_list_count = len( file_list )
  print_count = int( floor( 0.1 * float( file_list_count ) ) )
  i = 0
  for name in file_list:
    i = i + 1
    if print_info:
      if ( i % print_count == 0 ):
        fraction = 0.1 * float( i / print_count )
        print '... %d percent done' % ( int( 100. * fraction ) )
    if ( name == 'SPACE' ):
      continue
    if ( ( name[ 2 ] == 'D' ) and ( name[ 9 ] == '.' ) and ( name[ 13 ] == ';' ) ):
      user = ehex_to_decimal( name[ 10 : 13 ] )
      if ( user != userid ):
        continue
      if ( name[ 0 : 9 ] == 'CAD000000' ):
        continue
      tab = name[ 0 : 2 ]
      cno = ehex_to_decimal( name[ 3 : 6 ] )
      vers = ehex_to_decimal( name[ 6 : 9 ] )
#      if ( tab in [ 'MA', 'UV' ] ):# , 'CB' ] ):
      if ( tab in [ 'MA', 'UV', 'CB' ] ):
        if ( vers != 1 ):
          error_count = error_count + 1
          if print_info:
            print '... found file %s with version %d > 1' % ( name, vers )
            print '...... disk = %d, cno = %d, ext = %s' % ( disk, cno, tab )
          if delete:
            remove_file( file_path + '/' + name )
            if print_info:
              print '...... deleted file'
          continue
        uvim = get_aips_cno( disk, cno )
        if ( uvim is None ):
          error_count = error_count + 1
          if print_info:
            print '... found uncataloged file %s' % ( name )
            print '...... disk = %d, cno = %d, ext = %s' % ( disk, cno, tab )
          if delete:
            remove_file( file_path + '/' + name )
            if print_info:
              print '...... deleted file'
          continue
        if ( ( ( tab == 'MA' ) and ( not is_image( uvim ) ) ) or
            ( ( tab == 'UV' ) and ( not is_uv( uvim ) ) ) ):
          error_count = error_count + 1
          if print_info:
            print '... found file %s with type-match error' % ( name )
            print '...... disk = %d, cno = %d, ext = %s' % ( disk, cno, tab )
          if delete:
            remove_file( file_path + '/' + name )
            if print_info:
              print '...... deleted file'
          continue
        continue
      
      elif ( tab in [ 'AN', 'AT', 'BL', 'BP', 'CC', 'CD', 'CL', 'CQ', 'CT',
          'FG', 'FQ', 'GC', 'HI', 'IM', 'MC', 'MF', 'NI', 'NX', 'OB', 'OF',
          'OT', 'PC', 'PS', 'SN', 'SU', 'SY', 'TY', 'VL', 'VZ', 'WX',
          'CG', 'SC' ] ): # 'CB' ] ):
#          'CG', 'SC', 'CB' ] ):
        try:
          uvim = get_aips_cno( disk, cno )
        except:
          uvim = None
        if ( uvim is None ):
          error_count = error_count + 1
          if print_info:
            print '... found orphane table %s' % ( name )
            print '...... disk = %d, cno = %d, ext = %s' % ( disk, cno, tab )
          if delete:
            remove_file( file_path + '/' + name )
            if print_info:
              print '...... deleted file'
          continue
        table_ok = table_exists( uvim, tab, vers )
        if ( not table_ok ):
          error_count = error_count + 1
          if print_info:
            print '... found uncataloged table %s' % ( name )
            print '...... disk = %d, cno = %d, ext = %s' % ( disk, cno, tab )
          if delete:
            remove_file( file_path + '/' + name )
            if print_info:
              print '...... deleted file'
          continue
        continue
      
      else:
        error_count = error_count + 1
        if print_info:
          print '... found unknown file type %s' % ( name )
          print '...... disk = %d, cno = %d, ext = %s' % ( disk, cno, tab )
        if delete:
#          remove_file( file_path + '/' + name )
          if print_info:
            print '...... file NOT deleted'
  
  return error_count

###############################################################################

def time_average_uv_data( uv, factor, scan_break_min = 5. ):
  
  # determine new integration time
  ifactor = int( round( float( factor ) ) )
  if ( float( factor ) != float( ifactor ) ):
    print 'WARNING: rounding factor to the nearest integer'
  if ( ifactor < 1 ):
    raise error( 'cannot time-average data with factor < 1' )
  dtime = find_integration_time( uv ) / 86400.
  new_dtime = float( ifactor ) * dtime
  
  # detect scans
  time_list = get_time_list( uv )
  time_array = array( time_list, dtype = float64 )
  dtime_array = time_array[ 1 : ] - time_array[ : -1 ]
  sel = awhere( dtime_array > scan_break_min / 1440. )
  scan_list = [ 0 ]
  for [ s ] in sel:
    scan_list.append( s )
    scan_list.append( s + 1 )
  scan_list.append( len( time_array ) - 1 )
  scan_list = [ [ time_list[ scan_list[ i ] ], time_list[ scan_list[ i + 1 ] ] ]
      for i in range( 0, len( scan_list ), 2 ) ]
  
  # time average data using SPLAT
  temp_uv = get_aips_file( uv.disk, uv.name, 'TEMP', -1, 'UV' )
  call_aips_task( 'SPLAT', indata = uv, outdata = temp_uv, douvcomp = 0,
      solint = new_dtime * 1440. )
  time_list = get_time_list( temp_uv )
  time_array = array( time_list, dtype = float64 )
  
  # calculate new time grid
  new_time_array = time_array.copy()
  for scan in scan_list:
    interval = scan[ 1 ] - scan[ 0 ] + dtime
    count = int( floor( interval / new_dtime ) )
    if ( interval - float( count ) * new_dtime > 0.1 * dtime ):
      count = count + 1
    new_interval = float( count ) * new_dtime 
    start_time = 0.5 * ( scan[ 0 ] + scan[ 1 ] - new_interval )
    for i in range( count ):
      time_min = start_time + float( i ) * new_dtime + 0.01 * dtime
      time_max = start_time + float( i + 1 ) * new_dtime + 0.01 * dtime
      sel = awhere( ( time_array >= time_min ) & ( time_array < time_max ) )
      if ( len( sel ) > 0 ):
        time_mid = start_time + ( float( i ) + 0.5 ) * new_dtime
        new_time_array = aput( new_time_array, sel, time_mid )
  new_time_list = [ float32( x ) for x in new_time_array.tolist() ]
  
  # write new times
  wiz_temp_uv = wizardry( temp_uv )
  for group in wiz_temp_uv:
    index = time_list.index( group.time )
    group.time = new_time_list[ index ]
    group.update()
  # bug fix
  for group in wiz_temp_uv:
    break
  
  # sort data in TB
  tavg_uv = get_aips_file( uv.disk, uv.name, 'TAVG', -1, 'UV' )
  call_aips_task( 'UVSRT', indata = temp_uv, outdata = tavg_uv, sort = 'TB' )
  temp_uv.zap()
  del wiz_temp_uv
  
  return tavg_uv

###############################################################################

def get_reference_date( uvim ):
  if is_uv( uvim ):
    if ( uvim.table_highver( 'AN' ) > 1 ):
      raise error( 'currently only single AN tables supported' )
    an = uvim.table( 'AN', 1 )
    datum = an.keywords[ 'RDATE' ]
    if ( uvim.header.date_obs != datum ):
      raise error( 'reference date in AN table differs from UV header' )
  else:
    datum = uvim.header.date_obs
  return datum

###############################################################################

def calculate_offset_date( datum, offset ):
  ymd = datum.split( '-' )
  day = date( int( ymd[ 0 ] ), int( ymd[ 1 ] ), int( ymd[ 2 ] ) ).toordinal()
  new_day = date.fromordinal( day + offset )
  new_datum = '%4d-%2d-%2d' % ( new_day.year, new_day.month, new_day.day )
  return new_datum

###############################################################################

def change_uv_reference_date( uv, ref_uv, offset = 2 ):
  
  # get original time list and gst list
  time_list = get_time_list( uv )
  gst_list = get_gst_list( uv, time_list = time_list )
  
  # copy several AN header values
  an = ref_uv.table( 'AN', 0 )
  wiz_uv = wizardry( uv )
  wiz_uv.header[ 'date_obs' ] = ref_uv.header[ 'date_obs' ]
  wiz_uv.header.update()
  wiz_an = wiz_uv.table( 'AN', 0 )
  wiz_an.keywords[ 'RDATE' ] = an.keywords[ 'RDATE' ]
  wiz_an.keywords[ 'GSTIA0' ] = an.keywords[ 'GSTIA0' ]
  wiz_an.keywords[ 'DEGPDY' ] = an.keywords[ 'DEGPDY' ]
  wiz_an.keywords[ 'TIMSYS' ] = an.keywords[ 'TIMSYS' ]
  try:
    wiz_an.keywords[ 'DATUTC' ] = an.keywords[ 'DATUTC' ]
  except:
    pass
  try:
    wiz_an.keywords[ 'IATUTC' ] = an.keywords[ 'IATUTC' ]
  except:
    pass
  wiz_an.close()
  
  # calculate new time stamps
  temp_gst_list = get_gst_list( uv, time_list = time_list )
  time_array = array( time_list, dtype = float64 )
  gst_array = array( gst_list, dtype = float64 )
  temp_gst_array = array( temp_gst_list, dtype = float64 )
  deg_per_day = an.keywords[ 'DEGPDY' ]
  new_time_array = time_array * deg_per_day + gst_array - temp_gst_array
  new_time_array = amodulo( new_time_array, 360. ) + float( offset ) * 360.
  sel = awhere( new_time_array[ 1 : ] < new_time_array[ : -1 ] )
  for s in sel:
    new_time_array[ s[ 0 ] + 1 : ] = new_time_array[ s[ 0 ] + 1 : ] + 360.
  new_time_list = array( new_time_array / deg_per_day, dtype = float32 ).tolist()
  
  # write new time stamps to UV data
  n = 0
  for group in wiz_uv:
    while ( group.time != time_list[ n ] ):
      n = n + 1
      if ( n == len( time_list ) ):
        raise error( 'time stamp not found' )
    group.time = new_time_list[ n ]
    group.update()
  del wiz_uv
  
  return

###############################################################################

def get_units( uvim ):
  return uvim.header.bunit

###############################################################################

def set_units( uvim, units = 'JY/BEAM' ):
  wiz_uvim = wizardry( uvim )
  wiz_uvim.header.bunit = units
  wiz_uvim.header.update()
  wiz_uvim = wizardry( uvim )
  del wiz_uvim
  return

###############################################################################

def calculate_spectral_index_map( images, noises, snr_cutoff = 2.,
    make_error_map = False, error_cutoff = 0.4, minus_convention = False ):
  
  # perform minimal checks
  if ( len( images ) != len( noises ) ):
    raise error( 'number of images and noises provided do not match' )
  if ( len( images ) != 2 ):
    raise error( 'can only handle 2 images at this time' )
  [ im1, im2 ] = images
  size1 = get_image_size( im1 )
  size2 = get_image_size( im2 )
  if ( size1 != size2 ):
    raise error( 'image sizes do not match' )
  beam1 = around( get_beam_size( im1 ), decimals = 2 ).tolist()
  beam2 = around( get_beam_size( im2 ), decimals = 2 ).tolist()
  if ( beam1 != beam2 ):
    raise error( 'image beam sizes do not match' )
  
  # make copies for spectral index map and error map
  spix_map = get_aips_file( im1.disk, im1.name, 'SPIX', -1, 'MA' )
  call_aips_task( 'MOVE', indata = im1, outdata = spix_map,
      userid = get_aips_userid() )
  spix_pix = get_image_pixels( spix_map ) * 0. + get_aips_magic_value()
  if make_error_map:
    dspix_map = get_aips_file( im1.disk, im1.name, 'DSPIX', -1, 'MA' )
    call_aips_task( 'MOVE', indata = im1, outdata = dspix_map,
        userid = get_aips_userid() )
    dspix_pix = get_image_pixels( dspix_map ) * 0. + get_aips_magic_value()
  
  # get image pixels
  [ n1, n2 ] = noises
  f1 = get_frequency( im1 )
  f2 = get_frequency( im2 )
  pix1 = get_image_pixels( im1 )
  pix2 = get_image_pixels( im2 )
  
  # select pixels for SPIX calculation
  mask = ( pix1 != get_aips_magic_value() )
  mask = mask & ( pix2 != get_aips_magic_value() )
  mask = mask & ( pix1 > snr_cutoff * n1 )
  mask = mask & ( pix2 > snr_cutoff * n2 )
  sel = awhere( mask )
  
  # calculate SPIX and uncertainty
  pix1 = aget( pix1, sel )
  pix2 = aget( pix2, sel )
  spix = log( pix2 / pix1 ) / log( f2 / f1 )
  spix_pix = aput( spix_pix, sel, spix )
  if make_error_map:
    dspix2 = ( 1. / log( f2 / f1 ) )**2 * ( ( n1 / pix1 )**2 + ( n2 / pix2 )**2 )
    dspix = sqrt( dspix2 )
    dspix_pix = aput( dspix_pix, sel, dspix )
    sel = awhere( dspix_pix > error_cutoff )
    spix_pix = aput( spix_pix, sel, get_aips_magic_value() )
    dspix_pix = aput( dspix_pix, sel, get_aips_magic_value() )
  
  # write spectral index map
  if minus_convention:
    spix_pix = -spix_pix
  set_image_pixels( spix_map, spix_pix, units = '' )
  set_units( spix_map, units = '' )
  if make_error_map:
    set_image_pixels( dspix_map, dspix_pix, units = '' )
    set_units( dspix_map, units = '' )
    return [ spix_map, dspix_map ]
  else:
    return spix_map

###############################################################################

def baseline_to_index( uv, baseline, antenna_count = 0 ):
  if ( antenna_count == 0 ):
    ant_count = get_antenna_count( uv )
  else:
    ant_count = antenna_count
  if ( ( min( baseline ) < 1 ) or ( max( baseline ) > ant_count ) ):
    raise error( 'baseline indices are out of range' )
  [ i, j ] = baseline
  if ( i == j ):
    raise error( 'invalid baseline with equal indices' )
  if ( i > j ):
    [ j, i ] = baseline
  index = ( i - 1 ) * ant_count - ( i * ( i - 1 ) / 2 ) + ( j - i - 1 )
  return index

###############################################################################

def index_to_baseline( uv, index, antenna_count = 0 ):
  if ( antenna_count == 0 ):
    ant_count = get_antenna_count( uv )
  else:
    ant_count = antenna_count
  if ( ( index < 0 ) or ( index >= ant_count * ( ant_count - 1 ) / 2 ) ):
    raise error( 'baseline index is out of range' )
  for i in range( 1, ant_count ):
    # j = ant_count
    upper_index = i * ( ant_count - 1 ) - ( i * ( i - 1 ) / 2 ) - 1
    if ( upper_index >= index ):
      break
  j = ant_count - ( upper_index - index )
  return [ i, j ]

###############################################################################

def get_antenna_count( uv ):
  return len( uv.antennas )

###############################################################################

def get_baseline_count( uv ):
  antenna_count = len( uv.antennas )
  return antenna_count * ( antenna_count - 1 ) / 2

###############################################################################

def get_aips_message_file_name( userid = get_aips_userid() ):
  file_path = os.getenv( 'DA01' )
  user = zero_pad_ehex( decimal_to_ehex( userid ), 3 )
  file_name = file_path + '/' + 'MSD' + user + '000.' + user + ';'
  return file_name

###############################################################################

def get_beam_area( im ):
  beam =  convert_beam_size( im, to_pixel = True )
  area = pi * beam[ 0 ] * beam[ 1 ] / ( 4. * log( 2. ) )
  return area

###############################################################################

def get_image_flux( im ):
  pixels = get_image_pixels( im )
  pixels = aput( pixels, awhere( pixels == get_aips_magic_value() ), 0. )
  beam_area = get_beam_area( im )
  flux = pixels.sum() / beam_area
  return flux

###############################################################################

def write_fits_table( uvim, table, file_name, overwrite = True, version = 0,
    fits_table = False ):
  if ( not table_exists( uvim, table, version ) ):
    raise error( 'table does not exists' )
  if ( ( not overwrite ) and ( file_exists( file_name ) ) ):
    raise error( 'table FITS file %s already exists' % ( file_name ) )
  elif ( overwrite and ( file_exists( file_name ) ) ):
    remove_file( file_name )
  if is_uv( uvim ):
    uvim_type = 'UV'
  else:
    uvim_type = 'MA'
  temp_file = get_aips_file( uvim.disk, uvim.name, uvim.klass, -1, uvim_type )
  call_aips_task( 'TASAV', indata = uvim, outdata = temp_file )
  if ( version == 0 ):
    table_version = uvim.table_highver( table )
  else:
    table_version = version
  for [ tversion, tname ] in temp_file.tables:
    tname = ( tname.split() )[ -1 ]
    if ( ( tversion == table_version ) and ( tname == table ) ):
      continue
    temp_file.zap_table( tname, tversion )
  if fits_table:
    call_aips_task( 'FITAB', indata = temp_file, intype = uvim_type,
        dataout = file_name, douvcomp = 0 )
  else:
    call_aips_task( 'FITTP', indata = temp_file, intype = uvim_type,
        dataout = file_name )
  temp_file.zap()
  return

###############################################################################

def convert_xy_to_rl( uv ):
  stokes_index = uv.header.ctype.index( 'STOKES' )
  stokes_rpix = uv.header.crpix[ stokes_index ]
  stokes_rval = uv.header.crval[ stokes_index ]
  stokes_delta = uv.header.cdelt[ stokes_index ]
  stokes_count = uv.header.naxis[ stokes_index ]
  stokes = unique( array( [ stokes_rval + stokes_delta * ( float( i + 1 ) - stokes_rpix )
      for i in range( stokes_count ) ], dtype = int16 ) )
  stokes_xy = alltrue( array( [ s in [ -5,-6,-7,-8 ] for s in stokes ], dtype = bool ) )
  if stokes_xy:
    wiz_uv = wizardry( uv )
    wiz_uv.header.crval[ stokes_index ] += 4.
    wiz_uv.header.update()
  return

###############################################################################

def get_stokes( uv ):
  stokes_index = uv.header.ctype.index( 'STOKES' )
  stokes_rpix = uv.header.crpix[ stokes_index ]
  stokes_rval = uv.header.crval[ stokes_index ]
  stokes_delta = uv.header.cdelt[ stokes_index ]
  stokes_count = uv.header.naxis[ stokes_index ]
  stokes = [ int( stokes_rval + stokes_delta * ( float( i + 1 ) - stokes_rpix ) )
      for i in range( stokes_count ) ]
  stokes_list = []
  for s in stokes:
    if ( s > 0 ):
      stokes_strings = [ 'I', 'Q', 'U', 'V' ]
      stokes_list.append( stokes_strings[ s - 1 ] )
    elif ( s < 0 ):
      stokes_strings = [ 'RR', 'LL', 'RL', 'LR', 'XX', 'YY', 'XY', 'YX' ]
      stokes_list.append( stokes_strings[ - s - 1 ] )
    else:
      return False
  return stokes_list

###############################################################################

def find_integration_time( uv, regularity = [ 0.6, 0.4 ] ):
  time_list = get_time_list( uv )
  time_array = array( time_list, dtype = float64 )
  dtime_array = time_array[ 1 : ] - time_array[ : -1 ]
  dtime = median( dtime_array )
  sel = awhere( abs( dtime_array - dtime ) < regularity[ 1 ] * dtime )
  if ( float( len( sel ) ) < regularity[ 0 ] * float( len( time_list ) ) ):
    dtime_array[ : -1 ] = dtime_array[ : -1 ] + dtime_array[ 1 : ]
    dtime = median( dtime_array )
    sel = awhere( abs( dtime_array - dtime ) < regularity[ 1 ] * dtime )
    if ( float( len( sel ) ) < regularity[ 1 ] * float( len( time_list ) ) ):
      raise error( 'time grid is too irregular to determine integration time' )
  dtime = mean( aget( dtime_array, sel ) )
  return dtime * 86400.

###############################################################################

def find_scans( uv, max_gap = 2. ):
  dtime = find_integration_time( uv ) / 86400.
  time_list = get_time_list( uv )
  last_time = time_list[ 0 ]
  scan_list = []
  scan_entry = [ last_time - 0.5 * dtime ]
  for time in time_list[ 1 : ]:
    if ( time - last_time > max_gap / ( 24. * 60. ) ):
      scan_entry.append( last_time + 0.5 * dtime )
      scan_list.append( scan_entry )
      scan_entry = [ time - 0.5 * dtime ]
    last_time = time
  scan_entry.append( last_time + 0.5 * dtime )
  scan_list.append( scan_entry )
  return scan_list

###############################################################################

def convert_to_galactic( im ):
  radec = get_radec( im )
  lb = equatorial_to_galactic( radec, epoch = get_epoch( im ) )
  dlb = equatorial_to_galactic( [ radec[ 0 ], radec[ 1 ] + 1. ], epoch = get_epoch( im ) )
  [ r, p ] = calculate_angular_separation( lb, dlb )
  for ctype in im.header.ctype:
    if ( ctype.find( 'RA--' ) != - 1 ):
      ra_index = im.header.ctype.index( ctype )
    if ( ctype.find( 'DEC-' ) != - 1 ):
      dec_index = im.header.ctype.index( ctype )
  wiz_im = wizardry( im )
  wiz_im.header.ctype[ ra_index ] = 'GLON' + wiz_im.header.ctype[ ra_index ][ 4 : ]
  wiz_im.header.ctype[ dec_index ] = 'GLAT' + wiz_im.header.ctype[ dec_index ][ 4 : ]
  wiz_im.header.crval[ ra_index ] = lb[ 0 ]
  wiz_im.header.crval[ dec_index ] = lb[ 1 ]
  wiz_im.header.crota[ dec_index ] = wiz_im.header.crota[ dec_index ] - p
  wiz_im.header.update()
  wiz_im = wizardry( im )
  del wiz_im
  return

###############################################################################

def convert_to_equatorial( im, epoch = 2000. ):
  lb = get_glonlat( im )
  radec = galactic_to_equatorial( lb, epoch = epoch )
  dlb = equatorial_to_galactic( [ radec[ 0 ], radec[ 1 ] + 1. ], epoch = epoch )
  [ r, p ] = calculate_angular_separation( lb, dlb )
  for ctype in im.header.ctype:
    if ( ctype.find( 'GLON' ) != - 1 ):
      glon_index = im.header.ctype.index( ctype )
    if ( ctype.find( 'GLAT' ) != - 1 ):
      glat_index = im.header.ctype.index( ctype )
  wiz_im = wizardry( im )
  wiz_im.header.ctype[ glon_index ] = 'RA--' + wiz_im.header.ctype[ glon_index ][ 4 : ]
  wiz_im.header.ctype[ glat_index ] = 'DEC-' + wiz_im.header.ctype[ glat_index ][ 4 : ]
  wiz_im.header.crval[ glon_index ] = radec[ 0 ]
  wiz_im.header.crval[ glat_index ] = radec[ 1 ]
  wiz_im.header.crota[ glat_index ] = around( amodulo( 
      wiz_im.header.crota[ glat_index ] + p + 180., 360 ) - 180., decimals = 4 )
  wiz_im.header.epoch = epoch
  wiz_im.header.update()
  wiz_im = wizardry( im )
  del wiz_im
  return

###############################################################################

def set_array_name( uv, array_name, force = False ):
  wiz_uv = wizardry( uv )
  wiz_uv.header.telescop = array_name
  wiz_uv.header.instrume = array_name
  wiz_uv.header.update()
  wiz_uv = wizardry( uv )
  an_table = wiz_uv.table( 'AN', 0 )
  an_table.keywords[ 'ARRNAM' ] = array_name
  an_table.close()
  del wiz_uv

###############################################################################


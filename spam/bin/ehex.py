#! /usr/bin/python

from sys import *

if len( argv ) == 2:
  user_id = int( argv[ 1 ] )
  div = [ 0, 0, 0 ]
  div[ 2 ] = user_id % 36
  div[ 1 ] = ( ( user_id - div[ 2 ] ) / 36 ) % 36
  div[ 0 ] = ( ( ( ( user_id - div[ 2 ] ) / 36 ) - div[ 1 ] ) / 36 ) % 36
  ehex = [ '', '', '' ]
  for i in range( 3 ):
    if div[ i ] <= 9:
      ehex[ i ] = '%d' % ( div[ i ] )
    else:
      ehex[ i ] = '%c' % ( 55 + div[ i ] )
  print '%c%c%c' % ( ehex[ 0 ], ehex[ 1 ], ehex[ 2 ] )

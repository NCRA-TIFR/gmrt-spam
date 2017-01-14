//*****************************************************************************

// module name
//%module _interface

//*****************************************************************************

// typemap to convert an input list of floats to an array of doubles
%typemap(in) double [ANY] ( double temp[ $1_dim0 ] )
{
  if ( $input == NULL )
    return ( NULL );
  if ( !convert_float_list_to_double_array( $input, temp, $1_dim0 ) )
    return NULL;
  $1 = &temp[ 0 ];
}

// typemap to convert arrays of doubles to output lists of floats
%typemap(out) p_double_2
{  if ( $1 == NULL )
    return ( NULL );
  PyObject *output = Py_BuildValue( "[ff]", $1[ 0 ], $1[ 1 ] );
  if ( $1 != NULL )
    free( $1 );
  $result = output;
}
%typemap(out) p_double_3
{  if ( $1 == NULL )
    return ( NULL );
  PyObject *output = Py_BuildValue( "[fff]", $1[ 0 ], $1[ 1 ], $1[ 2 ] );
  if ( $1 != NULL )
    free( $1 );
  $result = output;
}
%typemap(out) p_double_4
{  if ( $1 == NULL )
    return ( NULL );
  PyObject *output = Py_BuildValue( "[ffff]", $1[ 0 ], $1[ 1 ], $1[ 2 ], $1[ 3 ] );
  if ( $1 != NULL )
    free( $1 );
  $result = output;
}
%typemap(out) p_double_5
{  if ( $1 == NULL )
    return ( NULL );
  PyObject *output = Py_BuildValue( "[fffff]", $1[ 0 ], $1[ 1 ], $1[ 2 ], $1[ 3 ], $1[ 4 ] );
  if ( $1 != NULL )
    free( $1 );
  $result = output;
}
%typemap(out) p_double_6
{  if ( $1 == NULL )
    return ( NULL );
  PyObject *output = Py_BuildValue( "[ffffff]", $1[ 0 ], $1[ 1 ], $1[ 2 ], $1[ 3 ], $1[ 4 ], $1[ 5 ] );
  if ( $1 != NULL )
    free( $1 );
  $result = output;
}

// typemap to convert a multidimensional input list of floats to an array of doubles
%typemap( in ) p_double_array
{
  if ( $input == NULL )
    return ( NULL );
  $1 = convert_multi_float_list_to_double_array( $input );
  if ( $1 == NULL )
    return ( NULL );
}
// typemap to convert a multidimensional input array of doubles to a list of floats
%typemap( out ) p_double_array
{
  if ( $1 == NULL )
    return ( NULL );
  $result = convert_multi_double_array_to_float_list( $1 );
  if ( $result == NULL )
    return ( NULL );
}

// typemap to convert a multidimensional input list of integers to an array of longs
%typemap( in ) p_long_array
{
  if ( $input == NULL )
    return ( NULL );
  $1 = convert_multi_integer_list_to_long_array( $input );
  if ( $1 == NULL )
    return ( NULL );
}
// typemap to convert a multidimensional input array of longs to a list of integers
%typemap( out ) p_long_array
{
  if ( $1 == NULL )
    return ( NULL );
  $result = convert_multi_long_array_to_integer_list( $1 );
  if ( $result == NULL )
    return ( NULL );
}

//*****************************************************************************

%include "__interface.h"

//*****************************************************************************

// included in C wrapper
%{
#include "__types.h"
#include "__interface.h"
%}

//*****************************************************************************

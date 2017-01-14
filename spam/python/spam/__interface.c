//*****************************************************************************

// Python header files
#include <Python.h>

// user header files
#include "__types.h"
#include "__interface.h"

//*****************************************************************************

int convert_float_list_to_double_array( PyObject *input, double *ptr, int size )
{
  int i, once, ok = 1;
  for ( once = 1; once; once-- )
  {
    if ( !PyList_Check( input ) )
    {
      PyErr_SetString( PyExc_TypeError, "input not a list" );
      ok = 0; break;
    }
    if ( PyObject_Length( input ) != size )
    {
      PyErr_SetString( PyExc_ValueError, "input list size mismatch" );
      ok = 0; break;
    }
    for ( i = 0; i < size; i++ )
    {
      PyObject *item = PyList_GetItem( input, i );
      if ( !PyFloat_Check( item ) )
      {
        PyErr_SetString( PyExc_TypeError, "input not a list of floats" );
        ok = 0; break;
      }
      ptr[ i ] = PyFloat_AsDouble( item );
    }
    if ( !ok ) break;
  }
//  Py_XDECREF( input );
  return ( ok );
}

//*****************************************************************************

int recursive_convert_float_list_to_double_array( PyObject *fl, p_double_array da, long dim_i, long *data_i )
{
  int once, ok = 1;
  long dim_size, i;
  PyObject *iter;
//  char info_string[ 100 ];
  for ( once = 1; once; once-- )
  {
    if ( dim_i < 0L  )
    {
      PyErr_SetString( PyExc_RuntimeError, "illegal dimension index lower than zero" );
      ok = 0; break;
    }
    if ( dim_i >= da->num_dims  )
    {
      PyErr_SetString( PyExc_ValueError, "input list has variable dimension sizes" );
      ok = 0; break;
    }
    if ( !PyList_Check( fl ) )
    {
      PyErr_SetString( PyExc_ValueError, "input list has variable content" );
      ok = 0; break;
    }
    dim_size = ( long )PyObject_Length( fl );
    if ( dim_size < 0 )
    {
      PyErr_SetString( PyExc_RuntimeError, "length of input list dimension could not be determined" );
      ok = 0; break;
    }
    if ( dim_size == 0 )
    {
      PyErr_SetString( PyExc_ValueError, "input list dimension is empty" );
      ok = 0; break;
    }
    if ( dim_size != da->dims[ dim_i ] )
    {
      PyErr_SetString( PyExc_ValueError, "input list has variable dimension sizes" );
      ok = 0; break;
    }
    for ( i = 0L; i < da->dims[ dim_i ]; i++ )
    {
      iter = PyList_GetItem( fl, i );
      if ( dim_i < da->num_dims - 1L )
      {
        ok = recursive_convert_float_list_to_double_array( iter, da, dim_i + 1L, data_i );
        if ( !ok ) break;
      }
      else // dim_i == da->num_dims - 1L
      {
        if ( !PyFloat_Check( iter ) )
        {
          if ( PyList_Check( fl ) )
            PyErr_SetString( PyExc_ValueError, "input list has variable dimension sizes" );
          else
            PyErr_SetString( PyExc_TypeError, "input is not a list of floats" );
          ok = 0; break;
        }
        da->data[ ( *data_i )++ ] = PyFloat_AsDouble( iter );
      }
    }
  }
  return ( ok );
}

//*****************************************************************************

p_double_array convert_multi_float_list_to_double_array( PyObject *fl )
{
  int once, ok = 1;
  long dim_size, dim_i, dim_sum, dims[ __TYPES_MAX_DIMS ], data_i;
  p_double_array da = NULL;
  PyObject *iter;
  for ( once = 1; once; once-- )
  {
    if ( fl == NULL )
    {
      PyErr_SetString( PyExc_TypeError, "input is not a list" );
      ok = 0; break;
    }
    if ( fl == Py_None )
    {
      da = ( p_double_array )malloc( sizeof( struct double_array ) );
      if ( da == NULL )
      {
        PyErr_SetString( PyExc_RuntimeError, "could not allocate memory for double_array" );
        ok = 0; break;
      }
      da->num_dims = - 1L;
      da->dims = NULL;
      da->data = NULL;
      break;
    }
    if ( !PyList_Check( fl ) )
    {
      PyErr_SetString( PyExc_TypeError, "input is not a list" );
      ok = 0; break;
    }
    dim_size = ( long )PyObject_Length( fl );
    if ( dim_size < 0L )
    {
      PyErr_SetString( PyExc_RuntimeError, "length of input list could not be determined" );
      ok = 0; break;
    }
    if ( dim_size == 0L )
    {
      PyErr_SetString( PyExc_ValueError, "input list is empty" );
      ok = 0; break;
    }
    dim_i = 0L;
    dims[ dim_i++ ] = dim_size;
    iter = PyList_GetItem( fl, 0 );
    if ( iter == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "could not get item from list" );
      ok = 0; break;
    }
    while ( ( dim_i < __TYPES_MAX_DIMS ) && ( PyList_Check( iter ) ) )
    {
      dim_size = ( long )PyObject_Length( iter );
      if ( dim_size < 0 )
      {
        PyErr_SetString( PyExc_RuntimeError, "length of input list could not be determined" );
        ok = 0; break;
      }
      if ( dim_size == 0 )
      {
        PyErr_SetString( PyExc_ValueError, "input list has empty dimension" );
        ok = 0; break;
      }
      dims[ dim_i++ ] = dim_size;
      iter = PyList_GetItem( iter, 0 );
      if ( iter == NULL )
      {
        PyErr_SetString( PyExc_RuntimeError, "could not get item from list" );
        ok = 0; break;
      }
    }
    if ( !ok ) break;
    if ( ( dim_i == __TYPES_MAX_DIMS ) && ( PyList_Check( iter ) ) )
    {
      PyErr_SetString( PyExc_ValueError, "input list has reached maximum number of dimensions" );
      ok = 0; break;
    }
    if ( !PyFloat_Check( iter ) )
    {
      PyErr_SetString( PyExc_TypeError, "input is not a list of floats" );
      ok = 0; break;
    }
    da = ( p_double_array )malloc( sizeof( struct double_array ) );
    if ( da == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "could not allocate memory for double_array" );
      ok = 0; break;
    }
    da->data = NULL;
    da->num_dims = dim_i;
    da->dims = ( long * )malloc( da->num_dims * sizeof( long ) );
    if ( da->dims == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "could not allocate memory for double_array dims" );
      ok = 0; break;
    }
    dim_sum = 1L;
    for ( dim_i = 0L; dim_i < da->num_dims; dim_i++ )
    {
      da->dims[ dim_i ] = dims[ dim_i ];
      dim_sum = dim_sum * dims[ dim_i ];
    }
    da->data = ( double * )malloc( dim_sum * sizeof( double ) );
    if ( da->data == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "could not allocate memory for double_array data" );
    }
    data_i = 0L;
    ok = recursive_convert_float_list_to_double_array( fl, da, 0L, &data_i );
    if ( !ok ) break;
  }
  Py_XDECREF( fl );
  if ( !ok )
  {
    free_double_array( da );
    return ( NULL );
  }
  else
    return ( da );
}

//*****************************************************************************

PyObject *recursive_convert_double_array_to_float_list( p_double_array da, long dim_i, long *data_i )
{
  int once, ok = 1;
  long i;
  PyObject *fl, *iter;
  for ( once = 1; once; once-- )
  {
    if ( ( dim_i < 0L  ) || ( dim_i >= da->num_dims ) )
    {
      PyErr_SetString( PyExc_RuntimeError, "double_array dimension index out of bounds" );
      ok = 0; break;
    }
    fl = PyList_New( da->dims[ dim_i ] );
    if ( fl == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "creation of new list dimension failed" );
      ok = 0; break;
    }
    for ( i = 0L; i < da->dims[ dim_i ]; i++ )
    {
      if ( dim_i < da->num_dims - 1L )
      {
        iter = recursive_convert_double_array_to_float_list( da, dim_i + 1L, data_i );
        if ( iter == NULL )
        {
          // error message has already been generated
          ok = 0; break;
        }
      }
      else // dim_i == da->num_dims - 1L
      {
        iter = PyFloat_FromDouble( da->data[ ( *data_i )++ ] );
      }
      if ( PyList_SetItem( fl, i, iter ) )
      {
        PyErr_SetString( PyExc_RuntimeError, "could not add new item to list" );
        ok = 0; break;
      }
    }
    if ( !ok ) break;
  }
  if ( !ok )
  {
    Py_XDECREF( fl );
    Py_XDECREF( iter );
    return ( NULL );
  }
  else
  {
    return ( fl );
  }
}

//*****************************************************************************

PyObject *convert_multi_double_array_to_float_list( p_double_array da )
{
  int once, ok = 1;
  long dim_i, da_i;
  PyObject *fl;
  for ( once = 1; once; once-- )
  {
    if ( da == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "pointer to double_array is invalid" );
      ok = 0; break;
    }
    if ( ( da->num_dims <= 0L  ) || ( da->num_dims > __TYPES_MAX_DIMS ) )
    {
      PyErr_SetString( PyExc_RuntimeError, "number of dimensions of double_array is out of bounds" );
      ok = 0; break;
    }
    for ( dim_i = 0L; dim_i < da->num_dims; dim_i++ )
    {
      if ( da->dims[ dim_i ] <= 0L  )
      {
        PyErr_SetString( PyExc_RuntimeError, "double_array dimension size out of bounds" );
        ok = 0; break;
      }
    }
    if ( !ok ) break;
    da_i = 0L;
    fl = recursive_convert_double_array_to_float_list( da, 0L, &da_i );
    if ( fl == NULL )
    {
      // error message has already been generated
      ok = 0; break;
    }
  }
  if ( da != NULL )
    free_double_array( da );
  if ( !ok )
  {
    Py_XDECREF( fl ); 
    return ( NULL );
  }
  else
    return ( fl );
}

//*****************************************************************************

int recursive_convert_integer_list_to_long_array( PyObject *il, p_long_array la, long dim_i, long *data_i )
{
  int once, ok = 1;
  long dim_size, i;
  PyObject *iter;
//  char info_string[ 100 ];
  for ( once = 1; once; once-- )
  {
    if ( dim_i < 0L  )
    {
      PyErr_SetString( PyExc_RuntimeError, "illegal dimension index lower than zero" );
      ok = 0; break;
    }
    if ( dim_i >= la->num_dims  )
    {
      PyErr_SetString( PyExc_ValueError, "input list has variable dimension sizes" );
      ok = 0; break;
    }
    if ( !PyList_Check( il ) )
    {
      PyErr_SetString( PyExc_ValueError, "input list has variable content" );
      ok = 0; break;
    }
    dim_size = ( long )PyObject_Length( il );
    if ( dim_size < 0 )
    {
      PyErr_SetString( PyExc_RuntimeError, "length of input list dimension could not be determined" );
      ok = 0; break;
    }
    if ( dim_size == 0 )
    {
      PyErr_SetString( PyExc_ValueError, "input list dimension is empty" );
      ok = 0; break;
    }
    if ( dim_size != la->dims[ dim_i ] )
    {
      PyErr_SetString( PyExc_ValueError, "input list has variable dimension sizes" );
      ok = 0; break;
    }
    for ( i = 0L; i < la->dims[ dim_i ]; i++ )
    {
      iter = PyList_GetItem( il, i );
      if ( dim_i < la->num_dims - 1L )
      {
        ok = recursive_convert_integer_list_to_long_array( iter, la, dim_i + 1L, data_i );
        if ( !ok ) break;
      }
      else // dim_i == la->num_dims - 1L
      {
        if ( !PyInt_Check( iter ) )
        {
          if ( PyList_Check( il ) )
            PyErr_SetString( PyExc_ValueError, "input list has variable dimension sizes" );
          else
            PyErr_SetString( PyExc_TypeError, "input is not a list of integers" );
          ok = 0; break;
        }
        la->data[ ( *data_i )++ ] = PyInt_AsLong( iter );
      }
    }
  }
  return ( ok );
}

//*****************************************************************************

p_long_array convert_multi_integer_list_to_long_array( PyObject *il )
{
  int once, ok = 1;
  long dim_size, dim_i, dim_sum, dims[ __TYPES_MAX_DIMS ], data_i;
  p_long_array la = NULL;
  PyObject *iter;
  for ( once = 1; once; once-- )
  {
    if ( il == NULL )
    {
      PyErr_SetString( PyExc_TypeError, "input is not a list" );
      ok = 0; break;
    }
    if ( il == Py_None )
    {
      la = ( p_long_array )malloc( sizeof( struct long_array ) );
      if ( la == NULL )
      {
        PyErr_SetString( PyExc_RuntimeError, "could not allocate memory for long_array" );
        ok = 0; break;
      }
      la->num_dims = - 1L;
      la->dims = NULL;
      la->data = NULL;
      break;
    }
    if ( !PyList_Check( il ) )
    {
      PyErr_SetString( PyExc_TypeError, "input is not a list" );
      ok = 0; break;
    }
    dim_size = ( long )PyObject_Length( il );
    if ( dim_size < 0L )
    {
      PyErr_SetString( PyExc_RuntimeError, "length of input list could not be determined" );
      ok = 0; break;
    }
    if ( dim_size == 0L )
    {
      PyErr_SetString( PyExc_ValueError, "input list is empty" );
      ok = 0; break;
    }
    dim_i = 0L;
    dims[ dim_i++ ] = dim_size;
    iter = PyList_GetItem( il, 0 );
    if ( iter == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "could not get item from list" );
      ok = 0; break;
    }
    while ( ( dim_i < __TYPES_MAX_DIMS ) && ( PyList_Check( iter ) ) )
    {
      dim_size = ( long )PyObject_Length( iter );
      if ( dim_size < 0 )
      {
        PyErr_SetString( PyExc_RuntimeError, "length of input list could not be determined" );
        ok = 0; break;
      }
      if ( dim_size == 0 )
      {
        PyErr_SetString( PyExc_ValueError, "input list has empty dimension" );
        ok = 0; break;
      }
      dims[ dim_i++ ] = dim_size;
      iter = PyList_GetItem( iter, 0 );
      if ( iter == NULL )
      {
        PyErr_SetString( PyExc_RuntimeError, "could not get item from list" );
        ok = 0; break;
      }
    }
    if ( !ok ) break;
    if ( ( dim_i == __TYPES_MAX_DIMS ) && ( PyList_Check( iter ) ) )
    {
      PyErr_SetString( PyExc_ValueError, "input list has reached maximum number of dimensions" );
      ok = 0; break;
    }
    if ( !PyInt_Check( iter ) )
    {
      PyErr_SetString( PyExc_TypeError, "input is not a list of integers" );
      ok = 0; break;
    }
    la = ( p_long_array )malloc( sizeof( struct long_array ) );
    if ( la == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "could not allocate memory for long_array" );
      ok = 0; break;
    }
    la->data = NULL;
    la->num_dims = dim_i;
    la->dims = ( long * )malloc( la->num_dims * sizeof( long ) );
    if ( la->dims == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "could not allocate memory for long_array dims" );
      ok = 0; break;
    }
    dim_sum = 1L;
    for ( dim_i = 0L; dim_i < la->num_dims; dim_i++ )
    {
      la->dims[ dim_i ] = dims[ dim_i ];
      dim_sum = dim_sum * dims[ dim_i ];
    }
    la->data = ( long * )malloc( dim_sum * sizeof( long ) );
    if ( la->data == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "could not allocate memory for long_array data" );
    }
    data_i = 0L;
    ok = recursive_convert_integer_list_to_long_array( il, la, 0L, &data_i );
    if ( !ok ) break;
  }
  Py_XDECREF( il ); 
  if ( !ok )
  {
    free_long_array( la );
    return ( NULL );
  }
  else
    return ( la );
}

//*****************************************************************************

PyObject *recursive_convert_long_array_to_integer_list( p_long_array la, long dim_i, long *data_i )
{
  int once, ok = 1;
  long i;
  PyObject *il, *iter;
  for ( once = 1; once; once-- )
  {
    if ( ( dim_i < 0L  ) || ( dim_i >= la->num_dims ) )
    {
      PyErr_SetString( PyExc_RuntimeError, "double_array dimension index out of bounds" );
      ok = 0; break;
    }
    il = PyList_New( la->dims[ dim_i ] );
    if ( il == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "creation of new list dimension failed" );
      ok = 0; break;
    }
    for ( i = 0L; i < la->dims[ dim_i ]; i++ )
    {
      if ( dim_i < la->num_dims - 1L )
      {
        iter = recursive_convert_long_array_to_integer_list( la, dim_i + 1L, data_i );
        if ( iter == NULL )
        {
          // error message has already been generated
          ok = 0; break;
        }
      }
      else // dim_i == da->num_dims - 1L
      {
        iter = PyInt_FromLong( la->data[ ( *data_i )++ ] );
      }
      if ( PyList_SetItem( il, i, iter ) )
      {
        PyErr_SetString( PyExc_RuntimeError, "could not add new item to list" );
        ok = 0; break;
      }
    }
    if ( !ok ) break;
  }
  if ( !ok )
  {
    Py_XDECREF( il );
    Py_XDECREF( iter );
    return ( NULL );
  }
  else
  {
    return ( il );
  }
}

//*****************************************************************************

PyObject *convert_multi_long_array_to_integer_list( p_long_array la )
{
  int once, ok = 1;
  long dim_i, la_i;
  PyObject *il;
  for ( once = 1; once; once-- )
  {
    if ( la == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "pointer to long_array is invalid" );
      ok = 0; break;
    }
    if ( ( la->num_dims <= 0L  ) || ( la->num_dims > __TYPES_MAX_DIMS ) )
    {
      PyErr_SetString( PyExc_RuntimeError, "number of dimensions of long_array is out of bounds" );
      ok = 0; break;
    }
    for ( dim_i = 0L; dim_i < la->num_dims; dim_i++ )
    {
      if ( la->dims[ dim_i ] <= 0L  )
      {
        PyErr_SetString( PyExc_RuntimeError, "long_array dimension size out of bounds" );
        ok = 0; break;
      }
    }
    if ( !ok ) break;
    la_i = 0L;
    il = recursive_convert_long_array_to_integer_list( la, 0L, &la_i );
    if ( il == NULL )
    {
      // error message has already been generated
      ok = 0; break;
    }
  }
  if ( la != NULL )
    free_long_array( la );
  if ( !ok )
  {
    Py_XDECREF( il ); 
    return ( NULL );
  }
  else
    return ( il );
}

//*****************************************************************************

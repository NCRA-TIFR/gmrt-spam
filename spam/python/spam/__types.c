//*****************************************************************************

// Python header files
#include <Python.h>

// standard C header files
#include <stdio.h>
#include <stdlib.h>

// user header files
#include "__types.h"

//*****************************************************************************

p_double_array create_double_array( long num_dims, long *dims, double *data )
{
  long dim_sum, dim_i, data_i;
  p_double_array da;
  int once, ok = 1;
  for ( once = 1; once; once-- )
  {
    da = ( p_double_array )malloc( sizeof( struct double_array ) );
    if ( da == NULL )
    {
      ok = 0; break;
    }
    da->data = NULL;
    da->num_dims = num_dims;
    da->dims = ( long * )malloc( da->num_dims * sizeof( long ) );
    if ( da->dims == NULL )
    {
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
      ok = 0; break;
    }
    if ( data != NULL )
    {
      for ( data_i = 0L; data_i < dim_sum; data_i++ )
        da->data[ data_i ] = data[ data_i ];
    }
  }
  if ( !ok )
  {
    free_double_array( da );
    return ( NULL );
  }
  else
    return ( da );
}

//*****************************************************************************

void free_double_array( p_double_array da )
{
  if ( da != NULL )
  {
    if ( da->dims != NULL )
      free( da->dims );
    if ( da->data != NULL )
      free( da->data );
    free( da );
  }
}

//*****************************************************************************

int check_double_array( p_double_array da, long num_dims, long *dims )
{
  int once, ok = 1;
  long dim_i;
  for ( once = 1; once; once-- )
  {
    if ( da == NULL )
    {
      ok = 0; break;
    }
    if ( da->num_dims != num_dims )
    {
      ok = 0; break;
    }
    for ( dim_i = 0L; dim_i < da->num_dims; dim_i++ )
    {
      if ( ( dims[ dim_i ] > 0L ) && ( da->dims[ dim_i ] != dims[ dim_i ] ) )
      {
        ok = 0; break;
      }
    }
    if ( !ok ) break;
  }
  return ( ok );
}

//*****************************************************************************

long double_array_dims_to_index( p_double_array da, long *dims )
{
  int once, ok = 1;
  long dim_i, index;
  for ( once = 1; once; once-- )
  {
    if ( da == NULL )
    {
      ok = 0; break;
    }
    index = 0L;
    for ( dim_i = 0L; dim_i < da->num_dims; dim_i++ )
    {
      if ( ( dims[ dim_i ] < 0L ) || ( dims[ dim_i ] >= da->dims[ dim_i ] ) )
      {
        ok = 0; break;
      }
      index = da->dims[ dim_i ] * index + dims[ dim_i ];
    }
    if ( !ok ) break;
  }
  if ( !ok )
    index = - 1L;
  return ( index );
}

//*****************************************************************************

p_long_array create_long_array( long num_dims, long *dims, long *data )
{
  long dim_sum, dim_i, data_i;
  p_long_array la;
  int once, ok = 1;
  for ( once = 1; once; once-- )
  {
    la = ( p_long_array )malloc( sizeof( struct long_array ) );
    if ( la == NULL )
    {
      ok = 0; break;
    }
    la->data = NULL;
    la->num_dims = num_dims;
    la->dims = ( long * )malloc( la->num_dims * sizeof( long ) );
    if ( la->dims == NULL )
    {
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
      ok = 0; break;
    }
    if ( data != NULL )
    {
      for ( data_i = 0L; data_i < dim_sum; data_i++ )
        la->data[ data_i ] = data[ data_i ];
    }
  }
  if ( !ok )
  {
    free_long_array( la );
    return ( NULL );
  }
  else
    return ( la );
}

//*****************************************************************************

void free_long_array( p_long_array la )
{
  if ( la != NULL )
  {
    if ( la->dims != NULL )
      free( la->dims );
    if ( la->data != NULL )
      free( la->data );
    free( la );
  }
}

//*****************************************************************************

int check_long_array( p_long_array la, long num_dims, long *dims )
{
  int once, ok = 1;
  long dim_i;
  for ( once = 1; once; once-- )
  {
    if ( la == NULL )
    {
      ok = 0; break;
    }
    if ( la->num_dims != num_dims )
    {
      ok = 0; break;
    }
    for ( dim_i = 0L; dim_i < la->num_dims; dim_i++ )
    {
      if ( ( dims[ dim_i ] > 0L ) && ( la->dims[ dim_i ] != dims[ dim_i ] ) )
      {
        ok = 0; break;
      }
    }
    if ( !ok ) break;
  }
  return ( ok );
}

//*****************************************************************************

long long_array_dims_to_index( p_long_array la, long *dims )
{
  int once, ok = 1;
  long dim_i, index;
  for ( once = 1; once; once-- )
  {
    if ( la == NULL )
    {
      ok = 0; break;
    }
    index = 0L;
    for ( dim_i = 0L; dim_i < la->num_dims; dim_i++ )
    {
      if ( ( dims[ dim_i ] < 0L ) || ( dims[ dim_i ] >= la->dims[ dim_i ] ) )
      {
        ok = 0; break;
      }
      index = la->dims[ dim_i ] * index + dims[ dim_i ];
    }
    if ( !ok ) break;
  }
  if ( !ok )
    index = - 1L;
  return ( index );
}

//*****************************************************************************

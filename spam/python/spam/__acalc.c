//*****************************************************************************

// Python header files
#include <Python.h>

// standard C header files
#include <stdio.h>
#include <string.h>
#include <math.h>

// user header files
#include "__types.h"
#include "__acalc.h"

//*****************************************************************************

double sign( double x )
{
  double s = ( x < 0. ? -1. : 1. ) ;
  return ( s );
}

//*****************************************************************************

double degrees( double rad )
{
  double deg = 180. * rad / M_PI;
  return ( deg );
}

//*****************************************************************************

double radians( double deg )
{
  double rad = M_PI * deg / 180.;
  return ( rad );
}

//*****************************************************************************

double modulo( double a, double b )
{
  double c = fmod( a, b );
  if ( b * c < 0. )
  {
    c = c + b;
  }
  return ( c );
}

//*****************************************************************************

long int factorial( int n )
{
  long int fac, nn;
  nn = ( long int )n;
  fac = 1L;
  while ( nn > 0 )
  {
    fac = fac * nn;
    nn = nn - 1L;
  }
  return ( fac );
}

//*****************************************************************************

long int binomial( int n, int k )
{
  double binom;
  binom = ( long int )factorial( n ) / ( long int )( factorial( k ) * factorial( n - k ) );
  return ( binom );
}

//*****************************************************************************

p_double_2 complex_to_r_phi( double c[ 2 ] )
{
  p_double_2 rp = malloc( 2 * sizeof( double ) );
  rp[ 0 ] = sqrt( c[ 0 ] * c[ 0 ] + c[ 1 ] * c[ 1 ] );
  rp[ 1 ] = degrees( atan2( c[ 1 ], c[ 0 ] ) );
  return ( rp );
}

//*****************************************************************************

p_double_2 r_phi_to_complex( double rp[ 2 ] )
{
  p_double_2 c = malloc( 2 * sizeof( double ) );
  c[ 0 ] = rp[ 0 ] * cos( radians( rp[ 1 ] ) );
  c[ 1 ] = rp[ 0 ] * sin( radians( rp[ 1 ] ) );
  return ( c );
}

//*****************************************************************************

double maximum( p_double_array da )
{
  int once, ok = 1;
  long dim_i, da_index, da_size;
  double maximum;
  for ( once = 1; once; once-- )
  {
    if ( da == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "da is invalid" );
      ok = 0; break;
    }
    da_size = 1L;
    for ( dim_i = 0L; dim_i < da->num_dims; dim_i++ )
      da_size = da_size * da->dims[ dim_i ];
    maximum = da->data[ 0 ];
    for ( da_index = 0L; da_index < da_size; da_index++ )
      if ( da->data[ da_index ] > maximum )
        maximum = da->data[ da_index ];
  }

  // cleanup
  if ( da != NULL )
  {
    free_double_array( da ); da = NULL;
  }
  if ( !ok )
    return ( 0. );
  else
    return ( maximum );
}

//*****************************************************************************

double minimum( p_double_array da )
{
  int once, ok = 1;
  long dim_i, da_index, da_size;
  double minimum;
  for ( once = 1; once; once-- )
  {
    if ( da == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "da is invalid" );
      ok = 0; break;
    }
    da_size = 1L;
    for ( dim_i = 0L; dim_i < da->num_dims; dim_i++ )
      da_size = da_size * da->dims[ dim_i ];
    minimum = da->data[ 0 ];
    for ( da_index = 0L; da_index < da_size; da_index++ )
      if ( da->data[ da_index ] < minimum )
        minimum = da->data[ da_index ];
  }

  // cleanup
  if ( da != NULL )
  {
    free_double_array( da ); da = NULL;
  }
  if ( !ok )
    return ( 0. );
  else
  return ( minimum );
}

//*****************************************************************************

double median( p_double_array da )
{
  int once, ok = 1;
  long dim_i, da_mid_size, da_size, i, j, median_index;
  double median, *da_mid = NULL, median_2;
  for ( once = 1; once; once-- )
  {
    if ( da == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "da is invalid" );
      ok = 0; break;
    }
    da_size = 1L;
    for ( dim_i = 0L; dim_i < da->num_dims; dim_i++ )
      da_size = da_size * da->dims[ dim_i ];
    if ( ( da_size & 1L ) != 0L )
      da_mid_size = ( da_size + 1L ) / 2L;
    else
      da_mid_size = ( da_size + 2L ) / 2L;
    da_mid = ( double * )malloc( da_mid_size * sizeof( double ) );
    if ( da_mid == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "allocation of memory for da_mid failed" );
      ok = 0; break;
    }
    median = da->data[ 0 ];
    median_index = 0L;
    for ( i = 0L; i < da_mid_size; i++ )
    {
      da_mid[ i ] = da->data[ i ];
      if ( da->data[ i ] > median )
      {
        median = da->data[ i ];
        median_index = i;
      }
    }
    for ( i = da_mid_size; i < da_size; i++ )
    {
      if ( da->data[ i ] < median )
      {
        median = da->data[ i ];
        da_mid[ median_index ] = median;
        for ( j = 0L; j < da_mid_size; j++ )
        {
          if ( da_mid[ j ] > median )
          {
            median = da_mid[ j ];
            median_index = j;
          }
        }
      }
    }
    if ( ( da_size & 1L ) == 0L )
    {
      if ( median_index != 0L )
        median_2 = da_mid[ 0 ];
      else
        median_2 = da_mid[ 1 ];
      for ( i = 0L; i < da_mid_size; i++ )
        if ( ( i != median_index ) && ( da_mid[ i ] > median_2 ) )
          median_2 = da_mid[ i ];
      median = ( median + median_2 ) / 2.;
    }
  }
  if ( da != NULL )
  {
    free_double_array( da ); da = NULL;
  }
  if ( da_mid != NULL )
  {
    free( da_mid ); da_mid = NULL;
  }
  if ( !ok )
    return ( 0. );
  else
    return ( median );
}

//*****************************************************************************

double sum( p_double_array da )
{
  int once, ok = 1;
  long dim_i, da_index, da_size;
  double sum;
  for ( once = 1; once; once-- )
  {
    if ( da == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "da is invalid" );
      ok = 0; break;
    }
    da_size = 1L;
    for ( dim_i = 0L; dim_i < da->num_dims; dim_i++ )
      da_size = da_size * da->dims[ dim_i ];
    sum = 0.;
    for ( da_index = 0L; da_index < da_size; da_index++ )
      sum = sum + da->data[ da_index ];
  }

  // cleanup
  if ( da != NULL )
  {
    free_double_array( da ); da = NULL;
  }
  if ( !ok )
    return ( 0. );
  else
  return ( sum );
}

//*****************************************************************************

p_double_array square( p_double_array da )
{
  int once, ok = 1;
  long sq_size, dim_i, sq_index;
  p_double_array sq = NULL;
  for ( once = 1; once; once-- )
  {
    if ( da == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "da is invalid" );
      ok = 0; break;
    }
    sq = create_double_array( da->num_dims, da->dims, da->data );
    if ( sq == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "allocation of memory for sq failed" );
      ok = 0; break;
    }
    sq_size = 1L;
    for ( dim_i = 0L; dim_i < sq->num_dims; dim_i++ )
      sq_size = sq_size * sq->dims[ dim_i ];
    for ( sq_index = 0L; sq_index < sq_size; sq_index++ )
      sq->data[ sq_index ] = sq->data[ sq_index ] * sq->data[ sq_index ];
  }

  // cleanup
  if ( da != NULL )
  {
    free_double_array( da ); da = NULL;
  }
  if ( !ok )
  {
    if ( sq != NULL )
    {
      free_double_array( sq ); sq = NULL;
    }
  }
  return ( sq );
}

//*****************************************************************************

p_double_array multiply( p_double_array da, double factor )
{
  int once, ok = 1;
  long mul_size, dim_i, mul_index;
  p_double_array mul = NULL;
  for ( once = 1; once; once-- )
  {
    if ( da == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "da is invalid" );
      ok = 0; break;
    }
    mul = create_double_array( da->num_dims, da->dims, da->data );
    if ( mul == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "allocation of memory for mul failed" );
      ok = 0; break;
    }
    mul_size = 1L;
    for ( dim_i = 0L; dim_i < mul->num_dims; dim_i++ )
      mul_size = mul_size * mul->dims[ dim_i ];
    for ( mul_index = 0L; mul_index < mul_size; mul_index++ )
      mul->data[ mul_index ] = mul->data[ mul_index ] * factor;
  }

  // cleanup
  if ( da != NULL )
  {
    free_double_array( da ); da = NULL;
  }
  if ( !ok )
  {
    if ( mul != NULL )
    {
      free_double_array( mul ); mul = NULL;
    }
  }
  return ( mul );
}


//*****************************************************************************

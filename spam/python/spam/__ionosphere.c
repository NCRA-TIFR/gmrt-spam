//*****************************************************************************

// Python header files
#include <Python.h>

// standard C header files
#include <stdio.h>
#include <string.h>
#include <math.h>

// user header files
#include "__types.h"
#include "__sphere.h"
#include "__acalc.h"
#include "__ionosphere.h"

//*****************************************************************************

p_double_array calculate_phase_gradients( p_double_array phase_table )
{
  int once, ok = 1;
  long phase_dims[ 2 ], gradient_dims[ 2 ], phase_count, gradient_count, i, j;
  long phase_index, gradient_index;
  double X_i[ 2 ], phase_i, error_i, X_j[ 2 ], phase_j, error_j;
  double phase_ij, error_ij, ra_ji[ 2 ], phase_ji, error_ji;
  p_double_2 ra_ij = NULL;
  p_double_array gradient_table = NULL;
  for ( once = 1; once; once-- )
  {
    phase_dims[ 0 ] = 0L; phase_dims[ 1 ] = 4L;
    if ( !check_double_array( phase_table, 2L, &phase_dims[ 0 ] ) )
    {
      PyErr_SetString( PyExc_RuntimeError, "phase_table has wrong dimensions" );
      ok = 0; break;
    }
    phase_count = phase_table->dims[ 0 ];
    gradient_count = phase_count * ( phase_count - 1L );
    gradient_dims[ 0 ] = gradient_count; gradient_dims[ 1 ] = 4L;
    gradient_table = create_double_array( 2L, &gradient_dims[ 0 ], NULL );
    if ( gradient_table == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "allocation of memory for gradient_table failed" );
      ok = 0; break;
    }
    gradient_index = 0L;
    for ( i = 0; i < phase_count; i++ )
    {
      phase_dims[ 0 ] = i; phase_dims[ 1 ] = 0L;
      phase_index = double_array_dims_to_index( phase_table, &phase_dims[ 0 ] );
      if ( phase_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for phase_table failed" );
        ok = 0; break;
      }
      X_i[ 0 ] = phase_table->data[ phase_index++ ];
      X_i[ 1 ] = phase_table->data[ phase_index++ ];
      phase_i = phase_table->data[ phase_index++ ];
      error_i = phase_table->data[ phase_index++ ];
      for ( j = i + 1; j < phase_count; j++ )
      {
        phase_dims[ 0 ] = j; phase_dims[ 1 ] = 0L;
        phase_index = double_array_dims_to_index( phase_table, &phase_dims[ 0 ] );
        if ( phase_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for phase_table failed" );
          ok = 0; break;
        }
        X_j[ 0 ] = phase_table->data[ phase_index++ ];
        X_j[ 1 ] = phase_table->data[ phase_index++ ];
        phase_j = phase_table->data[ phase_index++ ];
        error_j = phase_table->data[ phase_index++ ];
        ra_ij = calculate_angular_separation( &X_i[ 0 ], &X_j[ 0 ] );
        if ( ra_ij == NULL )
        {
          ok = 0; break;
        }
        phase_ij = modulo( ( phase_j - phase_i ) + 180., 360. ) - 180.;
        error_ij = sqrt( error_i * error_i + error_j * error_j );
        gradient_table->data[ gradient_index++ ] = ra_ij[ 0 ];
        gradient_table->data[ gradient_index++ ] = ra_ij[ 1 ];
        gradient_table->data[ gradient_index++ ] = phase_ij;
        gradient_table->data[ gradient_index++ ] = error_ij;
        ra_ji[ 0 ] = ra_ij[ 0 ];
        ra_ji[ 1 ] = modulo( ( ra_ij[ 1 ] + 180. ) + 180., 360. ) - 180.;
        phase_ji = modulo( ( phase_i - phase_j ) + 180., 360. ) - 180.;
        error_ji = error_ij;
        gradient_table->data[ gradient_index++ ] = ra_ji[ 0 ];
        gradient_table->data[ gradient_index++ ] = ra_ji[ 1 ];
        gradient_table->data[ gradient_index++ ] = phase_ji;
        gradient_table->data[ gradient_index++ ] = error_ji;
        free( ra_ij ); ra_ij = NULL;
      }
      if ( !ok ) break;
    }
    if ( !ok ) break;
  }

  // cleanup
  if ( phase_table != NULL )
  {
    free_double_array( phase_table ); phase_table = NULL;
  }
  if ( ra_ij != NULL )
  {
    free( ra_ij ); ra_ij = NULL;
  }
  if ( !ok )
  {
    if ( gradient_table != NULL )
    {
      free_double_array( gradient_table ); gradient_table = NULL;
    }
  }
  return ( gradient_table );
}

//*****************************************************************************

p_double_2 estimate_phase_gradient( p_double_array X_table, p_double_array ref_X_table,
    p_double_array phase_table, p_double_array error_table )
{
  int once, ok = 1;
  long X_dims[ 2 ], ref_X_dims[ 2 ], phase_dims[ 1 ], error_dims[ 1 ], phase_count;
  long local_phase_dims[ 2 ], local_phase_index, phase_index, i_offset, gradient_dims[ 2 ];
  long local_phase_count, i, j, ref_X_index, X_index, error_index, gradient_index;
  long local_gradient_index, gradient_count, xyz_dims[ 2 ], xyz_index, new_gradient_index;
  long rot_xaz_dims[ 2 ], sel_z_count, rot_xaz_index, sel_z_index, sel_count; 
  double ref_X[ 3 ], sum_xz, sum_yz, gradient_angle, gradient_amp, rot_y;
  p_double_2 gradient = NULL;
  p_double_array local_gradient_table = NULL, new_gradient_table = NULL, xyz_table = NULL;
  p_double_array local_phase_table = NULL, gradient_table = NULL, rot_xaz_table = NULL;
  p_double_array sel_z_table = NULL;

  for ( once = 1; once; once-- )
  {
    X_dims[ 0 ] = 0L; X_dims[ 1 ] = 2L;
    if ( !check_double_array( X_table, 2L, &X_dims[ 0 ] ) )
    {
      PyErr_SetString( PyExc_RuntimeError, "X_table has wrong dimensions" );
      ok = 0; break;
    }
    phase_count = X_table->dims[ 0 ];
    ref_X_dims[ 0 ] = phase_count; ref_X_dims[ 1 ] = 2L;
    if ( !check_double_array( ref_X_table, 2L, &ref_X_dims[ 0 ] ) )
    {
      PyErr_SetString( PyExc_RuntimeError, "ref_X_table has wrong dimensions" );
      ok = 0; break;
    }
    phase_dims[ 0 ] = phase_count;
    if ( !check_double_array( phase_table, 1L, &phase_dims[ 0 ] ) )
    {
      PyErr_SetString( PyExc_RuntimeError, "phase_table has wrong dimensions" );
      ok = 0; break;
    }
    error_dims[ 0 ] = phase_count;
    if ( !check_double_array( error_table, 1L, &error_dims[ 0 ] ) )
    {
      PyErr_SetString( PyExc_RuntimeError, "error_table has wrong dimensions" );
      ok = 0; break;
    }
    i_offset = 0L;
    ref_X[ 0 ] = ref_X_table->data[ 0 ];
    ref_X[ 1 ] = ref_X_table->data[ 1 ];
    local_phase_count = 0L;
    for ( j = 0L; j <= phase_count; j++ )
    {
      if ( j < phase_count )
      {
        ref_X_dims[ 0 ] = j; ref_X_dims[ 1 ] = 0L;
        ref_X_index = double_array_dims_to_index( ref_X_table, &ref_X_dims[ 0 ] );
        if ( ref_X_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for ref_X_table failed" );
          ok = 0; break;
        }
      }
      if( ( ref_X_table->data[ ref_X_index++ ] == ref_X[ 0 ] ) && 
          ( ref_X_table->data[ ref_X_index++ ] == ref_X[ 1 ] ) && ( j < phase_count ) )
        local_phase_count++;
      else
      {
        local_phase_dims[ 0 ] = local_phase_count + 1L; local_phase_dims[ 1 ] = 4L;
        local_phase_table = create_double_array( 2L, &local_phase_dims[ 0 ], NULL );
        if ( local_phase_table == NULL )
        {
          PyErr_SetString( PyExc_RuntimeError, "allocation of memory for local_phase_table failed" );
          ok = 0; break;
        }
        X_dims[ 0 ] = i_offset; X_dims[ 1 ] = 0L;
        X_index = double_array_dims_to_index( X_table, &X_dims[ 0 ] );
        if ( X_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for X_table failed" );
          ok = 0; break;
        }
        phase_dims[ 0 ] = i_offset; phase_dims[ 1 ] = 0L;
        phase_index = double_array_dims_to_index( phase_table, &phase_dims[ 0 ] );
        if ( phase_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for phase_table failed" );
          ok = 0; break;
        }
        error_dims[ 0 ] = i_offset; error_dims[ 1 ] = 0L;
        error_index = double_array_dims_to_index( error_table, &error_dims[ 0 ] );
        if ( error_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for error_table failed" );
          ok = 0; break;
        }
        local_phase_index = 0L;
        local_phase_table->data[ local_phase_index++ ] = ref_X[ 0 ];
        local_phase_table->data[ local_phase_index++ ] = ref_X[ 1 ];
        local_phase_table->data[ local_phase_index++ ] = 0.;
        local_phase_table->data[ local_phase_index++ ] = 0.;
        for ( i = 0L; i < local_phase_count; i++ )
        {
          local_phase_table->data[ local_phase_index++ ] = X_table->data[ X_index++ ];
          local_phase_table->data[ local_phase_index++ ] = X_table->data[ X_index++ ];
          local_phase_table->data[ local_phase_index++ ] = phase_table->data[ phase_index++ ];
          local_phase_table->data[ local_phase_index++ ] = error_table->data[ error_index++ ];
        }
        local_gradient_table = calculate_phase_gradients( local_phase_table ); // frees local_phase_table
        local_phase_table = NULL;
        if ( local_gradient_table == NULL )
        {
          ok = 0; break;
        }
        if ( gradient_table == NULL )
        {
          gradient_table = local_gradient_table;
          local_gradient_table = NULL;
        }
        else
        {
          gradient_dims[ 0 ] = gradient_table->dims[ 0 ] + local_gradient_table->dims[ 0 ];
          gradient_dims[ 1 ] = gradient_table->dims[ 1 ];
          new_gradient_table = create_double_array( 2L, &gradient_dims[ 0 ], NULL );
          if ( new_gradient_table == NULL )
          {
            PyErr_SetString( PyExc_RuntimeError, "allocation of memory for new_gradient_table failed" );
            ok = 0; break;
          }
          gradient_index = 0L;
          new_gradient_index = 0L;
          for ( i = 0; i < gradient_table->dims[ 0 ]; i++ )
          {
            new_gradient_table->data[ new_gradient_index++ ] = gradient_table->data[ gradient_index++ ];
            new_gradient_table->data[ new_gradient_index++ ] = gradient_table->data[ gradient_index++ ];
            new_gradient_table->data[ new_gradient_index++ ] = gradient_table->data[ gradient_index++ ];
            new_gradient_table->data[ new_gradient_index++ ] = gradient_table->data[ gradient_index++ ];
          }
          local_gradient_index = 0L;
          for ( i = 0; i < local_gradient_table->dims[ 0 ]; i++ )
          {
            new_gradient_table->data[ new_gradient_index++ ] = local_gradient_table->data[ local_gradient_index++ ];
            new_gradient_table->data[ new_gradient_index++ ] = local_gradient_table->data[ local_gradient_index++ ];
            new_gradient_table->data[ new_gradient_index++ ] = local_gradient_table->data[ local_gradient_index++ ];
            new_gradient_table->data[ new_gradient_index++ ] = local_gradient_table->data[ local_gradient_index++ ];
          }
          free_double_array( local_gradient_table );
          local_gradient_table = NULL;
          free_double_array( gradient_table );
          gradient_table = new_gradient_table;
          new_gradient_table = NULL;
        }

        if ( j < phase_count )
        {
          ref_X_dims[ 0 ] = j; ref_X_dims[ 1 ] = 0L;
          ref_X_index = double_array_dims_to_index( ref_X_table, &ref_X_dims[ 0 ] );
          if ( ref_X_index < 0L )
          {
            PyErr_SetString( PyExc_RuntimeError, "index calculation for ref_X_table failed" );
            ok = 0; break;
          }
          ref_X[ 0 ] = ref_X_table->data[ ref_X_index++ ];
          ref_X[ 1 ] = ref_X_table->data[ ref_X_index++ ];
          i_offset = j;
          local_phase_count = 1L;
        }
      }
    }
    if ( !ok ) break;

    // estimate gradient angle
    gradient_count = gradient_table->dims[ 0 ];
    xyz_dims[ 0 ] = gradient_count; xyz_dims[ 1 ] = 3L;
    xyz_table = create_double_array( 2L, &xyz_dims[ 0 ], NULL );
    if ( xyz_table == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "allocation of memory for xyz_table failed" );
      ok = 0; break;
    }
    xyz_index = 0L;
    sel_count = 0L;
    for ( i = 0L; i < gradient_count; i++ )
    {
      gradient_dims[ 0 ] = i; gradient_dims[ 1 ] = 0L;
      gradient_index = double_array_dims_to_index( gradient_table, &gradient_dims[ 0 ] );
      if ( gradient_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for gradient_table failed" );
        ok = 0; break;
      }
      if ( gradient_table->data[ gradient_index ] > 0. )
      {
        xyz_table->data[ xyz_index++ ] = gradient_table->data[ gradient_index ] * 
            sin( radians( gradient_table->data[ gradient_index + 1L ] ) );
        xyz_table->data[ xyz_index++ ] = gradient_table->data[ gradient_index ] * 
            cos( radians( gradient_table->data[ gradient_index + 1L ] ) );
        xyz_table->data[ xyz_index++ ] = gradient_table->data[ gradient_index + 2L ] /
            gradient_table->data[ gradient_index ];
        sel_count++;
      }
    }
    if ( !ok ) break;
    xyz_table->dims[ 0 ] = sel_count;
    xyz_index = 0L;
    sum_xz = 0.; sum_yz = 0.;
    for ( i = 0L; i < sel_count; i++ )
    {
      xyz_dims[ 0 ] = i; xyz_dims[ 1 ] = 0L;
      xyz_index = double_array_dims_to_index( xyz_table, &xyz_dims[ 0 ] );
      if ( xyz_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for xyz_table failed" );
        ok = 0; break;
      }
      sum_xz = sum_xz + xyz_table->data[ xyz_index ] * xyz_table->data[ xyz_index + 2L ];
      sum_yz = sum_yz + xyz_table->data[ xyz_index + 1L ] * xyz_table->data[ xyz_index + 2L ];
    }
    if ( !ok ) break;
    gradient_angle = degrees( atan2( sum_yz, sum_xz ) );

    // estimate gradient amplitude
    // correct measured phase gradients for baseline angle
    rot_xaz_dims[ 0 ] = xyz_table->dims[ 0 ]; rot_xaz_dims[ 1 ] = 3L;
    rot_xaz_table = create_double_array( 2L, &rot_xaz_dims[ 0 ], NULL );
    if ( rot_xaz_table == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "allocation of memory for rot_xaz_table failed" );
      ok = 0; break;
    }
    sel_z_count = 0L;
    for ( i = 0L; i < xyz_table->dims[ 0 ]; i++ )
    {
      xyz_dims[ 0 ] = i; xyz_dims[ 1 ] = 0L;
      xyz_index = double_array_dims_to_index( xyz_table, &xyz_dims[ 0 ] );
      if ( gradient_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for xyz_table failed" );
        ok = 0; break;
      }
      rot_xaz_dims[ 0 ] = i; rot_xaz_dims[ 1 ] = 0L;
      rot_xaz_index = double_array_dims_to_index( rot_xaz_table, &rot_xaz_dims[ 0 ] );
      if ( rot_xaz_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for rot_xaz_table failed" );
        ok = 0; break;
      }
      rot_xaz_table->data[ rot_xaz_index ] = cos( radians( gradient_angle ) ) * xyz_table->data[ xyz_index ] +
          sin( radians( gradient_angle ) ) * xyz_table->data[ xyz_index + 1L ];
      rot_y = - sin( radians( gradient_angle ) ) * xyz_table->data[ xyz_index ] +
          cos( radians( gradient_angle ) ) * xyz_table->data[ xyz_index + 1L ];
      rot_xaz_table->data[ rot_xaz_index + 1L ] = degrees( atan2( rot_y, rot_xaz_table->data[ rot_xaz_index ] ) );
      rot_xaz_table->data[ rot_xaz_index + 2L ] = xyz_table->data[ xyz_index + 2L ] /
          cos( radians( rot_xaz_table->data[ rot_xaz_index + 1L ] ) );
      if ( ( rot_xaz_table->data[ rot_xaz_index ] > 0. ) && ( rot_xaz_table->data[ rot_xaz_index + 1L ] < 60. ) )
        sel_z_count++;
    }
    if ( !ok ) break;

    // only use positive gradients
    sel_z_table = create_double_array( 1L, &sel_z_count, NULL );
    if ( sel_z_table == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "allocation of memory for sel_z_table failed" );
      ok = 0; break;
    }
    sel_z_index = 0L;
    for ( i = 0L; i < rot_xaz_table->dims[ 0 ]; i++ )
    {
      rot_xaz_dims[ 0 ] = i; rot_xaz_dims[ 1 ] = 0L;
      rot_xaz_index = double_array_dims_to_index( rot_xaz_table, &rot_xaz_dims[ 0 ] );
      if ( rot_xaz_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for rot_xaz_table failed" );
        ok = 0; break;
      }
      if ( ( rot_xaz_table->data[ rot_xaz_index ] > 0. ) && ( rot_xaz_table->data[ rot_xaz_index + 1L ] < 60. ) )
        sel_z_table->data[ sel_z_index++ ] = rot_xaz_table->data[ rot_xaz_index + 2L ];
    }
    if ( !ok ) break;
    gradient_amp = median( sel_z_table ); // this call frees sel_z_table
    sel_z_table = NULL;

    // calculate gradient coefficients
    gradient = ( p_double_2 )malloc( 2 * sizeof( double ) );
    if ( gradient == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "allocation of memory for gradient failed" );
      ok = 0; break;
    }
    gradient[ 0 ] = gradient_amp * cos( radians( gradient_angle ) );
    gradient[ 1 ] = gradient_amp * sin( radians( gradient_angle ) );
  }

  // cleanup
  if ( X_table != NULL )
  {
    free_double_array( X_table ); X_table = NULL;
  }
  if ( ref_X_table != NULL )
  {
    free_double_array( ref_X_table ); ref_X_table = NULL;
  }
  if ( phase_table != NULL )
  {
    free_double_array( phase_table ); phase_table = NULL;
  }
  if ( error_table != NULL )
  {
    free_double_array( error_table ); error_table = NULL;
  }
  if ( local_phase_table != NULL )
  {
    free_double_array( local_phase_table ); local_phase_table = NULL;
  }
  if ( local_gradient_table != NULL )
  {
    free_double_array( local_gradient_table ); local_gradient_table = NULL;
  }
  if ( gradient_table != NULL )
  {
    free_double_array( gradient_table ); gradient_table = NULL;
  }
  if ( new_gradient_table != NULL )
  {
    free_double_array( new_gradient_table ); new_gradient_table = NULL;
  }
  if ( xyz_table != NULL )
  {
    free_double_array( xyz_table ); xyz_table = NULL;
  }
  if ( rot_xaz_table != NULL )
  {
    free_double_array( rot_xaz_table ); rot_xaz_table = NULL;
  }
  if ( sel_z_table != NULL )
  {
    free_double_array( sel_z_table ); sel_z_table = NULL;
  }
  if ( !ok )
  {
    if ( gradient != NULL )
    {
      free( gradient ); gradient = NULL;
    }
  }
  return ( gradient );
}

//*****************************************************************************

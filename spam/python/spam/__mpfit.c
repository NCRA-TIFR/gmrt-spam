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
#include "__mpfit.h"

//*****************************************************************************

double enorm( p_double_array vec, double rgiant, double rdwarf )
{
  int once, ok = 1;
  long vec_dims[ 1 ];
  double agiant, adwarf, mx, mn, ans;
  p_double_array vec2 = NULL;
  for ( once = 1; once; once-- )
  {
    vec_dims[ 0 ] = 0L;
    if ( !check_double_array( vec, 1L, &vec_dims[ 0 ] ) )
    {
      PyErr_SetString( PyExc_ValueError, "vec must be a one-dimensional array of doubles" );
      ok = 0; break;
    }
    agiant = rgiant / ( double )vec->dims[ 0 ];
    adwarf = rdwarf / ( double )vec->dims[ 0 ];

    // This is hopefully a compromise between speed and robustness.
    // Need to do this because of the possibility of over- or underflow.
    vec2 = create_double_array( vec->num_dims, vec->dims, vec->data );
    if ( vec2 == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "allocation of memory for vec2 failed" );
      ok = 0; break;
    }
    mx = maximum( vec2 ); // frees vec2
    vec2 = NULL;
    vec2 = create_double_array( vec->num_dims, vec->dims, vec->data );
    if ( vec2 == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "allocation of memory for vec2 failed" );
      ok = 0; break;
    }
    mn = minimum( vec2 ); // frees vec2
    vec2 = NULL;
    mx = ( fabs( mx ) > fabs( mn ) ? fabs( mx ) : fabs( mn ) );
    if ( mx == 0. )
    {
      ans = 0.;
      break;
    }
    if ( ( mx > agiant ) || ( mx < adwarf ) )
    {
      ans = mx * sqrt( sum( square( multiply( vec, 1. / mx ) ) ) ); // frees vec
      // TODO: proper error handling in statement above
      vec = NULL;
    }
    else
    {
      ans = sqrt( sum( square( vec ) ) ); // frees vec
      vec = NULL;
    }
  }

  // cleanup
  if ( vec != NULL )
  {
    free_double_array( vec ); vec = NULL;
  }
  if ( vec2 != NULL )
  {
    free_double_array( vec2 ); vec2 = NULL;
  }
  if ( !ok )
    ans = - 1.;
  return ( ans );
}

//*****************************************************************************

p_double_array qrfac( p_double_array a, int pivot, double machep, double rgiant, double rdwarf )
// a = m x n matrix
{
  int once, ok = 1;
  long m, n, a_dims[ 2 ], acnorm_dims[ 1 ], temp_dims[ 1 ], i, j, acnorm_index, ipvt_index2;
  long ipvt_dims[ 1 ], ipvt_index, minmn, rdiag_dims[ 1 ], kmax, rdiag_index, temp3;
  long lj, wa_dims[ 1 ], wa_index, ajj_dims[ 1 ], ajj_index, ajk_dims[ 1 ], ajk_index;
  long a_index, temp_index, k, lk, result_dims[ 2 ], result_index;
  double temp2, rmax, ajnorm, temp4;
  p_double_array acnorm = NULL, temp = NULL, rdiag = NULL, wa = NULL, ajj = NULL, ajk = NULL;
  p_double_array result = NULL;
  p_long_array ipvt = NULL;

  for ( once = 1; once; once-- )
  {
    a_dims[ 0 ] = 0L; a_dims[ 1 ] = 0L;
    if ( !check_double_array( a, 2L, &a_dims[ 0 ] ) )
    {
      PyErr_SetString( PyExc_ValueError, "a must be a two-dimensional matrix" );
      ok = 0; break;
    }
    m = a->dims[ 0 ];
    n = a->dims[ 1 ];

    // Compute the initial column norms and initialize arrays
    acnorm_dims[ 0 ] = n;
    acnorm = create_double_array( 1L, &acnorm_dims[ 0 ], NULL );
    if ( acnorm == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "allocation of memory for a failed" );
      ok = 0; break;
    }
    for ( j = 0L; j < n; j++ )
    {
      acnorm_dims[ 0 ] = j;
      acnorm_index = double_array_dims_to_index( acnorm, &acnorm_dims[ 0 ] );
      if ( acnorm_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for acnorm failed" );
        ok = 0; break;
      }
      acnorm->data[ acnorm_index ] = 0.;
    }
    if ( !ok ) break;
    for ( j = 0L; j < n; j++ )
    {
      temp_dims[ 0 ] = m;
      temp = create_double_array( 1L, &temp_dims[ 0 ], NULL );
      if ( temp == NULL )
      {
        PyErr_SetString( PyExc_RuntimeError, "allocation of memory for temp failed 1" );
        ok = 0; break;
      }
      for ( i = 0L; i < m; i++ )
      {
        a_dims[ 0 ] = i; a_dims[ 1 ] = j;
        a_index = double_array_dims_to_index( a, &a_dims[ 0 ] );
        if ( a_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for a failed" );
          ok = 0; break;
        }
        temp_dims[ 0 ] = i;
        temp_index = double_array_dims_to_index( temp, &temp_dims[ 0 ] );
        if ( temp_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for temp failed" );
          ok = 0; break;
        }
        temp->data[ temp_index ] = a->data[ a_index ];
      }
      if ( !ok ) break;
      temp2 = enorm( temp, rgiant, rdwarf ); // frees temp
      temp = NULL;
      if ( temp2 < 0. )
      {
        ok = 0; break; // error message already generated
      }
      acnorm_dims[ 0 ] = j;
      acnorm_index = double_array_dims_to_index( acnorm, &acnorm_dims[ 0 ] );
      if ( acnorm_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for acnorm failed" );
        ok = 0; break;
      }
      acnorm->data[ acnorm_index ] = temp2;
    }
    if ( !ok ) break;
    rdiag = create_double_array( acnorm->num_dims, acnorm->dims, acnorm->data );
    if ( rdiag == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "allocation of memory for rdiag failed" );
      ok = 0; break;
    }
    wa = create_double_array( acnorm->num_dims, acnorm->dims, acnorm->data );
    if ( wa == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "allocation of memory for wa failed" );
      ok = 0; break;
    }
    ipvt_dims[ 0 ] = n;
    ipvt = create_long_array( 1L, &ipvt_dims[ 0 ], NULL );
    if ( ipvt == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "allocation of memory for ipvt_table failed" );
      ok = 0; break;
    }
    for ( i = 0L; i < n; i++ )
    {
      ipvt_dims[ 0 ] = i;
      ipvt_index = long_array_dims_to_index( ipvt, &ipvt_dims[ 0 ] );
      if ( ipvt_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for ipvt failed" );
        ok = 0; break;
      }
      ipvt->data[ ipvt_index ] = i;
    }

    // Reduce a to r with householder transformations
    minmn = ( m < n ? m : n );
    for ( j = 0L; j < minmn; j++ )
    {
      if ( pivot == 1 )
      {
        // Bring the column of largest norm into the pivot position
        for ( i = j; i < n; i++ )
        {
          rdiag_dims[ 0 ] = i;
          rdiag_index = double_array_dims_to_index( rdiag, &rdiag_dims[ 0 ] );
          if ( rdiag_index < 0L )
          {
            PyErr_SetString( PyExc_RuntimeError, "index calculation for rdiag failed" );
            ok = 0; break;
          }
          rmax = 0.;
          if ( ( i == j ) || ( rdiag->data[ rdiag_index ] > rmax ) )
          {
            rmax = rdiag->data[ rdiag_index ];
            kmax = i;
          }

          // Exchange rows via the pivot only.  Avoid actually exchanging
          // the rows, in case there is lots of memory transfer.  The
          // exchange occurs later, within the body of MPFIT, after the
          // extraneous columns of the matrix have been shed.
          if ( kmax != j )
          {
            ipvt_dims[ 0 ] = j;
            ipvt_index = long_array_dims_to_index( ipvt, &ipvt_dims[ 0 ] );
            if ( ipvt_index < 0L )
            {
              PyErr_SetString( PyExc_RuntimeError, "index calculation for ipvt failed" );
              ok = 0; break;
            }
            ipvt_dims[ 0 ] = kmax;
            ipvt_index2 = long_array_dims_to_index( ipvt, &ipvt_dims[ 0 ] );
            if ( ipvt_index2 < 0L )
            {
              PyErr_SetString( PyExc_RuntimeError, "index calculation for ipvt failed" );
              ok = 0; break;
            }
            temp3 = ipvt->data[ ipvt_index ];
            ipvt->data[ ipvt_index ] = ipvt->data[ ipvt_index2 ];
            ipvt->data[ ipvt_index2 ] = temp3;
            rdiag_dims[ 0 ] = j;
            rdiag_index = double_array_dims_to_index( rdiag, &rdiag_dims[ 0 ] );
            if ( rdiag_index < 0L )
            {
              PyErr_SetString( PyExc_RuntimeError, "index calculation for rdiag failed" );
              ok = 0; break;
            }
            temp2 = rdiag->data[ rdiag_index ];
            rdiag_dims[ 0 ] = kmax;
            rdiag_index = double_array_dims_to_index( rdiag, &rdiag_dims[ 0 ] );
            if ( rdiag_index < 0L )
            {
              PyErr_SetString( PyExc_RuntimeError, "index calculation for rdiag failed" );
              ok = 0; break;
            }
            rdiag->data[ rdiag_index ] = temp2;
            wa_dims[ 0 ] = j;
            wa_index = double_array_dims_to_index( wa, &wa_dims[ 0 ] );
            if ( wa_index < 0L )
            {
              PyErr_SetString( PyExc_RuntimeError, "index calculation for wa failed" );
              ok = 0; break;
            }
            temp2 = wa->data[ wa_index ];
            wa_dims[ 0 ] = kmax;
            wa_index = double_array_dims_to_index( wa, &wa_dims[ 0 ] );
            if ( wa_index < 0L )
            {
              PyErr_SetString( PyExc_RuntimeError, "index calculation for wa failed" );
              ok = 0; break;
            }
            wa->data[ wa_index ] = temp2;
          }
        }
        if ( !ok ) break;
      }

      // Compute the householder transformation to reduce the jth
      // column of A to a multiple of the jth unit vector
      ipvt_dims[ 0 ] = j;
      ipvt_index = long_array_dims_to_index( ipvt, &ipvt_dims[ 0 ] );
      if ( ipvt_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for ipvt failed" );
        ok = 0; break;
      }
      lj = ipvt->data[ ipvt_index ];
      ajj_dims[ 0 ] = m - j;
      ajj = create_double_array( 1L, &ajj_dims[ 0 ], NULL );
      if ( ajj == NULL )
      {
        PyErr_SetString( PyExc_RuntimeError, "allocation of memory for ajj failed" );
        ok = 0; break;
      }
      for ( i = 0L; i < m - j; i++ )
      {
        a_dims[ 0 ] = j + i; a_dims[ 1 ] = lj;
        a_index = double_array_dims_to_index( a, &a_dims[ 0 ] );
        if ( a_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for a failed" );
          ok = 0; break;
        }
        ajj_dims[ 0 ] = i;
        ajj_index = double_array_dims_to_index( ajj, &ajj_dims[ 0 ] );
        if ( ajj_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for ajj failed" );
          ok = 0; break;
        }
        ajj->data[ ajj_index ] = a->data[ a_index ];
      }
      if ( !ok ) break;
      temp = create_double_array( ajj->num_dims, ajj->dims, ajj->data );
      if ( temp == NULL )
      {
        PyErr_SetString( PyExc_RuntimeError, "allocation of memory for temp failed 2" );
        ok = 0; break;
      }
      ajnorm = enorm( temp, rgiant, rdwarf ); // frees temp
      temp = NULL;
      if ( ajnorm < 0. )
      {
        ok = 0; break; // error message already generated
      }
      if ( ajnorm == 0. )
        break;
      a_dims[ 0 ] = j; a_dims[ 1 ] = lj;
      a_index = double_array_dims_to_index( a, &a_dims[ 0 ] );
      if ( a_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for a failed" );
        ok = 0; break;
      }
      if ( a->data[ a_index ] < 0. )
        ajnorm = - ajnorm;
      ajj = multiply( ajj, 1. / ajnorm ); // both frees and allocates ajj
      if ( ajj == NULL )
      {
        ok = 0; break; // error message is already generated
      }
      ajj_dims[ 0 ] = 0L;
      ajj_index = double_array_dims_to_index( ajj, &ajj_dims[ 0 ] );
      if ( ajj_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for ajj failed" );
        ok = 0; break;
      }
      ajj->data[ ajj_index ] = ajj->data[ ajj_index ] + 1.;
      for ( i = 0L; i < m - j; i++ )
      {
        a_dims[ 0 ] = j + i; a_dims[ 1 ] = lj;
        a_index = double_array_dims_to_index( a, &a_dims[ 0 ] );
        if ( a_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for a failed" );
          ok = 0; break;
        }
        ajj_dims[ 0 ] = i;
        ajj_index = double_array_dims_to_index( ajj, &ajj_dims[ 0 ] );
        if ( ajj_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for ajj failed" );
          ok = 0; break;
        }
        a->data[ a_index ] = ajj->data[ ajj_index ];
      }
      if ( !ok ) break;

      // Apply the transformation to the remaining columns
      // and update the norms
      if ( j + 1L < n )
      {
        for ( k = j + 1L; k < n; k++ )
        {
          ipvt_dims[ 0 ] = k;
          ipvt_index = long_array_dims_to_index( ipvt, &ipvt_dims[ 0 ] );
          if ( ipvt_index < 0L )
          {
            PyErr_SetString( PyExc_RuntimeError, "index calculation for ipvt failed" );
            ok = 0; break;
          }
          lk = ipvt->data[ ipvt_index ];

          ajk_dims[ 0 ] = m - j;
          ajk = create_double_array( 1L, &ajk_dims[ 0 ], NULL );
          if ( ajk == NULL )
          {
            PyErr_SetString( PyExc_RuntimeError, "allocation of memory for ajk failed" );
            ok = 0; break;
          }
          for ( i = 0L; i < m - j; i++ )
          {
            a_dims[ 0 ] = j + i; a_dims[ 1 ] = lk;
            a_index = double_array_dims_to_index( a, &a_dims[ 0 ] );
            if ( a_index < 0L )
            {
              PyErr_SetString( PyExc_RuntimeError, "index calculation for a failed" );
              ok = 0; break;
            }
            ajk_dims[ 0 ] = i;
            ajk_index = double_array_dims_to_index( ajk, &ajk_dims[ 0 ] );
            if ( ajk_index < 0L )
            {
              PyErr_SetString( PyExc_RuntimeError, "index calculation for ajk failed" );
              ok = 0; break;
            }
            ajk->data[ ajk_index ] = a->data[ a_index ];
          }
          if ( !ok ) break;
          a_dims[ 0 ] = j; a_dims[ 1 ] = lj;
          a_index = double_array_dims_to_index( a, &a_dims[ 0 ] );
          if ( a_index < 0L )
          {
            PyErr_SetString( PyExc_RuntimeError, "index calculation for a failed" );
            ok = 0; break;
          }
          if ( a->data[ a_index ] != 0. )
          {
            temp = create_double_array( ajk->num_dims, ajk->dims, ajk->data );
            if ( temp == NULL )
            {
              PyErr_SetString( PyExc_RuntimeError, "allocation of memory for temp failed 3" );
              ok = 0; break;
            }
            for ( i = 0L; i < m - j; i++ )
            {
              temp_dims[ 0 ] = i;
              temp_index = double_array_dims_to_index( temp, &temp_dims[ 0 ] );
              if ( temp_index < 0L )
              {
                PyErr_SetString( PyExc_RuntimeError, "index calculation for temp failed" );
                ok = 0; break;
              }
              ajj_dims[ 0 ] = i;
              ajj_index = double_array_dims_to_index( ajj, &ajj_dims[ 0 ] );
              if ( ajj_index < 0L )
              {
                PyErr_SetString( PyExc_RuntimeError, "index calculation for ajj failed" );
                ok = 0; break;
              }
              temp->data[ temp_index ] = temp->data[ temp_index ] * ajj->data[ ajj_index ];
            }
            if ( !ok ) break;
            temp2 = sum( temp ); // frees temp
            temp = NULL;
            a_dims[ 0 ] = j; a_dims[ 1 ] = lj;
            a_index = double_array_dims_to_index( a, &a_dims[ 0 ] );
            if ( a_index < 0L )
            {
              PyErr_SetString( PyExc_RuntimeError, "index calculation for a failed" );
              ok = 0; break;
            }
            temp2 = temp2 / a->data[ a_index ];
            for ( i = 0L; i < m - j; i++ )
            {
              ajk_dims[ 0 ] = i;
              ajk_index = double_array_dims_to_index( ajk, &ajk_dims[ 0 ] );
              if ( ajk_index < 0L )
              {
                PyErr_SetString( PyExc_RuntimeError, "index calculation for ajk failed" );
                ok = 0; break;
              }
              ajj_dims[ 0 ] = i;
              ajj_index = double_array_dims_to_index( ajj, &ajj_dims[ 0 ] );
              if ( ajj_index < 0L )
              {
                PyErr_SetString( PyExc_RuntimeError, "index calculation for ajj failed" );
                ok = 0; break;
              }
              a_dims[ 0 ] = j + i; a_dims[ 1 ] = lk;
              a_index = double_array_dims_to_index( a, &a_dims[ 0 ] );
              if ( a_index < 0L )
              {
                PyErr_SetString( PyExc_RuntimeError, "index calculation for a failed" );
                ok = 0; break;
              }
              a->data[ a_index ] = ajk->data[ ajk_index ] - ajj->data[ ajj_index ] * temp2;
            }
            if ( !ok ) break;
            rdiag_dims[ 0 ] = k;
            rdiag_index = double_array_dims_to_index( rdiag, &rdiag_dims[ 0 ] );
            if ( rdiag_index < 0L )
            {
              PyErr_SetString( PyExc_RuntimeError, "index calculation for rdiag failed" );
              ok = 0; break;
            }
            if ( ( pivot == 1 ) && ( rdiag->data[ rdiag_index ] != 0. ) )
            {
              a_dims[ 0 ] = j; a_dims[ 1 ] = lk;
              a_index = double_array_dims_to_index( a, &a_dims[ 0 ] );
              if ( a_index < 0L )
              {
                PyErr_SetString( PyExc_RuntimeError, "index calculation for a failed" );
                ok = 0; break;
              }
              temp2 = a->data[ a_index ] / rdiag->data[ rdiag_index ];
              temp4 = 1. - temp2 * temp2;
              temp2 = ( temp4 > 0. ? temp4 : 0. );
              rdiag->data[ rdiag_index ] = rdiag->data[ rdiag_index ] * sqrt( temp2 );
              wa_dims[ 0 ] = k;
              wa_index = double_array_dims_to_index( wa, &wa_dims[ 0 ] );
              if ( wa_index < 0L )
              {
                PyErr_SetString( PyExc_RuntimeError, "index calculation for wa failed" );
                ok = 0; break;
              }
              temp2 = rdiag->data[ rdiag_index ] / wa->data[ wa_index ];
              if ( 0.05 * temp2 * temp2 <= machep )
              {
                temp_dims[ 0 ] = m - ( j + 1L );
                temp = create_double_array( 1L, &temp_dims[ 0 ], NULL );
                if ( temp == NULL )
                {
                  PyErr_SetString( PyExc_RuntimeError, "allocation of memory for temp failed 4" );
                  ok = 0; break;
                }
                for ( i = 0L; i < m - ( j + 1L ); i++ )
                {
                  a_dims[ 0 ] = ( j + 1L ) + i; a_dims[ 1 ] = lk;
                  a_index = double_array_dims_to_index( a, &a_dims[ 0 ] );
                  if ( a_index < 0L )
                  {
                    PyErr_SetString( PyExc_RuntimeError, "index calculation for a failed" );
                    ok = 0; break;
                  }
                  temp_dims[ 0 ] = i;
                  temp_index = double_array_dims_to_index( temp, &temp_dims[ 0 ] );
                  if ( temp_index < 0L )
                  {
                    PyErr_SetString( PyExc_RuntimeError, "index calculation for temp failed" );
                    ok = 0; break;
                  }
                  temp->data[ temp_index ] = a->data[ a_index ];
                }
                if ( !ok ) break;
                temp2 = enorm( temp, rgiant, rdwarf ); // frees temp
                temp = NULL;
                if ( temp2 < 0. )
                {
                  ok = 0; break; // error message already generated
                }
                rdiag->data[ rdiag_index ] = temp2;
                wa->data[ wa_index ] = temp2;
              }
            }
          }
        }
        if ( !ok ) break;
      }
      rdiag_dims[ 0 ] = j;
      rdiag_index = double_array_dims_to_index( rdiag, &rdiag_dims[ 0 ] );
      if ( rdiag_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for rdiag failed" );
        ok = 0; break;
      }
      rdiag->data[ rdiag_index ] = - ajnorm;

      free_double_array( ajj ); ajj = NULL;
      free_double_array( ajk ); ajk = NULL;
    }
    if ( !ok ) break;

    // combine a, ipvt, rdiag and acnorm into one array for return
    result_dims[ 0 ] = m + 3L; result_dims[ 1 ] = n;
    result = create_double_array( 2L, &result_dims[ 0 ], NULL );
    if ( result == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "allocation of memory for result failed" );
      ok = 0; break;
    }
    for ( j = 0L; j < m + 3L; j++ )
    {
      for ( k = 0L; k < n; k++ )
      {
        result_dims[ 0 ] = j; result_dims[ 1 ] = k;
        result_index = double_array_dims_to_index( result, &result_dims[ 0 ] );
        if ( result_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for result failed" );
          ok = 0; break;
        }
        if ( j < m )
        {
          a_dims[ 0 ] = j; a_dims[ 1 ] = k;
          a_index = double_array_dims_to_index( a, &a_dims[ 0 ] );
          if ( a_index < 0L )
          {
            PyErr_SetString( PyExc_RuntimeError, "index calculation for a failed" );
            ok = 0; break;
          }
          result->data[ result_index ] = a->data[ a_index ];
        }
        else if ( j == m )
        {
          ipvt_dims[ 0 ] = k;
          ipvt_index = long_array_dims_to_index( ipvt, &ipvt_dims[ 0 ] );
          if ( ipvt_index < 0L )
          {
            PyErr_SetString( PyExc_RuntimeError, "index calculation for ipvt failed" );
            ok = 0; break;
          }
          result->data[ result_index ] = ( double )ipvt->data[ ipvt_index ];
        }
        else if ( j == m + 1L )
        {
          rdiag_dims[ 0 ] = k;
          rdiag_index = double_array_dims_to_index( rdiag, &rdiag_dims[ 0 ] );
          if ( rdiag_index < 0L )
          {
            PyErr_SetString( PyExc_RuntimeError, "index calculation for rdiag failed" );
            ok = 0; break;
          }
          result->data[ result_index ] = rdiag->data[ rdiag_index ];
        }
        else if ( j == m + 2L )
        {
          acnorm_dims[ 0 ] = k;
          acnorm_index = double_array_dims_to_index( acnorm, &acnorm_dims[ 0 ] );
          if ( acnorm_index < 0L )
          {
            PyErr_SetString( PyExc_RuntimeError, "index calculation for acnorm failed" );
            ok = 0; break;
          }
          result->data[ result_index ] = acnorm->data[ acnorm_index ];
        }
      }
      if ( !ok ) break;
    }
    if ( !ok ) break;
  }

  // cleanup
  if ( a != NULL )
  {
    free_double_array( a ); a = NULL;
  }
  if ( acnorm != NULL )
  {
    free_double_array( acnorm ); acnorm = NULL;
  }
  if ( temp != NULL )
  {
    free_double_array( temp ); temp = NULL;
  }
  if ( rdiag != NULL )
  {
    free_double_array( rdiag ); rdiag = NULL;
  }
  if ( ipvt != NULL )
  {
    free_long_array( ipvt ); ipvt = NULL;
  }
  if ( wa != NULL )
  {
    free_double_array( wa ); wa = NULL;
  }
  if ( ajj != NULL )
  {
    free_double_array( ajj ); ajj = NULL;
  }
  if ( ajk != NULL )
  {
    free_double_array( ajk ); ajk = NULL;
  }
  if ( !ok )
  {
    if ( result != NULL )
    {
      free_double_array( result ); result = NULL;
    }
  }
  return ( result );
}
// result = ( m + 3 ) x n matrix, consisting of:
//   a = m x n matrix
//   ipvt = n vector (int!!)
//   rdiag = n vector
//   acnorm = n vector

//*****************************************************************************

p_double_array qrsolv( p_double_array r, p_long_array ipvt, p_double_array diag,
    p_double_array qtb, p_double_array sdiag )
// r = n x n double matrix    
// ipvt = n long vector
// diag = n double vector
// qtb = n double vector
// sdiag = n double vector
{
  int once, ok = 1;
  long m, n, r_dims[ 2 ], r_index, ipvt_dims[ 1 ], ipvt_index, diag_dims[ 1 ], diag_index;
  long qtb_dims[ 1 ], sdiag_dims[ 1 ], sdiag_index, x_dims[ 1 ], x_index;
  long temp2_dims[ 1 ], temp2_index, i, j, k, nsing, l, wa_dims[ 1 ], wa_index;
  long result_dims[ 2 ], result_index;
  double temp, summ, qtbpj, tang, cotan, sine, cosine;
  p_double_array x = NULL, temp2 = NULL, wa = NULL, result = NULL;
  for ( once = 1; once; once-- )
  {
    r_dims[ 0 ] = 0L; r_dims[ 1 ] = 0L;
    if ( !check_double_array( r, 2L, &r_dims[ 0 ] ) )
    {
      PyErr_SetString( PyExc_ValueError, "r must be a two-dimensional matrix" );
      ok = 0; break;
    }
    m = r->dims[ 0 ];
    n = r->dims[ 1 ];
    if ( m != n )
    {
      PyErr_SetString( PyExc_ValueError, "r must be an n x n two-dimensional matrix" );
      ok = 0; break;
    }
    ipvt_dims[ 0 ] = n;
    if ( !check_long_array( ipvt, 1L, &ipvt_dims[ 0 ] ) )
    {
      PyErr_SetString( PyExc_ValueError, "ipvt must be a one-dimensional vector of length n" );
      ok = 0; break;
    }
    diag_dims[ 0 ] = n;
    if ( !check_double_array( diag, 1L, &diag_dims[ 0 ] ) )
    {
      PyErr_SetString( PyExc_ValueError, "diag must be a one-dimensional vector of length n" );
      ok = 0; break;
    }
    qtb_dims[ 0 ] = n;
    if ( !check_double_array( qtb, 1L, &qtb_dims[ 0 ] ) )
    {
      PyErr_SetString( PyExc_ValueError, "qtb must be a one-dimensional vector of length n" );
      ok = 0; break;
    }
    sdiag_dims[ 0 ] = n;
    if ( !check_double_array( sdiag, 1L, &sdiag_dims[ 0 ] ) )
    {
      PyErr_SetString( PyExc_ValueError, "sdiag must be a one-dimensional vector of length n" );
      ok = 0; break;
    }

    // copy r and (q transpose)*b to preserve input and initialize s.
    // in particular, save the diagonal elements of r in x.
    for ( j = 0L; j < n; j++ )
    {
      for ( i = j; i < n; i++ )
      {
        r_dims[ 0 ] = j; r_dims[ 1 ] = i;
        r_index = double_array_dims_to_index( r, &r_dims[ 0 ] );
        if ( r_index < 0L )
        {
          PyErr_SetString( PyExc_ValueError, "error while calculating index for r" );
          ok = 0; break;
        }
        temp = r->data[ r_index ];
        r_dims[ 0 ] = i; r_dims[ 1 ] = j;
        r_index = double_array_dims_to_index( r, &r_dims[ 0 ] );
        if ( r_index < 0L )
        {
          PyErr_SetString( PyExc_ValueError, "error while calculating index for r" );
          ok = 0; break;
        }
        r->data[ r_index ] = temp;
      }
      if ( !ok ) break;
    }
    if ( !ok ) break;
    x = create_double_array( 1L, &n, NULL );
    if ( x == NULL )
    {
      PyErr_SetString( PyExc_ValueError, "allocation of memory for x failed" );
      ok = 0; break;
    }
    for ( i = 0L; i < n; i++ )
    {
      r_dims[ 0 ] = i; r_dims[ 1 ] = i;
      r_index = double_array_dims_to_index( r, &r_dims[ 0 ] );
      if ( r_index < 0L )
      {
        PyErr_SetString( PyExc_ValueError, "error while calculating index for r" );
        ok = 0; break;
      }
      x_dims[ 0 ] = i;
      x_index = double_array_dims_to_index( x, &x_dims[ 0 ] );
      if ( x_index < 0L )
      {
        PyErr_SetString( PyExc_ValueError, "error while calculating index for x" );
        ok = 0; break;
      }
      x->data[ x_index ] = r->data[ r_index ];
    }
    if ( !ok ) break;
    wa = create_double_array( qtb->num_dims, qtb->dims, qtb->data );
    if ( wa == NULL )
    {
      PyErr_SetString( PyExc_ValueError, "allocation of memory for wa failed" );
      ok = 0; break;
    }

    // Eliminate the diagonal matrix d using a given rotation
    for ( j = 0L; j < n; j++ )
    {
      ipvt_dims[ 0 ] = j;
      ipvt_index = long_array_dims_to_index( ipvt, &ipvt_dims[ 0 ] );
      if ( ipvt_index < 0L )
      {
        PyErr_SetString( PyExc_ValueError, "error while calculating index for ipvt" );
        ok = 0; break;
      }
      l = ipvt->data[ ipvt_index ];
      diag_dims[ 0 ] = l;
      diag_index = double_array_dims_to_index( diag, &diag_dims[ 0 ] );
      if ( diag_index < 0L )
      {
        PyErr_SetString( PyExc_ValueError, "error while calculating index for diag" );
        ok = 0; break;
      }
      if ( diag->data[ diag_index ] == 0. )
        break;
      for ( i = j; i < n; i++ )
      {
        sdiag_dims[ 0 ] = i;
        sdiag_index = double_array_dims_to_index( sdiag, &sdiag_dims[ 0 ] );
        if ( sdiag_index < 0L )
        {
          PyErr_SetString( PyExc_ValueError, "error while calculating index for sdiag" );
          ok = 0; break;
        }
        sdiag->data[ sdiag_index ] = 0.;
      }
      if ( !ok ) break;
      sdiag_dims[ 0 ] = j;
      sdiag_index = double_array_dims_to_index( sdiag, &sdiag_dims[ 0 ] );
      if ( sdiag_index < 0L )
      {
        PyErr_SetString( PyExc_ValueError, "error while calculating index for sdiag" );
        ok = 0; break;
      }
      sdiag->data[ sdiag_index ] = diag->data[ diag_index ];

      // The transformations to eliminate the row of d modify only a
      // single element of (q transpose)*b beyond the first n, which
      // is initially zero.
      qtbpj = 0.;
      for ( k = j; k < n; k++ )
      {
        sdiag_dims[ 0 ] = k;
        sdiag_index = double_array_dims_to_index( sdiag, &sdiag_dims[ 0 ] );
        if ( sdiag_index < 0L )
        {
          PyErr_SetString( PyExc_ValueError, "error while calculating index for sdiag" );
          ok = 0; break;
        }
        if ( sdiag->data[ sdiag_index ] == 0. )
          break;
        r_dims[ 0 ] = k; r_dims[ 1 ] = k;
        r_index = double_array_dims_to_index( r, &r_dims[ 0 ] );
        if ( r_index < 0L )
        {
          PyErr_SetString( PyExc_ValueError, "error while calculating index for r" );
          ok = 0; break;
        }
        if ( fabs( r->data[ r_index ] ) < fabs( sdiag->data[ sdiag_index ] ) )
        {
          cotan = r->data[ r_index ] / sdiag->data[ sdiag_index ];
          sine = 0.5 / sqrt( 0.25 + 0.25 * cotan * cotan );
          cosine = sine * cotan;
        }
        else
        {
          tang = sdiag->data[ sdiag_index ] / r->data[ r_index ];
          cosine = 0.5 / sqrt( 0.25 + 0.25 * tang * tang );
          sine = cosine * tang;
        }

        // Compute the modified diagonal element of r and the
        // modified element of ((q transpose)*b,0).
        r->data[ r_index ] = cosine * r->data[ r_index ] + sine * sdiag->data[ sdiag_index ];
        wa_dims[ 0 ] = k;
        wa_index = double_array_dims_to_index( wa, &wa_dims[ 0 ] );
        if ( wa_index < 0L )
        {
          PyErr_SetString( PyExc_ValueError, "error while calculating index for wa" );
          ok = 0; break;
        }
        temp = cosine * wa->data[ wa_index ] + sine * qtbpj;
        qtbpj = - sine * wa->data[ wa_index ] + cosine * qtbpj;
        wa->data[ wa_index ] = temp;

        // Accumulate the transformation in the row of s
        if ( n > k + 1L )
        {
          for ( i = 0L; i < n - ( k + 1L ); i++ )
          {
            r_dims[ 0 ] = ( k + 1L ) + i; r_dims[ 1 ] = k;
            r_index = double_array_dims_to_index( r, &r_dims[ 0 ] );
            if ( r_index < 0L )
            {
              PyErr_SetString( PyExc_ValueError, "error while calculating index for r" );
              ok = 0; break;
            }
            sdiag_dims[ 0 ] = ( k + 1L ) + i;
            sdiag_index = double_array_dims_to_index( sdiag, &sdiag_dims[ 0 ] );
            if ( sdiag_index < 0L )
            {
              PyErr_SetString( PyExc_ValueError, "error while calculating index for sdiag" );
              ok = 0; break;
            }
            temp = cosine * r->data[ r_index ] + sine * sdiag->data[ sdiag_index ];
            sdiag->data[ sdiag_index ] = - sine * r->data[ r_index ] + cosine * sdiag->data[ sdiag_index ];
            r->data[ r_index ] = temp;
          }
          if ( !ok ) break;
        }
      }
      if ( !ok ) break;
      r_dims[ 0 ] = j; r_dims[ 1 ] = j;
      r_index = double_array_dims_to_index( r, &r_dims[ 0 ] );
      if ( r_index < 0L )
      {
        PyErr_SetString( PyExc_ValueError, "error while calculating index for r" );
        ok = 0; break;
      }
      sdiag_dims[ 0 ] = j;
      sdiag_index = double_array_dims_to_index( sdiag, &sdiag_dims[ 0 ] );
      if ( sdiag_index < 0L )
      {
        PyErr_SetString( PyExc_ValueError, "error while calculating index for sdiag" );
        ok = 0; break;
      }
      sdiag->data[ sdiag_index ] = r->data[ r_index ];
      x_dims[ 0 ] = j;
      x_index = double_array_dims_to_index( x, &x_dims[ 0 ] );
      if ( x_index < 0L )
      {
        PyErr_SetString( PyExc_ValueError, "error while calculating index for x" );
        ok = 0; break;
      }
      r->data[ r_index ] = x->data[ x_index ];
    }
    if ( !ok ) break;

    // Solve the triangular system for z.  If the system is singular
    // then obtain a least squares solution
    nsing = n;
    for ( i = 0L; i < n; i++ )
    {
      sdiag_dims[ 0 ] = i;
      sdiag_index = double_array_dims_to_index( sdiag, &sdiag_dims[ 0 ] );
      if ( sdiag_index < 0L )
      {
        PyErr_SetString( PyExc_ValueError, "error while calculating index for sdiag" );
        ok = 0; break;
      }
      if ( sdiag->data[ sdiag_index ] == 0. )
        break;
    }
    if ( !ok ) break;
    if ( i < n )
    {
      nsing = i;
      for ( i = nsing; i < n; i++ )
      {
        wa_dims[ 0 ] = i;
        wa_index = double_array_dims_to_index( wa, &wa_dims[ 0 ] );
        if ( wa_index < 0L )
        {
          PyErr_SetString( PyExc_ValueError, "error while calculating index for wa" );
          ok = 0; break;
        }
        wa->data[ wa_index ] = 0.;
      }
      if ( !ok ) break;
    }

    if ( nsing >= 1L )
    {
      sdiag_dims[ 0 ] = nsing - 1L;
      sdiag_index = double_array_dims_to_index( sdiag, &sdiag_dims[ 0 ] );
      if ( sdiag_index < 0L )
      {
        PyErr_SetString( PyExc_ValueError, "error while calculating index for sdiag" );
        ok = 0; break;
      }
      wa_dims[ 0 ] = nsing - 1L;
      wa_index = double_array_dims_to_index( wa, &wa_dims[ 0 ] );
      if ( wa_index < 0L )
      {
        PyErr_SetString( PyExc_ValueError, "error while calculating index for wa" );
        ok = 0; break;
      }
      wa->data[ wa_index ] = wa->data[ wa_index ] / sdiag->data[ sdiag_index ]; // Degenerate case
      for ( j = nsing - 2L; j >= 0L; j-- )
      {
        temp2_dims[ 0 ] = nsing - ( j + 1L );
        temp2 = create_double_array( 1L, &temp2_dims[ 0 ], NULL );
        if ( temp2 == NULL )
        {
          PyErr_SetString( PyExc_ValueError, "allocation of memory for temp2 failed" );
          ok = 0; break;
        }
        for ( i = 0L; i < nsing - ( j + 1L ); i++ )
        {
          r_dims[ 0 ] = ( j + 1L ) + i; r_dims[ 1 ] = j;
          r_index = double_array_dims_to_index( r, &r_dims[ 0 ] );
          if ( r_index < 0L )
          {
            PyErr_SetString( PyExc_ValueError, "error while calculating index for r" );
            ok = 0; break;
          }
          wa_dims[ 0 ] = ( j + 1L ) + i;
          wa_index = double_array_dims_to_index( wa, &wa_dims[ 0 ] );
          if ( wa_index < 0L )
          {
            PyErr_SetString( PyExc_ValueError, "error while calculating index for wa" );
            ok = 0; break;
          }
          temp2_dims[ 0 ] = i;
          temp2_index = double_array_dims_to_index( temp2, &temp2_dims[ 0 ] );
          if ( temp2_index < 0L )
          {
            PyErr_SetString( PyExc_ValueError, "error while calculating index for temp2" );
            ok = 0; break;
          }
          temp2->data[ temp2_index ] = r->data[ r_index ] * wa->data[ wa_index ];
        }
        if ( !ok ) break;
        summ = sum( temp2 ); // frees temp2
        temp2 = NULL;
        sdiag_dims[ 0 ] = j;
        sdiag_index = double_array_dims_to_index( sdiag, &sdiag_dims[ 0 ] );
        if ( sdiag_index < 0L )
        {
          PyErr_SetString( PyExc_ValueError, "error while calculating index for sdiag" );
          ok = 0; break;
        }
        wa_dims[ 0 ] = j;
        wa_index = double_array_dims_to_index( wa, &wa_dims[ 0 ] );
        if ( wa_index < 0L )
        {
          PyErr_SetString( PyExc_ValueError, "error while calculating index for wa" );
          ok = 0; break;
        }
        wa->data[ wa_index ] = ( wa->data[ wa_index ] - summ ) / sdiag->data[ sdiag_index ];
      }
      if ( !ok ) break;
    }

    // Permute the components of z back to components of x
    for ( i = 0L; i < n; i++ )
    {
      wa_dims[ 0 ] = i;
      wa_index = double_array_dims_to_index( wa, &wa_dims[ 0 ] );
      if ( wa_index < 0L )
      {
        PyErr_SetString( PyExc_ValueError, "error while calculating index for wa" );
        ok = 0; break;
      }
      ipvt_dims[ 0 ] = i;
      ipvt_index = long_array_dims_to_index( ipvt, &ipvt_dims[ 0 ] );
      if ( ipvt_index < 0L )
      {
        PyErr_SetString( PyExc_ValueError, "error while calculating index for ipvt" );
        ok = 0; break;
      }
      x_dims[ 0 ] = ipvt->data[ ipvt_index ];
      x_index = double_array_dims_to_index( x, &x_dims[ 0 ] );
      if ( x_index < 0L )
      {
        PyErr_SetString( PyExc_ValueError, "error while calculating index for x" );
        ok = 0; break;
      }
      x->data[ x_index ] = wa->data[ wa_index ];
    }
    if ( !ok ) break;

    // build result to be returned
    result_dims[ 0 ] = n + 2L; result_dims[ 1 ] = n;
    result = create_double_array( 2L, &result_dims[ 0 ], NULL );
    if ( result == NULL )
    {
      PyErr_SetString( PyExc_ValueError, "allocation of memory for result failed" );
      ok = 0; break;
    }
    for ( j = 0L; j < n + 2L; j++ )
    {
      for ( k = 0L; k < n; k++ )
      {
        result_dims[ 0 ] = j; result_dims[ 1 ] = k;
        result_index = double_array_dims_to_index( result, &result_dims[ 0 ] );
        if ( result_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for result failed" );
          ok = 0; break;
        }
        if ( j < n )
        {
          r_dims[ 0 ] = j; r_dims[ 1 ] = k;
          r_index = double_array_dims_to_index( r, &r_dims[ 0 ] );
          if ( r_index < 0L )
          {
            PyErr_SetString( PyExc_RuntimeError, "index calculation for r failed" );
            ok = 0; break;
          }
          result->data[ result_index ] = r->data[ r_index ];
        }
        else if ( j == n )
        {
          x_dims[ 0 ] = k;
          x_index = double_array_dims_to_index( x, &x_dims[ 0 ] );
          if ( x_index < 0L )
          {
            PyErr_SetString( PyExc_RuntimeError, "index calculation for x failed" );
            ok = 0; break;
          }
          result->data[ result_index ] = x->data[ x_index ];
        }
        else if ( j == n + 1L )
        {
          sdiag_dims[ 0 ] = k;
          sdiag_index = double_array_dims_to_index( sdiag, &sdiag_dims[ 0 ] );
          if ( sdiag_index < 0L )
          {
            PyErr_SetString( PyExc_RuntimeError, "index calculation for sdiag failed" );
            ok = 0; break;
          }
          result->data[ result_index ] = sdiag->data[ sdiag_index ];
        }
      }
      if ( !ok ) break;
    }
    if ( !ok ) break;
  }

  // cleanup
  if ( r != NULL )
  {
    free_double_array( r ); r = NULL;
  }
  if ( ipvt != NULL )
  {
    free_long_array( ipvt ); ipvt = NULL;
  }
  if ( diag != NULL )
  {
    free_double_array( diag ); diag = NULL;
  }
  if ( qtb != NULL )
  {
    free_double_array( qtb ); qtb = NULL;
  }
  if ( sdiag != NULL )
  {
    free_double_array( sdiag ); sdiag = NULL;
  }
  if ( x != NULL )
  {
    free_double_array( x ); x = NULL;
  }
  if ( wa != NULL )
  {
    free_double_array( wa ); wa = NULL;
  }
  if ( temp2 != NULL )
  {
    free_double_array( temp2 ); temp2 = NULL;
  }
  if ( !ok )
  {
    if ( result != NULL )
    {
      free_double_array( result ); result = NULL;
    }
  }
  return ( result );
}
// result is ( n + 2 ) x n matrix, consisting of:
// r = n x n double matrix
// x = n double vector
// sdiag = n double vector

//*****************************************************************************

p_double_array lmpar( p_double_array r, p_long_array ipvt, p_double_array diag,
    p_double_array qtb, p_double_array x, p_double_array sdiag,
    double delta, double par, double dwarf, double rgiant, double rdwarf )
// r = n x n double matrix
// ipvt = n long array
// diag = n double array
// qtb = n double array
// x = n double array
// sdiag = n double array
{
  int once, ok = 1;
  long i, j, k, r_dims[ 2 ], r_index, ipvt_dims[ 1 ], ipvt_index, diag_dims[ 1 ], diag_index;
  long m, n, qtb_dims[ 1 ], qtb_index, x_dims[ 1 ], x_index, sdiag_dims[ 1 ], sdiag_index;
  long result_dims[ 2 ], result_index, nsing, wa1_dims[ 1 ], wa1_index, wa2_dims[ 1 ], wa2_index;
  long iter, temp2_dims[ 1 ], temp2_index;
  double dxnorm, fp, parl, paru, gnorm, temp, parc, summ;
  p_double_array result = NULL, wa1 = NULL, wa2 = NULL, temp2 = NULL;
  p_double_array t1 = NULL, t3 = NULL, t4 = NULL, t5 = NULL;
  p_long_array t2 = NULL;

//  char info_string[ 1000 ];
//  sprintf( info_string, "iter = %ld", iter );
//  PyErr_SetString( PyExc_RuntimeError, info_string );
//  ok = 0; break;

  for ( once = 1; once; once-- )
  {
    r_dims[ 0 ] = 0L; r_dims[ 1 ] = 0L;
    if ( !check_double_array( r, 2L, &r_dims[ 0 ] ) )
    {
      PyErr_SetString( PyExc_ValueError, "r must be a two-dimensional matrix" );
      ok = 0; break;
    }
    m = r->dims[ 0 ];
    n = r->dims[ 1 ];
    if ( m != n )
    {
      PyErr_SetString( PyExc_ValueError, "r must be an n x n two-dimensional matrix" );
      ok = 0; break;
    }
    ipvt_dims[ 0 ] = n;
    if ( !check_long_array( ipvt, 1L, &ipvt_dims[ 0 ] ) )
    {
      PyErr_SetString( PyExc_ValueError, "ipvt must be a one-dimensional vector of length n" );
      ok = 0; break;
    }
    diag_dims[ 0 ] = n;
    if ( !check_double_array( diag, 1L, &diag_dims[ 0 ] ) )
    {
      PyErr_SetString( PyExc_ValueError, "diag must be a one-dimensional vector of length n" );
      ok = 0; break;
    }
    qtb_dims[ 0 ] = n;
    if ( !check_double_array( qtb, 1L, &qtb_dims[ 0 ] ) )
    {
      PyErr_SetString( PyExc_ValueError, "qtb must be a one-dimensional vector of length n" );
      ok = 0; break;
    }
    x_dims[ 0 ] = n;
    if ( !check_double_array( x, 1L, &x_dims[ 0 ] ) )
    {
      PyErr_SetString( PyExc_ValueError, "x must be a one-dimensional vector of length n" );
      ok = 0; break;
    }
    sdiag_dims[ 0 ] = n;
    if ( !check_double_array( sdiag, 1L, &sdiag_dims[ 0 ] ) )
    {
      PyErr_SetString( PyExc_ValueError, "sdiag must be a one-dimensional vector of length n" );
      ok = 0; break;
    }

    // Compute and store in x the gauss-newton direction.  If the
    // jacobian is rank-deficient, obtain a least-squares solution
    nsing = n;
    wa1 = create_double_array( qtb->num_dims, qtb->dims, qtb->data );
    if ( wa1 == NULL )
    {
      PyErr_SetString( PyExc_ValueError, "allocation of memory for wa1 failed" );
      ok = 0; break;
    }
    for ( i = 0L; i < n; i++ )
    {
      r_dims[ 0 ] = i; r_dims[ 1 ] = i;
      r_index = double_array_dims_to_index( r, &r_dims[ 0 ] );
      if ( r_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for r failed" );
        ok = 0; break;
      }
      if ( r->data[ r_index ] == 0. )
        break;
    }
    if ( !ok ) break;
    if ( i < n )
    {
      nsing = i;
      for ( i = nsing; i < n; i++ )
      {
        wa1_dims[ 0 ] = i;
        wa1_index = double_array_dims_to_index( wa1, &wa1_dims[ 0 ] );
        if ( wa1_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for wa1 failed" );
          ok = 0; break;
        }
        wa1->data[ wa1_index ] = 0.;
      }
      if ( !ok ) break;
    }

    if ( nsing >= 1L )
    {
      // *** Reverse loop ***
      for ( j = nsing - 1L; j >= 0L; j-- )
      {
        r_dims[ 0 ] = j; r_dims[ 1 ] = j;
        r_index = double_array_dims_to_index( r, &r_dims[ 0 ] );
        if ( r_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for r failed" );
          ok = 0; break;
        }
        wa1_dims[ 0 ] = j;
        wa1_index = double_array_dims_to_index( wa1, &wa1_dims[ 0 ] );
        if ( wa1_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for wa1 failed" );
          ok = 0; break;
        }
        wa1->data[ wa1_index ] = wa1->data[ wa1_index ] / r->data[ r_index ];
        if ( j > 0L )
        {
          temp = wa1->data[ wa1_index ];
          for ( i = 0L; i < j; i++ )
          {
            r_dims[ 0 ] = i; r_dims[ 1 ] = j;
            r_index = double_array_dims_to_index( r, &r_dims[ 0 ] );
            if ( r_index < 0L )
            {
              PyErr_SetString( PyExc_RuntimeError, "index calculation for r failed" );
              ok = 0; break;
            }
            wa1_dims[ 0 ] = i;
            wa1_index = double_array_dims_to_index( wa1, &wa1_dims[ 0 ] );
            if ( wa1_index < 0L )
            {
              PyErr_SetString( PyExc_RuntimeError, "index calculation for wa1 failed" );
              ok = 0; break;
            }
            wa1->data[ wa1_index ] = wa1->data[ wa1_index ] - r->data[ r_index ]* temp;
          }
          if ( !ok ) break;
        }
      }
      if ( !ok ) break;
    }

    // Note: ipvt here is a permutation array
    for ( i = 0L; i < n; i++ )
    {
      wa1_dims[ 0 ] = i;
      wa1_index = double_array_dims_to_index( wa1, &wa1_dims[ 0 ] );
      if ( wa1_index < 0L )
      {
        PyErr_SetString( PyExc_ValueError, "error while calculating index for wa1" );
        ok = 0; break;
      }
      ipvt_dims[ 0 ] = i;
      ipvt_index = long_array_dims_to_index( ipvt, &ipvt_dims[ 0 ] );
      if ( ipvt_index < 0L )
      {
        PyErr_SetString( PyExc_ValueError, "error while calculating index for ipvt" );
        ok = 0; break;
      }
      x_dims[ 0 ] = ipvt->data[ ipvt_index ];
      x_index = double_array_dims_to_index( x, &x_dims[ 0 ] );
      if ( x_index < 0L )
      {
        PyErr_SetString( PyExc_ValueError, "error while calculating index for x" );
        ok = 0; break;
      }
      x->data[ x_index ] = wa1->data[ wa1_index ];
    }
    if ( !ok ) break;

    // Evaluate the function at the
    // origin, and test for acceptance of the gauss-newton direction
    wa2_dims[ 0 ] = n;
    wa2 = create_double_array( 1L, &wa2_dims[ 0 ], NULL );
    if ( wa2 == NULL )
    {
      PyErr_SetString( PyExc_ValueError, "allocation of memory for wa2 failed" );
      ok = 0; break;
    }
    for ( i = 0L; i < n; i++ )
    {
      diag_dims[ 0 ] = i;
      diag_index = double_array_dims_to_index( diag, &diag_dims[ 0 ] );
      if ( diag_index < 0L )
      {
        PyErr_SetString( PyExc_ValueError, "error while calculating index for diag" );
        ok = 0; break;
      }
      x_dims[ 0 ] = i;
      x_index = double_array_dims_to_index( x, &x_dims[ 0 ] );
      if ( x_index < 0L )
      {
        PyErr_SetString( PyExc_ValueError, "error while calculating index for x" );
        ok = 0; break;
      }
      wa2_dims[ 0 ] = i;
      wa2_index = double_array_dims_to_index( wa2, &wa2_dims[ 0 ] );
      if ( wa2_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for wa2 failed" );
        ok = 0; break;
      }
      wa2->data[ wa2_index ] = diag->data[ diag_index ] * x->data[ x_index ];
    }
    if ( !ok ) break;
    temp2 = create_double_array( wa2->num_dims, wa2->dims, wa2->data );
    if ( temp2 == NULL )
    {
      PyErr_SetString( PyExc_ValueError, "allocation of memory for temp2 failed" );
      ok = 0; break;
    }
    dxnorm = enorm( temp2, rgiant, rdwarf ); // frees temp2
    temp2 = NULL;
    if ( dxnorm < 0. )
    {
      ok = 0; break; // error message already generated
    }
    fp = dxnorm - delta;
    if ( fp <= 0.1 * delta )
    {
      par = 0.;
      break; // exit function
    }

    // If the jacobian is not rank deficient, the newton step provides a
    // lower bound, parl, for the zero of the function.  Otherwise set
    // this bound to zero.
    parl = 0.;
    if ( nsing == n )
    {
      for ( i = 0L; i < n; i++ )
      {
        ipvt_dims[ 0 ] = i;
        ipvt_index = long_array_dims_to_index( ipvt, &ipvt_dims[ 0 ] );
        if ( ipvt_index < 0L )
        {
          PyErr_SetString( PyExc_ValueError, "error while calculating index for ipvt" );
          ok = 0; break;
        }
        diag_dims[ 0 ] = ipvt->data[ ipvt_index ];
        diag_index = double_array_dims_to_index( diag, &diag_dims[ 0 ] );
        if ( diag_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for diag failed" );
          ok = 0; break;
        }
        wa2_dims[ 0 ] = ipvt->data[ ipvt_index ];
        wa2_index = double_array_dims_to_index( wa2, &wa2_dims[ 0 ] );
        if ( wa2_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for wa2 failed" );
          ok = 0; break;
        }
        wa1_dims[ 0 ] = i;
        wa1_index = double_array_dims_to_index( wa1, &wa1_dims[ 0 ] );
        if ( wa1_index < 0L )
        {
          PyErr_SetString( PyExc_ValueError, "error while calculating index for wa1" );
          ok = 0; break;
        }
        wa1->data[ wa1_index ] = diag->data[ diag_index ] * wa2->data[ wa2_index ] / dxnorm;
      }
      if ( !ok ) break;

      // Degenerate case
      r_dims[ 0 ] = 0L; r_dims[ 1 ] = 0L;
      r_index = double_array_dims_to_index( r, &r_dims[ 0 ] );
      if ( r_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for r failed" );
        ok = 0; break;
      }
      wa1_dims[ 0 ] = 0L;
      wa1_index = double_array_dims_to_index( wa1, &wa1_dims[ 0 ] );
      if ( wa1_index < 0L )
      {
        PyErr_SetString( PyExc_ValueError, "error while calculating index for wa1" );
        ok = 0; break;
      }
      wa1->data[ wa1_index ] = wa1->data[ wa1_index ] / r->data[ r_index ];
      for ( j = 1L; j < n; j++ )
      {
        temp2_dims[ 0 ] = j;
        temp2 = create_double_array( 1L, &temp2_dims[ 0 ], NULL );
        if ( temp2 == NULL )
        {
          PyErr_SetString( PyExc_ValueError, "allocation of memory for temp2 failed" );
          ok = 0; break;
        }
        for ( i = 0L; i < j; i++ )
        {
          r_dims[ 0 ] = i; r_dims[ 1 ] = j;
          r_index = double_array_dims_to_index( r, &r_dims[ 0 ] );
          if ( r_index < 0L )
          {
            PyErr_SetString( PyExc_RuntimeError, "index calculation for r failed" );
            ok = 0; break;
          }
          wa1_dims[ 0 ] = i;
          wa1_index = double_array_dims_to_index( wa1, &wa1_dims[ 0 ] );
          if ( wa1_index < 0L )
          {
            PyErr_SetString( PyExc_ValueError, "error while calculating index for wa1" );
            ok = 0; break;
          }
          temp2_dims[ 0 ] = i;
          temp2_index = double_array_dims_to_index( temp2, &temp2_dims[ 0 ] );
          if ( temp2_index < 0L )
          {
            PyErr_SetString( PyExc_ValueError, "error while calculating index for temp2" );
            ok = 0; break;
          }
          temp2->data[ temp2_index ] = r->data[ r_index ] * wa1->data[ wa1_index ];
        }
        if ( !ok ) break;
        summ = sum( temp2 ); // frees temp2
        temp2 = NULL;
        r_dims[ 0 ] = j; r_dims[ 1 ] = j;
        r_index = double_array_dims_to_index( r, &r_dims[ 0 ] );
        if ( r_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for r failed" );
          ok = 0; break;
        }
        wa1_dims[ 0 ] = j;
        wa1_index = double_array_dims_to_index( wa1, &wa1_dims[ 0 ] );
        if ( wa1_index < 0L )
        {
          PyErr_SetString( PyExc_ValueError, "error while calculating index for wa1" );
          ok = 0; break;
        }
        wa1->data[ wa1_index ] = ( wa1->data[ wa1_index ] - summ ) / r->data[ r_index ];
      }
      if ( !ok ) break;

      temp2 = create_double_array( wa1->num_dims, wa1->dims, wa1->data );
      if ( temp2 == NULL )
      {
        PyErr_SetString( PyExc_ValueError, "allocation of memory for temp2 failed" );
        ok = 0; break;
      }
      temp = enorm( temp2, rgiant, rdwarf ); // frees temp2
      temp2 = NULL;
      if ( temp < 0. )
      {
        ok = 0; break; // error message already generated
      }
      parl = ( fp / delta ) / ( temp * temp );
    }

    // Calculate an upper bound, paru, for the zero of the function
    for ( j = 0L; j < n; j++ )
    {
      temp2_dims[ 0 ] = j + 1L;
      temp2 = create_double_array( 1L, &temp2_dims[ 0 ], NULL );
      if ( temp2 == NULL )
      {
        PyErr_SetString( PyExc_ValueError, "allocation of memory for temp2 failed" );
        ok = 0; break;
      }
      for ( i = 0L; i <= j; i++ )
      {
        r_dims[ 0 ] = i; r_dims[ 1 ] = j;
        r_index = double_array_dims_to_index( r, &r_dims[ 0 ] );
        if ( r_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for r failed" );
          ok = 0; break;
        }
        qtb_dims[ 0 ] = i;
        qtb_index = double_array_dims_to_index( qtb, &qtb_dims[ 0 ] );
        if ( qtb_index < 0L )
        {
          PyErr_SetString( PyExc_ValueError, "error while calculating index for qtb" );
          ok = 0; break;
        }
        temp2_dims[ 0 ] = i;
        temp2_index = double_array_dims_to_index( temp2, &temp2_dims[ 0 ] );
        if ( temp2_index < 0L )
        {
          PyErr_SetString( PyExc_ValueError, "error while calculating index for temp2" );
          ok = 0; break;
        }
        temp2->data[ temp2_index ] = r->data[ r_index ] * qtb->data[ qtb_index ];
      }
      if ( !ok ) break;
      summ = sum( temp2 ); // frees temp2
      temp2 = NULL;
      ipvt_dims[ 0 ] = j;
      ipvt_index = long_array_dims_to_index( ipvt, &ipvt_dims[ 0 ] );
      if ( ipvt_index < 0L )
      {
        PyErr_SetString( PyExc_ValueError, "error while calculating index for ipvt" );
        ok = 0; break;
      }
      diag_dims[ 0 ] = ipvt->data[ ipvt_index ];
      diag_index = double_array_dims_to_index( diag, &diag_dims[ 0 ] );
      if ( diag_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for diag failed" );
        ok = 0; break;
      }
      wa1_dims[ 0 ] = j;
      wa1_index = double_array_dims_to_index( wa1, &wa1_dims[ 0 ] );
      if ( wa1_index < 0L )
      {
        PyErr_SetString( PyExc_ValueError, "error while calculating index for wa1" );
        ok = 0; break;
      }
      wa1->data[ wa1_index ] = summ / diag->data[ diag_index ];
    }
    if ( !ok ) break;
    temp2 = create_double_array( wa1->num_dims, wa1->dims, wa1->data );
    if ( temp2 == NULL )
    {
      PyErr_SetString( PyExc_ValueError, "allocation of memory for temp2 failed" );
      ok = 0; break;
    }
    gnorm = enorm( temp2, rgiant, rdwarf ); // frees temp2
    temp2 = NULL;
    if ( gnorm < 0. )
    {
      ok = 0; break; // error message already generated
    }
    paru = gnorm / delta;
    if ( paru == 0. )
    {
      temp = ( delta < 0.1 ? delta : 0.1 );
      paru = dwarf / temp;
    }

    // If the input par lies outside of the interval (parl,paru), set
    // par to the closer endpoint
    par = ( par > parl ? par : parl );
    par = ( par < paru ? par : paru );
    if ( par == 0. )
      par = gnorm / dxnorm;

    // Beginning of an interation
    for ( iter = 0L; iter <= 10L; iter++ )
    {
      // Evaluate the function at the current value of par
      if ( par == 0. )
        par = ( dwarf > paru * 0.001 ? dwarf : paru * 0.001 );
      temp = sqrt( par );
      for ( i = 0L; i < n; i++ )
      {
        diag_dims[ 0 ] = i;
        diag_index = double_array_dims_to_index( diag, &diag_dims[ 0 ] );
        if ( diag_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for diag failed" );
          ok = 0; break;
        }
        wa1_dims[ 0 ] = i;
        wa1_index = double_array_dims_to_index( wa1, &wa1_dims[ 0 ] );
        if ( wa1_index < 0L )
        {
          PyErr_SetString( PyExc_ValueError, "error while calculating index for wa1" );
          ok = 0; break;
        }
        wa1->data[ wa1_index ] = temp * diag->data[ diag_index ];
      }
      if ( !ok ) break;

      // call qrsolv()
      t1 = create_double_array( r->num_dims, r->dims, r->data );
      if ( t1 == NULL )
      {
        PyErr_SetString( PyExc_ValueError, "allocation of memory for t1 failed" );
        ok = 0; break;
      }
      t2 = create_long_array( ipvt->num_dims, ipvt->dims, ipvt->data );
      if ( t2 == NULL )
      {
        PyErr_SetString( PyExc_ValueError, "allocation of memory for t2 failed" );
        ok = 0; break;
      }
      t3 = create_double_array( wa1->num_dims, wa1->dims, wa1->data );
      if ( t3 == NULL )
      {
        PyErr_SetString( PyExc_ValueError, "allocation of memory for t3 failed" );
        ok = 0; break;
      }
      t4 = create_double_array( qtb->num_dims, qtb->dims, qtb->data );
      if ( t4 == NULL )
      {
        PyErr_SetString( PyExc_ValueError, "allocation of memory for t4 failed" );
        ok = 0; break;
      }
      t5 = create_double_array( sdiag->num_dims, sdiag->dims, sdiag->data );
      if ( t5 == NULL )
      {
        PyErr_SetString( PyExc_ValueError, "allocation of memory for t5 failed" );
        ok = 0; break;
      }
      result = qrsolv( t1, t2, t3, t4, t5 ); // frees t1, t2, t3, t4, t5
      t1 = NULL; t2 = NULL; t3 = NULL; t4 = NULL; t5 = NULL;
      result_dims[ 0 ] = n + 2L; result_dims[ 1 ] = n;
      if ( !check_double_array( result, 2L, &result_dims[ 0 ] ) )
      {
        PyErr_SetString( PyExc_RuntimeError, "result has unexpected dimensions" );
        ok = 0; break;
      }
      for ( j = 0L; j < n + 2L; j++ )
      {
        for ( k = 0L; k < n; k++ )
        {
          result_dims[ 0 ] = j; result_dims[ 1 ] = k;
          result_index = double_array_dims_to_index( result, &result_dims[ 0 ] );
          if ( result_index < 0L )
          {
            PyErr_SetString( PyExc_RuntimeError, "index calculation for result failed" );
            ok = 0; break;
          }
          if ( j < n )
          {
            r_dims[ 0 ] = j; r_dims[ 1 ] = k;
            r_index = double_array_dims_to_index( r, &r_dims[ 0 ] );
            if ( r_index < 0L )
            {
              PyErr_SetString( PyExc_RuntimeError, "index calculation for r failed" );
              ok = 0; break;
            }
            r->data[ r_index ] = result->data[ result_index ];
          }
          else if ( j == n )
          {
            x_dims[ 0 ] = k;
            x_index = double_array_dims_to_index( x, &x_dims[ 0 ] );
            if ( x_index < 0L )
            {
              PyErr_SetString( PyExc_RuntimeError, "index calculation for x failed" );
              ok = 0; break;
            }
            x->data[ x_index ] = result->data[ result_index ];
          }
          else if ( j == n + 1L )
          {
            sdiag_dims[ 0 ] = k;
            sdiag_index = double_array_dims_to_index( sdiag, &sdiag_dims[ 0 ] );
            if ( sdiag_index < 0L )
            {
              PyErr_SetString( PyExc_RuntimeError, "index calculation for sdiag failed" );
              ok = 0; break;
            }
            sdiag->data[ sdiag_index ] = result->data[ result_index ];
          }
        }
        if ( !ok ) break;
      }
      if ( !ok ) break;
      free_double_array( result ); result = NULL;

      for ( i = 0L; i < n; i++ )
      {
        diag_dims[ 0 ] = i;
        diag_index = double_array_dims_to_index( diag, &diag_dims[ 0 ] );
        if ( diag_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for diag failed" );
          ok = 0; break;
        }
        x_dims[ 0 ] = i;
        x_index = double_array_dims_to_index( x, &x_dims[ 0 ] );
        if ( x_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for x failed" );
          ok = 0; break;
        }
        wa2_dims[ 0 ] = i;
        wa2_index = double_array_dims_to_index( wa2, &wa2_dims[ 0 ] );
        if ( wa2_index < 0L )
        {
          PyErr_SetString( PyExc_ValueError, "error while calculating index for wa2" );
          ok = 0; break;
        }
        wa2->data[ wa2_index ] = diag->data[ diag_index ] * x->data[ x_index ];
      }
      if ( !ok ) break;
      temp2 = create_double_array( wa2->num_dims, wa2->dims, wa2->data );
      if ( temp2 == NULL )
      {
        PyErr_SetString( PyExc_ValueError, "allocation of memory for temp2 failed" );
        ok = 0; break;
      }
      dxnorm = enorm( temp2, rgiant, rdwarf ); // frees temp2
      temp2 = NULL;
      if ( dxnorm < 0. )
      {
        ok = 0; break; // error message already generated
      }
      temp = fp;
      fp = dxnorm - delta;

      if ( ( fabs( fp ) <= 0.1 * delta ) ||
           ( ( parl == 0. ) && ( fp <= temp ) && ( temp < 0. ) ) ||
           ( iter == 10L ) )
        break; // end iterations

      // Compute the newton correction
      for ( i = 0L; i < n; i++ )
      {
        ipvt_dims[ 0 ] = i;
        ipvt_index = long_array_dims_to_index( ipvt, &ipvt_dims[ 0 ] );
        if ( ipvt_index < 0L )
        {
          PyErr_SetString( PyExc_ValueError, "error while calculating index for ipvt" );
          ok = 0; break;
        }
        diag_dims[ 0 ] = ipvt->data[ ipvt_index ];
        diag_index = double_array_dims_to_index( diag, &diag_dims[ 0 ] );
        if ( diag_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for diag failed" );
          ok = 0; break;
        }
        wa2_dims[ 0 ] = ipvt->data[ ipvt_index ];
        wa2_index = double_array_dims_to_index( wa2, &wa2_dims[ 0 ] );
        if ( wa2_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for wa2 failed" );
          ok = 0; break;
        }
        wa1_dims[ 0 ] = i;
        wa1_index = double_array_dims_to_index( wa1, &wa1_dims[ 0 ] );
        if ( wa1_index < 0L )
        {
          PyErr_SetString( PyExc_ValueError, "error while calculating index for wa1" );
          ok = 0; break;
        }
        wa1->data[ wa1_index ] = diag->data[ diag_index ] * wa2->data[ wa2_index ] / dxnorm;
      }
      if ( !ok ) break;

      for ( j = 0L; j < n - 1L; j++ )
      {
        sdiag_dims[ 0 ] = j;
        sdiag_index = double_array_dims_to_index( sdiag, &sdiag_dims[ 0 ] );
        if ( sdiag_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for sdiag failed" );
          ok = 0; break;
        }
        wa1_dims[ 0 ] = j;
        wa1_index = double_array_dims_to_index( wa1, &wa1_dims[ 0 ] );
        if ( wa1_index < 0L )
        {
          PyErr_SetString( PyExc_ValueError, "error while calculating index for wa1" );
          ok = 0; break;
        }
        wa1->data[ wa1_index ] = wa1->data[ wa1_index ] / sdiag->data[ diag_index ];
        temp = wa1->data[ wa1_index ];
        for ( i = j + 1L; i < n; i++ )
        {
          r_dims[ 0 ] = i; r_dims[ 1 ] = j;
          r_index = double_array_dims_to_index( r, &r_dims[ 0 ] );
          if ( r_index < 0L )
          {
            PyErr_SetString( PyExc_RuntimeError, "index calculation for r failed" );
            ok = 0; break;
          }
          wa1_dims[ 0 ] = i;
          wa1_index = double_array_dims_to_index( wa1, &wa1_dims[ 0 ] );
          if ( wa1_index < 0L )
          {
            PyErr_SetString( PyExc_ValueError, "error while calculating index for wa1" );
            ok = 0; break;
          }
          wa1->data[ wa1_index ] = wa1->data[ wa1_index ] - r->data[ r_index ] * temp;
        }
        if ( !ok ) break;
      }
      if ( !ok ) break;
      // Degenerate case
      sdiag_dims[ 0 ] = n - 1L;
      sdiag_index = double_array_dims_to_index( sdiag, &sdiag_dims[ 0 ] );
      if ( sdiag_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for sdiag failed" );
        ok = 0; break;
      }
      wa1_dims[ 0 ] = n - 1L;
      wa1_index = double_array_dims_to_index( wa1, &wa1_dims[ 0 ] );
      if ( wa1_index < 0L )
      {
        PyErr_SetString( PyExc_ValueError, "error while calculating index for wa1" );
        ok = 0; break;
      }
      wa1->data[ wa1_index ] = wa1->data[ wa1_index ] / sdiag->data[ sdiag_index ];

      temp2 = create_double_array( wa1->num_dims, wa1->dims, wa1->data );
      if ( temp2 == NULL )
      {
        PyErr_SetString( PyExc_ValueError, "allocation of memory for temp2 failed" );
        ok = 0; break;
      }
      temp = enorm( temp2, rgiant, rdwarf ); // frees temp2
      temp2 = NULL;
      if ( temp < 0. )
      {
        ok = 0; break; // error message already generated
      }
      parc = ( fp / delta ) / ( temp * temp );

      // Depending on the sign of the function, update parl or paru
      if ( fp > 0. )
        parl = ( parl > par ? parl : par );
      if ( fp < 0. )
        paru = ( paru < par ? paru : par );

      // Compute an improved estimate for par
      par = ( parl > par + parc ? parl : par + parc );

      // End of an iteration
    }
    if ( !ok ) break;
  }

  for ( once = 1; ok && once; once-- )
  {
    // build result to be returned
    result_dims[ 0 ] = n + 3L; result_dims[ 1 ] = n;
    result = create_double_array( 2L, &result_dims[ 0 ], NULL );
    if ( result == NULL )
    {
      PyErr_SetString( PyExc_ValueError, "allocation of memory for result failed" );
      ok = 0; break;
    }
    for ( j = 0L; j < n + 3L; j++ )
    {
      for ( k = 0L; k < n; k++ )
      {
        result_dims[ 0 ] = j; result_dims[ 1 ] = k;
        result_index = double_array_dims_to_index( result, &result_dims[ 0 ] );
        if ( result_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for result failed" );
          ok = 0; break;
        }
        if ( j < n )
        {
          r_dims[ 0 ] = j; r_dims[ 1 ] = k;
          r_index = double_array_dims_to_index( r, &r_dims[ 0 ] );
          if ( r_index < 0L )
          {
            PyErr_SetString( PyExc_RuntimeError, "index calculation for r failed" );
            ok = 0; break;
          }
          result->data[ result_index ] = r->data[ r_index ];
        }
        else if ( j == n )
        {
          result->data[ result_index ] = par;
        }
        else if ( j == n + 1L )
        {
          x_dims[ 0 ] = k;
          x_index = double_array_dims_to_index( x, &x_dims[ 0 ] );
          if ( x_index < 0L )
          {
            PyErr_SetString( PyExc_RuntimeError, "index calculation for x failed" );
            ok = 0; break;
          }
          result->data[ result_index ] = x->data[ x_index ];
        }
        else if ( j == n + 2L )
        {
          sdiag_dims[ 0 ] = k;
          sdiag_index = double_array_dims_to_index( sdiag, &sdiag_dims[ 0 ] );
          if ( sdiag_index < 0L )
          {
            PyErr_SetString( PyExc_RuntimeError, "index calculation for sdiag failed" );
            ok = 0; break;
          }
          result->data[ result_index ] = sdiag->data[ sdiag_index ];
        }
      }
      if ( !ok ) break;
    }
  }

  // cleanup
  if ( r != NULL )
  {
    free_double_array( r ); r = NULL;
  }
  if ( ipvt != NULL )
  {
    free_long_array( ipvt ); ipvt = NULL;
  }
  if ( diag != NULL )
  {
    free_double_array( diag ); diag = NULL;
  }
  if ( qtb != NULL )
  {
    free_double_array( qtb ); qtb = NULL;
  }
  if ( x != NULL )
  {
    free_double_array( x ); x = NULL;
  }
  if ( sdiag != NULL )
  {
    free_double_array( sdiag ); sdiag = NULL;
  }
  if ( wa1 != NULL )
  {
    free_double_array( wa1 ); wa1 = NULL;
  }
  if ( wa2 != NULL )
  {
    free_double_array( wa2 ); wa2 = NULL;
  }
  if ( temp2 != NULL )
  {
    free_double_array( temp2 ); temp2 = NULL;
  }
  if ( t1 != NULL )
  {
    free_double_array( t1 ); t1 = NULL;
  }
  if ( t2 != NULL )
  {
    free_long_array( t2 ); t2 = NULL;
  }
  if ( t3 != NULL )
  {
    free_double_array( t3 ); t3 = NULL;
  }
  if ( t4 != NULL )
  {
    free_double_array( t4 ); t4 = NULL;
  }
  if ( t5 != NULL )
  {
    free_double_array( t5 ); t5 = NULL;
  }
  if ( !ok )
  {
    if ( result != NULL )
    {
      free_double_array( result ); result = NULL;
    }
  }
  return ( result );
}
// result is an ( n + 3 ) x n double matrix, consisting of:
// r = n x n double matrix
// par = n double vector ( single value repeated )
// x = n double vector
// sdiag = n double vector

//*****************************************************************************

p_double_array calc_covar( p_double_array rr_table, p_long_array ipvt_table, double tol )
{
  int once, ok = 1, sing;
  long rr_dims[ 2 ], m, n, i, l, r_dims[ 2 ], k, j, jj, ipvt_dims[ 1 ], ipvt_index, r_index;
  long wa_dims[ 1 ], wa_index, ii;
  double tolr, temp, temp2;
  p_double_array r_table = NULL, wa_table = NULL;
  for ( once = 1; once; once-- )
  {
    if ( tol < 0. )
      tol = 1.e-14;
    rr_dims[ 0 ] = 0L; rr_dims[ 1 ] = 0L;
    if ( !check_double_array( rr_table, 2L, &rr_dims[ 0 ] ) )
    {
      PyErr_SetString( PyExc_ValueError, "rr_table must be a two-dimensional matrix" );
      ok = 0; break;
    }
    m = rr_table->dims[ 0 ];
    n = rr_table->dims[ 1 ];
    if ( m != n )
    {
      PyErr_SetString( PyExc_ValueError, "rr_table must be a square matrix" );
      ok = 0; break;
    }
    if ( check_long_array( ipvt_table, - 1L, NULL ) )
    {
      free_long_array( ipvt_table );
      ipvt_table = create_long_array( 1L, &n, NULL );
      if ( ipvt_table == NULL )
      {
        PyErr_SetString( PyExc_RuntimeError, "allocation of memory for ipvt_table failed" );
        ok = 0; break;
      }
      for ( i = 0L; i < n; i++ )
      {
        ipvt_dims[ 0 ] = i;
        ipvt_index = long_array_dims_to_index( ipvt_table, &ipvt_dims[ 0 ] );
        if ( ipvt_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for ipvt_table failed" );
          ok = 0; break;
        }
        ipvt_table->data[ ipvt_index ] = i;
      }
    }
    ipvt_dims[ 0 ] = 0L;
    if ( !check_long_array( ipvt_table, 1L, &ipvt_dims[ 0 ] ) )
    {
      PyErr_SetString( PyExc_ValueError, "ipvt_table must be a one-dimensional array" );
      ok = 0; break;
    }
    r_table = create_double_array( rr_table->num_dims, rr_table->dims, rr_table->data );
    if ( r_table == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "allocation of memory for r_table failed" );
      ok = 0; break;
    }

    // Form the inverse of r in the full upper triangle of r
    l = - 1L;
    r_dims[ 0 ] = 0L; r_dims[ 1 ] = 0L;
    r_index = double_array_dims_to_index( r_table, &r_dims[ 0 ] );
    if ( r_index < 0L )
    {
      PyErr_SetString( PyExc_RuntimeError, "index calculation for r_table failed" );
      ok = 0; break;
    }
    tolr = tol * fabs( r_table->data[ r_index ] );
    for ( k = 0L; k < n; k++ )
    {
      r_dims[ 0 ] = k; r_dims[ 1 ] = k;
      r_index = double_array_dims_to_index( r_table, &r_dims[ 0 ] );
      if ( r_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for r_table failed" );
        ok = 0; break;
      }
      if ( fabs( r_table->data[ r_index ] <= tolr ) )
        break;
      r_table->data[ r_index ] = 1. / r_table->data[ r_index ];
      for ( j = 0L; j < k; j++ )
      {
        r_dims[ 0 ] = k; r_dims[ 1 ] = k;
        r_index = double_array_dims_to_index( r_table, &r_dims[ 0 ] );
        if ( r_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for r_table failed" );
          ok = 0; break;
        }
        temp = r_table->data[ r_index ];
        r_dims[ 0 ] = j; r_dims[ 1 ] = k;
        r_index = double_array_dims_to_index( r_table, &r_dims[ 0 ] );
        if ( r_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for r_table failed" );
          ok = 0; break;
        }
        temp = temp * r_table->data[ r_index ];
        r_table->data[ r_index ] = 0.;
        for ( i = 0L; i <= j; i++ )
        {
          r_dims[ 0 ] = i; r_dims[ 1 ] = j;
          r_index = double_array_dims_to_index( r_table, &r_dims[ 0 ] );
          if ( r_index < 0L )
          {
            PyErr_SetString( PyExc_RuntimeError, "index calculation for r_table failed" );
            ok = 0; break;
          }
          temp2 = temp * r_table->data[ r_index ];
          r_dims[ 0 ] = i; r_dims[ 1 ] = k;
          r_index = double_array_dims_to_index( r_table, &r_dims[ 0 ] );
          if ( r_index < 0L )
          {
            PyErr_SetString( PyExc_RuntimeError, "index calculation for r_table failed" );
            ok = 0; break;
          }
          r_table->data[ r_index ] = r_table->data[ r_index ] - temp2;
        }
        if ( !ok ) break;
      }
      if ( !ok ) break;
      l = k;
    }
    if ( !ok ) break;
    break;

    // Form the full upper triangle of the inverse of (r transpose)*r
    // in the full upper triangle of r
    if ( l >= 0L )
    {
      for ( k = 0L; k <= l; k++ )
      {
        for ( j = 0L; j < k; j++ )
        {
          r_dims[ 0 ] = j; r_dims[ 1 ] = k;
          r_index = double_array_dims_to_index( r_table, &r_dims[ 0 ] );
          if ( r_index < 0L )
          {
            PyErr_SetString( PyExc_RuntimeError, "index calculation for r_table failed" );
            ok = 0; break;
          }
          temp = r_table->data[ r_index ];
          for ( i = 0L; i <= j; i++ )
          {
            r_dims[ 0 ] = i; r_dims[ 1 ] = k;
            r_index = double_array_dims_to_index( r_table, &r_dims[ 0 ] );
            if ( r_index < 0L )
            {
              PyErr_SetString( PyExc_RuntimeError, "index calculation for r_table failed" );
              ok = 0; break;
            }
            temp2 = temp * r_table->data[ r_index ];
            r_dims[ 0 ] = i; r_dims[ 1 ] = j;
            r_index = double_array_dims_to_index( r_table, &r_dims[ 0 ] );
            if ( r_index < 0L )
            {
              PyErr_SetString( PyExc_RuntimeError, "index calculation for r_table failed" );
              ok = 0; break;
            }
            r_table->data[ r_index ] = r_table->data[ r_index ] + temp2;
          }
          if ( !ok ) break;
        }
        if ( !ok ) break;
        r_dims[ 0 ] = k; r_dims[ 1 ] = k;
        r_index = double_array_dims_to_index( r_table, &r_dims[ 0 ] );
        if ( r_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for r_table failed" );
          ok = 0; break;
        }
        temp = r_table->data[ r_index ];
        for ( i = 0L; i <= k; i++ )
        {
          r_dims[ 0 ] = i; r_dims[ 1 ] = k;
          r_index = double_array_dims_to_index( r_table, &r_dims[ 0 ] );
          if ( r_index < 0L )
          {
            PyErr_SetString( PyExc_RuntimeError, "index calculation for r_table failed" );
            ok = 0; break;
          }
          r_table->data[ r_index ] = temp * r_table->data[ r_index ];
        }
        if ( !ok ) break;
      }
      if ( !ok ) break;
    }
    if ( !ok ) break;

    // Form the full lower triangle of the covariance matrix
    // in the strict lower triangle or and in wa
    wa_table = create_double_array( 1L, &n, NULL );
    if ( wa_table == NULL )
    {
      PyErr_SetString( PyExc_RuntimeError, "allocation of memory for wa_table failed" );
      ok = 0; break;
    }
    r_dims[ 0 ] = 0L; r_dims[ 1 ] = 0L;
    r_index = double_array_dims_to_index( r_table, &r_dims[ 0 ] );
    if ( r_index < 0L )
    {
      PyErr_SetString( PyExc_RuntimeError, "index calculation for r_table failed" );
      ok = 0; break;
    }
    for ( i = 0L; i < n; i++ )
    {
      wa_dims[ 0 ] = i;
      wa_index = double_array_dims_to_index( wa_table, &wa_dims[ 0 ] );
      if ( wa_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for wa_table failed" );
        ok = 0; break;
      }
      wa_table->data[ wa_index ] = r_table->data[ r_index ];
    }
    for ( j = 0L; j < n; j++ )
    {
      ipvt_dims[ 0 ] = j;
      ipvt_index = long_array_dims_to_index( ipvt_table, &ipvt_dims[ 0 ] );
      if ( ipvt_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for ipvt_table failed" );
        ok = 0; break;
      }
      jj = ipvt_table->data[ ipvt_index ];
      sing = ( j > l ? 1 : 0 );
      for ( i = 0L; i <= j; i++ )
      {
        r_dims[ 0 ] = i; r_dims[ 1 ] = j;
        r_index = double_array_dims_to_index( r_table, &r_dims[ 0 ] );
        if ( r_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for r_table failed" );
          ok = 0; break;
        }
        if ( sing )
        {
          r_table->data[ r_index ] = 0.;
        }
        ipvt_dims[ 0 ] = i;
        ipvt_index = long_array_dims_to_index( ipvt_table, &ipvt_dims[ 0 ] );
        if ( ipvt_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for ipvt_table failed" );
          ok = 0; break;
        }
        ii = ipvt_table->data[ ipvt_index ];
        temp = r_table->data[ r_index ];
        if ( ii > jj )
        {
          r_dims[ 0 ] = ii; r_dims[ 1 ] = jj;
          r_index = double_array_dims_to_index( r_table, &r_dims[ 0 ] );
          if ( r_index < 0L )
          {
            PyErr_SetString( PyExc_RuntimeError, "index calculation for r_table failed" );
            ok = 0; break;
          }
          r_table->data[ r_index ] = temp;
        }
        if ( ii < jj )
        {
          r_dims[ 0 ] = jj; r_dims[ 1 ] = ii;
          r_index = double_array_dims_to_index( r_table, &r_dims[ 0 ] );
          if ( r_index < 0L )
          {
            PyErr_SetString( PyExc_RuntimeError, "index calculation for r_table failed" );
            ok = 0; break;
          }
          r_table->data[ r_index ] = temp;
        }
      }
      if ( !ok ) break;
      wa_dims[ 0 ] = jj;
      wa_index = double_array_dims_to_index( wa_table, &wa_dims[ 0 ] );
      if ( wa_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for wa_table failed" );
        ok = 0; break;
      }
      r_dims[ 0 ] = j; r_dims[ 1 ] = j;
      r_index = double_array_dims_to_index( r_table, &r_dims[ 0 ] );
      if ( r_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for r_table failed" );
        ok = 0; break;
      }
      wa_table->data[ wa_index ] = r_table->data[ r_index ];
    }
    if ( !ok ) break;

    // Symmetrize the covariance matrix in r
    for ( j = 0L; j < n; j++ )
    {
      for ( i = 0; i <= j; i++ )
      {
        r_dims[ 0 ] = j; r_dims[ 1 ] = i;
        r_index = double_array_dims_to_index( r_table, &r_dims[ 0 ] );
        if ( r_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for r_table failed" );
          ok = 0; break;
        }
        temp = r_table->data[ r_index ];
        r_dims[ 0 ] = i; r_dims[ 1 ] = j;
        r_index = double_array_dims_to_index( r_table, &r_dims[ 0 ] );
        if ( r_index < 0L )
        {
          PyErr_SetString( PyExc_RuntimeError, "index calculation for r_table failed" );
          ok = 0; break;
        }
        r_table->data[ r_index ] = temp;
      }
      if ( !ok ) break;


      r_dims[ 0 ] = j; r_dims[ 1 ] = j;
      r_index = double_array_dims_to_index( r_table, &r_dims[ 0 ] );
      if ( r_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for r_table failed" );
        ok = 0; break;
      }
      wa_dims[ 0 ] = j;
      wa_index = double_array_dims_to_index( wa_table, &wa_dims[ 0 ] );
      if ( wa_index < 0L )
      {
        PyErr_SetString( PyExc_RuntimeError, "index calculation for wa_table failed" );
        ok = 0; break;
      }
      r_table->data[ r_index ] = wa_table->data[ wa_index ];
    }
    if ( !ok ) break;
  }

  // cleanup
  if ( rr_table != NULL )
  {
    free_double_array( rr_table ); 
    rr_table = NULL;
  }
  if ( ipvt_table != NULL )
  {
    free_long_array( ipvt_table ); 
    ipvt_table = NULL;
  }
  if ( wa_table != NULL )
  {
    free_double_array( wa_table ); 
    wa_table = NULL;
  }
  if ( !ok )
  {
    if ( r_table != NULL )
    {
      free_double_array( r_table );
      r_table = NULL;
    }
  }
  return ( r_table );
}

//*****************************************************************************

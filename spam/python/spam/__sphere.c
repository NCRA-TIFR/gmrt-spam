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
#include "__sphere.h"

//*****************************************************************************

p_double_2 hmsdms_to_degdeg( double hmsdms[ 6 ] )
{
  double ra, dec;
  p_double_2 rd = malloc( 2 * sizeof( double ) );
  ra = modulo( hmsdms[ 0 ], 24. );
  rd[ 0 ] = 15. * ( ra + ( hmsdms[ 1 ] / 60. ) + ( hmsdms[ 2 ] / 3600. ) );
  dec = fabs( hmsdms[ 3 ] );
  if ( dec > 90. )
    rd[ 1 ] = 90.;
  else
    rd[ 1 ] = dec + ( hmsdms[ 4 ] / 60. ) + ( hmsdms[ 5 ] / 3600. );
  if ( sign( hmsdms[ 3 ] ) < 0. )
    rd[ 1 ] = - rd[ 1 ];
  return ( rd );
}

//*****************************************************************************

p_double_6 degdeg_to_hmsdms( double degdeg[ 2 ] )
{
  double ra, dec;
  p_double_6 rd = malloc( 6 * sizeof( double ) );
  ra = modulo( degdeg[ 0 ], 360. ) / 15.;
  rd[ 0 ] = floor( ra );
  ra = 60. * ( ra - rd[ 0 ] );
  rd[ 1 ] = floor( ra );
  rd[ 2 ] = 60. * ( ra - rd[ 1 ] );
  dec = fabs( degdeg[ 1 ] );
  if ( dec > 90. )
  {
    rd[ 3 ] = 90.;
    rd[ 4 ] = 0.;
    rd[ 5 ] = 0.;
  }
  else
  {
    rd[ 3 ] = floor( dec );
    dec = 60. * ( dec - rd[ 3 ] );
    rd[ 4 ] = floor( dec );
    rd[ 5 ] = 60. * ( dec - rd[ 4 ] );
  }
  if ( sign( degdeg[ 1 ] ) < 0. )
    rd[ 3 ] = - rd[ 3 ];
  return ( rd );
}

//*****************************************************************************

p_double_6 degdeg_to_dmsdms( double degdeg[ 2 ] )
{
  double lon, slon, lat;
  p_double_6 ll = malloc( 6 * sizeof( double ) );
  lon = modulo( degdeg[ 0 ] + 180., 360. ) - 180.;
  slon = sign( lon );
  lon = fabs( lon );
  ll[ 0 ] = floor( lon );
  lon = 60. * ( lon - ll[ 0 ] );
  if ( slon < 0. )
    ll[ 0 ] = - ll[ 0 ];
  ll[ 1 ] = floor( lon );
  ll[ 2 ] = 60. * ( lon - ll[ 1 ] );
  lat = fabs( degdeg[ 1 ] );
  if ( lat > 90. )
  {
    ll[ 3 ] = 90.;
    ll[ 4 ] = 0.;
    ll[ 5 ] = 0.;
  }
  else
  {
    ll[ 3 ] = floor( lat );
    lat = 60. * ( lat - ll[ 3 ] );
    ll[ 4 ] = floor( lat );
    ll[ 5 ] = 60. * ( lat - ll[ 4 ] );
  }
  if ( sign( degdeg[ 1 ] ) < 0. )
    ll[ 3 ] = - ll[ 3 ];
  return ( ll );
}

//*****************************************************************************

p_double_2 calculate_angular_separation( double radec1[ 2 ], double radec2[ 2 ] )
{
  double ra1, dec1, ra2, dec2, hav_r, hav_a, t, dra, ddec;
  p_double_2 ra = malloc( 2 * sizeof( double ) );
  ra1 = modulo( radec1[ 0 ], 360. );
  dec1 = radec1[ 1 ];
  if ( fabs( dec1 ) > 90. )
    dec1 = sign( dec1 ) * 90.;
  ra2 = modulo( radec2[ 0 ], 360. );
  dec2 = radec2[ 1 ];
  if ( fabs( dec2 ) > 90. )
    dec2 = sign( dec2 ) * 90.;
  dra = modulo( ( ra2 - ra1 ) + 180., 360. ) - 180.;
  ddec = dec2 - dec1;
  if ( ( dra == 0. ) || ( fabs( dec1 ) == 90. ) || ( fabs( dec2 ) == 90. ) )
  {
    ra[ 0 ] = fabs( ddec );
    if ( dec2 >= dec1 )
      ra[ 1 ] = 0.;
    else
      ra[ 1 ] = 180.;
  }
//  else if ( ddec == 0. )
//  {
//    ra[ 0 ] = dra * cos( radians( dec1 ) );
//    if ( ra[ 0 ] > 0. )
//      ra[ 1 ] = 90.;
//    else
//      ra[ 1 ] = 270.;
//    ra[ 0 ] = fabs( ra[ 0 ] );
//  }
  else
  {
    // use Haversine formula (adapted from http://www.plutoproject.com/dist_pa2.cpp)
    hav_r = ( sin( radians( ddec / 2. ) ) * sin( radians( ddec / 2. ) ) + 
        cos( radians( dec1 ) ) * cos( radians( dec2 ) ) * 
        sin( radians( dra / 2. ) ) * sin( radians( dra / 2. ) ) );
    if ( hav_r > 1. )
      hav_r = 1.;
    ra[ 0 ] = degrees( 2. * asin( sqrt( hav_r ) ) );
    t = cos( radians( dec1 ) ) * sin( radians( ra[ 0 ] ) );
    hav_a = ( sin( radians( dec1 + ra[ 0 ] ) ) - sin( radians( dec2 ) ) ) / ( 2. * t );
    if ( hav_a > 1. )
      hav_a = 1.;
    ra[ 1 ] = degrees( 2. * asin( sqrt( hav_a ) ) );
    if ( sin( radians( dra ) ) < 0. )
      ra[ 1 ] = 360. - ra[ 1 ];
  }
  return ( ra );
}

//*****************************************************************************

p_double_2 calculate_offset_position( double radec[ 2 ], double radius, double angle )
// 0. <= radius <= 180.
{
  double ra, dec, a, b, c, A, B, bb, AA;
  p_double_2 rd = malloc( 2 * sizeof( double ) );
  ra = radec[ 0 ];
  dec = radec[ 1 ];
  if ( fabs( dec ) > 90. )
    dec = sign( dec ) * 90.;
  if ( radius <= 0. )
  {
    rd[ 0 ] = ra;
    rd[ 1 ] = dec;
  }
  else
  {
    a = radians( radius );
    c = radians( 90. - dec );
    B = radians( - angle );
    bb = sin( a ) * cos( B ) * sin( c ) + cos( a ) * cos( c );
    if ( fabs( bb ) > 1. )
      bb = sign( bb ) * 1.;
    b = acos( bb );
    if ( b == 0. )
      A = 0.;
    else
    {
      AA = sin( a ) * sin( B ) / sin( b );
      if ( fabs( AA ) > 1. )
        AA = sign( AA ) * 1.;
      A = asin( AA );
      if ( ( ( cos( a ) * sin( c ) - sin( a ) * cos( B ) * cos( c ) ) / sin( b ) ) < 0. )
        A = radians( 180. ) - A;
    }
    rd[ 0 ] = modulo( ra - degrees( A ), 360. );
    rd[ 1 ] = 90. - degrees( b );
  }
  return ( rd );
}

//*****************************************************************************

p_double_3 xyz_to_llr( double xyz[ 3 ] )
{
  p_double_3 llr = malloc( 3 * sizeof( double ) );
  llr[ 0 ] = modulo( degrees( atan2( xyz[ 1 ], xyz[ 0 ] ) ) + 180., 360. ) - 180.;
  llr[ 1 ] = degrees( atan2( xyz[ 2 ], sqrt( xyz[ 0 ] * xyz[ 0 ] + xyz[ 1 ] * xyz[ 1 ] ) ) );
  llr[ 2 ] = sqrt( xyz[ 0 ] * xyz[ 0 ] + xyz[ 1 ] * xyz[ 1 ] + xyz[ 2 ] * xyz[ 2 ] );
  return ( llr );
}

//*****************************************************************************

p_double_3 xyz_to_geo_llh( double xyz[ 3 ], int iterations_, double a_, double f_, double e2_ )
// # default Earth ellipticity definition (a,f) is WGS (1984)
// Note that longitude is defined as positive towards east, just like RA
{
  int iterations, i;
  double a, f, e2, phi, n;
  p_double_3 gllh = malloc( 3 * sizeof( double ) );
  if ( iterations_ < 0 ) iterations = 4;
  else iterations = iterations_;
  if ( a_ == 0. ) a = 6378137.;
  else a = a_;
  if ( f_ == 0. ) f = 1. / 298.257;
  else f = f_;
  if ( e2_ == 0. ) e2 = 6.6943799013e-3;
  else e2 = e2_;
  gllh[ 0 ] = modulo( degrees( atan2( xyz[ 1 ], xyz[ 0 ] ) ) + 180., 360. ) - 180.;
  gllh[ 1 ] = degrees( atan2( xyz[ 2 ], sqrt( xyz[ 0 ] * xyz[ 0 ] + xyz[ 1 ] * xyz[ 1 ] ) ) );
  gllh[ 2 ] = sqrt( xyz[ 0 ] * xyz[ 0 ] + xyz[ 1 ] * xyz[ 1 ] + xyz[ 2 ] * xyz[ 1 ] ) - a * sqrt( 1. - f );
  if ( iterations > 0 )
  {
    phi = radians( gllh[ 1 ] );
    for ( i = 0; i < iterations; i++ )
    {
      n = a / sqrt( 1. - e2 * ( sin( phi ) * sin( phi ) ) );
      gllh[ 2 ] = ( sqrt( xyz[ 0 ] * xyz[ 0 ] + xyz[ 1 ] * xyz[ 1 ] ) / cos( phi ) ) - n;
      phi = atan( xyz[ 2 ] / ( sqrt( xyz[ 0 ] * xyz[ 0 ] + xyz[ 1 ] * xyz[ 1 ] ) *
          ( 1. - e2 * ( n / ( n + gllh[ 2 ] ) ) ) ) );
    }
    gllh[ 1 ] = degrees( phi );
  }
  return ( gllh );
}

//*****************************************************************************

p_double_3 geo_llh_to_xyz( double gllh[ 3 ], double a_, double f_, double e2_ )
// # default Earth ellipticity definition (a,f) is WGS (1984)
// Note that longitude is defined as positive towards east, just like RA
{
  double a, f, e2, lamda, phi, n;
  p_double_3 xyz = malloc( 3 * sizeof( double ) );
  if ( a_ == 0. ) a = 6378137.;
  else a = a_;
  if ( f_ == 0. ) f = 1. / 298.257;
  else f = f_;
  if ( e2_ == 0. ) e2 = 6.6943799013e-3;
  else e2 = e2_;
  lamda = radians( gllh[ 0 ] );
  phi = radians( gllh[ 1 ] );
  n = a / sqrt( 1. - e2 * ( sin( phi ) * sin( phi ) ) );
  xyz[ 0 ] = ( n + gllh[ 2 ] ) * cos( phi ) * cos( lamda );
  xyz[ 1 ] = ( n + gllh[ 2 ] ) * cos( phi ) * sin( lamda );
  xyz[ 2 ] = ( n * ( 1. - e2 ) + gllh[ 2 ] ) * sin( phi );
  return ( xyz );
}

//*****************************************************************************

p_double_2 calculate_hour_angles_at_elevation_limit( double lat, double dec, double elevation_limit )
{
  double a, b, c, A, AA;
  p_double_2 ha = malloc( 2 * sizeof( double ) );
  // check for circumpolar sources
  if ( ( dec + lat >= 90. ) || ( dec + lat <= - 90. ) )
    ha[ 1 ] = 180.;
  // check for non-visible sources
  else if ( ( dec - lat >= 90. ) || ( dec - lat <= - 90. ) )
    ha[ 1 ] = 0.;
  else
  {
    a = radians( 90. - elevation_limit );
    b = radians( 90. - dec );
    c = radians( 90. - lat );
    AA = ( cos( a ) - cos( b ) * cos( c ) ) / (  sin( b ) * sin( c ) );
    if ( fabs( AA ) > 1. )
      AA = sign( AA ) * 1.;
    A = acos( AA );
    ha[ 1 ] = degrees( A );
  }
  ha[ 0 ] = - ha[ 1 ];
  return ( ha );
}

//*****************************************************************************

p_double_4 time_to_dhms( double time )
{
  double res;
  p_double_4 dhms = malloc( 4 * sizeof( double ) );
  res = fabs( time );
  dhms[ 0 ] = floor( res );
  res = 24. * ( time - dhms[ 0 ] );
  if ( time < 0. )
    dhms[ 0 ] = - dhms[ 0 ];
  dhms[ 1 ] = floor( res );
  res = 60. * ( res - dhms[ 1 ] );
  dhms[ 2 ] = floor( res );
  dhms[ 3 ] = 60. * ( res - dhms[ 2 ] );
  return ( dhms );
}

//*****************************************************************************

p_double_2 calculate_local_sky_position( double gllh[ 3 ], double radec[ 2 ], double gst )
// gst in degrees
{
  double ha, a, b, c, A, B, aa, BB;
  p_double_2 za = malloc( 2 * sizeof( double ) );
  if ( ( fabs( gllh[ 1 ] ) == 90. ) || ( fabs( radec[ 1 ] ) == 90. ) )
  {
    za[ 0 ] = modulo( fabs( radec[ 1 ] - gllh[ 1 ] ), 180. );
    if ( ( radec[ 1 ] == 90. ) || ( gllh[ 1 ] == - 90 ) )
      za[ 1 ] = 0.;
    else
      za[ 1 ] = 180.;
  }
  else
  {
    ha = modulo( gst - radec[ 0 ] + gllh[ 0 ], 360. );
    A = radians( - ha );
    b = radians( 90. - radec[ 1 ] );
    c = radians( 90. - gllh[ 1 ] );
    aa = sin( b ) * cos( A ) * sin( c ) + cos( b ) * cos( c );
    if ( fabs( aa ) > 1. )
      aa = sign( aa ) * 1.;
    a = acos( aa );
    BB = ( cos( b ) - cos( a ) * cos( c ) ) / ( sin( a ) * sin( c ) );
    if ( fabs( BB ) > 1. )
      BB = sign( BB ) * 1.;
    B = acos( BB );
    za[ 0 ] = degrees( a );
    za[ 1 ] = degrees( B );
    if ( sin( A ) < 0. )
      za[ 1 ] = 360. - za[ 1 ];
  }
  return ( za );
}

//*****************************************************************************

p_double_4 calculate_pierce_point( double axyz[ 3 ], double radec[ 2 ], double gst,
    double height, int iterations )
// height in meters
{
  int i,j;
  double *agllh, alon, alat, rot[ 3 ][ 3 ], *azz, aza, aaz, l2axyz, lsdxyz[ 3 ];
  double B, C, l2pxyz, lsxyz, pxyz[ 3 ], *pgllh, *ra, sdxyz[ 3 ];
  p_double_4 pp = malloc( 4 * sizeof( double ) );

  // initialize some variables
  agllh = xyz_to_geo_llh( axyz, iterations, 0., 0., 0. );
  alon = radians( agllh[ 0 ] );
  alat = radians( agllh[ 1 ] );
  if ( agllh[ 2 ] > height )
    return ( NULL );
  rot[ 0 ][ 0 ] = - sin( alon );
  rot[ 0 ][ 1 ] = - cos( alon ) * sin( alat );
  rot[ 0 ][ 2 ] = cos( alon ) * cos( alat );
  rot[ 1 ][ 0 ] = cos( alon );
  rot[ 1 ][ 1 ] = - sin( alon ) * sin( alat );
  rot[ 1 ][ 2 ] = sin( alon ) * cos( alat );
  rot[ 2 ][ 0 ] = 0.;
  rot[ 2 ][ 1 ] = cos( alat );
  rot[ 2 ][ 2 ] = sin( alat );

  // determine source unit vector
  azz = calculate_local_sky_position( agllh, radec, gst );
  aza = radians( azz[ 0 ] );
  aaz = radians( azz[ 1 ] );
  l2axyz = 0.;
  for ( i = 0; i < 3; i++ )
    l2axyz = l2axyz + axyz[ i ] * axyz[ i ];
  lsdxyz[ 0 ] = sin( aza ) * sin( aaz );
  lsdxyz[ 1 ] = sin( aza ) * cos( aaz );
  lsdxyz[ 2 ] = cos( aza );
  for ( i = 0; i < 3; i++ )
  {
    sdxyz[ i ] = 0.;
    for ( j = 0; j < 3; j++ )
      sdxyz[ i ] = sdxyz[ i ] + rot[ i ][ j ] * lsdxyz[ j ];
  }

  // determine xyz coordinates of pierce point through vector algebra
  pgllh = ( double * )NULL;
  B = 0.;
  for ( i = 0; i < 3; i++ )
    B = B + 2. * axyz[ i ] * sdxyz[ i ];
  l2pxyz = ( sqrt( l2axyz ) + ( height - agllh[ 2 ] ) );
  l2pxyz = l2pxyz * l2pxyz;
  for ( i = 0; i < iterations; i++ )
  {
    C = l2axyz - l2pxyz;
    lsxyz = ( sqrt( B * B - 4. * C ) - B ) / 2.;
    l2pxyz = 0.;
    for ( j = 0; j < 3; j++ )
    {
      pxyz[ j ] = axyz[ j ] + lsxyz * sdxyz[ j ];
      l2pxyz = l2pxyz + pxyz[ j ] * pxyz[ j ];
    }
    if ( pgllh != NULL )
      free( pgllh );
    pgllh = xyz_to_geo_llh( pxyz, iterations, 0., 0., 0. );
    l2pxyz = ( sqrt( l2pxyz ) + ( height - pgllh[ 2 ] ) );
    l2pxyz = l2pxyz * l2pxyz;
  }
  C = l2axyz - l2pxyz;
  lsxyz = ( sqrt( B * B - 4. * C ) - B ) / 2.;
  for ( i = 0; i < 3; i++ )
  {
    pp[ i ] = axyz[ i ] + lsxyz * sdxyz[ i ];
  }

  // determine zenith angle at pierce point
  if ( pgllh != NULL )
    free( pgllh );
  pgllh = xyz_to_geo_llh( pxyz, iterations, 0., 0., 0. );
  ra = calculate_angular_separation( agllh, pgllh );
  pp[ 3 ] = azz[ 0 ] - ra[ 0 ];

  // cleanup
  free( agllh );
  free( azz );
  free( pgllh );
  free( ra );

  return ( pp );
}

//*****************************************************************************

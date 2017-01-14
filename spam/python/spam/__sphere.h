//*****************************************************************************

#ifndef __sphere__h
#define __sphere__h

//*****************************************************************************

// Python header files
#include <Python.h>

// user header files
#include "__types.h"
#include "__acalc.h"

//*****************************************************************************

extern p_double_2 hmsdms_to_degdeg( double hmsdms[ 6 ] );

extern p_double_6 degdeg_to_hmsdms( double degdeg[ 2 ] );

extern p_double_6 degdeg_to_dmsdms( double degdeg[ 2 ] );

extern p_double_2 calculate_angular_separation( double degdeg1[ 2 ], double degdeg2[ 2 ] );

extern p_double_2 calculate_offset_position( double radec[ 2 ], double radius, double angle );

extern p_double_3 xyz_to_llr( double xyz[ 3 ] );

extern p_double_3 xyz_to_geo_llh( double xyz[ 3 ], int iterations, double a, double f, double e2 );

extern p_double_3 geo_llh_to_xyz( double gllh[ 3 ], double a, double f, double e2 );

extern p_double_2 calculate_hour_angles_at_elevation_limit( double lat, double dec, double elevation_limit );

extern p_double_4 time_to_dhms( double time );

extern p_double_2 calculate_local_sky_position( double gllh[ 3 ], double radec[ 2 ], double gst );

extern p_double_4 calculate_pierce_point( double axyz[ 3 ], double radec[ 2 ], double gst,
    double height, int iterations );

//*****************************************************************************

#endif

//*****************************************************************************

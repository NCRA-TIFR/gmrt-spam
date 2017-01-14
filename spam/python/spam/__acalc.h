//*****************************************************************************

#ifndef __acalc__h
#define __acalc__h

//*****************************************************************************

// Python header files
#include <Python.h>

// standard C header files
//#include <stdio.h>
//#include <string.h>
//#include <math.h>

// user header files
#include "__types.h"

//*****************************************************************************

extern double sign( double x );

extern double degrees( double rad );

extern double radians( double deg );

extern double modulo( double a, double b );

extern long int factorial( int n );

extern long int binomial( int n, int k );

extern p_double_2 complex_to_r_phi( double c[ 2 ] );

extern p_double_2 r_phi_to_complex( double rp[ 2 ] );

extern double maximum( p_double_array da );

extern double minimum( p_double_array da );

extern double median( p_double_array da );

extern double sum( p_double_array da );

extern p_double_array square( p_double_array da );

extern p_double_array multiply( p_double_array da, double factor );

//*****************************************************************************

#endif

//*****************************************************************************

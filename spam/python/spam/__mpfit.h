//*****************************************************************************

#ifndef __mpfit__h
#define __mpfit__h

//*****************************************************************************

// Python header files
#include <Python.h>

// user header files
#include "__types.h"
#include "__acalc.h"

//*****************************************************************************

extern double enorm( p_double_array vec, double rgiant, double rdwarf );

extern p_double_array calc_covar( p_double_array rr_table, p_long_array ipvt_table, double tol );

extern p_double_array qrfac( p_double_array a, int pivot, double machep, double rgiant, double rdwarf );

extern p_double_array qrsolv( p_double_array r, p_long_array ipvt, p_double_array diag,
    p_double_array qtb, p_double_array sdiag );

extern p_double_array lmpar( p_double_array r, p_long_array ipvt, p_double_array diag,
    p_double_array qtb, p_double_array x, p_double_array sdiag,
    double delta, double par, double dwarf, double rgiant, double rdwarf );

//*****************************************************************************

#endif

//*****************************************************************************

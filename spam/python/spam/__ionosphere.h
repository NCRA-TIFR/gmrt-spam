//*****************************************************************************

#ifndef __ionosphere__h
#define __ionosphere__h

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

extern p_double_array calculate_phase_gradients( p_double_array phase_table );

extern p_double_2 estimate_phase_gradient( p_double_array X_table,
    p_double_array ref_X_table, p_double_array phase_table, p_double_array error_table );

//*****************************************************************************

#endif

//*****************************************************************************

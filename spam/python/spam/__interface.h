//*****************************************************************************

#ifndef __interface__h
#define __interface__h

//*****************************************************************************

// Python header files
#include <Python.h>

// user header files
#include "__types.h"

//*****************************************************************************

extern int convert_float_list_to_double_array( PyObject *input, double *ptr, int size );

extern int recursive_convert_float_list_to_double_array( PyObject *fl, p_double_array da, long dim_i, long *data_i );

extern p_double_array convert_multi_float_list_to_double_array( PyObject *fl );

extern PyObject *recursive_convert_double_array_to_float_list( p_double_array da, long dim_i, long *data_i );

extern PyObject *convert_multi_double_array_to_float_list( p_double_array da );

extern int recursive_convert_integer_list_to_long_array( PyObject *fl, p_long_array la, long dim_i, long *data_i );

extern p_long_array convert_multi_integer_list_to_long_array( PyObject *fl );

extern PyObject *recursive_convert_long_array_to_integer_list( p_long_array la, long dim_i, long *data_i );

extern PyObject *convert_multi_long_array_to_integer_list( p_long_array la );

//*****************************************************************************

#endif

//*****************************************************************************

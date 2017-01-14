//*****************************************************************************

#ifndef __types__h
#define __types__h

//*****************************************************************************

// Python header files
#include <Python.h>

#define __TYPES_MAX_DIMS ( 40L )

//*****************************************************************************

typedef double *p_double_2;
typedef double *p_double_3;
typedef double *p_double_4;
typedef double *p_double_5;
typedef double *p_double_6;

typedef struct double_array
{
  long num_dims;
  long *dims;
  double *data;
} *p_double_array;

typedef struct long_array
{
  long num_dims;
  long *dims;
  long *data;
} *p_long_array;

//*****************************************************************************

extern p_double_array create_double_array( long num_dims, long *dims, double *data );

extern void free_double_array( p_double_array da );

extern int check_double_array( p_double_array da, long num_dims, long *dims );

extern long double_array_dims_to_index( p_double_array da, long *dims );

extern p_long_array create_long_array( long num_dims, long *dims, long *data );

extern void free_long_array( p_long_array la );

extern int check_long_array( p_long_array la, long num_dims, long *dims );

extern long long_array_dims_to_index( p_long_array la, long *dims );

//*****************************************************************************

#endif

//*****************************************************************************

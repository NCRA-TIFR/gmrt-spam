//*****************************************************************************

// module name
%module _sphere

// include other interface files
%import "__interface.i"

// include files scanned for wrapping
%include "__sphere.h"

// included in C wrapper
%{
#include "__types.h"
#include "__interface.h"
#include "__sphere.h"
%}

//*****************************************************************************

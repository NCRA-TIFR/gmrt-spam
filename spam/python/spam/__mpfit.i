//*****************************************************************************

// module name
%module _mpfit

// include other interface files
%import "__interface.i"

// include files scanned for wrapping
%include "__mpfit.h"

// included in C wrapper
%{
#include "__types.h"
#include "__interface.h"
#include "__mpfit.h"
%}

//*****************************************************************************

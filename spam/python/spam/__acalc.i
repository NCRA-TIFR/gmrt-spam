//*****************************************************************************

// module name
%module _acalc

// import other interface files
%import "__interface.i"

// include files scanned for wrapping
%include "__acalc.h"

// included in C wrapper
%{
#include "__types.h"
#include "__interface.h"
#include "__acalc.h"
%}

//*****************************************************************************

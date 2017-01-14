//*****************************************************************************

// module name
%module _ionosphere

// import other interface files
%import "__interface.i"

// include files scanned for wrapping
%include "__ionosphere.h"

// included in C wrapper
%{
#include "__types.h"
#include "__interface.h"
#include "__ionosphere.h"
%}

//*****************************************************************************

%module weaklens
%include cstring.i
%include cpointer.i
%pointer_class(double,dp)
%apply double& INOUT { double& a };
%feature("autodoc", 1);
%{
    #define SWIG_FILE_WITH_INIT
    #include "weaklens.h"
%}

%include "weaklens.h"

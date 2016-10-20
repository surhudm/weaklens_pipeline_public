%module tomography
%include cstring.i
%include cpointer.i
%pointer_class(double,dp)
%apply double& INOUT { double& a };
%feature("autodoc", 1);
%{
    #define SWIG_FILE_WITH_INIT
    #include "tomography.h"
%}
%include numpy.i
%init %{
import_array();
%}
%apply (float* IN_ARRAY1, int DIM1) {(float* pofz, int zbins)};
%apply (double* IN_ARRAY1, int DIM1) {(double* pofz, int zbins)};

%include "tomography.h"

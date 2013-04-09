%module phasefit
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}

%apply (double * INPLACE_ARRAY1, int DIM1){(double *time, int lt), (double *mag, int lm)};
%apply (double * INPLACE_ARRAY1, int DIM1){(double *magbin, int lbin)};
%apply (double * INPLACE_ARRAY2, int DIM1, int DIM2){(double *color, int nbin, int ntran)};

%{
#include <iostream.h>
#include "phasefit.h"
%}

%include "phasefit.h"
%clear (double *time, int lt), (double *mag, int lm);
%clear (double *magbin, int lbin);
%clear (double *color, int nbin, int ntran);

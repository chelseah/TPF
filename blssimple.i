%module blssimple
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}

%apply (double * INPLACE_ARRAY1, int DIM1){(double *time, int lt), (double *mag, int lm)};

%{
#include <iostream>
#include <math.h>
#include "blssimple.h"
%}

%include "blssimple.h"
%clear (double *time, int lt), (double *mag, int lm);

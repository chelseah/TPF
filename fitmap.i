%module fitmap
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}

%apply (double * INPLACE_ARRAY2, int DIM1, int DIM2){(double *cmap, int nx, int ny)};

%{
#include <iostream>
#include <math.h>
#include "fitmap.h"
%}

%include "fitmap.h"
%clear (double *cmap, int nx, int ny);

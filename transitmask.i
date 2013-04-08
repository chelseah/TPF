%module transitmask
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}

%{
#include <iostream>
#include <math.h>
#include "maskarray.h"
#include "transitmask.h"
%}

%include "maskarray.h"
%include "transitmask.h"

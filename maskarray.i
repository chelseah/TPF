%module maskarray
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}


%{
#include <iostream>
#include "maskarray.h"
%}

%include "maskarray.h"

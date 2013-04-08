# SWIG
INPHASE = phasefit.i
INFIT = maskarray.i transitmask.i fitmap.i
SoPHASE = phasefit.cc phasefit.h
SoFIT = maskarray.cc maskarray.h transitmask.cc transitmask.h fitmap.cc fitmap.h  
INTERFACES = $(INPHASE) $(INFIT)
WRAPPERS   = $(INTERFACES:.i=_wrap.cxx)
PROXIES    = $(INTERFACES:.i=.py      )
Sources  = $(SoPHASE) $(SoPHASE)

# Default target: build the tests

.PHONY : all
all: $(WRAPPERS) $(Sources)	
	./setup.py build_ext -i
# Test target: run the tests

# Rule: %.i -> %_wrap.cc
%_wrap.cxx: %.i %.h ./numpy.i
	swig -c++ -python $<

# Clean target
.PHONY : clean
clean:
	$(RM) -r build
	$(RM) *.so *.pyc *_wrap.h 
	$(RM) -r *.dSYM  
	$(RM) $(WRAPPERS)
	$(RM) $(PROXIES)
	$(RM) *.o           
	$(RM) *.out        
	$(RM) *~                              
                       

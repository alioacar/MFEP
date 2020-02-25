FC = gfortran
SRCDIR=$(shell pwd)/src
SAMPLINGDIR=$(shell pwd)/sampling

MODULES = $(SRCDIR)/vec_solve.mod $(SRCDIR)/functions.mod $(SRCDIR)/circulardatafinder.mod \
$(SRCDIR)/veclengthfinder.mod $(SRCDIR)/steepestpoint.mod
OBJECTS = $(SRCDIR)/vec_solve.o $(SRCDIR)/functions.o $(SRCDIR)/circulardatafinder.o \
$(SRCDIR)/veclengthfinder.o $(SRCDIR)/steepestpoint.o $(SRCDIR)/main.o
.PHONY: clean

main.exe: $(MODULES) $(OBJECTS)
	$(FC) $(OBJECTS) -o $(SRCDIR)/main.exe
	ln -s $(SRCDIR)/main.exe $(shell pwd)/mfep
	rm $(shell pwd)/*.mod


$(SRCDIR)/%.o : $(SRCDIR)/%.f90
	$(FC) -c -o $@ $<
$(SRCDIR)/%.mod : $(SRCDIR)/%.f90
	$(FC) -c -o $@ $<
clean:
	rm -f $(OBJECTS) $(MODULES) main.exe mfep
	rm -f $(SRCDIR)/main.exe
	rm -rf $(SAMPLINGDIR)/*

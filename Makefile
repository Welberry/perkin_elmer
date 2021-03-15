include ../../modules/makefiles/definitions

INSTALL_DIR = $(TOOLS_DIR)/bin

# Additional fortran compilation flags
FFLAGS =  -O

# Additional linker flags
LFLAGS =  -lmodules -static-libcxa

# Targets ...

TARGETS = pe2mar pestats pe_dark_stats darkcleaner dudpixelfinder darkfilter

all: $(TARGETS)

darkfilter: darkfilter.o cluster_functions.o
	$(FC) $(FFLAGS) -o $@ $^ $(LFLAGS)

darkcleaner: darkcleaner.o rannum_module.o
	$(FC) $(FFLAGS) -o $@ $^ $(LFLAGS)

pe2mar: pe2mar.o
	$(FC) $(FFLAGS) -o $@ $^ $(LFLAGS)

dudpixelfinder: dudpixelfinder.o
	$(FC) $(FFLAGS) -o $@ $^ $(LFLAGS)

pestats: pestats.o
	$(FC) $(FFLAGS) -o $@ $^ $(LFLAGS)

pemax: pemax.o
	$(FC) $(FFLAGS) -o $@ $^ $(LFLAGS)

pe_dark_stats: pe_dark_stats.o
	$(FC) $(FFLAGS) -o $@ $^ $(LFLAGS)

install: $(TARGETS)
	cp $(TARGETS) $(INSTALL_DIR)

clean:
	rm -f *.$(OBJSUFFIX) *.$(MODSUFFIX) $(TARGETS)


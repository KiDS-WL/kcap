include ${COSMOSIS_SRC_DIR}/config/compilers.mk

all: halofit_module.so

test:
	@echo "Alas, halofit has no tests"

clean:
	rm -f halofit_module.so *.o *.mod
	rm -rf halofit_module.so.dSYM/

halofit_module.so: halofit.F90 halofit_module.F90
	$(FC) $(FFLAGS) -shared -o $@ $+ $(LDFLAGS) -lcosmosis_fortran

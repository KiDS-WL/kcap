include ${COSMOSIS_SRC_DIR}/cosmosis/compilers.mk

USER_CFLAGS= -I../limber -I${GSL_INC}
USER_LDFLAGS= -L../limber -llimber -L${GSL_LIB} -lgsl -lgslcblas -lcosmosis -lm

all: interface.so

interface.so: shear_shear.c interface.c
	$(CC) $(CFLAGS) -shared -o $@ $+ $(LDFLAGS)

clean:
	rm -rf interface.so interface.so.dSYM

test:
	@echo Alas, shear/spectra has no tests


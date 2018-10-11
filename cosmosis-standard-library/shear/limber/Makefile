include ${COSMOSIS_SRC_DIR}/config/compilers.mk

ifeq ($(GSL_INC)$(GSL_LIB),) 
$(error To compile the limber code you you must activate the UPS package gsl or set the environment variables GSL_INC and GSL_LIB)
endif

USER_CFLAGS=-I ${GSL_INC} 
USER_LDFLAGS=-L${GSL_LIB} -lgsl -lgslcblas -lcosmosis -lm

# AR = ar rc

OBJS = limber.o interp2d.o utils.o
LIB = liblimber.a
EXEC = test_limber

default:$(LIB) $(EXEC)
all: $(EXEC)

$(EXEC): $(LIB) test_limber.o
	$(CC) $(CFLAGS) -o $@ $+ $(LDFLAGS) -L. -llimber

$(LIB): $(OBJS)
	$(AR) rc $@ $^


.PHONY : clean

clean: 
	rm -rf *.o $(LIB) $(LIBNAME) $(LIBNAME).dSYM $(EXEC) $(EXEC).dSYM $(SHLIBNAME) $(SHLIBNAME).dSYM

test:
	@echo Alas, shear/limber has no tests


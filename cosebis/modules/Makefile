#-----------------------------------
#--- objects to build
#-----------------------------------
USER_LDFLAGS += -L${GSL_LIB} -lgslcblas -lgsl
OBJS     = matrix.o function_cosebis.o ticker.o tostring.o tolower.o timestamp.o \
		   globaldef.o Integrate.o COSEBIs.o WnLog.o TpnLogRadian.o nrutil.o

LIBNAME  = libpslib.a

#-----------------------------------
#--- libraries
#-----------------------------------

# make objects (components) of library
%.o : %.cc 
	$(CXX) $(CXXFLAGS) -c $<
#%.o: %.cc
#	$(CXX) -o $@ -c $(CXXFLAGS) $*.cc


# make library
all: $(LIBNAME)

$(LIBNAME): $(OBJS)
	$(AR) rc $@ $^

#%.a:
#	rm -f $@; ar -r $@ $(OBJECTS)

.PHONY: clean
clean:
	rm -f $(LIBNAME)
	rm -f *.o
#-----------------------------------
#--- objects to build
#-----------------------------------
USER_LDFLAGS += -L${GSL_LIB} -lgslcblas -lgsl
OBJS     = matrix.o function_cosebis.o ticker.o tostring.o tolower.o timestamp.o globaldef.o Integrate.o COSEBIs.o WnLog.o TpnLogRadian.o nrutil.o

LIBNAME  = libpslib.a

HEADER   = $(OBJS:.o=.h) 

#-----------------------------------
#--- libraries
#-----------------------------------

# make objects (components) of library
%.o : %.cc 
	$(CXX) $(USER_LDFLAGS) $(CXXFLAGS) -c $<
#%.o: %.cc
#	$(CXX) -o $@ -c $(CXXFLAGS) $*.cc


# make library
all: $(HEADER) $(OBJS) 
	rm -f $(LIBNAME)
	ar rv $(LIBNAME) $(OBJS)
	ranlib $(LIBNAME)


#%.a:
#	rm -f $@; ar -r $@ $(OBJECTS)


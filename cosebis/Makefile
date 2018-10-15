# If you already have your own Makefile you can
# replace all of this, but you need to keep this line
# at the top:
include ${COSMOSIS_SRC_DIR}/config/compilers.mk

MODULE = ./modules

USER_CXXFLAGS += -I${GSL_INC}
USER_LDFLAGS += -lcosmosis -L${GSL_LIB} -lgslcblas -lgsl -lm -L$(MODULE) -lpslib -lstdc++ -lz -Wno-deprecated -I$(MODULE)
#fgas
#USER_CXXFLAGS += -I${GSL_INC}
#USER_LDFLAGS += -lcosmosis -L${GSL_LIB}  -lgsl -lgslcblas
# RUN "make all" to make all binaries or "make [name]" to make [name]

#-----------------------------------
#--- linker
#-----------------------------------

MYLIB     = libpslib.a

#-----------------------------------
#--- rules
#-----------------------------------
all: libcosebis_cl.so libcosebis_2pcfs.so libcosebis_likelihood.so lib2pcfs_likelihood.so

libcosebis_cl.so: COSEBIs_Cl_cosmosis_interface.cc $(MYLIB)
	$(CXX) $(CXXFLAGS) COSEBIs_Cl_cosmosis_interface.cc -shared -o libcosebis_cl.so $(LDFLAGS) $(USER_LDFLAGS)

libcosebis_2pcfs.so: COSEBIs_2PCFs_cosmosis_interface.cc $(MYLIB)
	$(CXX) $(CXXFLAGS) COSEBIs_2PCFs_cosmosis_interface.cc -shared -o libcosebis_2pcfs.so $(LDFLAGS) $(USER_LDFLAGS)

libcosebis_likelihood.so: COSEBIs_likelihood_cosmosis_interface.cc $(MYLIB)
	$(CXX) $(CXXFLAGS) COSEBIs_likelihood_cosmosis_interface.cc -shared -o libcosebis_likelihood.so $(LDFLAGS) $(USER_LDFLAGS)

lib2pcfs_likelihood.so: 2pcfs_likelihood_cosmosis_interface.cc $(MYLIB)
	$(CXX) $(CXXFLAGS) 2pcfs_likelihood_cosmosis_interface.cc -shared -o lib2pcfs_likelihood.so $(LDFLAGS) $(USER_LDFLAGS)

#fgas.so: fgas_cosmosis.cpp src/libclusters.a
#	$(CXX) $(CXXFLAGS) fgas_cosmosis.cpp -shared -o fgas.so $(LDFLAGS) src/libclusters.a


$(MYLIB):
	$(MAKE) -C $(MODULE) all
#src/libclusters.a:
#	$(MAKE) -C src

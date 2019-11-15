# If you already have your own Makefile you can
# replace all of this, but you need to keep this line
# at the top:
include ${COSMOSIS_SRC_DIR}/config/compilers.mk

MODULE = ./modules

USER_CXXFLAGS += -I${GSL_INC}
USER_LDFLAGS += -lcosmosis -L${GSL_LIB}  -lgsl -lgslcblas -lm -L$(MODULE) -lcolib -lstdc++ -lz -Wno-deprecated -I$(MODULE)

#-----------------------------------
#--- linker
#-----------------------------------

COLIB     = libcolib.a

#-----------------------------------
#--- rules
#-----------------------------------
all: $(MODULE)/$(COLIB) libcosebis_cl.so libcosebis_2pcfs.so libcosebis_likelihood.so \
 lib2pcfs_likelihood.so libbandpower_likelihood.so
	# lib2pcfs_cl.so
# libbandpower_theory.so libbandpower_real_theory.so

# tells how to make my library (another makefile)

$(MODULE)/$(COLIB)::
	cd $(MODULE) && $(MAKE)

libcosebis_cl.so: COSEBIs_Cl_cosmosis_interface.cc
	$(CXX) $(CXXFLAGS) COSEBIs_Cl_cosmosis_interface.cc -shared -o libcosebis_cl.so $(LDFLAGS) $(USER_LDFLAGS)

libcosebis_2pcfs.so: COSEBIs_2PCFs_cosmosis_interface.cc
	$(CXX) $(CXXFLAGS) COSEBIs_2PCFs_cosmosis_interface.cc -shared -o libcosebis_2pcfs.so $(LDFLAGS) $(USER_LDFLAGS)

libcosebis_likelihood.so: COSEBIs_likelihood_cosmosis_interface.cc
	$(CXX) $(CXXFLAGS) COSEBIs_likelihood_cosmosis_interface.cc -shared -o libcosebis_likelihood.so $(LDFLAGS) $(USER_LDFLAGS)

lib2pcfs_likelihood.so: 2pcfs_likelihood_cosmosis_interface.cc
	$(CXX) $(CXXFLAGS) 2pcfs_likelihood_cosmosis_interface.cc -shared -o lib2pcfs_likelihood.so $(LDFLAGS) $(USER_LDFLAGS)

lib2pcfs_cl.so: 2PCFs_Cl_cosmosis_interface.cc
	$(CXX) $(CXXFLAGS) 2PCFs_Cl_cosmosis_interface.cc -shared -o lib2pcfs_cl.so $(LDFLAGS) $(USER_LDFLAGS)

libbandpower_theory.so: BandPower_interface.cc
	$(CXX) $(CXXFLAGS) BandPower_interface.cc -shared -o libbandpower_theory.so $(LDFLAGS) $(USER_LDFLAGS)

libbandpower_real_theory.so: BandPower_interface_real.cc
	$(CXX) $(CXXFLAGS) BandPower_interface_real.cc -shared -o libbandpower_real_theory.so $(LDFLAGS) $(USER_LDFLAGS)


libbandpower_likelihood.so: BandPower_likelihood_interface.cc
	$(CXX) $(CXXFLAGS) BandPower_likelihood_interface.cc -shared -o libbandpower_likelihood.so $(LDFLAGS) $(USER_LDFLAGS)


# Clean up
.PHONY: clean
clean:
	cd $(MODULE) && $(MAKE) clean
	rm -f libcosebis_cl.so
	rm -f libcosebis_2pcfs.so
	rm -f libcosebis_likelihood.so
	rm -f lib2pcfs_likelihood.so
	rm -f libbandpower_theory.so
	rm -f libbandpower_real_theory.so
	rm -f libbandpower_likelihood.so
	test -n "./" && rm -rf ./libcosebis_cl.so.dSYM/
	test -n "./" && rm -rf ./libcosebis_2pcfs.so.dSYM/
	test -n "./" && rm -rf ./libcosebis_likelihood.so.dSYM/
	test -n "./" && rm -rf ./lib2pcfs_likelihood.so.dSYM/
	test -n "./" && rm -rf ./libbandpower_theory.so.dSYM/
	test -n "./" && rm -rf ./libbandpower_real_theory.so.dSYM/
	test -n "./" && rm -rf ./libbandpower_likelihood.so.dSYM/
	

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
all: $(MODULE)/$(COLIB) libcosebis_2pcfs.so libcosebis.so libcosebis_cov.so libbandpower.so libxipm_binned.so 

# tells how to make my library (another makefile)

$(MODULE)/$(COLIB)::
	cd $(MODULE) && $(MAKE)


libcosebis_2pcfs.so: COSEBIs_2PCFs_cosmosis_interface.cc
	$(CXX) $(CXXFLAGS) COSEBIs_2PCFs_cosmosis_interface.cc -shared -o libcosebis_2pcfs.so $(LDFLAGS) $(USER_LDFLAGS)

libcosebis.so: COSEBIs_cosmosis_interface.cc
	$(CXX) $(CXXFLAGS) COSEBIs_cosmosis_interface.cc -shared -o libcosebis.so $(LDFLAGS) $(USER_LDFLAGS)

libcosebis_cov.so: COSEBIs_covariance_cosmosis_interface.cc
	$(CXX) $(CXXFLAGS) COSEBIs_covariance_cosmosis_interface.cc -shared -o libcosebis_cov.so $(LDFLAGS) $(USER_LDFLAGS)

libbandpower.so: BandPower_interface.cc
	$(CXX) $(CXXFLAGS) BandPower_interface.cc -shared -o libbandpower.so $(LDFLAGS) $(USER_LDFLAGS)

libxipm_binned.so: xipm_binned.cc
	$(CXX) $(CXXFLAGS) xipm_binned.cc -shared -o libxipm_binned.so $(LDFLAGS) $(USER_LDFLAGS)


# Clean up
.PHONY: clean
clean:
	cd $(MODULE) && $(MAKE) clean
	rm -f libcosebis.so
	rm -f libcosebis_cov.so
	rm -f libcosebis_2pcfs.so
	rm -f libbandpower.so
	rm -f libxipm_binned.so
	test -n "./" && rm -rf ./libcosebis_2pcfs.so.dSYM/
	test -n "./" && rm -rf ./libcosebis.so.dSYM/
	test -n "./" && rm -rf ./libcosebis_cov.so.dSYM/
	test -n "./" && rm -rf ./libbandpower.so.dSYM/
	test -n "./" && rm -rf ./libxipm_binned.so.dSYM/

	
	

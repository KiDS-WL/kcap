include ${COSMOSIS_SRC_DIR}/config/compilers.mk


USER_CXXFLAGS+= -I${GSL_INC}
USER_LDFLAGS+=  -lcosmosis ${LAPACK_LINK} -L${GSL_LIB} -lgslcblas -lgsl 


all: libjla.a jla.so

libjla.a: src/jla.o
	$(AR) $(ARFLAGS) $@ $^

jla.so: interface.cc
	$(CXX) $(CXXFLAGS) -shared -o $@ $^ -L. -ljla $(LDFLAGS) -lstdc++
#test_jla: src/test.o libjla.a
#	$(CXX) -o $@ $^ $(CXXFLAGS) $(LFLAGS) $(CLASS_LFLAGS)

#src/test.o: src/test.cc
#	$(CXX) -c -o $@ $^  $(CXXFLAGS) $(CLASS_CXXFLAGS)

clean:
	rm -f src/*.o
	rm -f libjla.a
	rm -rf jla.so jla.so.dSYM
test:
	@echo Alas, supernovae/jla_v3 has no tests


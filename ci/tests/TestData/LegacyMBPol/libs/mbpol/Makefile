CXX = g++
CXXFLAGS = -O2

MBPOL_CPP = mbpol/libmbpol.cpp \
            mbpol/poly-2b-v6x.cpp mbpol/poly-2b-v6x-nogrd.cpp \
            mbpol/poly-3b-v2x.cpp mbpol/poly-3b-v2x-nogrd.cpp \
            mbpol/ttm4-es.cpp mbpol/ttm4-smear.cpp \
						mbpol/gammq.cpp mbpol/tang-toennies.cpp mbpol/ps.cpp \
            mbpol/x2b-dispersion.cpp mbpol/x2b-v9x.cpp \
						mbpol/x3b-v2x.cpp mbpol/mbpol.cpp

MBPOL_OBJ = mbpol/libmbpol.o \
            mbpol/poly-2b-v6x.o mbpol/poly-2b-v6x-nogrd.o \
            mbpol/poly-3b-v2x.o mbpol/poly-3b-v2x-nogrd.o \
            mbpol/ttm4-es.o mbpol/ttm4-smear.o mbpol/gammq.o \
						mbpol/tang-toennies.o mbpol/ps.o \
            mbpol/x2b-dispersion.o mbpol/x2b-v9x.o mbpol/x3b-v2x.o mbpol/mbpol.o

all: libmbpol.so

$(MBPOL_OBJ) :
	$(CXX) $(CXXFLAGS) -fPIC -c $(MBPOL_CPP)
	mv *.o mbpol

libmbpol.so: $(MBPOL_OBJ)
	$(CXX) -shared -o libmbpol.so -fPIC $(MBPOL_OBJ)
	mv libmbpol.so ../libmbpol.so

clean:
	rm -f *.so *.o

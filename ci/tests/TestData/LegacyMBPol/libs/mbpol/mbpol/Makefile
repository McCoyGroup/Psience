CXX = icc
CXXFLAGS = -O0 -fopenmp

AR = /usr/bin/ar

MBPOL_CPP = libmbpol.cpp \
            poly-2b-v6x.cpp poly-2b-v6x-nogrd.cpp \
            poly-3b-v2x.cpp poly-3b-v2x-nogrd.cpp \
            ttm4-es.cpp ttm4-smear.cpp gammq.cpp tang-toennies.cpp ps.cpp \
            x2b-dispersion.cpp x2b-v9x.cpp x3b-v2x.cpp mbpol.cpp

MBPOL_OBJ = libmbpol.o \
            poly-2b-v6x.o poly-2b-v6x-nogrd.o \
            poly-3b-v2x.o poly-3b-v2x-nogrd.o \
            ttm4-es.o ttm4-smear.o gammq.o tang-toennies.o ps.o \
            x2b-dispersion.o x2b-v9x.o x3b-v2x.o mbpol.o

all: libmbpol.so

libmbpol.so: $(MBPOL_OBJ)
	$(CXX) $(CXXFLAGS) -fPIC -c $(MBPOL_CPP)
	$(CXX) -shared -o libmbpol.so -fPIC $(MBPOL_OBJ)

clean:
	rm -f *.so *.o

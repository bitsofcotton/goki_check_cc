CXX=	clang++

# compiler flags.
CXXFLAGS+=	-I/usr/local/include/eigen3
#CXXFLAGS+=	-fopenmp -pthread
CXXFLAGS+=	-std=c++11
CXXFLAGS+=	-O2 -g2 -mtune=native
LDFLAGS+=	-lstdc++

CLEANFILES= *.o tools

all:	tools

clean:
	@rm -rf ${CLEANFILES}

tools.o:        tools.cc edgedetect.hh enlarge.hh fisheye.hh obj2vector.hh ppm2eigen.hh scancontext.hh tilt.hh


CXX=	clang++

# compiler flags.
CXXFLAGS+=	-I/usr/local/include/eigen3
#CXXFLAGS+=	-fopenmp -lgomp
CXXFLAGS+=	-std=c++11
CXXFLAGS+=	-Ofast -g0 -mtune=native
#CXXFLAGS+=	-O2 -g2 -mtune=native
LDFLAGS+=	-lc++

CLEANFILES= *.o tools

all:	tools

clean:
	@rm -rf ${CLEANFILES}

tools.o:        tools.cc enlarge.hh fisheye.hh match.hh tilt.hh redig.hh fileio.hh


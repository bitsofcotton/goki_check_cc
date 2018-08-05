CXX=	clang++

# compiler flags.
CXXFLAGS+=	-I/usr/local/include/eigen3
#CXXFLAGS+=	-fopenmp -lgomp
CXXFLAGS+=	-std=c++11
# or use this for external libraries.
#CXXFLAGS+=	-std=c++17
CXXFLAGS+=	-Ofast -g0 -mtune=native
#CXXFLAGS+=	-O2 -g2 -mtune=native
# Do not use this because of the slowness, so this implementation is for what around licenses.
#CXXFLAGS+=	-D_WITHOUT_EIGEN_
#CXXFLAGS+=	-D_WITH_GLTF2_
LDFLAGS+=	-lc++

CLEANFILES= *.o tools

all:	tools

clean:
	@rm -rf ${CLEANFILES}

tools.o:        tools.cc enlarge.hh match.hh redig.hh fileio.hh


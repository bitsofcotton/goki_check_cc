CXX=	clang++
LD=	${CXX}

# compiler flags.
CXXFLAGS+=	-I/usr/local/include/eigen3
#CXXFLAGS+=	-fopenmp -lgomp
CXXFLAGS+=	-std=c++11
# or use this for external libraries.
#CXXFLAGS+=	-std=c++17
CXXFLAGS+=	-O2 -g2
# -Ofast is something bugly in some compiler implementations but don't inspected the bug.
# -mtune=native is not for porting.
# Do not use this because of the slowness, so this implementation is for what around licenses.
#CXXFLAGS+=	-D_WITHOUT_EIGEN_
# Please read the library page before to use.
#CXXFLAGS+=	-D_WITH_GLTF2_
LDFLAGS+=	-lc++

CLEANFILES= *.o gokicheck

all:	gokicheck
clean:
	@rm -rf ${CLEANFILES}
tools.o: tools.cc enlarge.hh match.hh redig.hh fileio.hh
gokicheck: tools.o
	${LD} ${LDFLAGS} -o gokicheck tools.o


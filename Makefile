CXX=	clang++
LD=	${CXX}

# compiler flags.
CXXFLAGS+=	-I/usr/local/include/eigen3
#CXXFLAGS+=	-fopenmp -lomp -pthread
CXXFLAGS+=	-std=c++11
CXXFLAGS+=	-Ofast -mtune=native -gfull
# Do not use these because of the slowness,
# so this implementation is for what around licenses.
#CXXFLAGS+=	-D_WITHOUT_EIGEN_
# Accuracy reason, not needed.
#CXXFLAGS+=	-D_WITH_NO_FLOAT_
#CXXFLAGS+=	-D_WITH_MPFR_=512
LDFLAGS+=	-lc++

CLEANFILES= *.o gokicheck

all:	gokicheck
clean:
	@rm -rf ${CLEANFILES}
tools.o: tools.cc enlarge.hh match.hh redig.hh fileio.hh p0.hh p1.hh simplelin.hh ifloat.hh
gokicheck: tools.o
	${LD} ${LDFLAGS} -o gokicheck tools.o


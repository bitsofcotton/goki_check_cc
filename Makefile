CXX=	clang++
#CXX=	/usr/local/bin/eg++

# compiler flags.
CXXFLAGS+=	-Ofast -mtune=native -gfull
#CXXFLAGS+=	-Oz -mtune=native -gfull
MPFLAGS=	-I/usr/local/include -L/usr/local/lib -lomp -fopenmp
CXXFLAGS+=	-std=c++11
LDFLAGS+=	-lc++ -L/usr/local/lib
#LDFLAGS+=	-lestdc++ -L/usr/local/lib

CLEANFILES= *.o gokibin gokibin32 gokibinmp gokibin32mp

clean:
	@rm -rf ${CLEANFILES}

all:	gokibin gokibin32 gokibinmp gokibin32mp
gokibin:	tools.cc
	${CXX} ${CXXFLAGS} -static -o gokibin tools.cc
gokibin32:	tools.cc
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o gokibin32 tools.cc
gokibinmp:	tools.cc
	${CXX} ${CXXFLAGS} ${MPFLAGS} -o gokibinmp tools.cc
gokibin32mp:	tools.cc
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=32 -o gokibin32mp tools.cc


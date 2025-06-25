CXX=	clang++
#CXX=	eg++
#CXX=	c++

# compiler flags.
CXXFLAGS+=	-Ofast -mtune=native -gfull
#CXXFLAGS+=	-O2 -mtune=native -g3
#CXXFLAGS+=	-O2 -g3
#CXXFLAGS+=	-Oz -mtune=native -gfull
MPFLAGS=	-I/usr/local/include -L/usr/local/lib -lomp -fopenmp
#CXXFLAGS+=	-pg
#CXXFLAGS+=	--analyze
CXXFLAGS+=	-std=c++11
#CXXFLAGS+=	-std=gnu++98
LDFLAGS+=	-lc++ -L/usr/local/lib
#LDFLAGS+=	-lestdc++ -L/usr/local/lib

# N.B. sed -e s/static\ inline//g | sed -e s/inline//g
#CXXFLAGS+=	-D_OLDCPP_ -ftemplate-depth-99
#LDFLAGS+=	-lm

CLEANFILES= *.o gokibin gokibin32 gokibinmp gokibin32mp

clean:
	@rm -rf ${CLEANFILES}

all:	gokibin gokibin32 gokibinmp gokibin32mp
gokibin:
	${CXX} ${CXXFLAGS} -static -o gokibin tools.cc
gokibin32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o gokibin32 tools.cc
gokibinp:
	${CXX} ${CXXFLAGS} -static -D_PERSISTENT_ -o gokibinp tools.cc
gokibinmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -o gokibinmp tools.cc
gokibin32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=32 -o gokibin32mp tools.cc
gokibinpmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PERSISTENT_ -o gokibinpmp tools.cc


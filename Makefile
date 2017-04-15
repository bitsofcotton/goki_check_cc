# gcc.
CXX=	g++

# compiler flags.
CXXFLAGS=	-I/Users/kazunobu/altroot/include
CXXFLAGS+=	-Ofast
#CXXFLAGS+=	-O2 -g2
LDFLAGS=	-lstdc++

CLEANFILES= *.o *.dSYM tools

all:	tools

clean:
	@rm -rf ${CLEANFILES}


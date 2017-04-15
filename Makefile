# gcc.
CXX=	g++

# compiler flags.
CXXFLAGS+=	-Ofast
#CXXFLAGS+=	-O2 -g2
LDFLAGS=	-lstdc++

CLEANFILES= *.o *.dSYM tools

all:	tools

clean:
	@rm -rf ${CLEANFILES}


#
# Makefile to compile REMETA library
#
AR         ?= ar
CXX        ?= g++
CXXFLAGS    = -O3 -Wall -std=c++17
CFLAGS      = 
INC         = ../eigen-3.4.0/

ifneq ($(strip $(HTSLIB_PATH)),) # defined externally
	CXXFLAGS   += -DWITH_HTSLIB
endif

UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)
ifeq ($(UNAME_S),Darwin)
 ifeq ($(UNAME_M),x86_64)
	 CXXFLAGS += -arch x86_64
 endif
endif

OBJECTS       = bgz_writer.o regenie_ld_matrix_writer.o

remeta.a: ${OBJECTS}
	${AR} rcs $@ $^

%.o: %.cpp
	${CXX} ${CXXFLAGS} ${CPPFLAGS} -o $@ -c $< -I${INC}

clean:
	rm -f *.o *.a
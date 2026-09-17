#
# Makefile to compile PGEN library
#
# Pass WITH_S3=1 to include AWS S3 streaming support (plink2_s3.cc).
# Pass AWS_INCLUDES="-I/path/to/aws/include" when headers are not in the
# default search path (required when WITH_S3=1).
#
AR         ?= ar
CXX        ?= g++
CXXFLAGS    = -O3 -Wall -std=c++17 -fPIC
CFLAGS      =

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
 CXXFLAGS += -arch arm64
endif

# WITH_S3=1: compile plink2_s3.cc and enable USE_S3 in pgenlib_read.cc
ifeq ($(WITH_S3),1)
 CXXFLAGS   += -DUSE_S3 $(AWS_INCLUDES)
 S3_OBJECTS  = plink2_s3.o
else
 S3_OBJECTS  =
endif

OBJECTS  = $(patsubst %.cc,%.o,$(wildcard ./include/*.cc)) \
           $(patsubst %.cpp,%.o,$(wildcard *.cpp)) \
           $(patsubst %.cc,%.o,$(wildcard *.cc))
# plink2_s3.cc is in the root; only add it when WITH_S3=1 (avoid double-add)
OBJECTS := $(filter-out plink2_s3.o,$(OBJECTS)) $(S3_OBJECTS)

INC      = -I. -I./simde/ -I./include/

pgenlib.a: ${OBJECTS}
	${AR} rcs $@ $^

%.o: %.cpp
	${CXX} ${CXXFLAGS} -o $@ -c $< ${INC}

%.o: %.cc
	${CXX} ${CXXFLAGS} -o $@ -c $< ${INC}

clean:
	rm -f *.o *.a ./include/*.o

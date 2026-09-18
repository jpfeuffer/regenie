#
# Makefile to compile PGEN library
#
# Pass WITH_S3=1 to include S3/HTTP streaming support. That builds the bundled
# s3stream/ (curl + SigV4, no AWS SDK) along with plink2_s3.cc, and requires
# libcurl >= 7.75.0 at link time.
#
AR         ?= ar
CXX        ?= g++
CXXFLAGS    = -O3 -Wall -std=c++17 -fPIC
CFLAGS      =

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
 CXXFLAGS += -arch arm64
endif

# WITH_S3=1: compile plink2_s3.cc plus s3stream/, and enable USE_S3 in
# pgenlib_read.cc
ifeq ($(WITH_S3),1)
 CXXFLAGS   += -DUSE_S3 -DS3STREAM_ENABLE
 S3_OBJECTS  = plink2_s3.o $(patsubst %.cc,%.o,$(wildcard ./s3stream/*.cc))
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
	rm -f *.o *.a ./include/*.o ./s3stream/*.o

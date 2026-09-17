FLAG  =-O3 -march=native -ffp-contract=fast -fno-math-errno -fno-trapping-math -DNDEBUG -fPIC
AR   ?= ar
CXX  ?= g++

# Optional: path to Eigen headers (auto-detected from build tree or system)
EIGEN_INC ?= $(shell \
  if [ -d "../../build/_deps/eigen3-src" ]; then echo "-I../../build/_deps/eigen3-src"; \
  elif [ -d "/usr/include/eigen3" ]; then echo "-I/usr/include/eigen3"; \
  elif [ -d "/opt/homebrew/include/eigen3" ]; then echo "-I/opt/homebrew/include/eigen3"; \
  fi)

# OpenMP flags (only if header + compiler support exist)
OPENMP_FLAG ?= $(shell \
  echo '\#include <omp.h>' | $(CXX) -fopenmp -x c++ -E - > /dev/null 2>&1 && echo "-fopenmp" || \
  (echo '\#include <omp.h>' | $(CXX) -Xpreprocessor -fopenmp -x c++ -E - > /dev/null 2>&1 && echo "-Xpreprocessor -fopenmp -lomp") || \
  echo "")

libMvtnorm.a: mvtnorm.o
	${AR} rcs $@ $^ 

mvtnorm.o: mvtnorm.cpp mvtnorm.h mvtdst.hpp
	${CXX} $(FLAG) -std=c++17 $(EIGEN_INC) $(OPENMP_FLAG) -c mvtnorm.cpp -o mvtnorm.o

clean:
	rm -f *.o *.a

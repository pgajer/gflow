# src/tests/Makefile
R_HOME := $(shell R RHOME)

CXX := g++
CC := gcc
CXXFLAGS := -std=c++17 -O2 -g -stdlib=libc++
CFLAGS := -O2 -g

# Include paths
INCLUDES := -I..
INCLUDES += -I../../inst/include/gflow
INCLUDES += -I../../inst/include
INCLUDES += -I$(R_HOME)/include
INCLUDES += -I/usr/local/include/eigen3

# Eigen configuration flags
INCLUDES += -DEIGEN_DONT_VECTORIZE -DEIGEN_DONT_PARALLELIZE
INCLUDES += -D_Alignof=__alignof__

# Force-include eigen_config.hpp (C++ only)
CXX_INCLUDES := $(INCLUDES)
CXX_INCLUDES += -include ../../inst/include/gflow/eigen_config.hpp

# Other flags
INCLUDES += -DR_NO_REMAP
INCLUDES += -DDLL_API=

# Linker flags
LDFLAGS := -L$(R_HOME)/lib -lR -lRblas -lRlapack -lstdc++

# C++ source files
SOURCES_CPP := main.cpp
SOURCES_CPP += ../riem_dcx.cpp
SOURCES_CPP += ../riem_dcx_regression.cpp
SOURCES_CPP += ../riem_dcx_initialization.cpp
SOURCES_CPP += ../kNN_eigen.cpp
SOURCES_CPP += ../kNN.cpp
SOURCES_CPP += ../geometric_edge_pruning.cpp
SOURCES_CPP += ../set_wgraph.cpp
SOURCES_CPP += ../graph_conn_components.cpp
SOURCES_CPP += ../SEXP_cpp_conversion_utils.cpp
SOURCES_CPP += ../prune_long_edges.cpp
SOURCES_CPP += ../graph_maximal_packing.cpp
SOURCES_CPP += ../shortest_path_utils.cpp

# C source files
SOURCES_C := ../kernels.c

# ANN sources
ANN_SOURCES := $(wildcard ../ANN/*.cpp)

# Compile everything
test_runner: $(SOURCES_CPP) $(SOURCES_C) $(ANN_SOURCES)
	$(CXX) $(CXXFLAGS) $(CXX_INCLUDES) -o test_runner \
		$(SOURCES_CPP) $(SOURCES_C) $(ANN_SOURCES) \
		$(LDFLAGS)

run: test_runner
	./test_runner

clean:
	rm -f test_runner

.PHONY: run clean

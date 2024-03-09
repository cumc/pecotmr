# Makefile for compiling C++ code into an R package

# Compiler options
CXX = g++
CXXFLAGS = -std=c++11 -fopenmp -O3 -Wall

# R-related flags
R_INC := $(shell Rscript -e "Rcpp:::CxxFlags()")
R_LIB := $(shell Rscript -e "Rcpp:::LdFlags()")

# Source files and output shared object file
SRC = dentist.cpp RcppExports.cpp
TARGET = dentist.so

# Rule to compile the code
all: $(TARGET)

$(TARGET): $(SRC)
    $(CXX) $(CXXFLAGS) $(R_INC) -shared -o $@ $^ $(R_LIB)

clean:
    rm -f $(TARGET)


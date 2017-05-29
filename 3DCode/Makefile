#My first mae file 
# Keep it as simple as possible
#

# Compiler
# g++

CXX = g++

# Series of useful tokens for compilation and linking
LIBTOKEN = -L
LIBPATH = /home/sdea/armadillo-7.800.2
ARMALIB = -larmadillo

# This is to use optimized version of LAPACK and OPENBLAS in MacOS
LIB_OS = -framework Accelerate 

# Optimization flac
OPT = -O3
OPTFLAG = -march=native

DYLIB = $(LIBTOKEN) $(LIBPATH) $(ARMALIB)

# Search directory where to find the armadillo library
ARMADIR = /home/sdea/armadillo-7.800.2/include
IFLAG = -I

IDIR =  $(IFLAG) $(ARMADIR)

# Name of the output file
OUTNAME = main

$(OUTNAME): main.cpp 
	$(CXX) -o $(OUTNAME) main.cpp $(OPT) $(OPTFLAG) $(DYLIB) $(IDIR) 


############################################################################
#
#  Program:         SuperLU
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   May 10, 2015
#
#  Modified:	    
#		    
#
############################################################################
#
#  The name of the libraries to be created/linked to
#
SuperLUroot	= /usr/local
SUPERLULIB   	= $(SuperLUroot)/lib/libsuperlu.a
COINLIBroot = /usr/local
COINLIB		= $(COINLIBroot)/lib/libCoinUtils.so

#TMGLIB       	= libtmglib.a
MATGENLIB	= $(SuperLUroot)/TESTING/MATGEN/libmatgen.a

XSDK_INDEX_SIZE = 
HAVE_METIS      = 

# BLASDEF 	= -DUSE_VENDOR_BLAS
BLASLIB		= /usr/local/lib/libblas.a
LIBS		= $(SUPERLULIB) $(BLASLIB) $(COINLIB)
LIBS  += 

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = /usr/bin/ar
ARCHFLAGS    = cr
RANLIB       = /usr/bin/ranlib

CXX           = /usr/bin/g++
CXXFLAGS 	     =  -std=c++11 -Wall -O2
INCLUDEDIR  = -I$(SuperLUroot)/include -I$(COINLIBroot)/include/coin-or
FORTRAN	     = 

LOADER       = $(CXX)
LOADOPTS     =

#
# The directory in which Matlab is installed
#
# MATLAB	     = /Applications/MATLAB_R2015b.app



#######################################################################
#  This makefile creates the example programs for the linear equation
#  routines in SuperLU.  The files are grouped as follows:
#
#       SLINEXM -- Single precision real example routines
#       DLINEXM -- Double precision real example routines
#       CLINEXM -- Double precision complex example routines
#       ZLINEXM -- Double precision complex example routines
#
#  Example programs can be generated for all or some of the four different
#  precisions.  Enter make followed by one or more of the data types
#  desired.  Some examples:
#       make single
#       make single double
#  Alternatively, the command
#       make
#  without any arguments creates all four example programs.
#  The executable files are called
#       slinsol		slinsolx 
#       dlinsol		dlinsolx
#       clinsol		clinsolx
#       zlinsol		zlinsolx
#
#  To remove the object files after the executable files have been
#  created, enter
#       make clean
#  On some systems, you can force the source files to be recompiled by
#  entering (for example)
#       make single FRC=FRC
#
#######################################################################


TARGET         = main
SRCS 		   = main.cpp 
OBJS           = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(LOADER) $(LOADOPTS) $(OBJS) $(LIBS) -lm -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDEDIR) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)





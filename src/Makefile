##                           -*- Mode: Makefile -*- 
## Makefile --- 
## Author          : Amelia Bellamy-Royds
## Created         : Summer 2006
## Last Modified By: Marcel Turcotte
## Last Modified On: Mon Dec 11 12:16:37 2006
##
## Contributions:
##
##   Amelia Bellamy-Royds developed this version of the program based
##   on the earlier work by Luke Cen.  Both were students under the
##   supervision of Marcel Turcotte.  The initial version of the source 
##   code was developed by David Mathews, see below.
##
## Amelia B. Bellamy-Royds and Marcel Turcotte (2006) Can Clustal-Style 
##   Progressive Pairwise Alignment of Multiple Sequences Be Used in RNA 
##   Secondary Structure Prediction?
##
## This copyrighted source code is freely distributed under the terms
## of the GNU General Public License. 
##
## See LICENSE for details.
##
## This work is based on the source code of Dynalign:
##
## RNA Secondary Structure Prediction, Using the Algorithm of Zuker
## C++ version by David Mathews
## Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004
## Programmed for Isis Pharmaceuticals and 
## the Turner Lab Chemistry Department, University of Rochester

# CC = g++
# CC_FLAGS = -O2 -Wall -Wno-deprecated 

CC = CC
CC_FLAGS = -O2

BINDIR = ../bin

all : pdynalign

pdynalign : algorithm.cpp datatable.cpp efunctions.cpp pd_interface.cpp pdynalign.cpp platform.cpp profile.cpp stackclass.cpp structure.cpp ASSERT.h defines.h pdynalign.h
	$(CC) $(CC_FLAGS) -o pdynalign pd_interface.cpp -lm

install:
	mv pdynalign $(BINDIR)

clean :
	-rm -f pdynalign



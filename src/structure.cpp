/*                               -*- Mode: C -*- 
 * datable.cpp ---
 * Author          : Amelia Bellamy-Royds
 * Created         : Summer 2006
 * Last Modified By: Marcel Turcotte
 * Last Modified On: Mon Dec 11 11:10:44 2006
 *
 * Contributions:
 *
 *   Amelia Bellamy-Royds developed this version of the program based
 *   on the earlier work by Luke Cen.  Both were students under the
 *   supervision of Marcel Turcotte.  The initial version of the source 
 *   code was developed by David Mathews, see below.
 *
 * Amelia B. Bellamy-Royds and Marcel Turcotte (2006) Can Clustal-Style 
 *   Progressive Pairwise Alignment of Multiple Sequences Be Used in RNA 
 *   Secondary Structure Prediction?
 *
 * This copyrighted source code is freely distributed under the terms
 * of the GNU General Public License. 
 *
 * See LICENSE for details.
 *
 * This work is based on the source code of Dynalign:
 *
 * RNA Secondary Structure Prediction, Using the Algorithm of Zuker
 * C++ version by David Mathews
 * Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004
 * Programmed for Isis Pharmaceuticals and 
 * the Turner Lab Chemistry Department, University of Rochester
 */

/* Definition and methods for struct "structure"
	which describes a single sequence and its associated base-pairing structures
	At present, the pdynalign program only calculates a single structure for each
	sequence, however the data structure is built to store multiple sub-optimal structures.
	(Just change the defined constant MAXSTRUCTURE)
	*/

#if !defined(STRUCTURE_CPP)
#define STRUCTURE_CPP


#include "defines.h"


/* Structure:
	contains all the info for an RNA structure based on a single sequence*/
struct structure
{
	short int	numofbases;
	short int	*numseq; //numerical sequence

	short int	*hnumber; //"historical" numbering of nucleotides
	char		*nucs; //character sequence

	//short int pair[MAXFORCE][2],npair; //forced base-pairs
		//NOT USED

	short int nopair[MAXFORCE]; //forced non-pairing bases
		// stores position of nucleotides that forced to be single
	short int nnopair;
		// counter for the number of nopairs

	//short int **basepr; //NOT USED

	int energy[MAXSTRUCTURES+1];
	//int inter[3]; //NOT USED

	char ctlabel[MAXSTRUCTURES+1][SEQNAME_LN+1];
	bool intermolecular, allocated;
	//bool templated;
	//bool **tem; //NOT USED

	structure();
	~structure();
	void allocate(int size = MAXBASES);

/****************************************************************************
	numofbases = number of bases in sequence
				that is held by structure

	numseq[i] = a numeric that stands for the base in the ith position
				of the sequence,
			X = 0 //or any other non-recognized character, such as N
			A = 1
			C = 2
			G = 3
			U = 4
			- = 5 // gap (or I, "intermolecular linker")

	basepr[i][j] = base to which the jth base is paired in the ith structure
					(only i = 1 calculated)

	energy[i] = 10 * the Gibb's free energy of the ith structure, this
				is done so that integer math can be performed
				and the nearest tenth of a kcal/mol can be
				followed
				(only i = 1 calculated)

	ctlabel = a string of information for each of the structures

	hnumber array stores the historical numbering of a sequence
   		-- i.e. the numbering without gaps from the profile

	nucs is a character array to store the sequence information --
   		this will allow the program to keep T and U from getting confused

  	allocated -- whether memory has been allocated for arrays
  	intermolecular -- whether sequence contains intermolecular links
****************************************************************************/

};

/* Constructor of Structure */
structure::structure()
{
	int i;

	for (i=1;i<=MAXSTRUCTURES;i++) {
		energy[i]=0;
	}

	allocated = false;
	nnopair=0;
	//npair=0;
	intermolecular = false;
	//templated = false;
}

/* Destructor of structure */
structure::~structure()
{
	if (allocated) {
		delete[] numseq;
		delete[] hnumber;
		delete[] nucs;
	}
}


/* Allocate memory for structure */
void structure::allocate(int size)
{
	int i;
	numseq = new short int [2*size+1];
	hnumber = new short int [size+1];
	for( i=0; i<size+1; i++) {
		hnumber[i] = 0;
	}
	nucs = new char [size+2];
	allocated = true;
}



#endif

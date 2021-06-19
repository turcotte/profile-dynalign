/*                               -*- Mode: C -*- 
 * profile.cpp ---
 * Author          : Amelia Bellamy-Royds
 * Created         : Summer 2006
 * Last Modified By: Marcel Turcotte
 * Last Modified On: Mon Dec 11 11:09:53 2006
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

/* Definition and methods for struct "profile", plus struct "alignment"*/

#include "defines.h"

struct profile
{
	char			id[ID_LN_MAX];
		// maximum length of id of the profile

	int 			nSeq;
		// number of sequences in the profile

	short int	numofbases;
		// number of bases in the profile

	bool			allocated;
		// indicates memory has been allocated

	structure	*ct;
		// an array of structures for each sequence

	short int	*basepr;
		// basepr[i] stores the idx. of the element paired with i

	profile();		// profile constructor
	~profile();		// profile destructor
	void allocate(int nSeq, int numofbases);	// allocate memory for profile
};


/* alignment pair info */
struct alignment
{
	short pos;					// column position in the alignment
	char	indicator[1];	// either '(' or ')'
};


/* Constructor of profile */
profile::profile()
{
	nSeq = 0;
	numofbases = 0;
	allocated = false;
	strcpy(id,"");
}


/* Destructor of profile */
profile::~profile()
{
	if(allocated)
	{
		delete[] ct;
		delete[] basepr;
	}
}


/* Allocate memory for profile */
void profile::allocate(int numSeq, int numbases)
{
	ct = new structure[numSeq];
	basepr = new short int [numbases+1];
	allocated = true;
}


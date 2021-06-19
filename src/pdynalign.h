/*                               -*- Mode: C -*- 
 * pdynalign.h ---
 * Author          : Amelia Bellamy-Royds
 * Created         : Summer 2006
 * Last Modified By: Marcel Turcotte
 * Last Modified On: Mon Dec 11 11:09:10 2006
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

/* This file contains prototypes for all methods in pdynalign.cpp,
	algorithm.cpp, and efunctions.cpp
   It also "includes" the defines.h header and the classes datatable
   structure, profile, and stackclass */

#if !defined(dynalignh)
#define dynalignh


/***** Define constants plus the struct HDException *****/

#include "defines.h"		//define constants

/***** Define essential classes *****/

#include "datatable.cpp"	//define datatable class & methods
#include "structure.cpp"	//define structure class & methods
#include "profile.cpp"		//define profile class & methods
#include "stackclass.cpp"	//define stackclass class & methods

/***** This method in file pdynalign.cpp *****/

// Provide major computing functionality of the H-Dynalign algorithm
void pdynalign(profile *pf1, profile *pf2, alignment *align,
	short int maxseparation, short int gapincrease, datatable *data, bool singleinsert,
	int *totalscore=NULL);

/***** These methods in file efunctions.cpp *****/

// energy funtions for a profile
energy_t einternal_pf(int i,int j,int ip,int jp,profile *pf, datatable *data);

inline energy_t ehairpin_pf(short i, short j, profile *pf, datatable *data);

inline energy_t edangle5_pf(int i,int j,int ip,profile* pf,datatable* data);

inline energy_t edangle3_pf(int i,int j,int ip,profile* pf,datatable* data);

inline energy_t ebp_pf(int i,int j,int ip,int jp,profile *pf, datatable *data);

//calculates end of helix penalty for a profile
inline int penalty_pf(int i,int j,profile *pf,datatable *data);

// Energy functions for a single sequence
energy_t einternal(int i,int j,int ip,int jp,structure *ct, datatable *data);

inline energy_t ehairpin(short i, short j, structure *ct, datatable *data);

inline energy_t edangle5(int i,int j,int ip,structure* ct,datatable* data);

inline energy_t edangle3(int i,int j,int ip,structure* ct,datatable* data);

inline energy_t ebp(int i,int j,int ip,int jp,structure *ct, datatable *data);

// Calculates end of helix penalty for a sequence
inline int penalty(int i,int j,structure *ct,datatable *data);


/***** These methods in file algorithm.cpp *****/

// Computes energy penalty for terminal pairs
int penalty2(int i, int j, datatable *data);

// Outputs a set of ct files
void ctout_pf(profile *pf, char *ctoutdir);

// Converts base chars to integer numbers, and put into structure->numseq[count].
void tonum(char *base,structure *ct,int count);

// Convert a numeric value for a base to the familiar character
char *tobase (int i);

// Open and read data from a profile file
int openProfile (profile *pf, char *fileName);

// Gets a substring from a string
char *copySubStr(char *substr, char *str, size_t size);

// Verify if a char is one of "AaCcGgTtUuXxN-"
int isValidChar(char *s);

// Output an alignment of 2 profiles to a file specified by 'aout'
void alignout(alignment* align,char *aout, profile *pf1, profile *pf2,
			int totalscore, int maxsep, int gappen, bool insert);

// Check if one element can paired with another in a profile
short int get_inc(int i, int j, profile *pf);

// Append a nucleotide characters to each line of a profile in the alignment
void	append_nucs(char **line, profile *pf, int pos);

// Append a string to each line of a profile in the alignment.
void	append_str(char **line, profile *pf, char *s);

// Initialize the lines for alignment output
void init_lines(char **line, profile *pf);


#endif


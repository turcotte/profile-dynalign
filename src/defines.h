/*                               -*- Mode: C -*- 
 * defines.cpp ---
 * Author          : Amelia Bellamy-Royds
 * Created         : Summer 2006
 * Last Modified By: Marcel Turcotte
 * Last Modified On: Mon Dec 11 11:07:21 2006
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

// defines.h

#if !defined(DEFINES_H)
#define DEFINES_H

/***** Define a struct for reporting errors *****/

struct HDException {
	int errorvalue;
	char *message;
	HDException(int e, char *m) {errorvalue = e; message = m;}
};


// to run the program in debug mode, remove the comment on "#define DEBUG".
// to run in regular run mode, comment out "#define DEBUG".
//#define	DEBUG
#include	"ASSERT.h"  //will add in debugging if DEBUG defined

//#define LARGE_MEMORY
// define "LARGE_MEMORY" to use (4-byte) ints for all energy values,
// comment it out to use (2-byte) shorts
#ifdef LARGE_MEMORY

typedef int energy_t;
#define INFINITE 1000000  //an arbitrary value given to INFINITE

#else

typedef short energy_t;
#define INFINITE 8000
//a large energy value to give to impossible events for use in minimization
//should be less than (1/4)*(maxvalue(energy_t)) = (1/4)*(2^(8*sizeof(energy_t) -1)) 
//to prevent wrap-around arithmetic
//but more than 10*gap-penalty*number of sequences in a profile * expected max number of gaps in a profile (= maxseparation);
//to prevent it from being incorrectly selected as a minimum value


#endif

//the environment variable which contains the path for the directory of data files:
#define DATA_DIR "DYNALIGN_DIR"

#define MAXFIL 250    //maximum length of file names
#define MAXTLOOP 100 //maximum tetraloops allowed (info read from tloop)
#define MAXSTRUCTURES 1 //maximum number of structures stored per sequence
#define MAXBASES 3000   //maximum number of bases in a structure
#define MAXFORCE 600 //maximum number of bases that can be forced single

#define CT_LINE_LN 200	//length of a line in the ct files
#define COL 80  //this is the number of columns in an output file
#define NUMLEN 8  //maximum digits in a number

#define MINLOOP 3 //The minimum substructure with a basepair possible
#define MAXLOOP 20 //maximum size of unpaired nucs on each side of an internal loop

#define VALID_SYM	"AaCcGgTtUuXxN-I"
#define WHITESPACE "\n\r\t "
#define LINEBREAK "\n\r"
#define PATH_SEPARATOR '/'

#define DONE			0
#define IO_ERR			5
#define BAD_INPUT		10
#define FILE_ERR		25
#define ERR				-1

#define MBLFreeBase	eparam[6]
#define MBLClose	eparam[5]
#define MBLHelix	eparam[10]

/* associates to profile file format */
#define ID_LN_MAX		100		// maximum length of the profile id

#define NUM_SEQ_MAX	1000		// maximum number of sequences in a profile
#define SEQNAME_LN	20		//(max) length of the sequence name (should be longer than any
							//column mark-up headers used, but not to long
							//--determines indent in alignment files)
#define SEQ_LN_MAX	10000		// maximum length of a sequence (10 000 nucleotides)

#define STOCKHOLM_HEADER	"# STOCKHOLM 1.0" //format header
#define	STOCKHOLM_END	"//"	// end of profile indicator
#define STOCKHOLM_MARKUP_START "#" //indicates start of a mark-up (non-sequence) line

//"magic" mark-up headers used in Stockholm format:
#define FILE_MARKUP "#=GF" //a mark-up line applying to the entire file
#define SEQ_MARKUP "#=GS" //a mark-up line applying to a specific sequence
#define COL_MARKUP "#=GC" //a column-by-column mark-up line for the entire alignment
#define SEQ_COL_MARKUP "#=GR" //a column-by-column mark-up line for a specific sequence

//The following are types of file mark-up comments used by this program in its output
#define COMMENT_MARKER "CC" //mark-up line contains descriptive comments
#define ID_MARKER "ID" //mark-up gives a profile ID
#define NSEQ_MARKER "SQ" // mark-up gives the number of sequences in the profile
#define TOTAL_E_MARKER "EN" //mark-up gives the total alignment energy (incl. gap penalty)
#define AVG_E_MARKER "EN_struct_avg" //mark-up gives the average structural energy per sequence
#define GP_MARKER "GAP_PEN" // mark-up gives the gap penalty used
#define TOTAL_GAP_MARKER "N_GAPS" //mark-up gives the total number of gaps added
#define MAXSEP_MARKER "MAX_SEP" //mark-up gives the maximum separation used
#define BP_MARKER "BP_INSERTS_ALLOWED" //indicates basepair inserts were allowed (no value needed)

//The following is the column mark-up used
#define SS_MARKER "SS_cons" //mark-up is the consensus secondary structure (bracket notation)



#endif


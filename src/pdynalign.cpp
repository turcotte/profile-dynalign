/*                               -*- Mode: C -*- 
 * pdynalign.cpp --- simultaneous alignment and structure prediction of 2 RNA profiles
 * Author          : Amelia Bellamy-Royds
 * Created         : Summer 2006
 * Last Modified By: Marcel Turcotte
 * Last Modified On: Mon Dec 11 11:08:50 2006
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

/*
This file contains one function: pdynalign().
This is the main calculation function for the program.
*/

#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <string.h>

/* NOTE: order is important in these include statements! */

#include "pdynalign.h"		//prototype all general methods
                                //includes classes and constants

#include "algorithm.cpp"	//define utility functions
#include "efunctions.cpp"	//define energy-calculation functions


/*-----------------------------------------------------------------------------
Function:	pdynalign()

Description: Provide major computing functionality of the H-Dynalign algorithm
-----------------------------------------------------------------------------*/
void pdynalign( profile *pf1, profile *pf2, alignment *align,
			 short int maxseparation, short int gapincrease, datatable *data,
			 bool singleinsert, int *totalscore)
{
/* ---- declare local variables ---- */

	register short int i,j,k,l, k_,l_, c,d,e,f, size,maxsep,maxrange,gap;
	register energy_t en1, en2, en3, lowest;
	/*	i: index of the start of a sub-segment in sequence 1
		j: index of the end (inclusive) of a sub-segment in sequence 1
		k: index of the start of a sub-segment in sequence 2
		l: index of the end (inclusive) of a sub-segment in sequence 2

		k_:adjusted value of k for accessing the v and w arrays
		l_:adjusted value of k for accessing the v and w arrays
			NOTE: k_ and l_ have values from 0 to (2*maxseparation+1)
			representing the values k = (i-maxseparation) to (i+maxseparation)
			and l = (j-maxseparation) to (j+maxseparation)
			(inclusive), repsectively.
			Therefore, to convert from k to k_ and l to l_:
			k_ = k  - i + maxseparation
			k  = k_ + i - maxseparation
			l_ = l  - j + maxseparation
			l  = l_ + j - maxseparation

		c, d, e, f:
			variables representing start or end indices of possible sub-structures
			used to build up the best structure when aligning i<->j with k<->l
			(exact definition varies by case)

		en1, en2, en3:  temporary storage of energy values during calculations

		size: the length of the sub-segment i<->j, used as a counter for a "for"-loop
		maxsep: a register variable to store the value of maxseparation
			(the maximum difference in index between two aligned nucleotides in
				the respective sequences)
		maxrange: = 2*maxsep+2, is the range of k_ and l_ values
		gap: a register variable to store the value of gapincrease
			(the penalty per gap in the alignment)

		lowest: the minimum energy value for the entire structure

	*/

	short int	nBases1,nBases2;
		// number of bases (in the alignment) for profile 1 and 2, respectively
	short int	nSeq1,nSeq2;
		// number of sequences for profile 1 and 2, respectively

		bool **pair;
	/*	pair[0][i] indicates whether or not nucleotide i of sequence 1 is
			permitted to form pair-bonds;
		similarly pair[1][k] indicates whether nucleotide k of sequence 2
			is permitted to form bonds.
		Values are based on pre-defined constraints within the profile/sequence
		objects.
		Note that if any of the sequences in an aligned profile has a no-pair
		constraint at a particular index, then the entire column is constrained.
	*/

	int ****v,****w;
	int  **mine;
	/*	v[j][i][k_][l_] is the sum of the lowest free energies for fragments i to j and k to l
		with i and j paired, k and l paired, and i aligned to k and j aligned to l

		w[j][i][k_][l_] is the sum of the lowest free energies for fragments i to j and k to l
		with the fragments interior to multibranch loops and i aligned to k and j aligned to l

		mine[i][k_] is the minumum energy of the fragments 1 to i and 1 to k
		with i aligned to k and no constraints on the configuration
	*/


/* ---- initialize variables and arrays -----*/
	nBases1 = pf1->numofbases;
	nBases2 = pf2->numofbases;
	nSeq1 = pf1->nSeq;
	nSeq2 = pf2->nSeq;

	maxsep = maxseparation; //place maxseparation into a register short
	maxrange = 2*maxsep+2;  //calculate and store the total range of k_ and l_ values
	gap = gapincrease;//place the gap penalty in a register short

	// Initialize en1 and en2
	en1 = 0;
	en2 = 0;
	en3 = 0;

	// Initialize the pair array based on values in the profile objects
	pair = new bool *[2];  //i.e. bool[2][]
	ASSERT(pair != NULL);
	pair[0] = new bool [nBases1+1];
	ASSERT(pair[0] != NULL);
	pair[1] = new bool [nBases2+1];

	ASSERT(pair[1] != NULL);

	for (i=0; i<=nBases1; i++) pair[0][i] = true;
	for (k=0; k<=nBases2; k++) pair[1][k] = true;

	ASSERT(pf1->nSeq >= 1);
	for(j=0; j<(pf1->nSeq); j++) {
		//for each sequence in the profile

		ASSERT((pf1->ct+j)->nnopair >= 0 && (pf1->ct+j) != NULL);
		for(i=1; i<((pf1->ct[j]).nnopair)+1; i++) {
			//for each nucleotide index stored in
			//the nopair array for that sequence

			pair[0][(pf1->ct[j]).nopair[i]]=false;
			//set the value at that index
			//in the pair array to false
		}
	}
	//repeat for second profile:
	for(j=0; j<(pf2->nSeq); j++) {
		ASSERT((pf2->ct+j)->nnopair >= 0 && (pf2->ct+j) != NULL);
		for(i=1; i<((pf2->ct[j]).nnopair)+1; i++) {
			pair[1][(pf2->ct[j]).nopair[i]]=false;
		}
	}

	//allocate space for the v and w arrays
	//****Note: dimensions are [j][i][k_][l_] with i<j, k<l *******
	v = new int ***[nBases1+1];
	w = new int ***[nBases1+1];

	for (j=0; j <= nBases1; j++) {
		v[j] = new int **[j];
		w[j] = new int **[j];

		for (i=0; i < j; i++) {
			v[j][i] = new int *[maxrange];
			w[j][i] = new int *[maxrange];

			for (k_=0; k_ < (maxrange); k_++) {
				v[j][i][k_] = new int [maxrange];
				w[j][i][k_] = new int [maxrange];

				for (l_=0; l_ < (maxrange); l_++) {
					//initialize energy values to infinite
					v[j][i][k_][l_]=INFINITE;
					w[j][i][k_][l_] = INFINITE;
				}
			}
		}
	}

	//allocate space for the mine array.
	//initialize values based on an unstructured alignment
	//(i.e. based on required gaps to align sequences of different length)
	mine = new int *[nBases1+1];
	for (i=0; i <= nBases1; i++) {
		mine[i] = new int[maxrange];

		for (k_=0; k_ < maxrange; k_++){
			mine[i][k_] = gap*((k_> maxsep)?(k_ - maxsep)*nSeq1:(maxsep - k_)*nSeq2);
			//RECALL: k_ = k - i + maxsep;
			//if k_ > maxsep, then k > i:
			//to align sequences add gaps to first profile
			//(number of gaps = k-i = k_ - maxsep)
			//else
			//i >= k, therefore add gaps to second profile
			//(number of gaps = i-k = maxsep - k)
		}
	}

	//initialize the structure portion of the ct files:
	for (i=1;i<=nBases1;i++) pf1->basepr[i]=0;
	for (i=1;i<=nBases2;i++) pf2->basepr[i]=0;

	//initialize the alignment
	for (i = 1; i <= nBases1; i++) {
		align[i].pos = 0;
		strcpy(align[i].indicator, "");
	}



/*---- end of initialization steps ---- */


/*---- fill w and v arrays ---- */

/* calculate best score for the alignment of every possible segment
	i to j in sequence 1 aligned with segment k to l in sequence 2 */

	// size = the structure(loop) size closed by i<->j, include i,j
	// try every possible loop sizes of i<->j

	cout << "Computing energy tables: up to max value " << nBases1 << "\n";
	for (size=MINLOOP;size<=nBases1;size++) {

		//report progress:
		cout <<size<<"..."<<flush;
		if (size % 10 == 0) cout <<"\n" <<flush;


		// for each loop size, try different positions of i<->j
		for (i=1;i+size<=nBases1+1;i++) {
			j = i+size-1;			// distance between i, j is size-1.
			// try different positions of k
			for (k=min(i+maxsep,nBases2); (k>=i-maxsep)&&(k>=1); k--) {
				// try different positions of l
				for (l=max(max(j-maxsep,1),k+MINLOOP); (l<=j+maxsep)&&(l<=nBases2); l++)
				{
					//use k_ and l_ when refering to the energy arrays
					k_ = k-i+maxsep;
					l_ = l-j+maxsep;

					//filter out isolated base pairs
					e = get_inc(i+1, j-1, pf1);
					if ( (i-1>0) && (j+1<=nBases1) )
						e = e + get_inc(i-1, j+1, pf1);
						//i.e. e false only if both get_inc() = 0

					f = get_inc(k+1, l-1, pf2);
					if ( (k-1>0) && (l+1<=nBases2) )
						f = f + get_inc(k-1, l+1, pf2);

/* V Cases: i-j paired, k-l paired -- determines minimum value for v[j][i][k_][l_]*/
/* only if both pairs possible and at least one adjacent pair possible for each */
/* otherwise, energy of arrangement assumed infinite (initialization value for "v" array)*/
					if (get_inc(i, j, pf1) && get_inc(k, l, pf2) &&
						pair[0][i] && pair[0][j] && pair[1][k] && pair[1][l]
						&& e && f)
					{
						en3=INFINITE;

/* Case V1: hairpin loops closed by base-pairs i:j and k:l */
						v[j][i][k_][l_] = ehairpin_pf(i,j,pf1,data)+ehairpin_pf(k,l,pf2,data)
								+gap*((j-i)<(l-k)?((l-k)-(j-i))*nSeq1:((j-i)-(l-k))*nSeq2);

/* Case V2:  pair one internal helix from each profile, and then extend (with ss-sequences)
			as necessary to reach the i:j and k:l basepairs.
			Depending on whether or not the last nucleotides of the internal helix
			are adjacent to i+k, or j+l, or both, the resulting structure is either
			a helix extension, bulged helix, or bounded internal loop.
			*/

						for (c=i+1; (c <= i+MAXLOOP) && (c < j-MINLOOP); c++) {
							for (d=j-1; (d >= j-MAXLOOP) && (d > c+MINLOOP); d--) {

							  	if (!get_inc(c,d,pf1)) break; //new!

								for (e = max(k+1,c-maxsep);
										(e <= k+MAXLOOP) && (e < l-MINLOOP) && ((e-c)<=maxsep);
										e++) {
									for (f = min(l-1,d+maxsep);
											(f >= l-MAXLOOP) && (f > e+MINLOOP) && ((d-f)<=maxsep);
											f--)
									{
									  	if (!get_inc(e,f, pf2)) break; //new!

										if ( (c==i+1) && (d==j-1) ) {
											//seq1 is helical stack
											en1 = ebp_pf(i,j,i+1,j-1,pf1, data);
										}
										else {
											//internal loop (includes bulge)
											en1 = einternal_pf(i,j,c,d,pf1,data);
										}

										if ( (e==k+1) && (f==l-1) ) {
											//seq2 is helical stack
											en2 = ebp_pf(k,l,k+1,l-1,pf2, data);

											//also allow single base pair insertions into one sequence only
											if (singleinsert
													&& (c==i+2) && (d==j-2) && (j-2 > i+2)
													&& get_inc(i+1,j-1,pf1)
													//&& get_inc(i+2,j-2,pf1) //== get_inc(c,d,pf1)
													)
													/*i.e. if it is possible to create three
													  stacked base pairs in pf1 from
													  (i:j, i+1:j-1, c:d == i+2: j+2)
													  and two pairs from (k:l, e:f = k+1:l-1)*/
											{
												/*score two stacked base pairs (k:l, e:f)
												  against three stacked pairs (i:j, i+1:j-1, c:d)*/
												en1 = ebp_pf(i,j,i+1,j-1,pf1, data) +
														ebp_pf(i+1,j-1,i+2,j-2,pf1,data);
											}
										}
										else {
											//internal loop
											en2 = einternal_pf(k,l,e,f,pf2,data);

											//also allow single base pair insertions into one sequence only
											if (singleinsert
													&& (c == i+1) && (d == j-1)
													&& (e == k+2) && (f == l-2) && (l-2 > k+2)
													&& get_inc(k+1,l-1,pf2)
													//&& get_inc(k+2,l-2,pf2) //== get_inc(e,f,pf2)
													//&& get_inc(i+1,j-1,pf1) //== get_inc(c,d,pf1)
													)
													/*i.e. if (1) c,d are adjacent to i,j (but e,f
													  are not immediately adjacent to k,l)
													  and (2) it is possible to create three
													  stacked base pairs in pf2 from
													  (k:l, k+1:l-1, e:f == k+2: l+2)
													  and two pairs from (i:j, c:d = i+1:j-1)*/
											{
												/*score three stacked base pairs (k:l, k+1:l-1, e:f)
												  against two stacked pairs (i:j, c:d)*/
												en2 = ebp_pf(k,l,k+1,l-1,pf2, data) +
														ebp_pf(k+1,l-1,k+2,l-2,pf2,data);
											}
										}

										// en1 = energy change for adding i:j to c:d structure
										// en2 = energy change for adding k:l to e:f structure
										// en3 =	sum of en1 + en2
										//			+ energy for aligned internal segements
										//			+ new gap penalties
										//		--> current minimum for all examined structures
										en3 = min(en3,
											en1 + en2 + v[d][c][e-c+maxsep][f-d+maxsep]
												+ gap*((c-i)<(e-k)?((e-k)-(c-i))*nSeq1:((c-i)-(e-k))*nSeq2)
												+ gap*((j-d)<(l-f)?((l-f)-(j-d))*nSeq1:((j-d)-(l-f))*nSeq2)
												);

									}	// end of for(f)
								}	// end of for(e)
							}	//end of for(d)
						}	// end of for(c)

						v[j][i][k_][l_]=min(en3, v[j][i][k_][l_]);
						//i.e. min of best case V2 value against V1 value

						en1=INFINITE;

/* case V3: multiloop junction: two segments from pf 1, closed by i:j, aligned
			against two segments from pf2, closed by k:l.
			Energy values for the sub-segments taken from the w array
			(i.e. internal helix termini already calculated, as necessary)
			*/
				/* Energy values are the sum of the energies for the two
				alignments of internal multi-branch loop segments,
				plus the offset for closing the multi-branch loop,
				plus penalty for adding an additional helix onto the
				multi-branch loop (Note that these two are added at the end,
				since they are the same for all cases), plus, as necessary,
				penalties for additional free bases in an interior loop, and gap
				penalties whenever a segment is extended (in a particular direction)
				in one sequence but not the other.*/

						for (c = i+MINLOOP+1; c < j-MINLOOP; c++) {
							//c is index within pf1 of split between two segments
							for (d = max(k+MINLOOP+1, c-maxsep);
									(d < l-MINLOOP) && (d <= c+maxsep); d++)
								//d is index within pf2 of split between two segments
							{
								e = d-c+maxsep;
								//e is adjusted value of d for accessing v/w arrays
								//(i.e. equivalent to d_)

								//There are 16 cases, depending on whether or not an
								//additional unpaired nucleotide (i.e. i+1, j-1,
								//k+1, or l-1) is added on the end of the pre-calculated
								//segments, stacked on to the terminus of the new
								//helix (the one starting with i:j, k:l) :
								//consider them in this order (0 unstacked, 1 stacked)

								//For every unpaired nucleotide, must add a MBLFreeBase penalty
								//For every nucleotide added on one end, on one sequence
								//but not the other, add a gap penalty to sequnce without add-on
								//						FreeBases	Gaps
								//						Seq	1	2	1	2
								//		i	j	k	l
								//1		0	0	0	0		0	0	0	0
								//2		0	0	0	1		0	1	1	0
								//3		0	0	1	0		0	1	1	0
								//4		0	0	1	1		0	2	2	0
								//5		0	1	0	0		1	0	0	1
								//6		0	1	0	1		1	1	0	0
								//7		0	1	1	0		1	1	1	1
								//8		0	1	1	1		1	2	1	0
								//9		1	0	0	0		1	0	0	1
								//10	1	0	0	1		1	1	1	1
								//11	1	0	1	0		1	1	0	0
								//12	1	0	1	1		1	2	1	0
								//13	1	1	0	0		2	0	0	2
								//14	1	1	0	1		2	1	0	1
								//15	1	1	1	0		2	1	0	1
								//16	1	1	1	1		2	2	0	0

					/*	Note:	MBLClose is the multi-branch loop offset (closure)
								MBLFreeBase is the multi-branch loop free base penalty
								MBLHelix is the multi-branch loop helix penalty
					*/

		//case 1 - no stacks, no gaps
								//1		0	0	0	0		0	0	0	0
								en1 = min(en1,
										w[c][i+1][k_][e] + w[j-1][c+1][e][l_] );

		//case 16 - all four stacked, no gaps
								//16	1	1	1	1		2	2	0	0
								en1 = min( en1, w[c][i+2][k_][e] + w[j-2][c+1][e][l_]
										+ 2*(nSeq1+nSeq2)*data->MBLFreeBase
										+ edangle5_pf(l,k,l-1,pf2,data)
										+ edangle3_pf(k,l,k+1,pf2,data)
										+ edangle5_pf(j,i,j-1,pf1,data)
										+ edangle3_pf(i,j,i+1,pf1,data) );

		//case 6 - j and l stacked:
								//6		0	1	0	1		1	1	0	0
								en1 = min(en1,
										w[c][i+1][k_][e] + w[j-2][c+1][e][l_]
										+ (nSeq1+nSeq2)*data->MBLFreeBase
										+ edangle5_pf(j,i,j-1,pf1,data)
										+ edangle5_pf(l,k,l-1,pf2,data));

		//case 11 - i and k stacked
									//11	1	0	1	0		1	1	0	0
							en1 = min( en1,
										w[c][i+2][k_][e] + w[j-1][c+1][e][l_]
										+ (nSeq1+nSeq2)*data->MBLFreeBase
										+ edangle3_pf(k,l,k+1,pf2,data)
										+ edangle3_pf(i,j,i+1,pf1,data) );

								if (l_ -1>=0) {//boundary check for the w array
		//case 2 - stack on l
								//2		0	0	0	1		0	1	1	0
									en1 = min( en1,
										w[c][i+1][k_][e] + w[j-1][c+1][e][l_ -1]
										+ nSeq2*data->MBLFreeBase
										+ edangle5_pf(l,k,l-1,pf2,data)
										+ nSeq1*gap );

		//case 12 - i, k, and l stacked
								//12	1	0	1	1		1	2	1	0
									en1 = min( en1,
										w[c][i+2][k_][e] + w[j-1][c+1][e][l_ -1]
										+ (nSeq1 + 2*nSeq2)*data->MBLFreeBase
										+ edangle5_pf(l,k,l-1,pf2,data)
										+ edangle3_pf(i,j,i+1,pf1,data)
										+ edangle3_pf(k,l,k+1,pf2,data)
										+ nSeq1*gap );

									if (k_ -1>0) {
		//case 10 - l and i stacked
								//10	1	0	0	1		1	1	1	1
										en1 = min( en1,
											w[c][i+2][k_ -1][e] + w[j-1][c+1][e][l_ -1]
											+ (nSeq1+nSeq2)*data->MBLFreeBase
											+ edangle5_pf(l,k,l-1,pf2,data)
											+ edangle3_pf(i,j,i+1,pf1,data)
											+ (nSeq1+nSeq2)*gap );
									}

									if (k_ +1<maxrange) {
		//case 4 - k and l stacked
								//4		0	0	1	1		0	2	2	0
										en1 = min( en1,
											w[c][i+1][k_ +1][e] + w[j-1][c+1][e][l_ -1]
											+ 2*nSeq2*data->MBLFreeBase
											+ edangle5_pf(l,k,l-1,pf2,data)
											+ edangle3_pf(k,l,k+1,pf2,data)
											+ 2*nSeq1*gap );
									}
								}

								if (k_ -1>0) {

		//case 14 - i, j, and l stacked
								//14	1	1	0	1		2	1	0	1
									en1 = min( en1,
										w[c][i+2][k_ -1][e] + w[j-2][c+1][e][l_]
										+ (2*nSeq1 + nSeq2)*data->MBLFreeBase
										+ edangle5_pf(l,k,l-1,pf2,data)
										+ edangle5_pf(j,i,j-1,pf1,data)
										+ edangle3_pf(i,j,i+1,pf1,data)
										+ nSeq2*gap );

		//case 9 - i stacked
								//9		1	0	0	0		1	0	0	1
									en1 = min( en1,
										w[c][i+2][k_ -1][e] + w[j-1][c+1][e][l_]
										+ nSeq1*data->MBLFreeBase
										+ edangle3_pf(i,j,i+1,pf1,data)
										+ nSeq2*gap );

									if (l_ +1<maxrange) {
		//case 13 - i and j stacked
								//13	1	1	0	0		2	0	0	2
										en1 = min( en1,
										w[c][i+2][k_ -1][e] + w[j-2][c+1][e][l_ +1]
										+ 2*nSeq1*data->MBLFreeBase
										+ edangle5_pf(j,i,j-1,pf1,data)
										+ edangle3_pf(i,j,i+1,pf1,data)
										+ 2*nSeq2*gap );
									}
								}

								if (k_ +1<maxrange) {
		//case 8 - j, k, and l stacked
								//8		0	1	1	1		1	2	1	0
									en1 = min( en1,
										w[c][i+1][k_ +1][e] + w[j-2][c+1][e][l_]
										+ (nSeq1 + 2*nSeq2)*data->MBLFreeBase
										+ edangle5_pf(j,i,j-1,pf1,data)
										+ edangle5_pf(l,k,l-1,pf2,data)
										+ edangle3_pf(k,l,k+1,pf2,data)
										+ nSeq1*gap );

		//case 3 - stack on k
								//3		0	0	1	0		0	1	1	0
									en1 = min( en1,
										w[c][i+1][k_ +1][e] + w[j-1][c+1][e][l_]
										+ nSeq2*data->MBLFreeBase
										+ edangle3_pf(k,l,k+1,pf2,data)
										+ nSeq1*gap );

									if (l_ +1<maxrange) {
		//case 7 - j and k stacked
								//7		0	1	1	0		1	1	1	1
										en1 = min( en1,
											w[c][i+1][k_ +1][e] + w[j-2][c+1][e][l_ +1]
											+ (nSeq1+nSeq2)*data->MBLFreeBase
											+ edangle5_pf(j,i,j-1,pf1,data)
											+ edangle3_pf(k,l,k+1,pf2,data)
											+ (nSeq1+nSeq2)*gap );
									}
								}

								if (l_ +1<maxrange) {
		//case 15 - i, j, and k stacked
								//15	1	1	1	0		2	1	0	1
									en1 = min( en1,
										w[c][i+2][k_][e] + w[j-2][c+1][e][l_ +1]
										+ (2*nSeq1+nSeq2)*data->MBLFreeBase
										+ edangle3_pf(k,l,k+1,pf2,data)
										+ edangle5_pf(j,i,j-1,pf1,data)
										+ edangle3_pf(i,j,i+1,pf1,data)
										+ nSeq2*gap );

		//case 5 - j stacked
								//5		0	1	0	0		1	0	0	1
									en1 = min( en1,
										w[c][i+1][k_][e] + w[j-2][c+1][e][l_ +1]
										+ nSeq1*data->MBLFreeBase
										+ edangle5_pf(j,i,j-1,pf1,data)
										+ nSeq2*gap );
								}
							}	//end of for(d)
						}	//end of for(c)

						en1 = en1 + penalty_pf(i,j,pf1,data)
								+ penalty_pf(k,l,pf2,data)
								+ (nSeq1+nSeq2)*data->MBLClose
								+ (nSeq1+nSeq2)*data->MBLHelix;
							//adjust en1 by a penalty value if either i,j,k,l is a "u"
							//to acknowledge starting a helix with a weak basepair;
							//and, add on the values for closing a multi-branch loop
							//while creating an additional helix within the MBLoop

						v[j][i][k_][l_]=min(en1, v[j][i][k_][l_]);
							//minimum of best V3 case against the best from V1 and V2

					}	// end of if

					else v[j][i][k_][l_] = INFINITE;
						//i:j and/or k:l could not be part of a helix


/* W Cases: best energy for a fragment in a multiloop --  w[j][i][k_][l_] */

/* Case W1: The alignment of i<->j and k<->l is created by extending a pre-calculated alignment
			from the w array (i.e. not in the middle of a helix), by adding at most one
			nucleotide in each direction
			*/
			/* 	The 15 cases are based on which end(s) of the two subsequences
				have been extended.
				The energy is based on the energy of alignment of the two subsequences,
				plus penalty for extra bases in a multi-branch loop, plus penalty
				for added gaps if one subsequence is extended in a particular direction,
				but not the other.

				Note: MBLFreeBase is the multibranch-loop free-base penalty

				Recall that when i or j is reduced/increased, then
					k_ or l_ is reduced/increased relative to k or l,
					and so must be adjusted to keep k/l constant.

				*/
						//1 indicates nucleotide (i,j,k,or l) has been "added on"
						//0 indicates that it is part of the pre-calculated fragment:
						//all added-on bases are extra free bases within loop
						//bases added at a particular end on one sequence
						//but not the other increase gap penalty

						//		i	j	k	l		bases	gaps
						//				Sequence:	1	2	1	2
						//(no case 1 -- must extend by at least one base)
						//2		0	0	0	1		0	1	1	0
						//3		0	0	1	0		0	1	1	0
						//4		0	0	1	1		0	2	2	0
						//5		0	1	0	0		1	0	0	1
						//6		0	1	0	1		1	1	0	0
						//7		0	1	1	0		1	1	1	1
						//8		0	1	1	1		1	2	1	0
						//9		1	0	0	0		1	0	0	1
						//10	1	0	0	1		1	1	1	1
						//11	1	0	1	0		1	1	0	0
						//12	1	0	1	1		1	2	1	0
						//13	1	1	0	0		2	0	0	2
						//14	1	1	0	1		2	1	0	1
						//15	1	1	1	0		2	1	0	1
						//16	1	1	1	1		2	2	0	0

		//case 6		//6		0	1	0	1		1	1	0	0
					en1 = w[j-1][i][k_][l_]
						+ (nSeq1+nSeq2)*data->MBLFreeBase;

		//case 11		//11	1	0	1	0		1	1	0	0
					en1 = min(en1,
								w[j][i+1][k_][l_]
									+ (nSeq1+nSeq2)*data->MBLFreeBase);

		//case 16		//16	1	1	1	1		2	2	0	0
					if (j-1 > i+1)
					{	en1 = min(en1,
								w[j-1][i+1][k_][l_]
									+ 2*(nSeq1+nSeq2)*data->MBLFreeBase);
					}

					if (k_ >=1) {
		//case 9		//9		1	0	0	0		1	0	0	1
						en1=min(en1,
									w[j][i+1][k_ -1][l_]
									+ nSeq1*data->MBLFreeBase
									+ nSeq2*gap);

		//case 14		//14	1	1	0	1		2	1	0	1
						if (j-1>i+1) {
							en1 = min(en1,
									w[j-1][i+1][k_ -1][l_]
									+ (2*nSeq1+nSeq2)*data->MBLFreeBase
									+ nSeq2*gap);

							if (l_ +1<(maxrange)) {
		//case 13		//13	1	1	0	0		2	0	0	2
								en1 = min(en1,
										w[j-1][i+1][k_ -1][l_ +1]
									+ 2*nSeq1*data->MBLFreeBase
									+ 2*nSeq2*gap);
							}
						}

						if (l_ >=1) {
		//case 10		//10	1	0	0	1		1	1	1	1
							en1 = min(en1,
									w[j][i+1][k_ -1][l_ -1]
									+ (nSeq1+nSeq2)*data->MBLFreeBase
									+ (nSeq1+nSeq2)*gap);
						}
					}

					if (l_ +1<(maxrange)) {
		//case 5		//5		0	1	0	0		1	0	0	1
						en1 = min(en1,
									w[j-1][i][k_][l_ +1]
									+ nSeq1*data->MBLFreeBase
									+ nSeq2*gap);

						if (k_+1<(maxrange)) {
		//case 7		//7		0	1	1	0		1	1	1	1
							en1=min(en1,
									w[j-1][i][k_ +1][l_ +1]
									+ (nSeq1+nSeq2)*data->MBLFreeBase
									+ (nSeq1+nSeq2)*gap);
						}

						if (j-1>i+1) {
		//case 15		//15	1	1	1	0		2	1	0	1
							en1 = min(en1,
									w[j-1][i+1][k_][l_ +1]
									+ (2*nSeq1+nSeq2)*data->MBLFreeBase
									+ nSeq2*gap);
						}

					}

					if (k_ +1<(maxrange)) {
		//case 3		//3		0	0	1	0		0	1	1	0
						en1 = min(en1,
									w[j][i][k_ +1][l_]
									+ nSeq2*data->MBLFreeBase
									+ nSeq1*gap);

		//case 8		//8		0	1	1	1		1	2	1	0
						en1 = min(en1,w[j-1][i][k_ +1][l_]
									+ (nSeq1+2*nSeq2)*data->MBLFreeBase
									+ nSeq1*gap);

						if (l_ >=1) {
		//case 4		//4		0	0	1	1		0	2	2	0
							en1 = min(en1,
									w[j][i][k_ +1][l_ -1]
									+ 2*nSeq2*data->MBLFreeBase
									+ 2*nSeq1*gap);
						}
					}

					if (l_ >=1) {
		//case 2		//2		0	0	0	1		0	1	1	0

						en1=min(en1,
									w[j][i][k_][l_ -1]
									+ nSeq2*data->MBLFreeBase
									+ nSeq1*gap);

		//case 12		//12	1	0	1	1		1	2	1	0
						en1= min(en1,
									w[j][i+1][k_][l_ -1]
									+ (nSeq1 + 2*nSeq2)*data->MBLFreeBase
									+ nSeq1*gap);
					}


/* Case W2: Create a multibranch loop by terminating a helix -- take a value from
			the v array, and add 0 or 1 nucleotide on each end of each sequence
			*/

					/*	Must calculate cost of ending helix /starting internal loop
						as well as costs/benefits relating to dangling ends, freebases
						within an interior loop, and any gaps induced by extending
						in one sequence but not the other
						The minimum for W2 sub-cases is calculated without the
						internal-loop helix penalty (which is the same for all sub-cases),
						and then this is added and the minimum compared against the W1 minimum
					*/
					//There are 16 cases, depending on which of the four nucleotides
					//(i,j,k,l) are unpaired, stacked on the end of the pre-calculated
					//helix:  0 unstacked, 1 stacked
					//if they are unpaired, that counts as an MBL Free base
					//if one end has an unpaired base in one sequence but not the other
					//than the other sequence has a gap
					//							Bases	Gaps
					//						Seq	1	2	1	2
					//		i	j	k	l
					//1		0	0	0	0		0	0	0	0
					//2		0	0	0	1		0	1	1	0
					//3		0	0	1	0		0	1	1	0
					//4		0	0	1	1		0	2	2	0
					//5		0	1	0	0		1	0	0	1
					//6		0	1	0	1		1	1	0	0
					//7		0	1	1	0		1	1	1	1
					//8		0	1	1	1		1	2	1	0
					//9		1	0	0	0		1	0	0	1
					//10	1	0	0	1		1	1	1	1
					//11	1	0	1	0		1	1	0	0
					//12	1	0	1	1		1	2	1	0
					//13	1	1	0	0		2	0	0	2
					//14	1	1	0	1		2	1	0	1
					//15	1	1	1	0		2	1	0	1
					//16	1	1	1	1		2	2	0	0


					//note:
					//	k_ = k-i+maxsep;
					//	l_ = l-j+maxsep;

					//	so when addressing i+1 => address k_ -1 to keep k unchanged
					//	and simlarly: j-1 => l_ +1 to keep l unchanged

					/* Note:	MBLClose is the multi-branch loop offset
								MBLFreeBase is the multi-branch loop free base penalty
								MBLHelix is the multi-branch loop helix penalty
								penalty_pf() calculates whether an "au/gu"
									helix-termini penalty is appropriate
					*/

	//case 1 - nothing stacked:
					//1		0	0	0	0		0	0	0	0
					en2 =
								v[j][i][k_][l_]
								+ penalty_pf(i,j,pf1,data)
								+ penalty_pf(k,l,pf2,data) ;

	//case 6 - j and l stacked
					//6		0	1	0	1		1	1	0	0
					en2 = min( en2,
								v[j-1][i][k_][l_]
								+ (nSeq1+nSeq2)*data->MBLFreeBase
								+ edangle3_pf(j-1,i,j,pf1,data)
								+ edangle3_pf(l-1,k,l,pf2,data)
								+ penalty_pf(i,j-1,pf1,data)
								+ penalty_pf(k,l-1,pf2,data) );

	//case 11 - i and k stacked
					//11	1	0	1	0		1	1	0	0
					en2 = min( en2,
								v[j][i+1][k_][l_]
								+ (nSeq1+nSeq2)*data->MBLFreeBase
								+ edangle5_pf(k+1,l,k,pf2,data)
								+ edangle5_pf(i+1,j,i,pf1,data)
								+ penalty_pf(i+1,j,pf1,data)
								+ penalty_pf(k+1,l,pf2,data) );

	//case 16 - i, j, k, and l stacked
					//16	1	1	1	1		2	2	0	0
					if (j-1>i+1)
						en2 = min( en2, v[j-1][i+1][k_][l_]
								+ 2*(nSeq1+nSeq2)*data->MBLFreeBase +
								+ edangle3_pf(l-1,k+1,l,pf2,data)
								+ edangle5_pf(k+1,l-1,k,pf2,data)
								+ edangle3_pf(j-1,i+1,j,pf1,data)
								+ edangle5_pf(i+1,j-1,i,pf1,data)
								+ penalty_pf(i+1,j-1,pf1,data)
								+ penalty_pf(k+1,l-1,pf2,data));

					if (l_ -1>=0) {
	//case 2 - l stacked
					//2		0	0	0	1		0	1	1	0
						en2 = min( en2,
								v[j][i][k_][l_ -1]
								+ nSeq2*data->MBLFreeBase
								+ edangle3_pf(l-1,k,l,pf2,data)
								+ nSeq1*gap
								+ penalty_pf(i,j,pf1,data)
								+ penalty_pf(k,l-1,pf2,data) );

						if (k_ +1<maxrange) {
	//case 4 - l and k stacked
					//4		0	0	1	1		0	2	2	0
							en2 = min( en2,
								v[j][i][k_ +1][l_ -1]
								+ 2*nSeq2*data->MBLFreeBase
								+ edangle3_pf(l-1,k+1,l,pf2,data)
								+ edangle5_pf(k+1,l-1,k,pf2,data)
								+ 2*nSeq1*gap
								+ penalty_pf(i,j,pf1,data)
								+ penalty_pf(k+1,l-1,pf2,data) );
						}
					}

					if (k_ +1<maxrange) {
	//case 8 - j, k, and l stacked:
					//8		0	1	1	1		1	2	1	0
						en2 = min( en2,
								v[j-1][i][k_ +1][l_]
								+ (nSeq1 + 2*nSeq2)*data->MBLFreeBase
								+ edangle3_pf(j-1,i,j,pf1,data)
								+ edangle3_pf(l-1,k+1,l,pf2,data)
								+ edangle5_pf(k+1,l-1,k,pf2,data)
								+ nSeq1*gap
								+ penalty_pf(i,j-1,pf1,data)
								+ penalty_pf(k+1,l-1,pf2,data) );

	//case 3 - k stacked
					//3		0	0	1	0		0	1	1	0
						en2 = min( en2,
								v[j][i][k_ +1][l_]
								+ nSeq2*data->MBLFreeBase
								+ edangle5_pf(k+1,l,k,pf2,data)
								+ nSeq1*gap
								+penalty_pf(i,j,pf1,data)
								+ penalty_pf(k+1,l,pf2,data) );

						if (l_ +1<maxrange) {
	//case 7 - j and k stacked
					//7		0	1	1	0		1	1	1	1
							en2 = min( en2,
								v[j-1][i][k_ +1][l_ +1]
								+ (nSeq1+nSeq2)*data->MBLFreeBase
								+ edangle3_pf(j-1,i,j,pf1,data)
								+ edangle5_pf(k+1,l,k,pf2,data)
								+ (nSeq1+nSeq2)*gap
								+ penalty_pf(i,j-1,pf1,data)
								+ penalty_pf(k+1,l,pf2,data) );
						}
					}

					if (l_ +1<maxrange) {
	//case 5 - j stacked
					//5		0	1	0	0		1	0	0	1
						en2 = min( en2,
								v[j-1][i][k_][l_ +1]
								+ nSeq1*data->MBLFreeBase
								+ edangle3_pf(j-1,i,j,pf1,data)
								+ nSeq2*gap
								+ penalty_pf(i,j-1,pf1,data)
								+ penalty_pf(k,l,pf2,data));

	//case 15 - i,j, and k stacked
					//15	1	1	1	0		2	1	0	1
						if (j-1>i+1) {
							en2 = min( en2,
								v[j-1][i+1][k_][l_ +1]
								+ (2*nSeq1+nSeq2)*data->MBLFreeBase
								+ edangle5_pf(k+1,l,k,pf2,data)
								+ edangle3_pf(j-1,i+1,j,pf1,data)
								+ edangle5_pf(i+1,j-1,i,pf1,data)
								+ nSeq2*gap
								+ penalty_pf(i+1,j-1,pf1,data)
								+ penalty_pf(k+1,l,pf2,data));

	//case 13 - i and j stacked
					//13	1	1	0	0		2	0	0	2
							if (k_>=1) {
								en2 = min( en2,
								v[j-1][i+1][k_ -1][l_ +1]
								+ 2*nSeq1*data->MBLFreeBase
								+ edangle3_pf(j-1,i+1,j,pf1,data)
								+ edangle5_pf(i+1,j-1,i,pf1,data)
								+ 2*nSeq2*gap
								+ penalty_pf(i+1,j-1,pf1,data)
								+ penalty_pf(k,l,pf2,data));
							}
						}
					}

					if (k_>=1) {
	//case 9 - i alone is stacked
					//9		1	0	0	0		1	0	0	1
						en2 = min( en2,
								v[j][i+1][k_ -1][l_]
								+ nSeq1*data->MBLFreeBase
								+ edangle5_pf(i+1,j,i,pf1,data)
								+ nSeq2*gap
								+ penalty_pf(i+1,j,pf1,data)
								+ penalty_pf(k,l,pf2,data) );

	//case 14 - i, j, and l stacked
					//14	1	1	0	1		2	1	0	1
						if (j-1>i+1) {
							en2 = min( en2,
								v[j-1][i+1][k_ -1][l_]
								+ (2*nSeq1+nSeq2)*data->MBLFreeBase
								+ edangle3_pf(l-1,k,l,pf2,data)
								+ edangle3_pf(j-1,i+1,j,pf1,data)
								+ edangle5_pf(i+1,j-1,i,pf1,data)
								+ nSeq2*gap
								+ penalty_pf(i+1,j-1,pf1,data)
								+ penalty_pf(k,l-1,pf2,data) );
						}
					}

					if (l_ -1>=0) {
	//case 12 - i, k, and l stacked:
					//12	1	0	1	1		1	2	1	0
							en2 = min( en2,
								v[j][i+1][k_][l_ -1]
								+ (nSeq1 + 2*nSeq2)*data->MBLFreeBase
								+ edangle3_pf(l-1,k+1,l,pf2,data)
								+ edangle5_pf(i+1,j,i,pf1,data)
								+ edangle5_pf(k+1,l-1,k,pf2,data)
								+ nSeq1*gap
								+ penalty_pf(i+1,j,pf1,data)
								+ penalty_pf(k+1,l-1,pf2,data) );

						if (k_>=1) {
	//case 10 - l and i stacked
					//10	1	0	0	1		1	1	1	1
							en2 = min( en2,
								v[j][i+1][k_ -1][l_ -1]
								+ (nSeq1+nSeq2)*data->MBLFreeBase
								+ edangle3_pf(l-1,k,l,pf2,data)
								+ edangle5_pf(i+1,j,i,pf1,data)
								+ (nSeq1+nSeq2)*gap
								+ penalty_pf(i+1,j,pf1,data)
								+ penalty_pf(k,l-1,pf2,data) );
						}
					}

					en1 = min(en1,
								en2 + (nSeq1+nSeq2)*data->MBLHelix
								);

/* Case W3:  merger of two pairs of aligned sub-fragments:
				i<->c aligned with k<->d plus (c+1)<->j aligned with (d+1)<->l
				*/
					for (c=i+MINLOOP; c<j-MINLOOP; c++) {
						for (d=max(k+MINLOOP, c-maxsep);
									(d < l-MINLOOP)&&(d <= c+maxsep); d++) {

							e = d-c+maxsep;
							//e is the adjusted value of d (equivalent to d_)
							en1 = min(en1, w[c][i][k_][e] + w[j][c+1][e][l_]);
						}
					}

					w[j][i][k_][l_] = en1; //overall minimum value

				}		// end of for(l)
			}		// end of for(k)
		}		// end of for(i)
	}		// end of for(size)



/*----end of w/v steps ---- */

/*----fill the mine array---- */

	lowest = INFINITE;
	//calculate exterior energy:
	//store this in mine[i][k_];
	//mine[i][k_] refers to the lowest free energy for fragment 1<->i and 1<->k with i and k aligned

	cout << "\n\n Identifying minimum energy for each partial alignment:\n";
	for (i=1;i<=nBases1;i++) {

		//Report progress:
		cout << i <<"..."<<flush;
		if (i % 10 == 0) cout <<"\n" <<flush;


		c = min(nBases2,i+maxsep);//c is the bound for k
		for (k = max(1, i-maxsep); k <= c; k++) {
			k_ = k - i + maxsep;

			//First, consider the best alignment of subsequences consisting of
			//one fewer nucleotides on each strand, then add one nucleotide
			//on each strand, unpaired and unstacked
			//(does not change the number of gaps between strand):
			en1 = mine[i-1][k_];
				//Note that since k = k_ + i - maxsep,
				//reducing i reduces k without changing k_.

			//Second, consider the best alignment of each subsequence with a
			//subsequence one nucleotide shorter on the other strand,
			//and then add a gap:
			if (k_ +1 < maxrange)
				en1 = min(en1, mine[i-1][k_ +1] + nSeq2*gap);
				//i-1, k constant -- i is aligned with a gap in profile2
			if (k_ >= 1)
				en1 = min(en1, mine[i][k_ -1] + nSeq1*gap);
				//i constant, k-1 -- k is aligned with a gap in profile1

			//Third, break the subsequence 1<->i into two pieces:
			//an unstructured region from 1<->j, plus a helix from j+1<->i
			//and consider the sum of their best alignments with pieces:
			//1<->l (unstructured), and l+1<->k (helix) in the second sequence

			for (j = 0; j+MINLOOP < i; j++) {

				//d is upper bound for l:
				d = min(min(j+maxsep-1,nBases2), k-MINLOOP-1);

				for (l = max(0,j-maxsep); l<=d; l++) {

					l_ = l - j + maxsep;

					//check all possible alignments

					//must consider whether any of the nucleotides (j+1, i, l+1, k)
					//is actually an unpaired nucleotide
					//stacked onto the end of the a helix.
					//If there is an unpaired (stacked) nucleotide on one end
					//in one sequence but not the other, add a gap penalty
					//for the other sequence
					//(Note: don't have to keep track of free bases since this is in
					//the unstructured start/end regions, not within a multi-branch loop)

					//There are 16 cases:
					//consider them in this order (0 unstacked, 1 stacked)
					//							Gaps
					//						Seq	1	2
					//		j+1	i	l+1	k
					//1		0	0	0	0		0	0
					//2		0	0	0	1		1	0
					//3		0	0	1	0		1	0
					//4		0	0	1	1		2	0
					//5		0	1	0	0		0	1
					//6		0	1	0	1		0	0
					//7		0	1	1	0		1	1
					//8		0	1	1	1		1	0
					//9		1	0	0	0		0	1
					//10	1	0	0	1		1	1
					//11	1	0	1	0		0	0
					//12	1	0	1	1		1	0
					//13	1	1	0	0		0	2
					//14	1	1	0	1		0	1
					//15	1	1	1	0		0	1
					//16	1	1	1	1		0	0

					/*****Note that for these exterior loops****
									j < i
									l < k
					 *******************************************/


					//no stacking
					//case 1:
					//1		0	0	0	0		0	0
					en1 = min(en1, mine[j][l_]
							+ v[i][j+1][l_][k_]
							+ penalty_pf(i,j+1,pf1,data)
							+ penalty_pf(k,l+1,pf2,data) );

					//case 6:
					//6		0	1	0	1		0	0
					en1 = min(en1, mine[j][l_]
							+ v[i-1][j+1][l_][k_]
							+ edangle3_pf(i-1,j+1,i,pf1,data)
							+ edangle3_pf(k-1,l+1,k,pf2,data)
							+ penalty_pf(i-1,j+1,pf1,data)
							+ penalty_pf(k-1,l+1,pf2,data) );

					//case 11
					//11	1	0	1	0		0	0
					en1 = min(en1, mine[j][l_]
							+ v[i][j+2][l_][k_]
							+ edangle5_pf(j+2,i,j+1,pf1,data)
							+ edangle5_pf(l+2,k,l+1,pf2,data)
							+ penalty_pf(i,j+2,pf1,data)
							+ penalty_pf(l+2,k,pf2,data) );

					//case 16
					//16	1	1	1	1		0	0
					en1 = min(en1, mine[j][l_]
							+ v[i-1][j+2][l_][k_]
							+ edangle5_pf(j+2,i-1,j+1,pf1,data)
							+ edangle3_pf(i-1,j+2,i,pf1,data)
							+ edangle5_pf(l+2,k-1,l+1,pf2,data)
							+ edangle3_pf(k-1,l+2,k,pf2,data)
							+ penalty_pf(i-1,j+2,pf1,data)
							+ penalty_pf(k-1,l+2,pf2,data) );

					if (k_>0) {
						//case 2:
					//2		0	0	0	1		1	0
						en1 = min(en1, mine[j][l_]
								+ v[i][j+1][l_][k_ -1]
								+ edangle3_pf(k-1,l+1,k,pf2,data)
								+ penalty_pf(i,j+1,pf1,data)
								+ penalty_pf(k-1,l+1,pf2,data)
								+ nSeq1*gap );

						//case 12
					//12	1	0	1	1		1	0
						en1 = min(en1, mine[j][l_]
								+ v[i][j+2][l_][k_ -1]
								+ edangle5_pf(j+2,i,j+1,pf1,data)
								+ edangle5_pf(l+2,k-1,l+1,pf2,data)
								+ edangle3_pf(k-1,l+2,k,pf2,data)
								+ penalty_pf(i,j+2,pf1,data)
								+ penalty_pf(k-1,l+2,pf2,data)
								+ nSeq1*gap );

						if (l_+1 < maxrange) {
							//case 4
					//4		0	0	1	1		2	0
							en1 = min(en1, mine[j][l_]
									+ v[i][j+1][l_ +1][k_ -1]
									+ edangle3_pf(k-1,l+2,k,pf2,data)
									+ edangle5_pf(l+2,k-1,l+1,pf2,data)
									+ penalty_pf(i,j+1,pf1,data)
									+ penalty_pf(k-1,l+2,pf2,data)
									+ 2*nSeq1*gap );
						}
					}

					if (l_ +1 < maxrange) {
						//case 3:
					//3		0	0	1	0		1	0
						en1 = min(en1, mine[j][l_]
								+ v[i][j+1][l_ +1][k_]
								+ edangle5_pf(l+2,k,l+1,pf2,data)
								+ penalty_pf(i,j+1,pf1,data)
								+ penalty_pf(l+2,k,pf2,data)
								+ nSeq1*gap );

						//case 8:
					//8		0	1	1	1		1	0
						en1 = min(en1, mine[j][l_]
								+ v[i-1][j+1][l_ +1][k_]
								+ edangle3_pf(i-1,j+1,i,pf1,data)
								+ edangle5_pf(l+2,k-1,l+1,pf2,data)
								+ edangle3_pf(k-1,l+2,k,pf2,data)
								+ penalty_pf(i-1,j+1,pf1,data)
								+ penalty_pf(k-1,l+2,pf2,data)
								+ nSeq1*gap );

						if (k_+1<maxrange) {
							//case 7
					//7		0	1	1	0		1	1
							en1 = min(en1, mine[j][l_]
									+ v[i-1][j+1][l_ +1][k_ +1]
									+ edangle5_pf(l+2,k,l+1,pf2,data)
									+ edangle3_pf(i-1,j+1,i,pf1,data)
									+ penalty_pf(i-1,j+1,pf1,data)
									+ penalty_pf(l+2,k,pf2,data)
									+ (nSeq1+nSeq2)*gap );
						}
					}

					if (k_+1<maxrange) {
						//case5:
					//5		0	1	0	0		0	1
						en1 = min(en1, mine[j][l_]
								+ v[i-1][j+1][l_][k_ +1]
								+ edangle3_pf(i-1,j+1,i,pf1,data)
								+ penalty_pf(i-1,j+1,pf1,data)
								+ penalty_pf(k,l+1,pf2,data)
								+ nSeq2*gap );

						//case 15
					//15	1	1	1	0		0	1
						en1 = min(en1, mine[j][l_]
								+ v[i-1][j+2][l_][k_ +1]
								+ edangle5_pf(j+2,i-1,j+1,pf1,data)
								+ edangle3_pf(i-1,j+2,i,pf1,data)
								+ edangle5_pf(l+2,k,l+1,pf2,data)
								+ penalty_pf(i-1,j+2,pf1,data)
								+ penalty_pf(l+2,k,pf2,data)
								+ nSeq2*gap );

						if (l_ >0) {
							//case 13
					//13	1	1	0	0		0	2
							en1 = min(en1, mine[j][l_]
									+ v[i-1][j+2][l_ -1][k_ +1]
									+ edangle5_pf(j+2,i-1,j+1,pf1,data)
									+ edangle3_pf(i-1,j+2,i,pf1,data)
									+ penalty_pf(i-1,j+2,pf1,data)
									+ penalty_pf(k,l+1,pf2,data)
									+ 2*nSeq2*gap );
						}
					}

					if (l_ >0) {
						//case 9
					//9		1	0	0	0		0	1
						en1 = min(en1, mine[j][l_]
								+ v[i][j+2][l_ -1][k_]
								+ edangle5_pf(j+2,i,j+1,pf1,data)
								+ penalty_pf(i,j+2,pf1,data)
								+ penalty_pf(k,l+1,pf2,data)
								+ nSeq2*gap );

						//case 14
					//14	1	1	0	1		0	1
						en1 = min(en1, mine[j][l_]
								+ v[i-1][j+2][l_ -1][k_]
								+ edangle5_pf(j+2,i-1,j+1,pf1,data)
								+ edangle3_pf(i-1,j+2,i,pf1,data)
								+ edangle3_pf(k-1,l+1,k,pf2,data)
								+ penalty_pf(i-1,j+2,pf1,data)
								+ penalty_pf(k-1,l+1,pf2,data)
								+ nSeq2*gap );

						if (k_>0) {
							//case 10
					//10	1	0	0	1		1	1
							en1 = min(en1, mine[j][l_]
									+ v[i][j+2][l_ -1][k_ -1]
									+ edangle5_pf(j+2,i,j+1,pf1,data)
									+ edangle3_pf(k-1,l+1,k,pf2,data)
									+ penalty_pf(i,j+2,pf1,data)
									+ penalty_pf(k-1,l+1,pf2,data)
									+ (nSeq1+nSeq2)*gap );
						}
					}

				}		// end of for(l)
			}		// end of for(j)

			mine[i][k_] = en1;

			//keep track of the best minimum energy for the entire sequence
			//by adjusting the minimum energy up to i,k by the number of gaps
			//required to align the rest of the sequences (as unstructed regions)
			lowest = min(lowest, en1+ gap*(
										(nBases1-i)>(nBases2-k) ?
										((nBases1-i)-(nBases2-k))*nSeq2 :
										((nBases2-k)-(nBases1-i))*nSeq1
										)
							);

		}		// end of for(k)
	}		// end of for(i)


/* ---- end of mine steps ---- */


	if (totalscore!=NULL)
		*totalscore = lowest;
		//store the minimum energy value where it can be accessed
		//by the calling function


	cout << "\nOverall minimum energy of global alignment: " << float(lowest)/10 << " kcal/mol \n\n";

	/*
	cout << "Minimum e array\n";
	for (i=1; i < nBases1; i++){
	  for(k_=0; k_ < maxrange; k_++){
	    cout << mine[i][k_] << "\t";
	  }
	  cout <<"\n";
	}
	cout <<"\n";
	*/

	if (lowest == gap*(
						(nBases1>nBases2) ?
						(nBases1-nBases2)*nSeq2 :
						(nBases2-nBases1)*nSeq1)
					) {
		//minimum energy is unstructured RNA
		//return with empty structures
		return;
	}


/* ---- traceback ---- */
	//declarations for traceback portion:
	stackclass stack;		// class stackclass defined in stackclass.cpp
	register bool found;
	register int constantclosure;

	//now traceback:
	//find all exterior fragments and place them on a stack

	found = false;
	c = d = -1; //initialize to a non-meaningful value

	cout << "Identifying secondary structure regions\n" ;

	for (i=1; (i<=nBases1) && (!found); i++) {
		for (k_=0; (k_<maxrange) && (!found); k_++) {
			k = k_ + i-maxsep;

			if ( (k > 0) && (k <= nBases2) ) {
				if ( (mine[i][k_] + gap*(
										(nBases1-i)>(nBases2-k) ?
										((nBases1-i)-(nBases2-k))*nSeq2 :
										((nBases2-k)-(nBases1-i))*nSeq1
										)
						) ==lowest) {
					/* If the lowest energy value is equal to the energy of
						the alignment of 1<->i with 1<->k, plus whatever gaps
						are necessary to align the rest of the sequences together,
						then the lowest energy structure is simply single stranded
						sequences after index i,k.
						*/
					c = i;//store the index of i with the lowest mine value
					d = k_;//store the index of k_ with the lowest mine value

					cout << "Structured region ends at " << i << " vs. " << k << "\n";

					found = true; //break the for loops
					en3 = mine[i][k_];
						//store the energy value for the alignments up to i,k
				}
			}
		}
	}

	if (!found) {
	   HDException ex (ERR, "Traceback error at start!");
	   throw ex;
	}

	i = c; //recall stored indices
	k_ = d;

/****Fill the stack with initial structures****/
/*	Determine the structural elements used to build up the best
	energy value by working backwards.  Identify aligned regions
	bound by helices and put the start and end indices of the helix
	(in both sequences) into the stack, along with the best energy
	score (from the v array) for those aligned fragments.
	(Note: does not examine structure bounded by helix).

	Once a helix has been identified, restart search upstream of it
	until the beginning of the sequence is reached.
	*/
	while(i > 0) {
		found = false;
		k = k_ + i - maxsep;
/*This section is the reverse of the steps used to fill the minE array*/

		/*	First, try reducing the length of the sequences by at most
			one nucleotide each.
			If the energy is the same (+ gap penalty if appropriate),
			then the nucleotides at i/k were unstructured.
			Reset indices and energy values for next round.
		*/
		if (en3 == mine[i-1][k_]) {
			//try reducing i and k by one:
			found = true;
			en3 = mine[i-1][k_];
			i--;
		}
		else if (k_ + 1 < maxrange) {
			if (en3 == (mine[i-1][k_+1] + nSeq2*gap) ){
				//try reducing i but not k -- i aligned with a gap in profile 2
				found = true;
				en3 = mine[i-1][k_+1];
				i--;
				k_++;
			}
		}
		if ((!found)&&( k_ >= 1)) {
			if(en3 == (mine[i][k_ -1] + nSeq1*gap) ) {
				//try reducing k but not i -- k aligned with a gap in profile 1
				found = true;
				en3 = mine[i][k_ -1];
				k_--;
			}
		}

		if (!found) cout << "Bifurcation \n";

		//Else (i.e. reduced structure not found),
		//try to find a bifurcation in the structure
		//break the sequences around, j,l, such that
		//energy value for alignment of sequences 1<->i vs 1<->k
		//can be made from a merger of the best alignment from
		//1<->j vs 1<->l, plus a fragment bound by a helix.
		//16 cases based on whether the nucleotides at
		//(j+1),i,(l+1),k are stacked on the end of the helix as
		//dangling nucleotides (as opposed to being part of the helix)
		for (j = 0; (j+MINLOOP < i)&&(!found); j++) {

		  //create a bound for l
		  d = min(min(j+maxsep - 1, nBases2), k-MINLOOP-1);

			for (l = max(0,j-maxsep);(l <= d)&&(!found);
					l++) {
				/****** NOTE: j < i, l < k ******/

				l_ = l - j+maxsep;

				//as usual 16 cases are
				//1 indicates the nucleotide is stacked
				//0 indicates the nucleotide is part of the helix termini
				//							Gaps
				//						Seq	1	2
				//		j+1	i	l+1	k
				//1		0	0	0	0		0	0
				//2		0	0	0	1		1	0
				//3		0	0	1	0		1	0
				//4		0	0	1	1		2	0
				//5		0	1	0	0		0	1
				//6		0	1	0	1		0	0
				//7		0	1	1	0		1	1
				//8		0	1	1	1		1	0
				//9		1	0	0	0		0	1
				//10	1	0	0	1		1	1
				//11	1	0	1	0		0	0
				//12	1	0	1	1		1	0
				//13	1	1	0	0		0	2
				//14	1	1	0	1		0	1
				//15	1	1	1	0		0	1
				//16	1	1	1	1		0	0


				//case 1:	bond between (j+1):i, and (l+1):k, no stacking
				//1		0	0	0	0		0	0
				if (en3 == mine[j][l_]
							+ v[i][j+1][l_][k_]
							+ penalty_pf(i,j+1,pf1,data)
							+ penalty_pf(k,l+1,pf2,data))
				{
					found = true;
					en3 = mine[j][l_];
					stack.push(j+1,i,l_,k_,v[i][j+1][l_][k_]);
					i = j;
					k_ = l_;
				}

				//case 6:	bond between (j+1):(i-1), and (l+1):(k-1),
				//			i, k stacked
				//			no gaps
				//6		0	1	0	1		0	0
				else if(en3 == mine[j][l_]
								+ v[i-1][j+1][l_][k_]
								+ edangle3_pf(i-1,j+1,i,pf1,data)
								+ edangle3_pf(k-1,l+1,k,pf2,data)
								+ penalty_pf(i-1,j+1,pf1,data)
								+ penalty_pf(k-1,l+1,pf2,data))
				{
					found = true;
					en3 = mine[j][l_];
					stack.push(j+1,i-1,l_,k_,v[i-1][j+1][l_][k_]);
					i = j;
					k_ = l_;
				}

				//case 11:	bond between (j+2):(i), and (l+2):(k),
				//			j+1, l+1 stacked
				//			no gaps
				//11	1	0	1	0		0	0
				else if (en3 == mine[j][l_]
								+ v[i][j+2][l_][k_]
								+ edangle5_pf(j+2,i,j+1,pf1,data)
								+ edangle5_pf(l+2,k,l+1,pf2,data)
								+ penalty_pf(i,j+2,pf1,data)
								+ penalty_pf(l+2,k,pf2,data))
				{
					stack.push(j+2,i,l_,k_,v[i][j+2][l_][k_]);
					found = true;
					en3=mine[j][l_];
					i = j;
					k_ = l_;
				}

				//case 16:	bond between (j+2):(i-1), and (l+2):(k-1),
				//			j+1, i, l+1, k stacked
				//			no gaps
				//16	1	1	1	1		0	0
				else if(en3 == mine[j][l_]
								+ v[i-1][j+2][l_][k_]
								+ edangle5_pf(j+2,i-1,j+1,pf1,data)
								+ edangle3_pf(i-1,j+2,i,pf1,data)
								+ edangle5_pf(l+2,k-1,l+1,pf2,data)
								+ edangle3_pf(k-1,l+2,k,pf2,data)
								+ penalty_pf(i-1,j+2,pf1,data)
								+ penalty_pf(k-1,l+2,pf2,data))
				{
					stack.push(j+2,i-1,l_,k_,v[i-1][j+2][l_][k_]);
					found = true;
					en3=mine[j][l_];
					i = j;
					k_ = l_;
				}

				else if (k_>0) {
					//case 2::	bond between (j+1):(i), and (l+1):(k-1),
					//			k stacked
					//			1 gap in prof 1
				//2		0	0	0	1		1	0
				  if (en3 == mine[j][l_]
								+ v[i][j+1][l_][k_ -1]
								+ edangle3_pf(k-1,l+1,k,pf2,data)
								+ penalty_pf(i,j+1,pf1,data)
								+ penalty_pf(k-1,l+1,pf2,data)
								+ nSeq1*gap)
					{
					  stack.push(j+1,i,l_,k_ -1,v[i][j+1][l_][k_ -1]);
						found = true;
						en3=mine[j][l_];
						i = j;
						k_ = l_;
					}

					//case 12:	bond between (j+2):(i), and (l+2):(k-1),
					//			j+1, l+1, k stacked
					//			1 gap in profile 1
				//12	1	0	1	1		1	0
					else if(en3 == mine[j][l_]
									+ v[i][j+2][l_][k_ -1]
									+ edangle5_pf(j+2,i,j+1,pf1,data)
									+ edangle5_pf(l+2,k-1,l+1,pf2,data)
									+ edangle3_pf(k-1,l+2,k,pf2,data)
									+ penalty_pf(i,j+2,pf1,data)
									+ penalty_pf(k-1,l+2,pf2,data)
									+ nSeq1*gap)
					{
						stack.push(j+2,i,l_,(k_-1),v[i][j+2][l_][k_ -1]);
						found = true;
						en3=mine[j][l_];
						i = j;
						k_ = l_;
					}

					else if (l_+1<maxrange) {
						//case 4:	bond between (j+1):(i), and (l+2):(k-1),
						//			l+1, k stacked
						//			2 gaps in prof 1
				//4		0	0	1	1		2	0
						if	(en3 == mine[j][l_]
										+ v[i][j+1][l_+1][k_ -1]
										+ edangle3_pf(k-1,l+2,k,pf2,data)
										+ edangle5_pf(l+2,k-1,l+1,pf2,data)
										+ penalty_pf(i,j+1,pf1,data)
										+ penalty_pf(k-1,l+2,pf2,data)
										+ 2*nSeq1*gap)
						{
							stack.push(j+1,i,l_+1,k_ -1,v[i][j+1][l_+1][k_ -1]);
							found = true;
							en3=mine[j][l_];
							i = j;
							k_ = l_;
						}
					}
				}		// end of else if(k_>0)

				if ((l_+1<maxrange)&&(!found)) {
					//case 3:	bond between (j+1):(i), and (l+2):(k),
					//			l+1 stacked
					//			1 gap in prof 1
				//3		0	0	1	0		1	0
					if(en3 == mine[j][l_]
								+ v[i][j+1][l_+1][k_]
								+ edangle5_pf(l+2,k,l+1,pf2,data)
								+ penalty_pf(i,j+1,pf1,data)
								+ penalty_pf(l+2,k,pf2,data)
								+ nSeq1*gap)
					{
						stack.push(j+1,i,l_+1,k_,v[i][j+1][l_+1][k_]);
						found = true;
						en3=mine[j][l_];
						i = j;
						k_ = l_;
					}

					//case 8:	bond between (j+1):(i-1), and (l+2):(k-1),
					//			i, l+1, k stacked
					//			1 gap in prof 1
				//8		0	1	1	1		1	0
					else if(en3 == mine[j][l_]
									+ v[i-1][j+1][l_+1][k_]
									+ edangle3_pf(i-1,j+1,i,pf1,data)
									+ edangle5_pf(l+2,k-1,l+1,pf2,data)
									+ edangle3_pf(k-1,l+2,k,pf2,data)
									+ penalty_pf(i-1,j+1,pf1,data)
									+ penalty_pf(k-1,l+2,pf2,data)
									+ nSeq1*gap)
					{
						stack.push(j+1,i-1,l_+1,k_,v[i-1][j+1][l_+1][k_]);
						found = true;
						en3=mine[j][l_];
						i = j;
						k_ = l_;
					}

					else if (k_+1<maxrange) {
						//case 7	bond between (j+1):(i-1), and (l+2):(k),
						//			i, l+1 stacked
						//			1 gap in each profile
				//7		0	1	1	0		1	1
						if(en3 == mine[j][l_]
									+ v[i-1][j+1][l_+1][k_+1]
									+ edangle5_pf(l+2,k,l+1,pf2,data)
									+ edangle3_pf(i-1,j+1,i,pf1,data)
									+ penalty_pf(i-1,j+1,pf1,data)
									+ penalty_pf(l+2,k,pf2,data)
									+ (nSeq1+nSeq2)*gap)
						{
							stack.push(j+1,i-1,l_+1,k_+1,v[i-1][j+1][l_+1][k_+1]);
							found = true;
							en3=mine[j][l_];
							i = j;
							k_ = l_;
						}
					}
				}		// end of if(l+1 < maxrange)

				if ((k_+1 < maxrange)&&(!found)) {
					//case5:	bond between (j+1):(i-1), and (l+1):(k),
					//			i stacked
					//			1 gap in prof 2
				//5		0	1	0	0		0	1
					if(en3 == mine[j][l_]
								+ v[i-1][j+1][l_][k_+1]
								+ edangle3_pf(i-1,j+1,i,pf1,data)
								+ penalty_pf(i-1,j+1,pf1,data)
								+ penalty_pf(k,l+1,pf2,data)
								+ nSeq2*gap)
					{
						stack.push(j+1,i-1,l_,k_+1,v[i-1][j+1][l_][k_+1]);
						found = true;
						en3=mine[j][l_];
						i = j;
						k_ = l_;
					}

					//case 15	bond between (j+2):(i-1), and (l+2):(k),
					//			j+1, i, l+1, stacked
					//			1 gap in prof 2
				//15	1	1	1	0		0	1
					else if(en3 == mine[j][l_]
									+ v[i-1][j+2][l_][k_+1]
									+ edangle5_pf(j+2,i-1,j+1,pf1,data)
									+ edangle3_pf(i-1,j+2,i,pf1,data)
									+ edangle5_pf(l+2,k,l+1,pf2,data)
									+ penalty_pf(i-1,j+2,pf1,data)
									+ penalty_pf(l+2,k,pf2,data)
									+ nSeq2*gap)
					{

						stack.push(j+2,i-1,l_,k_+1,v[i-1][j+2][l_][k_+1]);
						found = true;
						en3=mine[j][l_];
						i = j;
						k_ = l_;
					}

					else if (l_>0) {
						//case 13	bond between (j+2):(i-1), and (l+1):(k),
						//			j+1, i stacked
						//			2 gaps in prof 2
				//13	1	1	0	0		0	2
						if(en3 == mine[j][l_]
									+ v[i-1][j+2][l_-1][k_+1]
									+ edangle5_pf(j+2,i-1,j+1,pf1,data)
									+ edangle3_pf(i-1,j+2,i,pf1,data)
									+ penalty_pf(i-1,j+2,pf1,data)
									+ penalty_pf(k,l+1,pf2,data)
									+ 2*nSeq2*gap)
						{
							stack.push(j+2,i-1,l_-1,k_+1,v[i-1][j+2][l_-1][k_+1]);
							found = true;
							en3=mine[j][l_];
							i = j;
							k_ = l_;
						}
					}
				}		// end of if(k_+1<maxrange&&!found)

				if (l_>0&&!found) {
					//case 9	bond between (j+2):(i), and (l+1):(k),
					//			j+1 stacked
					//			1 gap in prof 2
				//9		1	0	0	0		0	1
					if(en3 == mine[j][l_]
								+ v[i][j+2][l_-1][k_]
								+ edangle5_pf(j+2,i,j+1,pf1,data)
								+ penalty_pf(i,j+2,pf1,data)
								+ penalty_pf(k,l+1,pf2,data)
								+ nSeq2*gap)
					{
						stack.push(j+2,i,l_-1,k_,v[i][j+2][l_-1][k_]);
						found = true;
						en3=mine[j][l_];
						i = j;
						k_ = l_;
					}

					//case 14	bond between (j+2):(i-1), and (l+1):(k-1),
					//			j+1, i, k stacked
					//			1 gap in prof 2
				//14	1	1	0	1		0	1
					else if(en3 == mine[j][l_]
									+ v[i-1][j+2][l_-1][k_]
									+ edangle5_pf(j+2,i-1,j+1,pf1,data)
									+ edangle3_pf(i-1,j+2,i,pf1,data)
									+ edangle3_pf(k-1,l+1,k,pf2,data)
									+ penalty_pf(i-1,j+2,pf1,data)
									+ penalty_pf(k-1,l+1,pf2,data)
									+ nSeq2*gap)
					{
						stack.push(j+2,i-1,l_-1,k_,v[i-1][j+2][l_-1][k_]);
						found = true;
						en3=mine[j][l_];
						i = j;
						k_ = l_;
					}

					else if (k_>0) {
						//case 10	bond between (j+2):(i), and (l+2):(k-1),
						//			j+1, l+1 stacked
						//			1 gap in each profile
				//10	1	0	0	1		1	1
						if(en3 == mine[j][l_]
									+ v[i][j+2][l_-1][k_-1]
									+ edangle5_pf(j+2,i,j+1,pf1,data)
									+ edangle3_pf(k-1,l+1,k,pf2,data)
									+ penalty_pf(i,j+2,pf1,data)
									+ penalty_pf(k-1,l+1,pf2,data)
									+ (nSeq1+nSeq2)*gap)
						{
							stack.push(j+2,i,l_-1,k_-1, v[i][j+2][l_-1][k_-1]);
							found = true;
							en3=mine[j][l_];
							i = j;
							k_ = l_;
						}
					}
				}		// end of if(l_>0&&!found)

			}		// end of for(l)
		}		// end of for(j)

		if(!found) {//bifurcation was never found--issue error
			HDException ex (ERR, "Bifurcation not found!");
			throw ex;
		}

		cout << "Start of structure: " << i << " vs. " << (k_ + i - maxsep) << "\n";
		//cout << " -- energy equal to " << en3 << "\n";


	}// end of while(i>0)

/*Retrieve structure elements from the stack, and characterize them in detail*/

/*	The initial value pulled from the stack will be the most upstream helix.
	This is then characterized by comparing the score for the aligned fragments
	against scores for the various cases of building that fragment from
	shorter fragments.  The characterization is made recursive by adding the inner
	fragments, and their energy values, to the stack.  Recursion stops each
	time a hairpin loop is reached, at which time the next element pulled from
	the stack will be the closest downstream element.
	*/
	cout << "\nIdentifying aligned basepairs:\n";
	while (stack.pull(&i,&j, &k_, &l_, &en3)) {
		ASSERT(&i != NULL && i<=pf1->numofbases);
		ASSERT(&j != NULL && j<=pf1->numofbases);
		ASSERT(&k_ != NULL && k_<=maxrange);
		ASSERT(&l_ != NULL && l_<=maxrange);
		ASSERT(stack.getSize() >= 0);

		if (en3 >= INFINITE ) {
			HDException ex (ERR, "Value of 'INFINITE' too small for calculation parameters.");
			throw ex;

		}

		k = k_ +i-maxsep;
		l = l_ +j-maxsep;


		if (en3 == v[j][i][k_][l_]) {
			cout << i << " : "<<j<<" vs. "<<k<<" : "<<l<<"\n"<<flush;
			//i and j are paired and aligned to k and l
			pf1->basepr[i]=j;
			pf1->basepr[j]=i;
			pf2->basepr[k]=l;
			pf2->basepr[l]=k;
			align[i].pos = k;
			align[j].pos =l;
			if( i < j ) {
				strcpy(align[i].indicator, "(");
				strcpy(align[j].indicator, ")");
			}
			else {
				strcpy(align[i].indicator, ")");
				strcpy(align[j].indicator, "(");
			}
/*Reverse V cases:*/
			//now find the next pair, by recreating all the cases used
			//to fill in the v array.

/*Reverse V1 (hairpin -- no internal structure):*/
			if (en3 != ehairpin_pf(i,j,pf1,data) + ehairpin_pf(k,l,pf2,data)
					+ gap*((j-i)<(l-k)?((l-k)-(j-i))*nSeq1:((j-i)-(l-k))*nSeq2) )
			{
				//the fragment does not close hairpins, therefore, the internal
				// fragment does need to be characterized
				found= false;

/*Reverse V2: single internal helix (internal loop/helix extension/bulge)*/

				for (c=i+1; (c <= i+MAXLOOP)&&(c < j-MINLOOP)&&(!found); c++) {
					for (d=j-1; (d >= j-MAXLOOP)&&(d > c+MINLOOP)&&(!found); d--) {

						//break if c:d not a valid basepair
					    if (!get_inc(c,d,pf1)) break; //new!

						for (e=max(k+1,c-maxsep);
								(e <= k+MAXLOOP)&&(e < l-MINLOOP)&&((e-c)<=maxsep)&&(!found);
								e++) {

							for (f=min(l-1,d+maxsep);
									(f >= l-MAXLOOP)&&(f > e+MINLOOP)&&((d-f)<=maxsep)&&(!found);
									f--) {

								//break if e:f not a valid basepair
							  	if (!get_inc(e,f, pf2)) break; //new!

								if ((c==i+1)&&(d==j-1)) {
									//seq1 is helical stack
									en1 = ebp_pf(i,j,i+1,j-1,pf1, data);
								}
								else {
									//internal loop/bulge
									en1 = einternal_pf(i,j,c,d,pf1,data);
								}

								if ((e==k+1)&&(f==l-1)) {
									//seq2 is helical stack
									en2 = ebp_pf(k,l,k+1,l-1,pf2, data);
								}
								else {
									//internal loop
									en2 = einternal_pf(k,l,e,f,pf2,data);
								}

								if(en3 ==  en1 + en2 + v[d][c][e-c+maxsep][f-d+maxsep] +
										+ gap*((c-i)<(e-k)?((e-k)-(c-i))*nSeq1:((c-i)-(e-k))*nSeq2)
										+ gap*((j-d)<(l-f)?((l-f)-(j-d))*nSeq1:((j-d)-(l-f))*nSeq2))
								{
									found = true;
									stack.push(c,d,e-c+maxsep,f-d+maxsep,v[d][c][e-c+maxsep][f-d+maxsep]);
									//push the internal helix into the stack
								}
								else if (singleinsert) {
									if ( (c == i+2)&&(d == j-2)&&(j-2 > i+2)&&
											get_inc(i+1, j-1, pf1) && get_inc(i+2, j-2, pf1) &&
											(e == k+1)&&(f == l-1)&& get_inc(k+1, l-1, pf2) )
									{
										en1 = ebp_pf(i,j,i+1,j-1,pf1,data) +
												ebp_pf(i+1,j-1,i+2,j-2,pf1,data);

										if (en3 == en1+en2+v[d][c][e-c+maxsep][f-d+maxsep]+
												2*nSeq2*gap )
										{
											//base pairs with single base pair insertion in ct1
											pf1->basepr[i+1] = j-1;
											pf1->basepr[j-1] = i+1;
											stack.push(c,d,e-c+maxsep,f-d+maxsep,v[d][c][e-c+maxsep][f-d+maxsep]);
											found = true;
										}
									}
									else if ((e==k+2)&&(f==l-2)&&(l-2>k+2)&&(c==i+1)&&(d==j-1)&&
											get_inc(k+1, l-1, pf2) && get_inc(k+2, l-2, pf2) &&
											get_inc(i+1, j-1, pf1))
									{
										en2= ebp_pf(k,l,k+1,l-1,pf2,data) +
												ebp_pf(k+1,l-1,k+2,l-2,pf2,data);

										if (en3==en1+en2+v[d][c][e-c+maxsep][f-d+maxsep]+
												2*nSeq1*gap )
										{
											//base pairs with single base pair insertion in ct2
											pf2->basepr[k+1]=l-1;
											pf2->basepr[l-1]=k+1;
											stack.push(c,d,e-c+maxsep,f-d+maxsep,v[d][c][e-c+maxsep][f-d+maxsep]);
											found = true;
										}
									}
								}

							}		// end of for(f)
						}		// end of for(e)
					}		// end of for(d)
				}		// end of for(c)

/*Reverse V3 (multiloop junction):*/
				//If still not found, then innner structure must be a multiloop
				//junction closed by i-j pair aligned with k and l.
				//calculate the free energy of 2 fragments merged:
				for (c=i+MINLOOP+1; (c<j-MINLOOP)&&(!found); c++) {
					for (d=max(k+MINLOOP+1,c-maxsep); (d<l-MINLOOP)&&(d<=c+maxsep)&&(!found); d++) {

						constantclosure = penalty_pf(i,j,pf1,data)
											+ penalty_pf(k,l,pf2,data)
											+ (nSeq1+nSeq2)*data->MBLClose
											+ (nSeq1+nSeq2)*data->MBLHelix;
							//The net energy change for closing the multi-branch loop,
							//adding an extra helix on the multi-branch loop,
							//and starting a helix with the pairs (i:j) and (k:l).
							//This value will be part of every energy value in this set

						e = d-c+maxsep; //adjusted value of d for accessing arrays

						//must consider whether an unpaired nucleotide
						//is stacked onto each of the nucleotides in a helix
						//There are 16 cases:
						//consider them in this order (0 unstacked, 1 stacked)
						//						FreeBases	Gaps
						//						Seq	1	2	1	2
						//		i	j	k	l
						//1		0	0	0	0		0	0	0	0
						//2		0	0	0	1		0	1	1	0
						//3		0	0	1	0		0	1	1	0
						//4		0	0	1	1		0	2	2	0
						//5		0	1	0	0		1	0	0	1
						//6		0	1	0	1		1	1	0	0
						//7		0	1	1	0		1	1	1	1
						//8		0	1	1	1		1	2	1	0
						//9		1	0	0	0		1	0	0	1
						//10	1	0	0	1		1	1	1	1
						//11	1	0	1	0		1	1	0	0
						//12	1	0	1	1		1	2	1	0
						//13	1	1	0	0		2	0	0	2
						//14	1	1	0	1		2	1	0	1
						//15	1	1	1	0		2	1	0	1
						//16	1	1	1	1		2	2	0	0

						//case 1 - no stacks
						//1		0	0	0	0		0	0	0	0
						if(en3== w[c][i+1][k_][e] + w[j-1][c+1][e][l_]
								+ constantclosure)
						{
							stack.push(i+1,c,k_,e,w[c][i+1][k_][e]);
							stack.push(c+1,j-1,e,l_,w[j-1][c+1][e][l_]);
							found = true;
						}

						//case 16 - all four stacked
						//16	1	1	1	1		2	2	0	0
						else if(en3 ==w[c][i+2][k_][e] + w[j-2][c+1][e][l_]
								+ 2*(nSeq1+nSeq2)*data->MBLFreeBase
								+ edangle5_pf(l,k,l-1,pf2,data)
								+ edangle3_pf(k,l,k+1,pf2,data)
								+ edangle5_pf(j,i,j-1,pf1,data)
								+ edangle3_pf(i,j,i+1,pf1,data)
								+ constantclosure)
						{
							found = true;
							stack.push(i+2,c,k_,e,w[c][i+2][k_][e]);
							stack.push(c+1,j-2,e,l_,w[j-2][c+1][e][l_]);
						}

						//case 6 - j and l stacked:
						//6		0	1	0	1		1	1	0	0
						else if(en3 == w[c][i+1][k_][e] + w[j-2][c+1][e][l_]
								+ (nSeq1+nSeq2)*data->MBLFreeBase
								+ edangle5_pf(j,i,j-1,pf1,data)
								+ edangle5_pf(l,k,l-1,pf2,data)
								+ constantclosure)
						{
							found = true;
							stack.push(i+1,c,k_,e,w[c][i+1][k_][e]);
							stack.push(c+1,j-2,e,l_,w[j-2][c+1][e][l_]);
						}

						//case 11 - i and k stacked
						//11	1	0	1	0		1	1	0	0
						else if(en3 == w[c][i+2][k_][e] + w[j-1][c+1][e][l_]
								+ (nSeq1+nSeq2)*data->MBLFreeBase
								+ edangle3_pf(k,l,k+1,pf2,data)
								+ edangle3_pf(i,j,i+1,pf1,data)
								+ constantclosure)
						{

							found = true;
							stack.push(i+2,c,k_,e,w[c][i+2][k_][e]);
							stack.push(c+1,j-1,e,l_,w[j-1][c+1][e][l_]);
						}

						if ((l_-1 >= 0)&&(!found)) {
							//case 2 - stack on l
						//2		0	0	0	1		0	1	1	0
							if(en3 == w[c][i+1][k_][e] + w[j-1][c+1][e][l_-1]
								+ nSeq2*data->MBLFreeBase
								+ edangle5_pf(l,k,l-1,pf2,data)
								+ nSeq1*gap
								+ constantclosure)
							{
								found = true;
								stack.push(i+1,c,k_,e,w[c][i+1][k_][e]);
								stack.push(c+1,j-1,e,l_-1,w[j-1][c+1][e][l_-1]);
							}

							//case 12 - i, k, and l stacked
						//12	1	0	1	1		1	2	1	0
							else if(en3 == w[c][i+2][k_][e] + w[j-1][c+1][e][l_-1]
										+ (nSeq1 + 2*nSeq2)*data->MBLFreeBase
										+ edangle5_pf(l,k,l-1,pf2,data)
										+ edangle3_pf(i,j,i+1,pf1,data)
										+ edangle3_pf(k,l,k+1,pf2,data)
										+ nSeq1*gap
										+ constantclosure)
							{
								found = true;
								stack.push(i+2,c,k_,e,w[c][i+2][k_][e]);
								stack.push(c+1,j-1,e,l_-1,w[j-1][c+1][e][l_-1]);
							}

							if ((k_-1 > 0)&&(!found)) {
								//case 10 - l and i stacked
						//10	1	0	0	1		1	1	1	1
								if(en3 == w[c][i+2][k_-1][e] + w[j-1][c+1][e][l_-1]
										+ (nSeq1+nSeq2)*data->MBLFreeBase
										+ edangle5_pf(l,k,l-1,pf2,data)
										+ edangle3_pf(i,j,i+1,pf1,data)
										+ (nSeq1+nSeq2)*gap
										+ constantclosure)
								{
									found = true;
									stack.push(i+2,c,k_-1,e,w[c][i+2][k_-1][e]);
									stack.push(c+1,j-1,e,l_-1,w[j-1][c+1][e][l_-1]);
								}
							}

							if ((k_+1 < maxrange)&&(!found)) {
								//case 4 - k and l stacked
						//4		0	0	1	1		0	2	2	0
								if(en3 == w[c][i+1][k_+1][e] + w[j-1][c+1][e][l_-1]
										+ 2*nSeq2*data->MBLFreeBase
										+ edangle5_pf(l,k,l-1,pf2,data)
										+ edangle3_pf(k,l,k+1,pf2,data)
										+ 2*nSeq1*gap
										+ constantclosure)
								{
									found = true;
									stack.push(i+1,c,k_+1,e,w[c][i+1][k_+1][e]);
									stack.push(c+1,j-1,e,l_-1,w[j-1][c+1][e][l_-1]);
								}
							}

						}		// end of if(l_-1>=0&&!found)

						if ((k_-1 > 0)&&(!found)) {
							//case 14 - i, j, and l stacked
						//14	1	1	0	1		2	1	0	1
							if(en3 == w[c][i+2][k_-1][e] + w[j-2][c+1][e][l_]
									+ (2*nSeq1+nSeq2)*data->MBLFreeBase
									+ edangle5_pf(l,k,l-1,pf2,data)
									+ edangle5_pf(j,i,j-1,pf1,data)
									+ edangle3_pf(i,j,i+1,pf1,data)
									+ nSeq2*gap
									+ constantclosure)
							{
								found = true;
								stack.push(i+2,c,k_-1,e,w[c][i+2][k_-1][e]);
								stack.push(c+1,j-2,e,l_,w[j-2][c+1][e][l_]);
							}

							//case 9 - i stacked
						//9		1	0	0	0		1	0	0	1
							else if(en3 == w[c][i+2][k_-1][e] + w[j-1][c+1][e][l_]
										+ nSeq1*data->MBLFreeBase
										+ edangle3_pf(i,j,i+1,pf1,data)
										+ nSeq2*gap
										+ constantclosure)
							{
								found = true;
								stack.push(i+2,c,k_-1,e,w[c][i+2][k_-1][e]);
								stack.push(c+1,j-1,e,l_,w[j-1][c+1][e][l_]);
							}

							else if (l_+1 < maxrange) {
								//case 13 - i and j stacked
						//13	1	1	0	0		2	0	0	2
								if(en3 == w[c][i+2][k_-1][e] + w[j-2][c+1][e][l_+1]
										+ 2*nSeq1*data->MBLFreeBase
										+ edangle5_pf(j,i,j-1,pf1,data)
										+ edangle3_pf(i,j,i+1,pf1,data)
										+ 2*nSeq2*gap
										+ constantclosure)
								{
									found = true;
									stack.push(i+2,c,k_-1,e,w[c][i+2][k_-1][e]);
									stack.push(c+1,j-2,e,l_+1,w[j-2][c+1][e][l_+1]);
								}
							}

						}		// end of if(k_-1>0&&!found)

						if ((k_+1 < maxrange)&&(!found)) {
							//case 8 - j, k, and l stacked
						//8		0	1	1	1		1	2	1	0
							if(en3 == w[c][i+1][k_+1][e] + w[j-2][c+1][e][l_]
									+ (nSeq1 + 2*nSeq2)*data->MBLFreeBase
									+ edangle5_pf(j,i,j-1,pf1,data)
									+ edangle5_pf(l,k,l-1,pf2,data)
									+ edangle3_pf(k,l,k+1,pf2,data)
									+ nSeq1*gap
									+ constantclosure)
							{
								found = true;
								stack.push(i+1,c,k_+1,e,w[c][i+1][k_+1][e]);
								stack.push(c+1,j-2,e,l_,w[j-2][c+1][e][l_]);
							}

							//case 3 - stack on k
						//3		0	0	1	0		0	1	1	0
							else if(en3 == w[c][i+1][k_+1][e] + w[j-1][c+1][e][l_]
										+ nSeq2*data->MBLFreeBase
										+ edangle3_pf(k,l,k+1,pf2,data)
										+ nSeq1*gap
										+ constantclosure)
							{
								found = true;
								stack.push(i+1,c,k_+1,e,w[c][i+1][k_+1][e]);
								stack.push(c+1,j-1,e,l_,w[j-1][c+1][e][l_]);
							}

							else if (l_+1 < maxrange) {
								//case 7 - j and k stacked
						//7		0	1	1	0		1	1	1	1
								if ( en3 == w[c][i+1][k_+1][e] + w[j-2][c+1][e][l_+1]
										+ (nSeq1+nSeq2)*data->MBLFreeBase
										+ edangle5_pf(j,i,j-1,pf1,data)
										+ edangle3_pf(k,l,k+1,pf2,data)
										+ (nSeq1+nSeq2)*gap
										+ constantclosure)
								{
									found = true;
									stack.push(i+1,c,k_+1,e,w[c][i+1][k_+1][e]);
									stack.push(c+1,j-2,e,l_+1,w[j-2][c+1][e][l_+1]);
								}
							}

						}		// end of if(k_+1<maxrange&&!found)

						if ((l_+1 < maxrange)&&(!found)) {
							//case 15 - i, j, and k stacked
						//15	1	1	1	0		2	1	0	1
							if(en3 == w[c][i+2][k_][e] + w[j-2][c+1][e][l_+1]
									+ (2*nSeq1+nSeq2)*data->MBLFreeBase
									+ edangle3_pf(k,l,k+1,pf2,data)
									+ edangle5_pf(j,i,j-1,pf1,data)
									+ edangle3_pf(i,j,i+1,pf1,data)
									+ nSeq2*gap
									+ constantclosure)
							{
								found = true;
								stack.push(i+2,c,k_,e,w[c][i+2][k_][e]);
								stack.push(c+1,j-2,e,l_+1,w[j-2][c+1][e][l_+1]);
							}

							//case 5 - j stacked
						//5		0	1	0	0		1	0	0	1
							else if(en3 == w[c][i+1][k_][e] + w[j-2][c+1][e][l_+1]
										+ nSeq1*data->MBLFreeBase
										+ edangle5_pf(j,i,j-1,pf1,data)
										+ nSeq2*gap
										+ constantclosure)
							{
								found = true;
								stack.push(i+1,c,k_,e,w[c][i+1][k_][e]);
								stack.push(c+1,j-2,e,l_+1,w[j-2][c+1][e][l_+1]);
							}

						}		// end of if(l_+1<maxrange&&!found)

					}		// end of for(d)
				}		// end of for(c)

			}		// end of if(en3!= hairpin energy)

			if (!found) {
				HDException ex (ERR, "Traceback error -- cannot trace value in V array!");
				throw ex;
			}

		}		// end of if(en3== v[i][j][k_][l_])

/*Reverse W cases:*/
		else {	// the subfragments i<->j and k<->l are within a multi-branched loop
			found = false;
			cout << i << " , "<<j<<" vs. "<<k<<" , "<<l<<"\n"<<flush;

/*Reverse W1 (adding on to a shorter loop fragment):*/

			//1 indicates nucleotide (i,j,k,or l) has been "added on"
			//0 indicates that it is part of the pre-calculated fragment:
			//all added-on bases are extra free bases within loop
			//bases added at a particular end on one sequence
			//but not the other increase gap penalty

			//		i	j	k	l		bases	gaps
			//				Sequence:	1	2	1	2
			//(no case 1 -- must extend by at least one base)
			//2		0	0	0	1		0	1	1	0
			//3		0	0	1	0		0	1	1	0
			//4		0	0	1	1		0	2	2	0
			//5		0	1	0	0		1	0	0	1
			//6		0	1	0	1		1	1	0	0
			//7		0	1	1	0		1	1	1	1
			//8		0	1	1	1		1	2	1	0
			//9		1	0	0	0		1	0	0	1
			//10	1	0	0	1		1	1	1	1
			//11	1	0	1	0		1	1	0	0
			//12	1	0	1	1		1	2	1	0
			//13	1	1	0	0		2	0	0	2
			//14	1	1	0	1		2	1	0	1
			//15	1	1	1	0		2	1	0	1
			//16	1	1	1	1		2	2	0	0

			//case 6
			//6		0	1	0	1		1	1	0	0
			if (en3 == w[j-1][i][k_][l_]
						+ (nSeq1+nSeq2)*data->MBLFreeBase){
				found = true;
				stack.push(i,j-1,k_,l_,w[j-1][i][k_][l_]);
			}
			//case 11
			//11	1	0	1	0		1	1	0	0
			else if (en3 == w[j][i+1][k_][l_]
							+ (nSeq1+nSeq2)*data->MBLFreeBase) {
				found = true;
				stack.push(i+1,j,k_,l_,w[j][i+1][k_][l_]);
			}
			//case 16
			//16	1	1	1	1		2	2	0	0
			else if (j-1>i+1) {
				if (en3 == w[j-1][i+1][k_][l_]
							+ 2*(nSeq1+nSeq2)*data->MBLFreeBase) {
					found = true;
					stack.push(i+1,j-1,k_,l_,w[j-1][i+1][k_][l_]);
				}
			}

			if ((k_>=1)&&(!found)) {
				//case 9
			//9		1	0	0	0		1	0	0	1
				if(en3 == w[j][i+1][k_-1][l_]
							+ nSeq1*data->MBLFreeBase
							+ nSeq2*gap){
					found = true;
					stack.push(i+1,j,k_-1,l_,w[j][i+1][k_-1][l_]);
				}

				//case 14
			//14	1	1	0	1		2	1	0	1
				else if(j-1 > i+1) {
					if(en3 == w[j-1][i+1][k_-1][l_]
								+ (2*nSeq1+nSeq2)*data->MBLFreeBase
								+ nSeq2*gap){
						found = true;
						stack.push(i+1,j-1,k_-1,l_,w[j-1][i+1][k_-1][l_]);
					}
				}

				if ((l_+1 < maxrange)&&(!found)&&(j-1 > i+1)) {
					//case 13
			//13	1	1	0	0		2	0	0	2
					if(en3 == w[j-1][i+1][k_-1][l_+1]
								+ 2*nSeq1*data->MBLFreeBase
								+ 2*nSeq2*gap){
						found = true;
						stack.push(i+1,j-1,k_-1,l_+1,w[j-1][i+1][k_-1][l_+1]);
					}
				}

				if ((l_>=1)&&(!found)) {
					//case 10
			//10	1	0	0	1		1	1	1	1
					if(en3 == w[j][i+1][k_-1][l_-1]
								+ (nSeq1+nSeq2)*data->MBLFreeBase
								+ (nSeq1+nSeq2)*gap){
						found = true;
						stack.push(i+1,j,k_-1,l_-1,w[j][i+1][k_-1][l_-1]);
					}
				}
			}		// end of if(k_>=1&&!found)

			if ((l_+1 < maxrange)&&(!found)) {
				//case 5
			//5		0	1	0	0		1	0	0	1
				if(en3 == w[j-1][i][k_][l_+1]
							+ nSeq1*data->MBLFreeBase
							+ nSeq2*gap){
					found = true;
					stack.push(i,j-1,k_,l_+1,w[j-1][i][k_][l_+1]);
				}
				else if ((k_+1 < maxrange)) {
					//case 7
			//7		0	1	1	0		1	1	1	1
					if(en3 == w[j-1][i][k_+1][l_+1]
								+ (nSeq1+nSeq2)*data->MBLFreeBase
								+ (nSeq1+nSeq2)*gap){
						found = true;
						stack.push(i,j-1,k_+1,l_+1,w[j-1][i][k_+1][l_+1]);
					}
				}

				//case 15
			//15	1	1	1	0		2	1	0	1
				else if(j-1 > i+1) {
					if(en3 == w[j-1][i+1][k_][l_+1]
								+ (2*nSeq1+nSeq2)*data->MBLFreeBase
								+ nSeq2*gap) {
						found = true;
						stack.push(i+1,j-1,k_,l_+1,w[j-1][i+1][k_][l_+1]);
					}
				}
			}		// end of if(l_+1<(maxrange)&&!found)

			if ((k_+1 < maxrange)&&(!found)) {
				//case 3
			//3		0	0	1	0		0	1	1	0
				if(en3 == w[j][i][k_+1][l_]
							+ nSeq2*data->MBLFreeBase
							+ nSeq1*gap){
					found = true;
					stack.push(i,j,k_+1,l_,w[j][i][k_+1][l_]);
				}

				//case 8
			//8		0	1	1	1		1	2	1	0
				else if (en3 == w[j-1][i][k_+1][l_]
								+ (nSeq1 + 2*nSeq2)*data->MBLFreeBase
								+ nSeq1*gap) {
					found = true;
					stack.push(i,j,k_+1,l_,w[j][i][k_+1][l_]);
				}

				else if (l_ >= 1) {
					//case 4
			//4		0	0	1	1		0	2	2	0
					if(en3 == w[j][i][k_+1][l_-1]
								+ 2*nSeq2*data->MBLFreeBase
								+ 2*nSeq1*gap){
						found = true;
						stack.push(i,j,k_+1,l_-1,w[j][i][k_+1][l_-1]);
					}
				}
			}		// end of if(k_+1<(maxrange)&&!found)

			if ((l_ >= 1)&&(!found)) {
				//case 2
			//2		0	0	0	1		0	1	1	0
				if(en3 == w[j][i][k_][l_-1]
							+ nSeq2*data->MBLFreeBase
							+ nSeq1*gap) {
					found = true;
					stack.push(i,j,k_,l_-1,w[j][i][k_][l_-1]);
				}
				//case 12
			//12	1	0	1	1		1	2	1	0
				else if(en3==w[j][i+1][k_][l_-1]
							+ (nSeq1 + 2*nSeq2)*data->MBLFreeBase
							+ nSeq1*gap){
					found = true;
					stack.push(i+1,j,k_,l_-1,w[j][i+1][k_][l_-1]);
				}
			}		// end of if(l_>=1&&!found)

/*Reverse W2: (helix termini):*/
			//16 cases, depending on whether or not the ends of the fragments (i,j,k,l)
			//are stacked on the end of the helix (as opposed to being part of the
			//final basepair of the helix)
				//							Bases	Gaps
				//						Seq	1	2	1	2
				//		i	j	k	l
				//1		0	0	0	0		0	0	0	0
				//2		0	0	0	1		0	1	1	0
				//3		0	0	1	0		0	1	1	0
				//4		0	0	1	1		0	2	2	0
				//5		0	1	0	0		1	0	0	1
				//6		0	1	0	1		1	1	0	0
				//7		0	1	1	0		1	1	1	1
				//8		0	1	1	1		1	2	1	0
				//9		1	0	0	0		1	0	0	1
				//10	1	0	0	1		1	1	1	1
				//11	1	0	1	0		1	1	0	0
				//12	1	0	1	1		1	2	1	0
				//13	1	1	0	0		2	0	0	2
				//14	1	1	0	1		2	1	0	1
				//15	1	1	1	0		2	1	0	1
				//16	1	1	1	1		2	2	0	0

			constantclosure = (nSeq1+nSeq2)*data->MBLHelix;
			//Energy of an additional helix within this multi-branch loop,
			//part of every sub-case within W2

			//case 1 - nothing stacked:
				//1		0	0	0	0		0	0	0	0
			if((!found) && en3 ==
						v[j][i][k_][l_]
						+ constantclosure
						+ penalty_pf(i,j,pf1,data)
						+ penalty_pf(k,l,pf2,data) )
			{
				found = true;
				stack.push(i,j,k_,l_,v[j][i][k_][l_]);
			}

			//case 6 - j and l stacked
				//6		0	1	0	1		1	1	0	0
			if((!found) && en3 ==
					v[j-1][i][k_][l_]
					+ (nSeq1+nSeq2)*data->MBLFreeBase
					+ edangle3_pf(j-1,i,j,pf1,data)
					+ edangle3_pf(l-1,k,l,pf2,data)
					+ constantclosure
					+ penalty_pf(i,j-1,pf1,data)
					+ penalty_pf(l-1,k,pf2,data) )
			{
				found = true;
				stack.push(i,j-1,k_,l_,v[j-1][i][k_][l_]);
			}

			//case 11 - i and k stacked
				//11	1	0	1	0		1	1	0	0
			if((!found) && en3 ==
					v[j][i+1][k_][l_]
					+ (nSeq1+nSeq2)*data->MBLFreeBase
					+ edangle5_pf(k+1,l,k,pf2,data)
					+ edangle5_pf(i+1,j,i,pf1,data)
					+ constantclosure
					+ penalty_pf(i+1,j,pf1,data)
					+ penalty_pf(k+1,l,pf2,data))
			{
				found = true;
				stack.push(i+1,j,k_,l_,v[j][i+1][k_][l_]);
			}

			//case 16 - i, j, k, and l stacked
				//16	1	1	1	1		2	2	0	0
			if((!found)&&(j-1>i+1)) {
				if (en3 == v[j-1][i+1][k_][l_]
						+ 2*(nSeq1+nSeq2)*data->MBLFreeBase
						+ edangle3_pf(l-1,k+1,l,pf2,data)
						+ edangle5_pf(k+1,l-1,k,pf2,data)
						+ edangle3_pf(j-1,i+1,j,pf1,data)
						+ edangle5_pf(i+1,j-1,i,pf1,data)
						+ constantclosure
						+ penalty_pf(i+1,j-1,pf1,data)
						+ penalty_pf(k+1,l-1,pf2,data))
				{
					found = true;
					stack.push(i+1,j-1,k_,l_,v[j-1][i+1][k_][l_]);
				}
			}

			if ((!found)&&(l_-1 >= 0)) {
				//case 2 - l stacked
				//2		0	0	0	1		0	1	1	0
				if(en3 == v[j][i][k_][l_-1]
						+ nSeq2*data->MBLFreeBase
						+ edangle3_pf(l-1,k,l,pf2,data)
						+ constantclosure
						+ nSeq1*gap
						+ penalty_pf(i,j,pf1,data)
						+ penalty_pf(l-1,k,pf2,data))
				{
					found = true;
					stack.push(i,j,k_,l_-1,v[j][i][k_][l_-1]);
				}

				else if (k_+1 < maxrange) {
					//case 4 - l and k stacked
				//4		0	0	1	1		0	2	2	0
					if(en3 == v[j][i][k_+1][l_-1]
							+ 2*nSeq2*data->MBLFreeBase
							+ edangle3_pf(l-1,k+1,l,pf2,data)
							+ edangle5_pf(k+1,l-1,k,pf2,data)
							+ constantclosure
							+ 2*nSeq1*gap
							+ penalty_pf(i,j,pf1,data)
							+ penalty_pf(l-1,k+1,pf2,data))
					{
						found = true;
						stack.push(i,j,k_+1,l_-1,v[j][i][k_+1][l_-1]);
					}
				}
			}		// end of if(!found&&l_-1>=0)

			if ((!found)&&(k_+1 < maxrange)) {
				//case 8 - j, k, and l stacked:
				//8		0	1	1	1		1	2	1	0
				if(en3 == v[j-1][i][k_+1][l_]
						+ (nSeq1 + 2*nSeq2)*data->MBLFreeBase
						+ edangle3_pf(j-1,i,j,pf1,data)
						+ edangle3_pf(l-1,k+1,l,pf2,data)
						+ edangle5_pf(k+1,l-1,k,pf2,data)
						+ constantclosure
						+ nSeq1*gap
						+ penalty_pf(i,j-1,pf1,data)
						+ penalty_pf(k+1,l-1,pf2,data))
				{
					found = true;
					stack.push(i,j-1,k_+1,l_,v[j-1][i][k_+1][l_]);
				}

				//case 3 - k stacked
				//3		0	0	1	0		0	1	1	0
				else if(en3 == v[j][i][k_+1][l_]
							+ nSeq2*data->MBLFreeBase
							+ edangle5_pf(k+1,l,k,pf2,data)
							+ constantclosure
							+ nSeq1*gap
							+ penalty_pf(i,j,pf1,data)
							+ penalty_pf(k+1,l,pf2,data))
				{
					found = true;
					stack.push(i,j,k_+1,l_,v[j][i][k_+1][l_]);
				}

				else if (l_+1 < maxrange) {
					//case 7 - j and k stacked
				//7		0	1	1	0		1	1	1	1
					if(en3 == v[j-1][i][k_+1][l_+1]
							+ (nSeq1+nSeq2)*data->MBLFreeBase
							+ edangle3_pf(j-1,i,j,pf1,data)
							+ edangle5_pf(k+1,l,k,pf2,data)
							+ constantclosure
							+ (nSeq1+nSeq2)*gap
							+ penalty_pf(i,j-1,pf1,data)
							+ penalty_pf(k+1,l,pf2,data))
					{
						found = true;
						stack.push(i,j-1,k_+1,l_+1,v[j-1][i][k_+1][l_+1]);
					}
				}
			}		// end of if(!found&&k_+1<maxrange)

			if ((l_+1 < maxrange)&&(!found)) {
				//case 5 - j stacked
				//5		0	1	0	0		1	0	0	1
				if(en3 == v[j-1][i][k_][l_+1]
							+ nSeq1*data->MBLFreeBase
							+ edangle3_pf(j-1,i,j,pf1,data)
							+ constantclosure
							+ nSeq2*gap
							+ penalty_pf(i,j-1,pf1,data)
							+ penalty_pf(k,l,pf2,data))
				{
					found = true;
					stack.push(i,j-1,k_,l_+1,v[j-1][i][k_][l_+1]);
				}

				//case 15 - i,j, and k stacked
				//15	1	1	1	0		2	1	0	1
				else if(en3 == v[j-1][i+1][k_][l_+1]
							+ (2*nSeq1+nSeq2)*data->MBLFreeBase
							+ edangle5_pf(k+1,l,k,pf2,data)
							+ edangle3_pf(j-1,i+1,j,pf1,data)
							+ edangle5_pf(i+1,j-1,i,pf1,data)
							+ constantclosure
							+ nSeq2*gap
							+ penalty_pf(i+1,j-1,pf1,data)
							+ penalty_pf(k+1,l,pf2,data))
				{
					found = true;
					stack.push(i+1,j-1,k_,l_+1,v[j-1][i+1][k_][l_+1]);
				}

				else if (k_>=1&&(j-1>i+1)) {
					//case 13 - i and j stacked
				//13	1	1	0	0		2	0	0	2
					if(en3 == v[j-1][i+1][k_-1][l_+1]
							+ 2*nSeq1*data->MBLFreeBase
							+ edangle3_pf(j-1,i+1,j,pf1,data)
							+ edangle5_pf(i+1,j-1,i,pf1,data)
							+ constantclosure
							+ 2*nSeq2*gap
							+ penalty_pf(i+1,j-1,pf1,data)
							+ penalty_pf(k,l,pf2,data))
					{
						found = true;
						stack.push(i+1,j-1,k_-1,l_+1,v[j-1][i+1][k_-1][l_+1]);
					}
				}
			}		// end of if(l_+1<maxrange&&!found)

			if ((!found)&&(k_>=1)) {
				//case 9 - i alone is stacked
				//9		1	0	0	0		1	0	0	1
				if(en3 == v[j][i+1][k_-1][l_]
							+ nSeq1*data->MBLFreeBase
							+ edangle5_pf(i+1,j,i,pf1,data)
							+ constantclosure
							+ nSeq2*gap
							+ penalty_pf(i+1,j,pf1,data)
							+ penalty_pf(k,l,pf2,data))
				{
					found = true;
					stack.push(i+1,j,k_-1,l_,v[j][i+1][k_-1][l_]);
				}

				//case 14 - i, j, and l stacked
				//14	1	1	0	1		2	1	0	1
				else if(j-1 > i+1) {
					if(en3 == v[j-1][i+1][k_-1][l_]
							+ (2*nSeq1+nSeq2)*data->MBLFreeBase
							+ edangle3_pf(l-1,k,l,pf2,data)
							+ edangle3_pf(j-1,i+1,j,pf1,data)
							+ edangle5_pf(i+1,j-1,i,pf1,data)
							+ constantclosure
							+ nSeq2*gap
							+ penalty_pf(i+1,j-1,pf1,data)
							+ penalty_pf(k,l-1,pf2,data))
					{
						found = true;
						stack.push(i+1,j-1,k_-1,l_,v[j-1][i+1][k_-1][l_]);
					}
				}
			}		// end of if(!found&&k_>=1)

			if ((!found)&&(l_-1 >= 0)) {
				//case 12 - i, k, and l stacked:
				//12	1	0	1	1		1	2	1	0
				if(en3 == v[j][i+1][k_][l_-1]
							+ (nSeq1 + 2*nSeq2)*data->MBLFreeBase
							+ edangle3_pf(l-1,k+1,l,pf2,data)
							+ edangle5_pf(i+1,j,i,pf1,data)
							+ edangle5_pf(k+1,l-1,k,pf2,data)
							+ constantclosure
							+ nSeq1*gap
							+ penalty_pf(i+1,j,pf1,data)
							+ penalty_pf(l-1,k+1,pf2,data))
				{
					found = true;
					stack.push(i+1,j,k_,l_-1,v[j][i+1][k_][l_-1]);
				}
				else if (k_ >= 1) {
					//case 10 - l and i stacked
				//10	1	0	0	1		1	1	1	1
					if(en3 == v[j][i+1][k_-1][l_-1]
							+ (nSeq1+nSeq2)*data->MBLFreeBase
							+ edangle3_pf(l-1,k,l,pf2,data)
							+ edangle5_pf(i+1,j,i,pf1,data)
							+ constantclosure
							+ (nSeq1+nSeq2)*gap
							+ penalty_pf(i+1,j,pf1,data)
							+ penalty_pf(k,l-1,pf2,data))
					{
						found = true;
						stack.push(i+1,j,k_-1,l_-1,v[j][i+1][k_-1][l_-1]);
					}
				}
			}		// end of if(!found&&l_-1>=0)

/*Reverse W3 (two merged fragments)*/
			//calculate the free energy of 2 fragments merged:
			for (c = i+MINLOOP; (c < j-MINLOOP)&&(!found); c++) {
				for (d = max(k+MINLOOP,c-maxsep);
						(d < l-MINLOOP)&&(d <= c+maxsep)&&(!found); d++) {
					e = d-c+maxsep;

					if(en3 ==w[c][i][k_][e]+w[j][c+1][e][l_]) {
						found = true;
						stack.push(i,c,k_,e,w[c][i][k_][e]);
						stack.push(c+1,j,e,l_,w[j][c+1][e][l_]);
					}
				}
			}

			if (!found) {
				HDException ex (ERR, "Traceback error -- cannot trace value in W array!");
				throw ex;
			}

		}		// end of else

	}		// end of while(stack.pull...)

	cout << "Finished calculations.\n\n";

/****Free memory****/
	for (i=0; i<=nBases1; i++) {
		for (j=0; j<i; j++) {
			for (k=0; k<maxrange; k++) {
				ASSERT(v[i][j][k] != NULL);
				delete[] v[i][j][k];
				ASSERT(w[i][j][k] != NULL);
				delete[] w[i][j][k];
			}
			ASSERT(v[i][j] != NULL);
			delete[] v[i][j];
			ASSERT(w[i][j] != NULL);
			delete[] w[i][j];
		}
		ASSERT(v[i] != NULL);
		delete[] v[i];
		ASSERT(w[i] != NULL);
		delete[] w[i];
	}
	ASSERT(v != NULL);
	delete[] v;
	ASSERT(w != NULL);
	delete[] w;

	for (i=0;i<=nBases1;i++) {
		ASSERT(mine[i] != NULL);
		delete[] mine[i];
	}
	ASSERT(mine != NULL);
	delete[] mine;

	ASSERT(pair[0] != NULL);
	delete[] pair[0];
	ASSERT(pair[0] != NULL);
	delete[] pair[1];
	ASSERT(pair != NULL);
	delete[] pair;

}	// end of pdynalign()


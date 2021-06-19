/*                               -*- Mode: C -*- 
 * efunctions.cpp ---
 * Author          : Amelia Bellamy-Royds
 * Created         : Summer 2006
 * Last Modified By: Marcel Turcotte
 * Last Modified On: Mon Dec 11 11:07:34 2006
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
This file contains various energy-calculation functions which are called by
the main pdynalign function.  The functions have two versions: one for a single
sequnce and one for a profile, which sums the scores for all sequences in
the profile.

einternal()	and einternal_pf();
ehairpin()	and ehairping_pf();
edangle5()	and edangle5_pf();
edangle3()	and edangle3_pf();
ebp()		and ebp_pf();
penalty()	and penalty_pf();

Also contains
get_inc() //determines whether a particular basepair is possible (for a profile)

*/
#include "pdynalign.h"
#include <math.h>
#include <string.h>


/*--------------------------------------------------------------------------------------
Function:	einternal_pf()

Description:	Calculate the energy of bulge/internal loops for a profile

Return value:	Sum of the energy
--------------------------------------------------------------------------------------*/
energy_t einternal_pf(int i,int j,int ip,int jp,profile *pf, datatable *data)
{
	int 			k;
	energy_t	sum;	// sum energy for sequences

	ASSERT( pf != NULL && data != NULL );
	sum = 0;
	for( k=0; k < (pf->nSeq); k++) {
		ASSERT( (pf->ct+k) != NULL );
		sum = sum + einternal(i, j, ip, jp, (pf->ct+k), data);
		if (sum >= INFINITE) return sum;
	}

	return sum;
}


/*--------------------------------------------------------------------------------------
Function:	einternal()

Description:	Calculate the energy of a bulge/internal loop for a sequence
		where i is paired to j; ip is paired to jp; ip > i; j > jp
--------------------------------------------------------------------------------------*/
energy_t einternal(int i,int j,int ip,int jp,structure *ct, datatable *data)
{
int energy,size,size1,size2,loginc, lopsid;//tlink,count,key,e[4]
	/* size,size1,size2 = size of a loop
		energy = energy calculated
		loginc = the value of a log used in large hairpin loops
	*/
		size1 = ip-i-1;
		size2 = j - jp - 1;

    //a typical internal or bulge loop:
		if (size1==0||size2==0) {//bulge loop
			size = size1+size2;
			if (size==1) {
				energy = data->stack[ct->numseq[i]][ct->numseq[j]]
						[ct->numseq[ip]][ct->numseq[jp]]
						+ data->bulge[size] + data->eparam[2];
			}
			else if (size>30) {
				loginc = int((data->prelog)*log(double ((size)/30.0)));
				energy = data->bulge[30] + loginc + data->eparam[2];
				energy = energy + penalty(i,j,ct,data) + penalty(jp,ip,ct,data);
			}
			else {
				energy = data->bulge[size] + data->eparam[2];
				energy = energy + penalty(i,j,ct,data) + penalty(jp,ip,ct,data);
      }
		}		// end of if
		else {										//internal loop
			size = size1 + size2;
			lopsid = abs(size1-size2);

			if (size>30) {
				loginc = int((data->prelog)*log((double ((size))/30.0)));
				if ((size1==1||size2==1)&&data->gail) {
					energy = data->tstki[ct->numseq[i]][ct->numseq[j]][1][1] +
					data->tstki[ct->numseq[jp]][ct->numseq[ip]][1][1] +
					data->inter[30] + loginc + data->eparam[3] +
					min(data->maxpen,(lopsid*
					data->poppen[min(2,min(size1,size2))]));
				}
				else {
					energy=data->tstki[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]]
					+	data->tstki[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]]
					+	data->inter[30] + loginc + data->eparam[3] +
					min(data->maxpen,(lopsid * data->poppen[min(2,min(size1,size2))]));
				}
			}
			else if ((size1==2)&&(size2==2))					//2x2 internal loop
				energy = data->iloop22[ct->numseq[i]][ct->numseq[ip]][ct->numseq[j]][ct->numseq[jp]]
						[ct->numseq[i+1]][ct->numseq[i+2]][ct->numseq[j-1]][ct->numseq[j-2]];

			else if ((size1==1)&&(size2==2)) {				//2x1 internal loop
				energy = data->iloop21[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]]
						[ct->numseq[j-1]][ct->numseq[jp+1]][ct->numseq[ip]][ct->numseq[jp]];
			}

			else if ((size1==2)&&(size2==1)) {				//1x2 internal loop
				energy = data->iloop21[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]]
						[ct->numseq[ip-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]];
			}

			else if (size==2) //a single mismatch
				energy = data->iloop11[ct->numseq[i]][ct->numseq[i+1]][ct->numseq[ip]]
						[ct->numseq[j]][ct->numseq[j-1]][ct->numseq[jp]];

			else if ((size1==1||size2==1)&&data->gail) { //this loop is lopsided
			//note, we treat this case as if we had a loop composed of all As
			//if and only if the gail rule is set to 1 in miscloop.dat
				energy = data->tstki[ct->numseq[i]][ct->numseq[j]][1][1] +
				data->tstki[ct->numseq[jp]][ct->numseq[ip]][1][1] +
				data->inter[size] + data->eparam[3] +
				min(data->maxpen,(lopsid *	data->poppen[min(2,min(size1,size2))]));
			}

			else {
				energy = data->tstki[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]]
				+	data->tstki[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]]
				+	data->inter[size] + data->eparam[3]
				+	min(data->maxpen,(lopsid*data->poppen[min(2,min(size1,size2))]));
			}

		}		// end of else

		return energy;
}		// end of einternal()


/*--------------------------------------------------------------------------------------
Function:	ehairpin_pf()

Description:	Calculate the energy of hairpins for a profile

Return value:	Sum of the energy
--------------------------------------------------------------------------------------*/
inline energy_t ehairpin_pf(short i,short j,profile *pf, datatable *data)
{
	int 			k;
	energy_t	sum;	// sum energy for sequences

	ASSERT( pf != NULL && data != NULL );
	sum = 0;
	for( k=0; k<pf->nSeq; k++) {
		ASSERT( (pf->ct+k) != NULL );
		sum = sum + ehairpin(i,j,pf->ct+k,data);
		if (sum >= INFINITE) return sum;
	}

	return sum;
}


/*--------------------------------------------------------------------------------------
Function:	ehairpin()

Description:	Calculate the energy of hairpins for a sequence

Return value:	Energy of the sequence
--------------------------------------------------------------------------------------*/
inline energy_t ehairpin(short i, short j, structure *ct, datatable *data) {
	int energy,size,loginc,tlink,count,key,k;
	size = j-i-1;

		if (size>30) {

			loginc = int((data->prelog)*log((double ((size))/30.0)));

			energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]]
				+ data->hairpin[30]+loginc+data->eparam[4];
		}
		else if (size<3) {
     		energy = data->hairpin[size] + data->eparam[4];
			if (ct->numseq[i]==4||ct->numseq[j]==4) energy = energy+6;
		}
		else if (size==4) {
			tlink = 0;
			key = (ct->numseq[j])*3125 + (ct->numseq[i+4])*625 +
					(ct->numseq[i+3])*125 + (ct->numseq[i+2])*25+(ct->numseq[i+1])*5+(ct->numseq[i]);
			for (count=1;count<=data->numoftloops&&tlink==0;count++) {
				if (key==data->tloop[count][0]) tlink = data->tloop[count][1];
			}
			energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
					[ct->numseq[i+1]][ct->numseq[j-1]]
					+ data->hairpin[size] + data->eparam[4] + tlink;
		}
		else if (size==3) {
			tlink = 0;
			key = (ct->numseq[j])*625 + (ct->numseq[i+3])*125 + (ct->numseq[i+2])*25
					+(ct->numseq[i+1])*5+(ct->numseq[i]);
			for (count=1;count<=data->numoftriloops&&tlink==0;count++) {
				if (key==data->triloop[count][0]) tlink = data->triloop[count][1];
			}

			energy =	data->hairpin[size] + data->eparam[4] + tlink
         	+penalty(i,j,ct,data);
		}

		else {
			energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
					[ct->numseq[i+1]][ct->numseq[j-1]]
					+ data->hairpin[size] + data->eparam[4];
		}

		//check for GU closeure preceded by GG
			if (ct->numseq[i]==3&&ct->numseq[j]==4) {
				if ((i>2&&i<ct->numofbases)||(i>ct->numofbases+2)) {
					if (ct->numseq[i-1]==3&&ct->numseq[i-2]==3) {
						energy = energy + data->gubonus;
					}
				}
			}

      //check for a poly-c loop
			tlink = 1;
			for (k=1;(k<=size)&&(tlink==1);k++) {
				if (ct->numseq[i+k] != 2) tlink = 0;
			}
			if (tlink==1) {  //this is a poly c loop so penalize
				if (size==3) energy = energy + data->c3;
				else energy = energy + data->cint + size*data->cslope;
			}

			return energy;
}


/*--------------------------------------------------------------------------------------
Function:	edangle5_pf()

Description:	Calculate the energy of dangle5's for a profile

Return value:	Sum of the energy
--------------------------------------------------------------------------------------*/
inline energy_t edangle5_pf(int i,int j,int ip,profile *pf, datatable *data)
{
	int 			k;
	energy_t	sum;	// sum energy for sequences

	ASSERT( pf != NULL && data != NULL );
	sum = 0;
	for( k=0; k<pf->nSeq; k++) {
		ASSERT( (pf->ct+k) != NULL );
		sum = sum + edangle5(i,j,ip,pf->ct+k,data);
		if (sum >= INFINITE) return sum;
	}

	return sum;
}


/*--------------------------------------------------------------------------------------
Function:	edangle5()

Description:	Calculate the energy of dangle5's for a sequence

Return value:	Energy of the sequence
--------------------------------------------------------------------------------------*/
inline energy_t edangle5(int i,int j,int ip,structure* ct,datatable* data)
{
	return data->dangle[ct->numseq[j]][ct->numseq[i]][ct->numseq[ip]][2];
}

/*--------------------------------------------------------------------------------------
Function:	edangle3_pf()

Description:	Calculate the energy of dangle3's for a profile

Return value:	Sum of the energy
--------------------------------------------------------------------------------------*/
inline energy_t edangle3_pf(int i,int j,int ip,profile *pf, datatable *data)
{
	int 			k;
	energy_t	sum;	// sum energy for sequences

	ASSERT( pf != NULL && data != NULL );
	sum = 0;
	for( k=0; k<pf->nSeq; k++) {
		ASSERT( (pf->ct+k) != NULL );
		sum = sum + edangle3(i,j,ip,pf->ct+k,data);
		if (sum >= INFINITE) return sum;
	}

	return sum;
}


/*--------------------------------------------------------------------------------------
Function:	edangle3()

Description:	Calculate the energy of dangle3's for a sequence

Return value:	Energy of the sequence
--------------------------------------------------------------------------------------*/
inline energy_t edangle3(int i,int j,int ip,structure* ct,datatable* data)
{
	return data->dangle[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][1];
}


/*--------------------------------------------------------------------------------------
Function:	ebp_pf()

Description:	Calculate the energy of bp's for a profile

Return value:	Sum of the energy
--------------------------------------------------------------------------------------*/
inline energy_t ebp_pf(int i,int j,int ip,int jp,profile *pf, datatable *data)
{
	int 			k;
	energy_t	sum;	// sum energy for sequences

	ASSERT( pf != NULL && data != NULL );
	sum = 0;
	for( k=0; k<pf->nSeq; k++) {
		ASSERT( (pf->ct+k) != NULL );
		sum = sum + ebp(i,j,ip,jp,pf->ct+k,data);
		if (sum >= INFINITE) return sum;
	}

	return sum;
}


/*--------------------------------------------------------------------------------------
Function:	ebp()

Description:	Calculate the energy of bp's for a sequence

Return value:	Energy of the sequence
--------------------------------------------------------------------------------------*/
inline energy_t ebp(int i,int j,int ip,int jp,structure *ct, datatable *data)
{
	return data->stack[(ct->numseq[i])][(ct->numseq[j])]
				[(ct->numseq[ip])][(ct->numseq[jp])];

}


/*--------------------------------------------------------------------------------
Function:	penalty_pf()

Description: 	Calculates the end penalty for a profile

Return Value:	Sum of the penalty
--------------------------------------------------------------------------------*/
inline int penalty_pf(int i,int j,profile* pf, datatable *data) {
	int k,sum;

	ASSERT( pf != NULL && data != NULL );
	sum = 0;
	for( k=0; k < pf->nSeq; k++) {
		ASSERT( (pf->ct+k) != NULL );
		sum += penalty(i,j,pf->ct+k,data);
		if (sum >= INFINITE) return sum;
	}

	return sum;
}


/*--------------------------------------------------------------------------------
Function:	penalty()

Description: 	Calculates whether a terminal pair i,j requires the end penalty

Return Value: Penalty
--------------------------------------------------------------------------------*/
inline int penalty(int i,int j,structure* ct, datatable *data) {
	if (ct->numseq[i]==4||ct->numseq[j]==4)
		return data->auend;
	else return 0;//no end penalty
}



/*-----------------------------------------------------------------------------
Function:			get_inc()

Description:	Check if one element can paired with another in a profile

Return value:	If yes, return 1, otherwise return 0
-----------------------------------------------------------------------------*/
short int get_inc(int i, int j, profile *pf)
{
	int 			k;
	short int	e;
	short int inc[6][6] =
			{
				{0,0,0,0,0,0},		//	\	X	A	C	G	U	I
				{0,0,0,0,1,0},		//	X	0	0	0	0	0	0
				{0,0,0,1,0,0},		//	A	0	0	0	0	1	0
				{0,0,1,0,1,0},		//	C	0	0	0	1	0	0
				{0,1,0,1,0,0},		//	G	0	0	1	0	1	0
				{0,0,0,0,0,0}		//	U	0	1	0	1	0	0
			};						//	I	0	0	0	0	0	0

	ASSERT( pf != NULL );
	e = 1;

	// go through all sequences; e = 1 only when all i,j's in each sequences can
	// be paired with each other
	for( k=0; k<pf->nSeq; k++) {
		ASSERT( (pf->ct+k) != NULL );
		e = e & inc[(pf->ct+k)->numseq[i]][(pf->ct+k)->numseq[j]];
	}

	ASSERT( (e == 1) || (e == 0) );

	return e;
}


/*                               -*- Mode: C -*- 
 * datable.cpp ---
 * Author          : Amelia Bellamy-Royds
 * Created         : Summer 2006
 * Last Modified By: Marcel Turcotte
 * Last Modified On: Mon Dec 11 11:10:56 2006
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

#include <stdio.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <stdlib.h>

#include "defines.h"

struct datatable //this structure contains all the info read from thermodynamic
							//data files
{
	float prelog;
		//read from the miscloop file
		//"Extrapolation for large loops based on polymer theory
		// internal, bulge or hairpin loops > 30:
		// dS(T)=dS(30)+param*ln(n/30)"

 	energy_t

		//read from the miscloop file:
		maxpen,
			// maximum penalty for the Ninio correction
		poppen [5],
			//array of values for the Ninio f(m) correction
		eparam[11],
			//Various energetic parameters:
			//eparam[1...4] hard-coded to 0
			//eparam[5] multi-branch loop offset
			//eparam[6] multi-branch loop free base penalty
			//eparam[7] hard-coded to 30
			//eparam[8] hard-coded to 30
			//eparam[9] hard-coded to -500
			//eparam[10] multi-branch loop helix penalty
		efn2a,
			// efn2 multi-branch loop offset
		efn2b,
			// efn2 multi-branch loop free base penalty
		efn2c,
			// efn2 multi-branch loop helix penalty
		auend,
			// penalty for a terminal au (or gu???) pair on a helix
		gubonus,
			// bonus for GGG hairpin
		cslope,
			// c-hairpin slope
		cint,
			// c-hairpin intercept
		c3,
			// c-hairpin of 3
		init,
			// intermolecular initiation free energy
		gail,
			// GAIL Rule (Grossly Asymmetric Interior Loop Rule)
			// (on/off <-> 1/0)


		//read from the dangle file:
		dangle[6][6][6][3],
			//energies for a single strand dangling off
			//the end of a helix.

			//for dangle[i][j][k][l]:
			//if l = 1 dangle nucleotide is on 3' strand of helix
			//else if l=2 dangle nucleotide is on 5' strand of helix
			//i indicates nucleotide on 3' strand of helix
			//j indicates nucleotide on 5' strand of helix
			//k indicates the dangle nucleotide
			/*
				e.g. dangle [1][4][3][1] == dangle [a][u][g][1]:

				5' --> 3' AG
				3' <-- 5' U
			*/

		//read from the loop file:
		inter[31],
			//destabilizing energy of internal loop, based on size of loop
			//(from 1 to 30 nucleotides,  <=3 infinite)
		bulge[31],
			//destabilizing energy of bulge loop, based on size
			//(from 1 to 30 nucleotides)
		hairpin[31],
			//destabilizing energy of hairpin loop, based on size
			//(from 1 to 30 nucleotides, <=2 infinite)

		//read from the stack file:
		stack[6][6][6][6],
			//energies for one base pair stacked against another
			//in a helix.

			//for stack[i][j][k][l]:
			//i indicates nucleotide on 5' end of first strand
			//j indicates nucleotide on 3' end of second strand
			//k indicates stacked nucleotide on first strand
			//l indicates stacked nucleotide on second strand
			/*
				e.g. stack [1][4][3][2] == stack [a][u][g][c]:

				5' --> 3' AG
				3' <-- 5' UC
			*/

		//read from the tstackh file:
		tstkh[6][6][6][6],
			//energies for a pair of non-bonded nuclei stacked against a
			//basepair at the termini of a helix.

			//for stack[i][j][k][l]:
			//i indicates nucleotide on 5' end of first strand (helix)
			//j indicates nucleotide on 3' end of second strand (helix)
			//k indicates stacked nucleotide on first strand
			//l indicates stacked nucleotide on second strand
			/*
				e.g. tstkh [1][4][3][1] == tstkh [a][u][g][a]:

				5' --> 3' AG
				3' <-- 5' UA
			*/

		//read from the tstacki file:
		tstki[6][6][6][6],
			//energies for a pair of non-bonded nuclei stacked against a
			//basepair at the termini of a helix.
			//??? not sure what the difference is between tstackh and tstacki
			//but the energy values are more positive in tstacki ???

			//for stack[i][j][k][l]:
			//i indicates nucleotide on 5' end of first strand (helix)
			//j indicates nucleotide on 3' end of second strand (helix)
			//k indicates stacked nucleotide on first strand
			//l indicates stacked nucleotide on second strand
			/*
				e.g. tstki [1][4][3][1] == tstki [a][u][g][a]:

				5' --> 3' AG
				3' <-- 5' UA
			*/

		//read from the tloop file:
		tloop[MAXTLOOP+1][2],
			//stabilizing energies of particular 6-nucleotide loop sequences:
			//first dimension is the coded sequence
			//code: sum for n=1 to 6 of
			//	(numeric value of nucleotide n in sequence)*5^(n-1)
			//second dimension is the value
		numoftloops,
			//counter for the max. index of a valid value in the tloop array

		//read from the int22 file:
		iloop22[6][6][6][6][6][6][6][6],
			//energy of 2+2 symmetric internal loops:
			//iloop22[a][b][c][d][j][l][k][m] =
			//5' --> 3'
			//a j l b
			//|     |
			//c k m d
			//3' <-- 5'
			//only values for valid basepairs stored

		//read from the int21 file:
		iloop21[6][6][6][6][6][6][6],
			//energy of 2+1 assymetric internal loops:
			//iloop22[a][b][c][d][e][f][g] =
			//5'-->3'
			//a  c  f
			//|     |
			//b d e g
			//3'<--5'
			//only values for valid basepairs stored

		//read from the triloop file (NOTE file is currently empty!)
		triloop[MAXTLOOP+1][2],
			//same as for tloop, except only reads a 5-nucleotide sequence
		numoftriloops,
			//counter for the max index of valid values in triloop array

		//read from the coaxial file
		coax[6][6][6][6],
			//"this is the array that keeps track of
			//flush coaxial stacking (no intervening nucs)
			//data arrangement of coax: data->coax[a][b][c][d]
			//5'b-c3'
			//3'a d5'
			//this means the helix backbone is continuous between nucs b and c"

		//read from the tstackcoax file
		tstackcoax[6][6][6][6],
			//"data arrangement of tstackcoax:
			//5'a-c -> strand continues into stack
			//3'b-d -> strand does not continue to stack
			//pair between a-b, c-d is a mismatch"

		//read from the coaxstack file
		coaxstack[6][6][6][6],
			//"data arrangement of coaxstack:
			//5'a-c ->strand contnues into stack
			//3'b d ->strand does not continue to stack
			//pair between a-b, mismatch between c-d
			//backbone is discontinuous between b and d"

		//read from the tstack file
		tstack[6][6][6][6],
			//"this is the terminal mismatch data used in intermolecular folding
			//add to the tstack table the case where X (represented as 0) is looked up.
			//also add the case where 5 (the intermolecular linker) is looked up,
			//this is actually a dangling end, not a terminal mismatch."

		//read from the tstackm file
		tstkm[6][6][6][6],
			//terminal mismatches and pairs (difference with tstackh,
			//tstacki not explained, but tstackm has finite values if
			//either base pair is valid
			//tstkm[a][b][c][d]:
			//5'ac
			//3'bd

		//read from the int11 file:
		iloop11[6][6][6][6][6][6];
			//energy of 1+1 symmetric internal loops:
			//iloop22[a][b][c][d][e][f] =
			//a b c
			//|   |
			//d e f
			//only values for valid basepairs stored


		//NOTE:  for all tables indexed by nucleotide identity, the code is
		// 0=X, 1=A, 2=C, 3=G, 4=U, 5=I
		//NOTE(2):  all energy values are stored as 10x the energy in
		//	kcal/mol, so that precision to 1/10th of a kcal can be
		//	maintained as an integer

	datatable();

};


/* Constructor of datatable: initializes values of multi-dimensional arrays*/
datatable::datatable()
{
	int a,b,c,d,e,f,g,h;

	for (a=0;a<=5;a++) {
		for (b=0;b<=5;b++) {
			for (c=0;c<=5;c++) {
				for (d=0;d<=5;d++) {
					for (e=0;e<=5;e++) {
						for (f=0;f<=5;f++) {
							iloop11[a][b][c][d][e][f] = 0;

							for (g=0;g<=5;g++) {
								iloop21[a][b][c][d][e][f][g] = INFINITE;

								for (h=0;h<=5;h++) {
									iloop22[a][b][c][d][e][f][g][h] = INFINITE;
								}
							}
						}
					}
				}
			}
		}
	}

};



/*--------------------------------------------------------------------------------
Function:	penalty2()

Description: Calculates whether a terminal pair i,j requires the end penalty
where i,j represent numerical codes for the bases
--------------------------------------------------------------------------------*/
int penalty2(int i,int j, datatable *data) {
	if (i==4||j==4)
		return data->auend;
	else return 0;//no end penalty
}


// Converts a base to a numeric
short int tonumi(char *base);

/*--------------------------------------------------------------------------------
Function:	tonumi()

Description: 	Simply Converts base chars to integer numbers.

Return Value:	Converted interger
--------------------------------------------------------------------------------*/
short int tonumi(char *base)	{
	short int	a;
	if (!strcmp(base,"A")||!strcmp(base,"B")) (a = 1);
	else if(!strcmp(base,"C")||!strcmp(base,"Z")) (a = 2);
	else if(!strcmp(base,"G")||!strcmp(base,"H")) (a = 3);
	else if(!strcmp(base,"U")||!strcmp(base,"V")) (a = 4);
	else if(!strcmp(base,"T")||!strcmp(base,"W")) (a = 4);
	else (a=0);  //this is for others, like X

	return a;
}



// Gets thermodynamic data from data files
int opendat (char *loop2,char *stackf,char *tstackh,char *tstacki,
		char *tloop,char *miscloop, char *danglef, char *int22, char *int21,
      char *coax,char *tstackcoax,char *coaxstack,
      char *tstack,char *tstackm, char *triloop, char *int11, datatable* data);


/*****************************************************************************
 Function: opendat()

 Description:	Function opens data files to read thermodynamic data
*****************************************************************************/
int opendat (char *loop2,char *stackf,char *tstackh,char *tstacki,
		char *tloop,char *miscloop, char *danglef, char *int22, char *int21,
		char *coax,char *tstackcoax,char *coaxstack,char *tstack,char *tstackm,
		char *triloop, char *int11, datatable* data)
{
	char lineoftext[100], base[110];
	int count,i,j,k,l, m, a, b, c, d,e,f,g;
	int	numFiles = 16;		// number of data files
	char *fileName[16/*numfiles*/]={loop2, stackf, tstackh, tstacki, tloop, miscloop, danglef,
			int22, int21, coax, tstackcoax, coaxstack, tstack, tstackm,	triloop, int11};
	float temp;
	FILE *check;

	ifstream ml1;
	ifstream lo1;
	ifstream st1;
	ifstream th1;
	ifstream ti1;
	ifstream tl1;
	ifstream da1;
	ifstream in1;
	ifstream in2;
	ifstream tri;
	ifstream co1;
	ifstream co2;
	ifstream co3;
	ifstream st2;
	ifstream tsm;
	ifstream i11;

	for( i=0; i<numFiles; i++ ) {
		if( (check = fopen(fileName[i], "r")) == NULL ) {
			fclose(check);
			char *message = new char[100];
			strcpy(message, "Cannot find data file ");
			strcat(message, fileName[i]);
			throw new HDException (IO_ERR, message );

       		}
	}

	//open data files for reading
	ml1.open(miscloop);
	lo1.open(loop2);
	st1.open(stackf);
	th1.open(tstackh);
	ti1.open(tstacki);
	tl1.open(tloop);
	da1.open(danglef);
	in1.open(int22);
	in2.open(int21);
	tri.open(triloop);
	co1.open(coax);
	co2.open(tstackcoax);
	co3.open(coaxstack);
	st2.open(tstack);
	tsm.open(tstackm);
	i11.open(int11);

/* Read information from miscloop */

	// the key sequence "-->" now indicates a record
	ml1>>lineoftext;
	while(strcmp(lineoftext,"-->")) {
		ml1>>lineoftext;
	}

	ml1 >> (data->prelog); //extrapolation for large data loops
	data->prelog = (data->prelog)*10;

	ml1>>lineoftext;
	while(strcmp(lineoftext,"-->")) {
		ml1>>lineoftext;
	}

	ml1 >> temp; //ninio equation maximum correction
	data->maxpen = int (temp*10 + .5);

	ml1>>lineoftext;
	while(strcmp(lineoftext,"-->")) {
		ml1>>lineoftext;
	}

	/* reads float values, converts to int and assigns into array 'poppen' */
	for (count=1;count<= 4;count ++) { //ninio f(m) array
		ml1 >> temp;
		(data->poppen[count])= (int) (temp*10 + .5);
	}

	ml1>>lineoftext;
	while(strcmp(lineoftext,"-->")) {
		ml1>>lineoftext;
	}

	data->eparam[1] = 0;						 // assign some variables that are
	data->eparam[2] = 0;						 //	"hard-wired" into code
	data->eparam[3] = 0;
	data->eparam[4] = 0;

	ml1 >> temp;
	data->eparam[5] = short (floor (temp*10+.5));  //constant multi-loop penalty (offset)

	ml1 >> temp;
	data->eparam[6] = short (floor (temp*10+.5));  //multi-loop free-base penalty

	data->eparam[7] = 30;
	data->eparam[8] = 30;
	data->eparam[9] = -500;

	ml1 >> temp;
	data->eparam[10] = short (floor (temp*10+.5)); //multi-loop helix penalty

	ml1>>lineoftext;
	while(strcmp(lineoftext,"-->")) {
		ml1>>lineoftext;
	}

	ml1 >> temp;								//efn2 multi-branch loop offset

	if (ml1.peek()==EOF) {
		//these are old energy rules -- treat the other constants properly
		data->efn2a = data->eparam[5];
		data->efn2b = data->eparam[6];
		data->efn2c = data->eparam[10];
		data->auend=0;
		data->gubonus=0;
		data->cslope = 0;
		data->cint=0;
		data->c3=0;
		data->init=0;
		data->gail = 0;
	}
	else {
		data->efn2a = short (floor (temp*10+.5));  //constant multi-loop penalty for efn2 (offset)

		ml1 >> temp;
		data->efn2b= short (floor(temp*10+.5));  	//efn2 multi-loop free-base penalty

		ml1 >> temp;
		data->efn2c= short (floor(temp*10+.5)); 	//efn2 multi-loop helix penalty

		//now read the terminal AU penalty:
		ml1>>lineoftext;
		while(strcmp(lineoftext,"-->")) {
			ml1>>lineoftext;
		}

		ml1>> temp;
		data->auend = short (floor (temp*10+.5));	//terminal au penalty

		//now read the GGG hairpin bonus:
		ml1>>lineoftext;
		while(strcmp(lineoftext,"-->")) {
			ml1>>lineoftext;
		}

		ml1 >> temp;
		data->gubonus = short (floor (temp*10+.5));  //bonus for triple G hairpin

		//now read the poly c hairpin penalty slope:
		ml1>>lineoftext;
		while(strcmp(lineoftext,"-->")) {
			ml1>>lineoftext;
		}

		ml1 >> temp;
		data->cslope = short (floor (temp*10+.5));	//poly-c hairpin slope

		//now read the poly c hairpin penalty intercept:
		ml1>>lineoftext;
		while(strcmp(lineoftext,"-->")) {
			ml1>>lineoftext;
		}

		ml1 >> temp;
		data->cint = short (floor (temp*10+.5));	//poly-c hairpin intercept

		//now read the poly c penalty for a loop of 3:
		ml1>>lineoftext;
		while(strcmp(lineoftext,"-->")) {
			ml1>>lineoftext;
		}

		ml1 >> temp;
		data->c3 = short (floor (temp*10+.5));		//c hairpin of 3

		ml1>>lineoftext;
		while(strcmp(lineoftext,"-->")) {
			ml1>>lineoftext;
		}

		ml1 >> temp;
		data->init = short (floor (temp*10+.5));	//intermolecular initiation free energy

		//now read the GAIL rule indicator
		ml1>>lineoftext;
		while(strcmp(lineoftext,"-->")) {
			ml1>>lineoftext;
		}

		ml1>> temp;
		data->gail = short (floor (temp+.5));		//GAIL rule on or off (1 or 0)
	}


/*	read info from dangle */

	//add to dangle the case where X (represented as 0) is looked up
	for (l = 1;l <=2; l++){
		for (i = 0;i <=5; i++){
			if ((i!=0)&&(i!=5)) {
				for (count=1;count <=60;count++) da1 >> lineoftext;
				//read in text headers for each row of tables
			}
			for (j=0;j<=5; j++) {
				for (k=0;k<=5; k++) {
					if ((i==0)||(j==0)||(k==0)) {
							data->dangle[i][j][k][l] = 0;
					}
					else if ((i==5)||(j==5)||(k==5)) {
						data->dangle[i][j][k][l] = 0;
					}
					else {
						da1 >> lineoftext; //read in next data value
						if (strcmp(lineoftext,".")){
							data->dangle[i][j][k][l] = short (floor (10*(atof(lineoftext))+.5));
						}
						else data->dangle[i][j][k][l] = INFINITE;
					}
				}
			}
		}
	}


/*	read info from loop for internal loops, hairpin loops, and bulge loops */

	for (count = 1; count <=26; count++) lo1 >> lineoftext; //get past text in file

	for (i=1;i <= 30; i++) {
		lo1 >> lineoftext;//get past the size column in table
		lo1 >> lineoftext;

		if (strcmp(lineoftext,".")){
			data->inter[i] = short (floor (10*(atof(lineoftext))+.5));
		}
		else data->inter[i] = INFINITE;

		lo1 >> lineoftext;

		if (strcmp(lineoftext,"."))
			data->bulge[i] = short (floor(10*(atof(lineoftext))+.5));
		else data->bulge[i] = INFINITE;

		lo1 >> lineoftext;

		if (strcmp(lineoftext,".")){
			data->hairpin[i] = short (floor(10*(atof(lineoftext))+.5));
		}
		else data->hairpin[i] = INFINITE;
	}


/* Read info from stack */

	//add to the stack table the case where X (represented as 0) is looked up:
	for (count=1;count<=42;count++) st1 >> lineoftext;//get past text in file

	for (i=0;i<=5;i++) {
		if ((i!=0) && (i!=5))
		  for (count=1;count<=60;count++){
		    st1 >> lineoftext;
			//read in headers for each row of tables
		  }
		for (k=0; k<=5; k++) {
			for (j=0; j<=5; j++) {
				for (l=0; l<=5; l++) {
					if ((i==0)||(j==0)||(k==0)||(l==0)) {
						data->stack[i][j][k][l]=0;
					}
					else if ((i==5)||(j==5)||(k==5)||(l==5)) {
						data->stack[i][j][k][l] = INFINITE;
					}
					else {
						st1 >> lineoftext;
						if (strcmp(lineoftext,".")){
							data->stack[i][j][k][l] =short (floor(10*(atof(lineoftext))+.5));
						}
						else {
							data->stack[i][j][k][l] = INFINITE;
						}
					}
				}
			}
		}
	}


/* Read info from tstackh */

	//add to the tstackh table the case where X (represented as 0) is looked up:
	for (count=1;count<=46;count++) th1 >> lineoftext;//get past text in file
	for (i=0;i<=5;i++) {
		if ((i!=0)&&(i!=5)) for (count=1;count<=60;count++) th1 >> lineoftext;
		for (k=0;k<=5;k++) {
			for (j=0;j<=5;j++) {
				for (l=0;l<=5;l++) {
				  if ((i==0)||(j==0)||(k==0)||(l==0)) { //0 is an unknown nucleotide
						data->tstkh[i][j][k][l]=0;
					}
					else if ((i==5)||(j==5)||(k==5)||(l==5)) { //5 is a gap
						data->tstkh[i][j][k][l] = INFINITE;
					}
					else {
						th1 >> lineoftext;
						if (strcmp(lineoftext,".")){
							data->tstkh[i][j][k][l] =short (floor(10*(atof(lineoftext))+.5));
						}
						else data->tstkh[i][j][k][l] = INFINITE;
					}
				}
			}
		}
	}

/* Read info from tstacki */

	//add to the tstacki table the case where X (represented as 0) is looked up:
	for (count=1;count<=46;count++) ti1 >> lineoftext;//get past text in file
	for (i=0;i<=5;i++) {
		if ((i!=0)&&(i!=5)) for (count=1;count<=60;count++) ti1 >> lineoftext;
		for (k=0;k<=5;k++) {
			for (j=0;j<=5;j++) {
				for (l=0;l<=5;l++) {
					if ((i==0)||(j==0)||(k==0)||(l==0)) {
						data->tstki[i][j][k][l]=0;
					}
					else if ((i==5)||(j==5)||(k==5)||(l==5)) {
						data->tstki[i][j][k][l] = INFINITE;
					}
					else {
						ti1 >> lineoftext;
						if (strcmp(lineoftext,".")){
							data->tstki[i][j][k][l] =short (floor (10*(atof(lineoftext))+.5));
						}
						else data->tstki[i][j][k][l] = INFINITE;
					}
				}
			}
		}
	}


/*	Read info from tloops */
	for (count=1;count<=3;count++)	tl1 >> lineoftext;//get past text in file
	data->numoftloops=0;
	tl1>>lineoftext;

	for (count=1;count<=MAXTLOOP&&!tl1.eof();count++){
		(data->numoftloops)++;

		*(base+1) = 0;

		*base = *lineoftext;
		data->tloop[data->numoftloops][0] = tonumi(base);

		*base = *(lineoftext+1);
		data->tloop[data->numoftloops][0] = data->tloop[data->numoftloops][0]+
				5*tonumi(base);

		*base = *(lineoftext+2);
		data->tloop[data->numoftloops][0] = data->tloop[data->numoftloops][0]+
				25*tonumi(base);

		*base = *(lineoftext+3);
		data->tloop[data->numoftloops][0] = data->tloop[data->numoftloops][0]+
				125*tonumi(base);

		*base = *(lineoftext+4);
		data->tloop[data->numoftloops][0] = data->tloop[data->numoftloops][0]+
				625*tonumi(base);

		*base = *(lineoftext+5);
		data->tloop[data->numoftloops][0] = data->tloop[data->numoftloops][0]+
				3125*tonumi(base);

		tl1 >> temp;
		data->tloop[data->numoftloops][1] = short (floor (10*temp+0.5));
		tl1 >> lineoftext;

/*

		strcpy(base,lineoftext);
		strcpy(base+1,"\0");
		data->tloop[data->numoftloops][0] = tonumi(base);
		strcpy(base,lineoftext+1);
		strcpy(base+1,"\0");
		data->tloop[data->numoftloops][0] = data->tloop[data->numoftloops][0]+
				5*tonumi(base);
		strcpy(base,lineoftext+2);
		strcpy(base+1,"\0");
		data->tloop[data->numoftloops][0] = data->tloop[data->numoftloops][0]+
				25*tonumi(base);
		strcpy(base,lineoftext+3);
		strcpy(base+1,"\0");
		data->tloop[data->numoftloops][0] = data->tloop[data->numoftloops][0]+
				125*tonumi(base);
		strcpy(base,lineoftext+4);
		strcpy(base+1,"\0");
		data->tloop[data->numoftloops][0] = data->tloop[data->numoftloops][0]+
				625*tonumi(base);
		strcpy(base,lineoftext+5);
		strcpy(base+1,"\0");
		data->tloop[data->numoftloops][0] = data->tloop[data->numoftloops][0]+
				3125*tonumi(base);
		tl1 >> temp;
		data->tloop[data->numoftloops][1] = short (floor (10*temp+0.5));
		tl1 >> lineoftext;
*/
	}


	//Read the 2x2 internal loops
	//key iloop22[a][b][c][d][j][l][k][m] =
	//a j l b
	//c k m d

	for (count=1;count<=340;count++) in1 >> lineoftext;//get past text in file

	for (i=1;i<=36;i++) {//read each of 21 tables
		for (j=1;j<=39;j++) in1 >> lineoftext;//get past text in file

		*(base +1) = 0;

		*base = *(lineoftext);
		a = tonumi(base);

		for (j=1;j<=3;j++) in1 >> lineoftext;

		*base = *(lineoftext);
		b = tonumi(base);

		in1>>lineoftext;
		*base =*(lineoftext);
		c = tonumi(base);

		for (j=1;j<=3;j++) in1 >> lineoftext;

		*base = *(lineoftext);
		d = tonumi(base);
/*
  strcpy(base,lineoftext);
		strcpy(base+1, "\0");
		a = tonumi(base);

		for (j=1;j<=3;j++) in1 >> lineoftext;

		strcpy(base, lineoftext);
		strcpy(base+1, "\0");
		b = tonumi(base);
		in1>>lineoftext;
		strcpy(base, lineoftext);
		strcpy(base+1, "\0");
		c = tonumi(base);

		for (j=1;j<=3;j++) in1 >> lineoftext;

		strcpy(base, lineoftext);
		strcpy(base+1, "\0");
		d = tonumi(base);
*/

		for (j=1;j<=3;j++) in1 >> lineoftext;//get past text in file

		for (j=1;j<=4;j++) {
			for (k=1;k<=4;k++) {
				for (l=1;l<=4;l++) {
					for (m=1;m<=4;m++) {
						in1 >> temp;
						data->iloop22[a][b][c][d][j][l][k][m] = short (floor(10*temp+0.5));
						//no longer need to store the reverse order at same time because
						//the tables contain redundancy:
						//	data->iloop22[d][c][b][a][m][k][l][j] = floor(10*temp+0.5);
					}
				}
			}
		}
	}

	//Read the 2x1 internal loop data
	for (i=1;i<=58;i++) in2 >> lineoftext; //get past text at top of file
	for (i=1;i<=6;i++) { //read each row of tables
		for (e=1;e<=4;e++) {
			for (j=1;j<=66;j++) in2 >> lineoftext; //get past text in file

			*(base+1) = 0;
			in2 >> lineoftext;
			*base = *lineoftext;
			a = tonumi(base);

			for (j=1;j<=11;j++) in2 >> lineoftext; //get past text in file

			in2 >> lineoftext;
			*base = *lineoftext;
			b = tonumi(base);
/*

			in2 >> lineoftext;
			strcpy(base,lineoftext);
			strcpy(base+1,"\0");
			a = tonumi(base);

			for (j=1;j<=11;j++) in2 >> lineoftext; //get past text in file

			in2 >> lineoftext;
			strcpy(base,lineoftext);
			strcpy(base+1,"\0");
			b = tonumi(base);
*/

			for (j=1;j<=35;j++) in2 >> lineoftext; //get past text in file

			for (c=1;c<=4;c++) {
				for (j=1;j<=6;j++) {
					switch (j) {
						case 1:
							f = 1;
							g = 4;
							break;
						case 2:
							f = 2;
							g = 3;
							break;
						case 3:
							f = 3;
							g = 2;
							break;
						case 4:
							f = 4;
							g = 1;
							break;
						case 5:
							f = 3;
							g = 4;
							break;
						case 6:
							f = 4;
							g = 3;
							break;
					}
					for (d=1;d<=4;d++) {
						in2 >> temp;
						data->iloop21[a][b][c][d][e][f][g]=short (floor(10*temp+0.5));
					}
				}
			}
		}
	}

/*	Read info from triloops */
	for (count=1;count<=3;count++)	tri >> lineoftext;//get past text in file
	data->numoftriloops=0;
	tri>>lineoftext;

	for (count=1;count<=MAXTLOOP&&!tri.eof();count++){
		(data->numoftriloops)++;

		*(base+1) = 0;

		*base = *(lineoftext);
		data->triloop[data->numoftriloops][0] = tonumi(base);

		*base = *(lineoftext+1);
		data->triloop[data->numoftriloops][0] = data->triloop[data->numoftriloops][0]+
				5*tonumi(base);

		*base = *(lineoftext+2);
		data->triloop[data->numoftriloops][0] = data->triloop[data->numoftriloops][0]+
				25*tonumi(base);

		*base = *(lineoftext+3);
		data->triloop[data->numoftriloops][0] = data->triloop[data->numoftriloops][0]+
				125*tonumi(base);

		*base = *(lineoftext+4);
		data->triloop[data->numoftriloops][0] = data->triloop[data->numoftriloops][0]+
				625*tonumi(base);

		tri >> temp;
		data->triloop[data->numoftriloops][1] = short (floor (10*temp+0.5));
		tri >> lineoftext;
/*
strcpy(base,lineoftext);
		strcpy(base+1,"\0");
		data->triloop[data->numoftriloops][0] = tonumi(base);
		strcpy(base,lineoftext+1);
		strcpy(base+1,"\0");
		data->triloop[data->numoftriloops][0] = data->triloop[data->numoftriloops][0]+
				5*tonumi(base);
		strcpy(base,lineoftext+2);
		strcpy(base+1,"\0");
		data->triloop[data->numoftriloops][0] = data->triloop[data->numoftriloops][0]+
				25*tonumi(base);
		strcpy(base,lineoftext+3);
		strcpy(base+1,"\0");
		data->triloop[data->numoftriloops][0] = data->triloop[data->numoftriloops][0]+
				125*tonumi(base);
		strcpy(base,lineoftext+4);
		strcpy(base+1,"\0");
		data->triloop[data->numoftriloops][0] = data->triloop[data->numoftriloops][0]+
				625*tonumi(base);
		tri >> temp;
		data->triloop[data->numoftriloops][1] = short (floor (10*temp+0.5));
		tri >> lineoftext;
*/
	}


/* Read info from coax */

	//add to the stack table the case where X (represented as 0) is looked up:
	//this is the array that keeps track of flush coaxial stacking (no intervening nucs)
	//data arrangement of coax: data->coax[a][b][c][d]
	//5'b-c3'
	//3'a d5'
	//this means the helix backbone is continuous between nucs b and c

	for (count=1;count<=42;count++) co1 >> lineoftext;//get past text in file
	for (i=0;i<=5;i++) {
		if ((i!=0)&&(i!=5)) for (count=1;count<=60;count++) co1 >> lineoftext;
		for (k=0;k<=5;k++) {
			for (j=0;j<=5;j++) {
				for (l=0;l<=5;l++) {
					if ((i==0)||(j==0)||(k==0)||(l==0)) {
						data->coax[j][i][k][l]=0;
					}
					else if ((i==5)||(j==5)||(k==5)||(l==5)) {
						data->coax[j][i][k][l] = INFINITE;
					}
					else {
						co1 >> lineoftext;

						if (strcmp(lineoftext,".")){
							data->coax[j][i][k][l] =short (floor(10*(atof(lineoftext))+.5));
						}
						else data->coax[j][i][k][l] = INFINITE;
					}
				}
			}
		}
	}

/* Read info from tstackcoax */

	//add to the tstackh table the case where X (represented as 0) is looked up:
	//data arrangement of tstackcoax:
	//5'a-c -> strand continues into stack
	//3'b-d -> strand does not continue to stack
	//pair between a-b, c-d is a mismatch

	for (count=1;count<=46;count++) co2 >> lineoftext;//get past text in file
	for (i=0;i<=5;i++) {
		if (!(i==0||i==5)) for (count=1;count<=60;count++) co2 >> lineoftext;
		for (k=0;k<=5;k++) {
			for (j=0;j<=5;j++) {
				for (l=0;l<=5;l++) {
					if ((i==0)||(j==0)||(k==0)||(l==0)||(i==5)||(j==5)||(k==5)||(l==5)) {
						data->tstackcoax[i][j][k][l]=0;
					}
					else {
						co2 >> lineoftext;
						if (strcmp(lineoftext,".")){
							data->tstackcoax[i][j][k][l] =short (floor(10*(atof(lineoftext))+.5));
						}
						else data->tstackcoax[i][j][k][l] = INFINITE;
					}
				}
			}
		}
	}

/* Read info from coaxstack */

	//add to the tstackh table the case where X (represented as 0) is looked up:
	//data arrangement of coaxstack:
	//5'a-c ->strand contnues into stack
	//3'b d ->strand does not continue to stack
	//pair between a-b, mismatch between c-d
	//backbone is discontinuous between b and d

	for (count=1;count<=46;count++) co3 >> lineoftext;//get past text in file
	for (i=0;i<=4;i++) {
		if (!(i==0||i==5)) for (count=1;count<=60;count++) co3 >> lineoftext;
		for (k=0;k<=4;k++) {
			for (j=0;j<=4;j++) {
				for (l=0;l<=4;l++) {
					if ((i==0)||(j==0)||(k==0)||(l==0)||(i==5)||(j==5)||(k==5)||(l==5)) {
						data->coaxstack[i][j][k][l]=0;
					}
					else {
						co3 >> lineoftext;
						if (strcmp(lineoftext,".")){
							data->coaxstack[i][j][k][l] =short (floor(10*(atof(lineoftext))+.5));
						}
						else data->coaxstack[i][j][k][l] = INFINITE;
					}
				}
			}
		}
	}

/* Read info from tstack */

	//this is the terminal mismatch data used in intermolecular folding
	//add to the tstack table the case where X (represented as 0) is looked up.
	//also add the case where 5 (the intermolecular linker) is looked up,
	//this is actually a dangling end, not a terminal mismatch.

	for (count=1;count<=46;count++) st2 >> lineoftext;//get past text in file
	for (i=0;i<=5;i++) {
		if ((i!=0)&&(i!=5)) for (count=1;count<=60;count++) st2 >> lineoftext;
		for (k=0;k<=5;k++) {
			for (j=0;j<=5;j++) {
				for (l=0;l<=5;l++) {
					if ((i==0)||(j==0)||(k==0)||(l==0)) {
						data->tstack[i][j][k][l]=0;
					}
					else if ((i==5)||(j==5)) {
							data->tstack[i][j][k][l] = INFINITE;
					}
					else if ((k==5)||(l==5)) {
					//include "5", linker for intermolecular for case of flush ends
						if ((k==5)&&(l==5)) {//flush end
							data->tstack[i][j][k][l]=0;
						}
						else if (k==5) {//5' dangling end
							//look up number for dangling end
							data->tstack[i][j][k][l] = data->dangle[i][j][l][2]+penalty2(i,j,data);
						}
						else if (l==5) {//3' dangling end
							data->tstack[i][j][k][l] = data->dangle[i][j][k][1]+penalty2(i,j,data);
						}
					}
					else {
						st2 >> lineoftext;
						if (strcmp(lineoftext,".")){
							data->tstack[i][j][k][l] =short (floor (10*(atof(lineoftext))+.5));
						}
						else data->tstack[i][j][k][l] = INFINITE;
					}
				}
			}
		}
	}

/* Read info from tstackm */

	//add to the tstackm table the case where X (represented as 0) is looked up:

	for (count=1;count<=46;count++) tsm >> lineoftext;//get past text in file
	for (i=0;i<=4;i++) {
		if (i!=0) for (count=1;count<=60;count++) tsm >> lineoftext;
		for (k=0;k<=4;k++) {
			for (j=0;j<=4;j++) {
				for (l=0;l<=4;l++) {
					if ((i==0)||(j==0)||(k==0)||(l==0)) {
						data->tstkm[i][j][k][l]=0;
					}
					else {
						tsm >> lineoftext;
						if (strcmp(lineoftext,".")){
							data->tstkm[i][j][k][l] =short (floor(10*(atof(lineoftext))+.5));
						}
						else data->tstkm[i][j][k][l] = INFINITE;
					}
				}
			}
		}
	}


	//data arrangement for 1x1 loop tables iloop11[a][b][c][d][e][f]:
	//abc
	//def

	//Read the 1x1 internal loop data
	//encode the data like:  abc
	//                       def where b-e is a mismatch
	for (i=1;i<=58;i++) i11 >> lineoftext; //get past text at top of file
	for (i=1;i<=6;i++) { //read each row of table
		if (i==1) {
				a = 1;
				d = 4;
		}
		else if (i==2) {
				a = 2;
				d = 3;
		}
		else if (i==3) {
				a = 3;
				d = 2;
		}
		else if (i==4) {
				a = 4;
				d = 1;
		}
		else if (i==5) {
				a = 3;
				d = 4;
		}
		else {
			a = 4;
				d = 3;
		}
		for (j=1;j<=114;j++) i11 >> lineoftext;//get past text
		for (b=1;b<=4;b++) {
			for (j=1;j<=6;j++) {
				if (j==1) {
					c = 1;
					f = 4;
				}
				else if (j==2) {
					c = 2;
					f = 3;
				}
				else if (j==3) {
					c = 3;
					f = 2;
				}
				else if (j==4) {
					c = 4;
					f = 1;
				}
				else if (j==5) {
					c = 3;
					f = 4;
				}
				else {
					c = 4;
					f = 3;
				}
				for (e=1;e<=4;e++) {
					i11 >> temp;
					data->iloop11[a][b][c][d][e][f]=short (floor(10*temp+0.5));
				}
			}
		}
	}

	return 1;
}	// end of opendat()


/*                               -*- Mode: C -*- 
 * pd_interface.cpp ---
 * Author          : Amelia Bellamy-Royds
 * Created         : Summer 2006
 * Last Modified By: Marcel Turcotte
 * Last Modified On: Mon Dec 11 12:17:45 2006
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

/* main method, plus a method with hard-coded data file names*/

#include "pdynalign.cpp"


/*--------------------------------------------------------------------------------
Function:	getdat()

Description:	Gets the names of all the data files to be openned. All data files
		should located in a directory identified by an environment variable
		specifed by the defined constant DATA_DIR, if this environment
		variable is not defined, then the data files must be in the current
		working directory.
--------------------------------------------------------------------------------*/
void getdat(char *loop, char *stackf, char *tstackh, char *tstacki,
		char *tloop, char *miscloop, char *danglef, char *int22,
    char *int21,char *coax, char *tstackcoax,
    char *coaxstack, char *tstack, char *tstackm,char *triloop, char *int11)
{

  //look for an environment variable to identify the directory of data files,
  //otherwise use current dir

  char * env = getenv(DATA_DIR);
  char * dir;
  if (env == NULL) {
    dir = "./";
  }
  else {
    dir = new char[strlen(env) + 2];
    strcpy(dir, env);
    strcat(dir, "/");
  }
  cout << "looking for data files in directory " << dir << "\n";

	strcpy (loop, dir);
	strcpy (stackf, dir);
	strcpy (tstackh, dir);
	strcpy (tstacki, dir);
	strcpy (tloop, dir);
	strcpy (miscloop, dir);
	strcpy (danglef, dir);
	strcpy (int22, dir);
	strcpy (int21, dir);
	strcpy (coax, dir);
	strcpy (tstackcoax, dir);
	strcpy (coaxstack, dir);
	strcpy (tstack, dir);
	strcpy (tstackm, dir);
	strcpy (triloop, dir);
	strcpy (int11, dir);

	strcat (loop,"loop.dat");
	strcat (stackf,"stack.dat");
	strcat (tstackh,"tstackh.dat");
	strcat (tstacki,"tstacki.dat");
	strcat (tloop,"tloop.dat");
	strcat (miscloop,"miscloop.dat");
	strcat (danglef,"dangle.dat");
	strcat (int22,"int22.dat");
	strcat (int21,"int21.dat");
	strcat (coax,"coaxial.dat");
	strcat (tstackcoax,"tstackcoax.dat");
	strcat (coaxstack,"coaxstack.dat");
	strcat (tstack,"tstack.dat");
	strcat (tstackm,"tstackm.dat");
	strcat (triloop,"triloop.dat");
	strcat (int11,"int11.dat");
}


int main(int argc, char* argv[])
{
	/* MAXFIL - maximum length of file names, defined in 'defines.h'.  */
	char inpf1[MAXFIL],inpf2[MAXFIL],aout[MAXFIL];
	char outct1[MAXFIL],outct2[MAXFIL];
	structure ct1,ct2;
	profile pf1,pf2;
	int i,totalscore,imaxseparation,igapincrease;
	datatable data;
	char loop[MAXFIL],stackf[MAXFIL],tstackh[MAXFIL],tstacki[MAXFIL],
			tloop[MAXFIL],miscloop[MAXFIL],danglef[MAXFIL],int22[MAXFIL],
      int21[MAXFIL],coax[MAXFIL],tstackcoax[MAXFIL],
      coaxstack[MAXFIL],tstack[MAXFIL],tstackm[MAXFIL],triloop[MAXFIL],int11[MAXFIL];
	alignment	*align = NULL;

	bool insert;
	float fgap;

	cout << "Profile-Dynalign 1.0 - alignment and structure prediction of two RNA profiles\n\n";

	cout << "Copyright (C) 2006 University of Ottawa\n";
        cout << "All Rights Reserved\n\n";

	cout <<  "This program is distributed under the terms of the\n";
	cout <<  "GNU General Public License. See the source code\n";
	cout <<  "for details.\n\n";

	try {

	  cout << "Loading data files...\n";

		/* Gets names of all data files. Implemented in this file. */
		getdat (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,
			int21,coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11);

	  /* Collect info. from data files into 'data'.
		 opendat is implemented in 'algorithm.cpp'.
		 data is a datatable which defined in 'datatable.cpp' */
		opendat (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,int21,
			coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11, &data) ;
		cout << "ready\n\n";

		cout << "Enter the name of the first profile:\n";
		cin >> inpf1;
		cout << ">>Opening profile file " << inpf1 << "\n";
		openProfile(&pf1, inpf1);
		cout << "\n";

		cout << "Enter the name of the second profile:\n";
		cin >> inpf2;
		cout << ">>Opening profile file " << inpf2 << "\n";
		openProfile(&pf2, inpf2);
		cout << "\n";

		cout << "Enter the max separation:\n";
		cin >> imaxseparation;
		if (imaxseparation <= 0) throw HDException(BAD_INPUT, "Maximum separation value must be a positive integer");
		cout << ">>Using maximum separation " << imaxseparation << "\n\n";

		cout << "Enter the gap penalty (in kcal/mol, to the nearest tenth of a kcal) :\n";
		cin >> fgap;
		igapincrease = (int) (fgap*10.0);		/* transfer gap penalty to integer. */
		cout << ">>Using gap penalty " << (float(igapincrease)/10) << " kcal/mol\n\n";

		cout << "Allow single BP inserts into one sequence? (1/0)  ";
		cin >> i;
		if (i) insert=true;	/* 1 -> true, trun on insert mode. */
		else insert=false;	/* 0 -> false, turn off insert mode. */
		cout << ">>" << (i? "  yes ": "  no  ") << "\n\n";

		cout << "Enter the directory for output ct files of the first profile:\n";
		cin >> outct1;
		cout << ">>Output directory 1:  " << outct1 << "\n\n" ;

		cout << "Enter the directory for output ct files of the second profile:\n";
		cin >> outct2;
		cout << ">>Output directory 2:  " << outct2 << "\n\n";

		cout << "Enter the file name for the output of the alignment:\n";
		cin >> aout;
		cout << ">>Output alignment file:  " << aout << "\n\n";

		align = new alignment [pf1.numofbases+1];

		pdynalign(&pf1,&pf2,align,imaxseparation,igapincrease,&data,insert,&totalscore);

		ctout_pf(&pf1,outct1);

		ctout_pf(&pf2,outct2);

		alignout( align, aout, &pf1, &pf2, totalscore, imaxseparation, igapincrease, insert);
		//cout << "!";

	}
	catch (HDException& ex) {
		//any errors during pdynalign
		cout << "!!Improper Termination -- see error stream for more!!\n" << flush;
		cerr << ex.message << "\n" << flush;
		if (align != NULL) delete[] align;
		exit (ex.errorvalue);
	}
	catch (...) {
		cerr << "Unknown exception occurred";
	}

	ASSERT( align != NULL );
	delete[] align;

	cout << "************************************************************\n\n";

	return 0;
}




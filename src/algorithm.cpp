/*                               -*- Mode: C -*- 
 * algorithm.cpp ---
 * Author          : Amelia Bellamy-Royds
 * Created         : Summer 2006
 * Last Modified By: Marcel Turcotte
 * Last Modified On: Mon Dec 11 12:40:40 2006
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
  This file contains various utility functions for the program.

  tonum() converts a character-based sequence into numerical code
  tonumi() converts a single character to a number code
  push() and pull() for stackstructs
  ctout_pf() outputs individual ct files for each sequence in a profile
  tobase() reverses tonumi()
  openProfile() opens a profile file and creates a profile object
  copySubStr()
  isValidChar()
  alignout() creates output alignment profile file -- uses next 3 methods:
  append_nucs()
  append_str()
  init_lines()

*/

#include <stdio.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <stdlib.h>

#include "platform.cpp"
#include "pdynalign.h"

/*--------------------------------------------------------------------------------
  Function:	tonum()

  Description: 	Converts base chars to integer numbers, and put into
  structure->numseq[count]. Once encounter "a", "c", "g", "u", or "t"(they
  are forced to be single), append the position of these chars to the array
  nopair[] and increase counter nnopair

  A - 1
  C - 2
  G - 3
  U - 4
  I - 5
  '-'- 5
  others -	0
  --------------------------------------------------------------------------------*/

void tonum(char *base,structure *ct,int count)	{
  if (!strcmp(base,"A")) (ct->numseq[count] = 1);
  else if(!strcmp(base,"B")) {
    (ct->numseq[count] = 1);
  }
  else if(!strcmp(base,"a")) {
    ct->numseq[count]=1;
    ct->nnopair++;
    ct->nopair[ct->nnopair] = count;
  }
  else if(!strcmp(base,"C")) (ct->numseq[count] = 2);
  else if(!strcmp(base,"Z")) {
    (ct->numseq[count] = 2);
  }
  else if(!strcmp(base,"c")) {
    ct->numseq[count] = 2;
    ct->nnopair++;
    ct->nopair[ct->nnopair] = count;
  }
  else if(!strcmp(base,"G")) (ct->numseq[count] = 3);
  else if(!strcmp(base,"H")) {
    (ct->numseq[count] = 3);
  }
  else if(!strcmp(base,"g")) {
    ct->numseq[count] = 3;
    ct->nnopair++;
    ct->nopair[ct->nnopair] = count;
  }

  else if(!strcmp(base,"U")||!strcmp(base,"T")) (ct->numseq[count] = 4);
  else if(!strcmp(base,"V")||!strcmp(base,"W")) {
    (ct->numseq[count] = 4);
  }
  else if(!strcmp(base,"u")||!strcmp(base,"t")) {
    ct->numseq[count] = 4;
    ct->nnopair++;
    ct->nopair[ct->nnopair] = count;
  }
  else if(!strcmp(base,"I")) {
    ct->numseq[count] = 5;
    ct->intermolecular= true;
  }
  else if(!strcmp(base,"-")) {
    ct->numseq[count] = 5; //
    //ct->intermolecular= true;
  }
  else (ct->numseq[count]=0);  //this is for others, like X

  return;
}	// end of tonum()


/*"stackstruct" is replaced by "stackclass"
 *--------------------------------------------------------------------------------
 Function:	push()

 Description: Push data into a stack
 --------------------------------------------------------------------------------*
 void push(stackstruct *stack,int a,int b,int c,int d)
 {
 (stack->sp)++;
 stack->stk[stack->sp][0]= a;
 stack->stk[stack->sp][1]= b;
 stack->stk[stack->sp][2]= c;
 stack->stk[stack->sp][3]= d;
 }

 *--------------------------------------------------------------------------------
 Function:	pull()

 Description: Pull data from a stack
 --------------------------------------------------------------------------------*
 void pull(stackstruct *stack,int *i,int *j,int *open,int *null,int *stz)
 {
 if (stack->sp==0) {
 *stz = 1;
 return;
 }
 else {
 *stz = 0;
 *i = stack->stk[stack->sp][0];
 *j = stack->stk[stack->sp][1];
 *open= stack->stk[stack->sp][2];
 *null= stack->stk[stack->sp][3];
 stack->sp--;
 }
 }

*/

/*-----------------------------------------------------------------------------
  Function: ctout_pf()

  Description:	Outputs a set of ct files from a profile to a directory. The
  profile is specified by the	pointer 'pf' and the directory is specified
  by 'ctoutdir'
  -----------------------------------------------------------------------------*/
void ctout_pf (profile *pf, char *ctoutdir) {
  int count,i;	//length
  char line[CT_LINE_LN]; // number[2*NUMLEN], base[2];
  char ctoutfile[MAXFIL + 2*SEQNAME_LN];
  structure *ct;

  cout << "Creating output ct files for profile " << pf ->id << ":\n" ;

  for (count = 0; count < pf->nSeq; count++) {
    ct = pf->ct + count; //i.e. &(pf->ct[count])
    strcpy(ctoutfile, ctoutdir);
    strcat(ctoutfile, "/");
    strcat(ctoutfile, ct->ctlabel[1]);
    strcat(ctoutfile, ".ct");

    FILE *ctfile;
    ctfile=fopen(ctoutfile,"w");	// Open ctfile for writing
    if (ctfile == NULL)
      {
	//warn of non-fatal error and go onto next file
	cerr << "Warning: Could not create output file "<< ctoutfile << "\n" ;
	break;
      }


    // Output the header
    strcpy(line,"");
    sprintf(line,"%5i",ct->numofbases);

    sgifix  // This corrects a difference on the sgi computers

      //NOTE that the program does not currently calculate energies for individual structures:
      strcat(line,"   dG = ");
    //gcvt((float (ct->energy[1]))/10,6,number);
    strcat(line,"0.0");
    strcat(line,"  [initially  ");
    strcat(line,"0.0");
    strcat(line,"]    ");

    strcat(line,ct->ctlabel[1]);
    strcat(line, "\n");
    fputs (line,ctfile);

    // Output nucleotides and their associated info
    for (i=1;i<pf->numofbases;i++) {
      if ( ct->numseq[i] != 5 ) { //i.e. not a gap

	ASSERT( ct->hnumber[i] != 0 );
	ASSERT( pf->basepr[i]==0 || ct->hnumber[pf->basepr[i]] != 0 );

	sprintf(line,"%5i%2c%8i%5i%5i%5i\n", ct->hnumber[i], ct->nucs[i],
		ct->hnumber[i]-1,ct->hnumber[i]+1,ct->hnumber[pf->basepr[i]],
		ct->hnumber[i]);
	fputs(line,ctfile);
      }
    }

    i = pf->numofbases;
    if ( ct->numseq[i] != 5 ) { //i.e. not a gap
      sprintf(line,"%5i %1c%8i%5i%5i%5i\n", ct->hnumber[i],ct->nucs[i],
	      ct->hnumber[i]-1,0,ct->hnumber[pf->basepr[i]],ct->hnumber[i]);
      fputs(line,ctfile);
    }

    fclose (ctfile);	// Close ctfile
    cout << "Wrote file " << ctoutfile << "\n";
  }

  return;
} // end of ctout_pf()


/*-------------------------------------------------------------------------
  Function: tobase()

  Description:	Converts an integer to a base symbol
  //NOT USED (original sequence stored as character array nucs[], instead)
  -------------------------------------------------------------------------* /
  char *tobase (int i)
  {
  if (i==1) return "A";
  else if (i==2) return "C";
  else if (i==3) return "G";
  else if (i==4) return "U";
  else if (i==0) return "X";
  else if (i==5) return "I";
  else return "?";   //to track down mistakes
  }
*/


/*---------------------------------------------------------------------------------------
  Function: openProfile()

  Description:	Opens a profile file and collects data for the profile structure 'pf'

  Throws:			HDException if encounters error opening or parsing file

  Return value:	DONE
  ---------------------------------------------------------------------------------------*/
int openProfile (profile *pf, char *fileName)
{
  int		i, j, k, nSeq, seqLength, linecount;
  bool done;

  FILE	*pfFile;  //pointer to the file object
  char line[SEQ_LN_MAX]; //a single line from the file
  //may contain up to a full sequence
  char	base[2];	// stores one base character in the sequence
  //as a null-terminated string
  int seqlines[NUM_SEQ_MAX]; //line numbers within the file containing sequences
  char * cursor;
  structure	*ctptr;		// a pointer to a structure 'ct'

  //char	sq[sizeof(int)];		// string stores the number of sequences
  //char	*buf[NUM_SEQ_MAX+5];		// buffer for reading the file
  //int		bufSize;			// actual size of the buffer

  // open file
  if( (pfFile = fopen(fileName, "r")) == NULL)
    {
      char *message = new char[100];
      strcpy(message,"Error: Cannot open ");
      strcat(message, fileName);
      HDException ex(IO_ERR, message);
      throw ex;
    }

  linecount = 0;
  if (fgets(line, SEQ_LN_MAX, pfFile) != NULL) {
    //check that first line matches the STOCKHOLM 1.0 header:
    if (strncmp(line, STOCKHOLM_HEADER, strlen(STOCKHOLM_HEADER))) {
      //Issue warning -- not expected file header
      //(but try to read file anyways, in case header is all that is wrong)
      cerr << "Input file " << fileName << " does not have expected format header\n";
      cerr << " -- expected format: " << STOCKHOLM_HEADER <<".\n";
      rewind(pfFile);
    }
    else {
      //okay, continue at next line
      linecount = 1;
    }
  }
  else if (feof(pfFile)) {
    //file was empty -- throw exception
    char *message = new char[100];
    strcpy(message,"Error: Empty input file ");
    strcat(message, fileName);
    HDException ex(BAD_INPUT, message);
    throw ex;
  }
  else {
    //there was an error when reading the file -- throw exception
    char *message = new char[100];
    strcpy(message,"Error: Cannot read from ");
    strcat(message, fileName);
    HDException ex(IO_ERR, message);
    throw ex;
  }

  nSeq = 0;
  seqLength = 0;
  done = false;
  //scan through file and identify all sequences
  while ((!done)&&fgets(line, SEQ_LN_MAX, pfFile) != NULL) {
    linecount++;

    //check that a full line was copied
    if (strlen(line) >= SEQ_LN_MAX - 1) {
      //file line too long -- throw exception
      char *message = new char[200];
      strcpy(message,"Error: Line too long to parse in file ");
      strcat(message, fileName);
      HDException ex(BAD_INPUT, message);
      throw ex;
    }
    if (!strncmp(line, STOCKHOLM_END, strlen(STOCKHOLM_END))){
      //then have reached the "end of alignment" indicator
      done = true;
      break;
    }


    //NOTE: strcspn(cursor, WHITESPACE) returns the number of characters
    //			from cursor until the start of whitespace
    //			(space, tab, carriage return or linebreak)
    //		strspn(cursor, WHITESPACE) returns the number of characters
    //			from cursor to the first non-whitespace character
    //		LINEBREAK is only carriage return or newline

    if (!strncmp(line, STOCKHOLM_MARKUP_START, strlen(STOCKHOLM_MARKUP_START) )){
      //line is a markup line, so look for useful markup:

      if (!strncmp(line, FILE_MARKUP, strcspn(line, WHITESPACE) )){

	//line is a file markup (#=GF)
	cursor = line + strlen(FILE_MARKUP);
	cursor += strspn(cursor, WHITESPACE);
	if (!strncmp(cursor, ID_MARKER, strcspn(cursor, WHITESPACE) )) {
	  //line is an ID markup (#=GF ID)
	  //save profile ID
	  cursor += strlen(ID_MARKER);
	  cursor += strspn(cursor, WHITESPACE);
	  copySubStr(pf->id, cursor, min(strcspn(cursor, LINEBREAK), ID_LN_MAX-1) );
	  //copySubStr copies a set-sized substring, and adds null-terminator
	}
	//parse other FILE_MARKUP lines here
	//(cursor is set to beginning of the markup-type token)
      }
      //parse other MARKUP lines here
    }
    else {
      //then assume this is a sequence line,
      //record the line number and increase the sequence count:
      seqlines[nSeq++] = linecount;
      //determine the length of the sequence
      cursor = line + strcspn(line, WHITESPACE);
      cursor += strspn(cursor, WHITESPACE);
      seqLength = max(seqLength, strcspn(cursor, LINEBREAK));
    }

  }
  if (!done) {
    if (feof(pfFile)) {
      //finished because reached end of file
      //issue warning -- no end of alignment indicator
      cerr << "Input file " << fileName << " does not have expected end of alignment indicator \n";
      cerr << " -- expected indicator: " << STOCKHOLM_END <<".\n";
    }
    else {
      //there was an error when reading the file -- throw exception
      char *message = new char[100];
      strcpy(message,"Error: While reading from ");
      strcat(message, fileName);
      HDException ex(IO_ERR, message);
      throw ex;
    }
  }
  if (nSeq == 0) {
    //throw exception -- no sequences found
    char *message = new char[100];
    strcpy(message,"Error: No sequences in file ");
    strcat(message, fileName);
    HDException ex(BAD_INPUT, message);
    throw ex;
  }
  if (strlen(pf->id) == 0 ) {
    //then an id markup was not found
    //use the filename instead
    cursor = strrchr(fileName, PATH_SEPARATOR);
    if ((cursor != NULL)&&(strlen(cursor)<ID_LN_MAX-1))
      strcpy(pf->id, cursor + 1);
    else {
      copySubStr(pf->id, fileName, min(strlen(fileName), ID_LN_MAX-1)  );
    }
  }

  //initialize sequence array of profile
  //and read through file a second time, parsing the sequence lines
  pf->nSeq = nSeq;
  pf->numofbases = seqLength;
  pf->allocate(nSeq, seqLength);
  for (i=0; i < nSeq; i++) {
    pf->ct[i].allocate(seqLength+1); //to allow for indexing starting at 1
  }

  rewind(pfFile);
  linecount = 0;
  for (i=0; i<nSeq; i++) {
    //read to the appropriate line in the file
    do {
      fgets(line, SEQ_LN_MAX, pfFile);
    }while (++linecount < seqlines[i]) ;

    ctptr = &pf->ct[i];		// point to the ith sequence in profile

    // get label for a sequence and store as the first (and only) label
    copySubStr(ctptr->ctlabel[1], line, min(strcspn(line, WHITESPACE), SEQNAME_LN -1) );
    cursor = line +strcspn(line, WHITESPACE);
    cursor += strspn(cursor, WHITESPACE);
    //cursor now points to the start of the sequence

    if ((seqLength = strcspn(cursor, LINEBREAK)) != pf->numofbases) {
      //issue warning -- sequences in alignment not all same length
      cerr << "Sequence " << ctptr->ctlabel[i] << " in file " << fileName;
      cerr << " is not the full alignment length -- Gaps added at end.";
    }

    // get bases for a sequence
    k = 1;	// k stores history number
    base[1] = 0; //null-terminator for single-character string
    for (j=0; j < seqLength; j++)
      {

	base[0] = cursor[j];
	if( !isValidChar(base) )
	  {
	    char *message = new char[100];
	    strcpy(message, "Error: Invalid nucleotide in profile ");
	    strcat(message, fileName);
	    HDException ex (FILE_ERR, message);
	    throw ex;
	  }

	//Note: j counts from 0 to seqLength - 1;
	//arrays go from 1 to seqLength
	tonum(base, ctptr, j+1);	// convert base to int, and put in ct->numseq[j]
	ctptr->nucs[j+1]=base[0];	// stores the base
	if(ctptr->nucs[j+1]!='-' && ctptr->nucs[j+1]!='I') {
	  ctptr->hnumber[j+1] = k++;	// stores original position of the base
	}

      }
    while (j < pf-> numofbases) {
      //add gaps at the end if seqLength != numofbases
      ctptr->numseq[++j] = 5;
      ctptr->numseq[j] = '-';
    }

    // set number of (non-gap) bases for the sequence
    ctptr->numofbases = k-1;

  }

  /*
  // read all the lines and store in buf
  i = 0;
  do {
  ASSERT(i<=NUM_SEQ_MAX);

  buf[i] = new char[SEQ_LN_MAX];

  ASSERT(buf[i] != NULL);

  } while(fgets(buf[i++], SEQ_LN_MAX, pfFile) != NULL);	// note: last char in buf is "\n"

  bufSize = i;

  // construct the profile structure

  // 1st line, verify format of the file by comparing the header
  ASSERT((buf[0] != NULL) && (strlen(buf[0])>=1));
  if( strncmp(buf[0], PROFILE_HEADER, strlen(buf[0])-1) != 0 )
  {
  char *message = new char[150];
  strcpy(message,"Error: Incorrect sequence format for ");
  strcat(message, fileName);
  strcat(message, ".  File should be in STOCKHOLM 1.0 format.");
  HDException ex (FILE_ERR, message);
  throw ex;
  }

  // 2nd line, get profile id and set it to the structure profile
  ASSERT( strlen(buf[1])-ID_POS-1 <= ID_LN_MAX );
  copySubStr( pf->id, buf[1]+ID_POS, strlen(buf[1])-ID_POS-1 );

  // 3rd line, get number of sequences in profile
  ASSERT( strlen(buf[1])-SQ_POS-1 <= NUM_SEQ_MAX );
  copySubStr( sq, buf[2]+SQ_POS, strlen(buf[2])-SQ_POS-1 );
  pf->nSeq = atoi(sq);	// set nSeq in structure profile

  // following lines are sequences
  int	lineNum = 4;		// the line number of the first sequence in the file
  ASSERT(buf[lineNum] != NULL);
  pf->numofbases = strlen(buf[lineNum]) - SEQNAME_LN - 1;	// get number of bases

  // allocate memory for ct's in the profile
  pf->allocate(pf->nSeq,pf->numofbases);
  ASSERT(pf->ct != NULL);
  for (i=0; i<pf->nSeq; i++)
  {
  pf->ct[i].allocate(pf->numofbases+1);	//'+1' since 1st base put in idx. 1
  }

  // construct each ct in the profile
  int j = 0;
  for (i=0; i<pf->nSeq; i++)
  {
  ctptr = &pf->ct[i];		// point to the ith sequence in profile

  // get label for a sequence
  copySubStr(ctptr->ctlabel[1], buf[i+lineNum], SEQNAME_LN);
  *(strstr(ctptr->ctlabel[1], " ")) = '\0';

  // get bases for a sequence
  k = 1;	// k stores history number
  for (j=0; j<pf->numofbases; j++)
  {
  copySubStr(base, buf[i+lineNum]+j+SEQNAME_LN, 1);
  strcpy(base+1, "\0");
  if( !isValidChar(base) )
  {
  char *message = new char[100];
  strcpy(message, "Error: Invalid nucleotide in profile ");
  strcat(message, fileName);
  HDException ex (FILE_ERR, message);
  throw ex;
  }

  tonum(base, ctptr, j+1);	// convert base to int, and put in ct->numseq[j]
  ctptr->nucs[j+1]=base[0];	// stores the base
  if(ctptr->nucs[j+1]!='-' && ctptr->nucs[j+1]!='I') {
  ctptr->hnumber[j+1] = k++;	// stores original position of the base
  }
  }

  // set number of bases for a sequence
  ctptr->numofbases = k-1;
  }

  // free memory
  for(i=0; i<bufSize; i++) {
  ASSERT(buf[i]!= NULL);
  delete[] buf[i];
  }
  */
  // Show profile info read from the file
  structure *tmpct;

  cout << "\n************* open profile ********************\n";
  cout << "profile id: " << pf->id<<"   ";
  cout << "num of seq: " << pf->nSeq << "   ";
  cout << "total length: " << pf->numofbases << "\n";
  for (i=0; i<pf->nSeq; i++) {
    tmpct = &pf->ct[i];

    // print label and numofbases
    cout << tmpct->ctlabel[1] << ": ( 1-" << tmpct->numofbases << " )\n";

    // print sequence
    for (j=0; j<pf->numofbases; j++) {
      cout << tmpct->nucs[j+1];
    }
    cout << "\n";
  }

  fclose(pfFile);

  return DONE;
}	// end of openProfile()


/*---------------------------------------------------------------------------------------
  Function: copySubStr()

  Description:
  Copies a substring from the string 'str' to the string 'substr' given the size to
  copy. Add null terminator '\0' to the end of the 'substr'.

  Return value:
  A pointer points to the char array 'substr'

  Note:
  Before calling the function, make sure 'substr' is large enough to store the
  substring.
  ---------------------------------------------------------------------------------------*/
char *copySubStr(char *substr, char *str, size_t size)
{
  ASSERT( substr != NULL && str != NULL );
  ASSERT( size <= strlen(str) );

  char	*strStart = substr;

  while( size-- > 0 ) {
    *substr++ = *str++;
  }
  *substr = '\0';

  return	strStart;
}


/*---------------------------------------------------------------------------
  Function: isValidChar()

  Description:	Verify that if the symbol read from the sequence is a valid
  symbol.

  Return Value:	Return 0 if the char is not found, otherwise return 1
  ---------------------------------------------------------------------------*/
int isValidChar(char *s)
{
  if(strstr(VALID_SYM, s) == NULL)
    return 0;		// not a valid symbol
  else
    return 1;
}

/*--------------------------------------------------------------------------------------
  Function:	alignout()

  Description:	Output the alignment
  aout - name of the output file
  --------------------------------------------------------------------------------------*/
// void alignout(short* align,char *aout, profile *pf1, profile *pf2,int totalscore) {
void alignout(alignment* align, char *aout, profile *pf1, profile *pf2,
	      int totalscore, int maxsep, int gappen, bool insert)
{
  ostream * out;
  ofstream outfile;

  short	i,j,last,next, nGaps1, nGaps2;
  char	**line1; //[pf1->nSeq][(pf1->numofbases)+(pf2->numofbases)+100];
  char	**line2; //[pf2->nSeq][(pf1->numofbases)+(pf2->numofbases)+100];
  char	*line3; //[(pf1->numofbases)+(pf2->numofbases)+100];
  short int	nBases1,nBases2;

  cout << "Creating alignment file " << aout << "...";

  nBases1 = pf1->numofbases;
  nBases2 = pf2->numofbases;

  // allocate memory for lines
  line1 = new char *[pf1->nSeq];
  for( i=0; i<pf1->nSeq; i++) {
    line1[i] = new char[nBases1 + nBases2 + 2*SEQNAME_LN];
  }

  line2 = new char *[pf2->nSeq];
  for( i=0; i<pf2->nSeq; i++) {
    line2[i] = new char[nBases1 + nBases2 + 2*SEQNAME_LN];
  }

  line3 = new char[nBases1 + nBases2 + 2*SEQNAME_LN];

  // initial lines
  init_lines(line1, pf1);
  init_lines(line2, pf2);
  strcpy(line3, COL_MARKUP);
  strcat(line3, " ");
  strcat(line3, SS_MARKER);
  while( strlen(line3) < SEQNAME_LN ) {
    strcat(line3, " ");
  }

  // setup lines for output
  last = 0;
  nGaps1 = nGaps2 = 0;
  for (i=1; i<=nBases1; i++) {
    if (last==nBases2) {
      //nothing more can be put down in line 2 for the second sequence
      append_nucs(line1, pf1, i);
      append_str(line2, pf2, "-");
      strcat(line3, ".");
      nGaps2++;
    }
    else if (align[i].pos > 0) {
      //this nucleotide is aligned to something
      while (align[i].pos > last+1) {
	//need to lay nucleotides down in second sequence
	append_str(line1, pf1, "-");
	last++;
	append_nucs(line2, pf2, last);
	strcat(line3, ".");
	nGaps1++;
      }
      //now put i in alignment with last+1
      append_nucs(line1, pf1, i);
      last++;
      append_nucs(line2, pf2, last);
      strcat(line3, align[i].indicator);
    }
    else {//align[i]=0
      next=0;
      for (j=i+1;(j<=nBases1)&&(next==0);j++) next = align[j].pos;
      if (next==last+1) {
	//no nucs from seq2 can be put down next to seq1
	append_nucs(line1, pf1, i);
	append_str(line2, pf2, "-");
	strcat(line3, ".");
	nGaps2++;
      }
      else {
	append_nucs(line1, pf1, i);
	last++;
	append_nucs(line2, pf2, last);
	strcat(line3, ".");
      }
    }

  }		// end of for(i)

  //put down any nucs from seq2 that did not go down:
  for (i=last+1;i<=nBases2;i++) {
    append_str(line1, pf1, "-");
    append_nucs(line2, pf2, i);
    strcat(line3, ".");
    nGaps1++;
  }

  outfile.open( aout );

  if ( outfile.is_open() ) {

    out = &outfile;

  } else {

      //Issue warning and set out to point to stdout

      cerr << "Could not create alignment output file, file is not open, " << aout << "\n";
      cerr << "Alignment printed to standard output instead. \n";
      cout << "\nFinal alignment:\n";
      out = &cout;

  }
  
  // set profile id to be the basename of the output file
  char	*pf_id = new char[ID_LN_MAX];
  strcpy(pf_id, aout);
  if( strrchr(pf_id, '.') != NULL ) {
    *(strrchr(pf_id, '.')) = '\0';
  }
  if( strrchr(pf_id, '/') != NULL ) {
    pf_id = strrchr(pf_id, '/')+1;
  }

  //output headers:
  *out << STOCKHOLM_HEADER << "\n";
  *out << FILE_MARKUP << " " << COMMENT_MARKER << " ";
  *out << "Profile-Dynalign alignment of " ;
  *out << pf1->id << " (" << pf1->nSeq << " sequences) and ";
  *out << pf2->id << " (" << pf2->nSeq << " sequences)\n";
  *out << FILE_MARKUP << " " << ID_MARKER << " "<< pf_id << "\n";
  *out << FILE_MARKUP << " " << NSEQ_MARKER << " "<< pf1->nSeq + pf2->nSeq << "\n";

  //output energy values and parameters:

  char number[2*NUMLEN];
  int totalgaps;
  gcvt((float (totalscore))/10, 6, number);
  *out << FILE_MARKUP << " " << TOTAL_E_MARKER << " " << number << " kcal/mol \n";

  *out << FILE_MARKUP << " " << MAXSEP_MARKER << " " << maxsep << "\n";
  if (insert)
    *out << FILE_MARKUP << " " << BP_MARKER << "\n";

  gcvt((float (gappen))/10, 6, number);
  *out << FILE_MARKUP << " " << GP_MARKER << " " << number << " kcal/mol \n";

  totalgaps = nGaps1*(pf1->nSeq) + nGaps2*(pf2->nSeq);
  *out << FILE_MARKUP << " " << TOTAL_GAP_MARKER << " " << totalgaps << "\n";

  //calculate average energy once gaps removed:
  gcvt((
	(float (totalscore - gappen*totalgaps)/(pf1->nSeq + pf2->nSeq))/10
	), 6, number );
  *out << FILE_MARKUP << " " << AVG_E_MARKER << " " << number << " kcal/mol \n";

  //output each sequence line:
  for( i=0; i<pf1->nSeq; i++) {
    *out << line1[i] << "\n";
    //output any sequence-specific mark-up here
  }

  for( i=0; i<pf2->nSeq; i++) {
    *out << line2[i] << "\n";
    //output any sequence-specific mark-up here
  }

  //output alignment mark-up line
  *out << line3 << "\n";
  *out << STOCKHOLM_END << "\n";	// Put down 'end of profile' indicator

  if (outfile.is_open()) {
    outfile.close();
  }

  // free memory of lines
  for( i=0; i< pf1->nSeq; i++) {
    ASSERT( line1[i] != NULL);
    delete[] line1[i];
  }
  ASSERT( line1 != NULL);
  delete[] line1;

  for( i=0; i<pf2->nSeq; i++) {
    ASSERT( line2[i] != NULL);
    delete[] line2[i];
  }
  ASSERT( line2 != NULL);
  delete[] line2;

  ASSERT( line3 != NULL);
  delete[] line3;

  ASSERT( pf_id != NULL);
  delete[] pf_id;

  cout << "done\n";
}

/*-----------------------------------------------------------------------------
  Function:			append_nucs()

  Description: Append a nucleotide character to a line in the alignment. The
  nucleotide is specified by its 'pos' in the profile. Repeat this for all
  sequences in the profile, each sequence according to a line.
  -----------------------------------------------------------------------------*/
void	append_nucs(char **line, profile *pf, int pos)
{
  ASSERT( line != NULL && pf != NULL && pos <= pf->numofbases);

  int	i;
  int	length;

  length = strlen(line[0]);

  ASSERT( length <= 2*pf->numofbases+100 );

  for( i=0; i<pf->nSeq; i++) {
    line[i][length+1] = '\0';
    line[i][length] = (pf->ct+i)->nucs[pos];
  }
}


/*-----------------------------------------------------------------------------
  Function:			append_str()

  Description: Append a string to each line of a profile in the alignment.
  -----------------------------------------------------------------------------*/
void	append_str(char **line, profile *pf, char *s)
{
  ASSERT( line != NULL && s != NULL );

  int	i;

  //ASSERT( strlen(*line) <= (unsigned short int)(2*pf->numofbases+100)  );

  for( i=0; i<pf->nSeq; i++) {
    strcat(line[i], s);
  }
}


/*-----------------------------------------------------------------------------
  Function:			init_lines()

  Description: Initialize the lines for alignment output
  -----------------------------------------------------------------------------*/
void init_lines(char **line, profile *pf)
{
  ASSERT( line != NULL && pf != NULL );

  int i;

  for( i=0; i<pf->nSeq; i++) {
    if( strlen((pf->ct+i)->ctlabel[1]) > SEQNAME_LN ) {
      HDException ex(FILE_ERR, "Sequence name too long");
      throw ex;
    }

    strcpy(line[i], (pf->ct+i)->ctlabel[1]);

    while( strlen(line[i]) < SEQNAME_LN ) {
      strcat(line[i], " ");
    }
  }
}

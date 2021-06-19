/*                               -*- Mode: C -*- 
 * platform.cpp ---
 * Author          : Amelia Bellamy-Royds
 * Created         : Summer 2006
 * Last Modified By: Marcel Turcotte
 * Last Modified On: Mon Dec 11 11:09:26 2006
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

#if !defined (platform)
#define platfrom

//SGI platform code for RNAstructure
//# define bool int
//# define true 1
//# define false 0
//# define TRUE 1
//# define FALSE 0
//# define try /* not supported on SGI */
//# define catch /* not supported on SGI */
//# define (xalloc) /* not defined on SGI */
//# define floor ffloor
//# define floor /*floor*/
#define sgifix strcat(line," ");
#define binary in


#include <math.h>


int min(int i, int j) {

	if (i>j) return j;

	return i;

}

int max (int i,int j) {

	if (i<j) return j;

	return i;

}

int pow10(int i) {
	int j,n;

	if (i==0) return 1;
	j = 1;

	for (n=1;n<=i;n++) {
		j = 10*j;
	}

	return j;
}



void itoa (int x,char *ch,int i)
{
	float y;

	y = float (x);
	gcvt(y,6,ch);

}




#endif

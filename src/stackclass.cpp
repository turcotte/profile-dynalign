/*                               -*- Mode: C -*- 
 * stackclass.cpp ---
 * Author          : Amelia Bellamy-Royds
 * Created         : Summer 2006
 * Last Modified By: Marcel Turcotte
 * Last Modified On: Mon Dec 11 11:10:32 2006
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

/* Definition and methods for class "stackclass" */

class stackclass
{
  short size,**stack,max;
  energy_t *estack;

	void allocate_stack() {
		short i;

		stack=new short int *[max];
		estack = new energy_t [max];
		for (i=0;i<max;i++) stack[i] = new short int [4];
	}

public:

	// Constructor
	stackclass(short int stacksize = 50) {
		max = stacksize;
		size = 0;
		allocate_stack();
	}

	bool pull(short int *i,short int *j, short int *k_,
		short int *l_, energy_t *energy) {

		if (size==0) return false;
		else {
			size--;
			*i = stack[size][0];
			*j = stack[size][1];
			*k_ = stack[size][2];
			*l_ = stack[size][3];
			*energy = estack[size];
			return true;
		}
	}		// end of pull

	void push(short int i,short int j, short int k_,
		short int l_, energy_t energy){

		short x;


		if (size == max) {

			// allocate more space:
			stackclass *temp;
			temp = new stackclass(max);
			for (x=0;x<max;x++) {
				temp->push(stack[x][0],stack[x][1],stack[x][2],stack[x][3],estack[x]);
			}
			delete_array();
			max = 2*max;

			allocate_stack();
			for (x=0;x<(max/2);x++) {
				temp->pull(&stack[x][0],&stack[x][1],&stack[x][2],&stack[x][3],&estack[x]);
			}
		}

		stack[size][0] = i;
		stack[size][1] = j;
		stack[size][2] = k_;
		stack[size][3] = l_;
		estack[size] = energy;
		size++;
	}		// end of push

	short getSize() {
		return size;
	}

	void delete_array() {
		short i;

		for (i=0;i<max;i++) delete[] stack[i];
		delete[] stack;
		delete[] estack;
	}

	// Destructor
	~stackclass() {
		delete_array();
	}

};

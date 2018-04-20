/*
tools.cpp
Eric Ng (100446517)
April 20, 2018
MCSC6040G Final Exam
External file containing the procedures used in the two main files.

Procedures:
 - Uniform RNG
 - Box Muller Transform
 - Rule 73 for Cellular Automaton
*/

#include "stdlib.h"
#include "stdio.h"
#include "math.h"

#define PI 3.14159265359

// RNG for arbitrary interval (uniform distribution)
double rnd()
{
	return (double)rand() / RAND_MAX;
}

// Box-Muller Transform
double BMT()
{
	double r1, r2, bm;
	do {
		r1 = rnd();
	} while (r1 == 0);
	r2 = rnd();
	bm = (double)sqrt(-2 * log(r1))*cos(2 * PI*r2);
	return bm;
}

// Rule 73
int rule73(int nbhd[3])
{
	int state, sum;
	sum = nbhd[0] + nbhd[1] + nbhd[2];
	if (sum == 1) // 100, 010, 001 update
	{
		state = 0;
	}
	else if (nbhd[1] == 1)
	{
		if (nbhd[0] != nbhd[2]) // 110 and 011 update
		{
			state = 1;
		}
		else // 111 update
		{
			state = 0;
		}
	}
	else
	{
		state = -nbhd[0] + 1; // 101 and 000 update
	}
	return state;
}
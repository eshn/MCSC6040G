// A1IVS.c for Assignment 1 Question 1 of MCSC6040
// Eric Ng (100446517)

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#define N 100000

int i;
float rand_num;

int main(){
	// Pointer for I/O
	FILE *outfile;
	outfile = fopen("rand_num.dat", "w");
	for (i = 0; i < N ; i++)
	{
		// Generates random number in [0,2]
		rand_num = 2 * (float)rand() / RAND_MAX;
		// Computes corresponding x*
		rand_num = acos(1 - rand_num);
		// Write to file
		fprintf(outfile, "%f\n", rand_num);
	}
	fclose(outfile);
}

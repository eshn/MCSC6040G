/*
FinalQ2-CA.cpp
Eric Ng (100446517)
April 20, 2018
MCSC6040G Final Exam
Code for Question 2 (Cellular Automata)
Requires tools.cpp and tools.h to run
*/

#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "tools.h" // Header for library containing update rules. See tools.cpp

#define N 100
#define TMAX 100
#define SEED 0 // Flag for seed. 0 - Off, 1 - On
#define Q 0 // Flag for different initialization. 0 - all 1s except central cell (question 2b). 1 - randomized (question 2c)

int main()
{
	if (SEED == 1)
	{
		srand(time(NULL));
	}

	int prev[N];
	int t = 0;
	FILE *outfile;

	// Initialization
	int state[N];
	if (Q == 0) // Initialization for 2b.
	{
		outfile = fopen("Final2b-CA.dat", "w");
		for (int i = 0; i < N; i++)
		{
			if (i == int(N / 2.0))
			{
				state[i] = 0;
			}
			else
			{
				state[i] = 1;
			}
		}
		state[int(N / 2.0)] = 0; // Sets central cell to 0
	}

	if (Q == 1) // Initialization for 2c.
	{
		outfile = fopen("Final2c-CA.dat", "w");
		double r;
		for (int i = 0; i < N; i++)
		{
			r = rnd(); // External procedure for uniform RNG. See rnd() in tools.h and tools.cpp.
			if (r > 0.5)
			{
				state[i] = 1;
			}
			else
			{
				state[i] = 0;
			}
		}
	}
	

	while(t < TMAX)
	{
		// Copies previous state
		memcpy(prev, state, N * sizeof(int));

		if (t == 0) // Writes initial state
		{
			for (int i = 0; i < N; i++)
			{
				fprintf(outfile, "%d ", state[i]);
			}
			fprintf(outfile, "\n");
		}
		// Update steps
		for (int i = 0; i < N; i++)
		{
			// Gets neighborhood values
			int nbhd[3], ind;
			for (int j = 0; j < 3; j++)
			{
				ind = i + j - 1;
				if (ind == -1 || ind == N) // Periodic Boundary COnditions
				{
					if (ind == -1)
					{
						nbhd[j] = prev[N - 1];
					}
					else
					{
						nbhd[j] = prev[0];
					}
				}
				else
				{
					nbhd[j] = prev[ind];
				}
			}
			state[i] = rule73(nbhd); // Calls rule from external procedure. See tools.cpp

			fprintf(outfile, "%d ", state[i]);
		}
		fprintf(outfile, "\n");
		t++;
	}
	fclose(outfile);
}
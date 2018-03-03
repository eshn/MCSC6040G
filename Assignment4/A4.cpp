// A4.cpp for Assignment 4 of MCSC6040
// Eric Ng (100446517)

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"

#define SEED 0 // Flag to enable psuedo-randomized seed
#define N 100000 // Number of particles
#define N_BIN 1000 // Number of bins for plotting

// Procedure to generate random numbers (so I don't have to type so much)
double rnd()
{
	return (double)rand() / RAND_MAX;
}

int main() {
	if (SEED == 1)
	{
		srand(time(NULL));
	}

	double x, KE; // Position and KE variable

	int t = 0, dt = 1;
	double prob;
	double m = 1.0;
	double x_bin = 100; // Bin width. Higher for smoother curve. Smaller for coarser curve.

	double* lost_E = new double[N_BIN]; // Array of variable length to store lost energy by bin


	FILE *outx;
	outx = fopen("A4b.dat", "w");

	// Bin initialization
	for (int i = 0; i < N_BIN; i++)
	{
		lost_E[i] = 0;
	}

	// Simulate particle by particle (since they are fully independent, saves RAM)
	for (int i = 0; i < N; i++)
	{
		if (i % 1000 == 0) // Prints progress
		{
			printf("Particle %d of %d.\n", i, N);
		}
		// Particle initialization
		x = 0.0;
		KE = rnd() + 4.0;
		t = 0;
		// Velocity update
		while (KE > 0.1)
		{
			x += sqrt(2 * KE / m);
			prob = rnd();
			if (prob < 0.05)
			{
				KE -= 0.2;
				lost_E[int(x / x_bin)] += 2.0;
			}
			t += dt;
		}
	}
	
	for (int i = 0; i < N_BIN; i++)
	{
		fprintf(outx, "%f %f\n", i*x_bin, lost_E[i]);
	}
	fclose(outx);
}
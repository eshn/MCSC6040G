#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"


#define N 50 // Population size
#define INF_INIT 5 // Percentage of population initially infected.
#define A 3 // Minimum particle distance for interaction

#define REC_MEAN 10 // Mean of recovery time
#define REC_VAR 10 // Variance of recovery time
#define L 25 // Domain size
#define SEED 0 // Flag for seed
#define TMAX 5000
#define PI 3.14159265359


// RNG Procedure
double rnd()
{
	return (double)rand() / RAND_MAX;
}

// Probability of infection (by distance)
double disease_prob(double x, double a)
{
	double prob;
	prob = exp(1 / (a*a))*exp(-1 / (a*a - x * x))*exp(-4 * x*x / (a*a));
	return prob;
}

// Box-Muller Transform with arbitrary mean and variance
double BMT(double m, double v)
{
	double r1, r2, bm;
	do{
		r1 = rnd();
	} while (r1 == 0);
		r2 = rnd();
	bm = (double)sqrt(-2 * v * log(r1))*cos(2 * PI*r2) + m;
	return bm;
}

int main()
{
	if (SEED == 1)
	{
		srand(time(NULL));
	}

	double prob = 0.01;

	double x[N], y[N], vx[N], vy[N];
	int state[N]; // State of the particle. 1 - Susceptible. 2 - Infected. 3 - Recovered.
	int TtR[N]; // Time to Recover (after infection)
	int t = 0;
	double r, dist;
	
	// Initialization
	for (int i = 0; i < N; i++)
	{
		x[i] = L * rnd();
		y[i] = L * rnd();

		vx[i] = L / 50.0 * (rnd() - 0.5);
		vy[i] = L / 50.0 * (rnd() - 0.5);

		TtR[i] = 0;

		state[i] = 1;
	}

	// Randomly selects infected individuals
	for (int i = 0; i < int(N * INF_INIT / 100.0); i++)
	{
		int ind;
		ind = int(N*rnd());
		state[ind] = 2;
		TtR[ind] = int(BMT((double)REC_MEAN, (double)REC_VAR)); // Recovery time based on Gaussian
	}

	FILE *outfile;
	outfile = fopen("out.xyz", "w");

	while (t < TMAX)
	{
		fprintf(outfile, "%i\n", N);
		fprintf(outfile, "title\n");

		for (int i = 0; i < N; i++)
		{
			// Checks probability of changing velocity
			r = rnd();
			if (r < prob)
			{
				vx[i] = (L / 100.0)*(2 * rnd() - 1);
				vy[i] = (L / 100.0)*(2 * rnd() - 1);
			}

			// Position update
			x[i] += vx[i];
			y[i] += vy[i];

			// Wall Reflection
			if (x[i] < 0)
			{
				x[i] = -x[i];
				vx[i] = -vx[i];
			}
			else if (x[i] > L)
			{
				x[i] = 2.0 * L - x[i];
				vx[i] = -vx[i];
			}
			if (y[i] < 0)
			{
				y[i] = -y[i];
				vy[i] = -vy[i];
			}
			else if (y[i] > L)
			{
				y[i] = 2.0 * L - y[i];
				vy[i] = -vy[i];
			}
			if (state[i] == 2) // Subtracts 1 timestep from recovery time. If 0, change status to recovered.
			{
				if (TtR[i] == 0)
				{
					state[i] = 3;
				}
				else {
					TtR[i] -= 1;
				}
			}
		}
		for (int i = 0; i < N; i++)
		{
			// Disease interaction
			if (state[i] == 2)
			{
				for (int j = 0; j < N; j++)
				{
					if (state[j] == 1)
					{
						dist = sqrt(pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2));
						if (dist < A)
						{
							prob = disease_prob(dist, A);
							r = rnd();
							if (r < prob)
							{
								state[j] = 2;
								TtR[j] = int(BMT((double)REC_MEAN, (double)REC_VAR)); // Adds recovery time
							}
						}
					}
				}
			}
		}
		for (int i = 0; i < N; i++)
		{
			if (state[i] == 1)
			{
				fprintf(outfile, "s%i %f %f 0\n", i, x[i], y[i]);
			}
			else if (state[i] == 2)
			{
				fprintf(outfile, "i%i %f %f 0\n", i, x[i], y[i]);
			}
			else if (state[i] == 3)
			{
				fprintf(outfile, "r%i %f %f 0\n", i, x[i], y[i]);
			}
		}
		
		if (t % 1000 == 0)
		{
			printf("Iteration %d of %d.\n", t, TMAX);
		}
		t++;
	}

	fclose(outfile);
}
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "tools.h"


#define N 500 // Population size
#define INFECTED_INIT 10 // Percentage of population initially infected.
#define A 1 // Minimum particle distance for interaction

//////////// TIME IN SECONDS //////////////
#define REC_MEAN 0.5 // Mean of recovery time
#define SUS_MEAN 0.5 // Mean time to become susceptible

#define L 25 // Domain size
#define SEED 0 // Flag for seed
#define TMAX 50

int main()
{
	if (SEED == 1)
	{
		srand(time(NULL));
	}

	double x[N], y[N], vx[N], vy[N];
	int state[N]; // State of the particle. 1 - Susceptible. 2 - Infected. 3 - Recovered.
	double t = 0.0, dt = 0.05;
	int total_S, total_I, total_R, intsteps = 0;
	double r, dist;
	
	// Initialization
	for (int i = 0; i < N; i++)
	{
		x[i] = L * rnd();
		y[i] = L * rnd();

		vx[i] = 2.0 * rnd() - 1.0;
		vy[i] = 2.0 * rnd() - 1.0;
		
		state[i] = 1;
	}

	// Randomly selects infected individuals
	total_I = 0;
	for (int i = 0; i < int(N * INFECTED_INIT / 100.0); i++)
	{
		int ind;
		ind = int(N*rnd());
		do {
			ind = int(N*rnd());
		} while (state[ind] == 2);
		state[ind] = 2;
		total_I += 1;
	}
	total_S = N - total_I;
	total_R = 0;


	FILE *position, *SIRstat, *contact_i, *ItoR, *RtoS;
	position = fopen("pos.xyz", "w");
	SIRstat = fopen("SIRstat.dat", "w");
	contact_i = fopen("InfectedContact.dat", "w");
	ItoR = fopen("ItoR.dat", "w");
	RtoS = fopen("RtoS.dat", "w");

	fprintf(SIRstat, "%d %d %d\n", total_S, total_I, total_R);
	while (t < TMAX)
	{
		fprintf(position, "%i\n", N);
		fprintf(position, "title\n");

		int contact = 0; // tracks contact rate
		int ItoRcount = 0, RtoScount = 0;
		double prob = 0.01;
		
		// Position Updates only
		for (int i = 0; i < N; i++)
		{
			// Checks probability of changing velocity
			/*r = rnd();
			if (r < prob)
			{
				vx[i] = (L / 100.0)*(2 * rnd() - 1);
				vy[i] = (L / 100.0)*(2 * rnd() - 1);
			}*/

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
		}
		for (int i = 0; i < N; i++)
		{
			/*if (state[i] == 3) // Remaining time of recovered. If 0, change to susceptible.
			{

			}*/
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
							contact += 1;
							prob = disease_prob(dist, A, dt);
							r = rnd();
							if (r < prob)
							{
								state[j] = 2;
								total_I += 1;
								total_S -= 1;
							}
						}
					}
				}
				r = rnd() * 100.0;
				if (r < REC_MEAN)
				{
					state[i] = 3;
					total_I -= 1;
					total_R += 1;
					ItoRcount += 1;
				}
			}
		}
		for (int i = 0; i < N; i++)
		{
			if (state[i] == 1)
			{
				fprintf(position, "s%i %f %f %d\n", i, x[i], y[i], state[i]);
			}
			else if (state[i] == 2)
			{
				fprintf(position, "i%i %f %f %d\n", i, x[i], y[i], state[i]);
			}
			else if (state[i] == 3)
			{
				fprintf(position, "r%i %f %f %d\n", i, x[i], y[i], state[i]);
			}
		}
		fprintf(SIRstat, "%d %d %d\n", total_S, total_I, total_R);
		fprintf(contact_i, "%d %d %d\n", contact, total_S, total_I);
		fprintf(ItoR, "%d %d\n", ItoRcount, total_I);
		fprintf(RtoS, "%d %d\n", RtoScount, total_R);
		intsteps++;
		t += dt;
	}
	fclose(contact_i);
	fclose(ItoR);
	fclose(RtoS);
	fclose(position);
	fclose(SIRstat);
}
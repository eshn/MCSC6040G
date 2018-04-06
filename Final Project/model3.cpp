////// MODEL 3 - SIRS WITH WANING IMMUNITIES ///////

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "tools.h"


#define N 500 // Population size
#define INFECTED_INIT 10 // Percentage of population initially infected.
#define A 1 // Minimum particle distance for interaction

//////////// RECOVERY AND IMMUNITY PROBABILITIES //////////////
#define REC_MEAN 0.5 // Recovery percentage per unit time
#define SUS_MEAN 0.5 // Immunity loss percentage per unit time

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
	int state[N]; // State of the particle. 1 - Susceptible. 2 - Infected. 3 - Recovered. 4 - Waning Immunity.
	double t = 0.0, dt = 0.05;
	int total_S, total_I, total_R, total_W, intsteps = 0;
	double r, dist;

	int total_infected = 0, total_secondary = 0;
	int secondary[N];
	
	// Initialization
	for (int i = 0; i < N; i++)
	{
		x[i] = L * rnd();
		y[i] = L * rnd();

		vx[i] = 2.0 * rnd() - 1.0;
		vy[i] = 2.0 * rnd() - 1.0;

		secondary[i] = 0;
		
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
	total_W = 0;


	FILE *position, *SIRstat, *contact_i, *ItoR, *RtoS;
	position = fopen("pos.xyz", "w");
	SIRstat = fopen("SIRstat.dat", "w");
	contact_i = fopen("InfectedContact.dat", "w");
	ItoR = fopen("ItoR.dat", "w");
	RtoS = fopen("RtoS.dat", "w");

	fprintf(SIRstat, "%d %d %d %d\n", total_S, total_I, total_R, total_W);
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
			// Immunity loss mechanic
			if (state[i] == 3)
			{
				r = rnd() * 100.0;
				if (r < SUS_MEAN)
				{
					state[i] = 4;
					total_R -= 1;
					total_W += 1;
					RtoScount += 1;
				}
			}
			// Disease interaction
			if (state[i] == 2)
			{
				for (int j = 0; j < N; j++)
				{
					if (state[j] == 1 || state[j] == 4)
					{
						dist = sqrt(pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2));
						if (dist < A)
						{
							contact += 1;
							prob = disease_prob(dist, A, dt);
							r = rnd();
							if (r < prob)
							{
								total_I += 1;
								if (state[j] == 1)
								{
									total_S -= 1;
								}
								else if (state[j] == 4)
								{
									total_W -= 1;
								}
								state[j] = 2;
								secondary[i] += 1;
							}
						}
					}
				}
				r = rnd() * 100.0;
				// Recovery Mechanic
				if (r < REC_MEAN)
				{
					state[i] = 3;
					total_I -= 1;
					total_R += 1;
					ItoRcount += 1;
					total_infected += 1;
					total_secondary += secondary[i];
					secondary[i] = 0;
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
		fprintf(SIRstat, "%d %d %d %d\n", total_S, total_I, total_R, total_W);
		fprintf(contact_i, "%d %d %d\n", contact, total_S + total_W, total_I);
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
	printf("Average number of secondary infections: %f \n", (double)total_secondary / (double)total_infected);
	printf("Total number of secondary infections: %d \n", total_secondary);
	printf("Total number of infections: %d \n", total_infected);
}
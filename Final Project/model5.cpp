#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "tools.h"


#define N 100 // Population size
#define INFECTED_INIT 5 // Percentage of population initially infected.
#define A 1 // Minimum particle distance for interaction

//////////// TIME IN SECONDS //////////////
#define REC_MEAN 50 // Mean of recovery time
#define REC_VAR 100 // Variance of recovery time
#define SUS_MEAN 20 // Mean time to become susceptible
#define SUS_VAR 100 // Variance of time to become susceptible

#define L 25 // Domain size
#define SEED 0 // Flag for seed
#define TMAX 2000

int main()
{
	if (SEED == 1)
	{
		srand(time(NULL));
	}

	double x[N], y[N], vx[N], vy[N];
	int state[N]; // State of the particle. 1 - Susceptible. 2 - Infected. 3 - Recovered.
	int movestate[N]; // Move state of the partcile. 1 - Run. 0 - Tumble.
	int TtR[N], TtS[N]; // Time to Recover (after infection), Susceptable (after recovery)
	double t = 0.0, dt = 0.05, vmax = 0.5;
	int total_S, total_I, total_R, intsteps = 0;
	double r, dist, prob;
	double RTinit_prob = 0.5, RtoTprob = 0.25, TtoRprob = 0.25, RTrand;
	
	// Initialization
	for (int i = 0; i < N; i++)
	{
		x[i] = L * rnd();
		y[i] = L * rnd();

		if (rnd() < RTinit_prob) // Run and Tumble state initialization
		{
			movestate[i] = 1;
			vx[i] = 2.0*vmax * rnd() - vmax;
			vy[i] = 2.0*vmax * rnd() - vmax;
		}
		else
		{
			movestate[i] = 0;
			vx[i] = 0.0;
			vy[i] = 0.0;
		}

		TtR[i] = 0;
		TtS[i] = 0;

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
		TtR[ind] = int(BMT(0.0, (double)REC_VAR) / dt); // Recovery time based on Gaussian
		do { // Prevents negative recovery time
			TtR[ind] = int(BMT(0.0, (double)REC_VAR) / dt);
		} while (TtR[ind] < 0);
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
		
		// Position Updates only
		for (int i = 0; i < N; i++)
		{
			
			RTrand = rnd(); // Checks for changes in move state and position updates.
			if (movestate[i] == 0) // Tumble phase
			{
				x[i] += sqrt(2 * dt) * (2.0*vmax * rnd() - vmax); // Brownian
				y[i] += sqrt(2 * dt) * (2.0*vmax * rnd() - vmax);
				if (RTrand < TtoRprob) // Checks for change in move state
				{
					movestate[i] = 1; // Changes state, assign new velocity
					vx[i] = 2.0*vmax * rnd() - vmax;
					vy[i] = 2.0*vmax * rnd() - vmax;
				}
			}
			else if (movestate[i] == 1) // Run phase
			{
				if (state[i] == 2)
				{
					x[i] += sqrt(2 * dt) * (2.0*vmax * rnd() - vmax);	// position update
					y[i] += sqrt(2 * dt) * (2.0*vmax * rnd() - vmax);
				}
				else
				{
					x[i] += vx[i];
					y[i] += vy[i];
				}
				if (RTrand < RtoTprob) // checks for change in move state
				{
					movestate[i] = 0; // changes state
				}
			}
			
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
			if (state[i] == 2) // Remaining time of infection. If 0, change to recovered.
			{
				if (TtR[i] == 0)
				{
					state[i] = 3;
					total_R += 1;
					total_I -= 1;
					TtS[i] = int(BMT((double)SUS_MEAN, (double)SUS_VAR) / dt); // Adds time to become susceptible again
					do {
						TtS[i] = int(BMT((double)SUS_MEAN, (double)SUS_VAR) / dt);
					} while (TtS[i] < 0);
					ItoRcount += 1;
				}
				else {
					TtR[i] -= 1;
				}
			}
			else if (state[i] == 3) // Remaining time of recovered. If 0, change to susceptible.
			{
				if (TtS[i] == 0)
				{
					state[i] = 1;
					total_S += 1;
					total_R -= 1;
					RtoScount += 1;
				}
				else
				{
					TtS[i] -= 1;
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
							contact += 1;
							prob = disease_prob(dist, A, dt);
							r = rnd();
							if (r < prob)
							{
								state[j] = 2;
								total_I += 1;
								total_S -= 1;
								TtR[j] = int(BMT((double)REC_MEAN, (double)REC_VAR) / dt); // Adds recovery time
								do {
									TtR[j] = int(BMT((double)REC_MEAN, (double)REC_VAR) / dt);
								} while (TtR[j] < 0);
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
// A3Brownian.cpp for Assignment 3 of MCSC6040
// Eric Ng (100446517)

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"

#define WALL 0 // Bottom BC. 1 for reflecting. 0 for absorbing (way faster).
#define OBS 1 // Set to 1 if there are obstacles. Set to 0 for no obstacle.
#define SEED 0 // Flag to enable psuedo-randomized seed
#define N 100
#define n 20 // Number of obstacles in each coordinate direction. Update line 14 as well if needed.

int main() {
	if (SEED == 1)
	{
		srand(time(NULL));
	}
	////////////////////// CHANGE KT HERE////////////////////
	double kT = 0.4;
	double dt = 0.001;

	int L = 50;
	double x[N], y[N];
	double x_eff, y_eff;
	double x_obs[20][20], y_obs[20][20];	// positions of obstacles. max value of array is k = ceil(19.5/sigma)
	double offset = 1.25;
	double m = 1.0, zeta = 1.0, eps = 1.0, sigma = 1.0;
	double dx, dy;
	double t = 0.0;
	double r1, r2, xi_x, xi_y;
	double PI = 3.14159265359;
	double r, theta, F, Fx, Fy;
	double flag[N];
	double t_total = 0, FPT = 0;
	double MFPT, varFPT, sdFPT, seFPT;
	int count = N; //Counts how many particles are in the domain

	// Initializing obstacle positions
	if (OBS == 1)
	{
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				x_obs[i][j] = offset + 2.5*i*sigma;
				y_obs[i][j] = offset + 2.5*j*sigma;
			}

		}
	}

	// Initializing particle positions
	for (int i = 0; i < N; i++)
	{
		x[i] = L / 2.0;
		y[i] = L + 5.0;
	}

	// File I/O
	FILE *outfile;
	outfile = fopen("A3Plinko.xyz", "w");

	int intsteps = 0;
	while (count > 0)
		{
		if (intsteps % 10 == 0) {
			fprintf(outfile, "%i\n", N+n*n);
			fprintf(outfile, "title\n");
		}
		// Position update with Brownian motion
		for (int i = 0; i < N; i++) {
			if (y[i] > 0)
			{
				// Random Motion (Box-Muller transform)
				r1 = (double)rand() / RAND_MAX;
				r2 = (double)rand() / RAND_MAX;
				if (r1 == 0)
				{
					r1 = (double)rand() / RAND_MAX;
				}
				xi_x = (double)sqrt(-2 * log(r1))*cos(2 * PI*r2);
				xi_y = (double)sqrt(-2 * log(r1))*sin(2 * PI*r2);

				// External Force
				Fx = 0;
				Fy = -1;

				// Obstacle Interaction
				if (OBS == 1)
				{
					if (y[i] < L)
					{
						// Translation to "lower left" block
						x_eff = x[i] - int(x[i] / 2.5)*2.5;
						y_eff = y[i] - int(y[i] / 2.5)*2.5;
						// Checks the 4 obstacles in the block
						for (int j = 0; j < 2; j++) {
							for (int k = 0; k < 2; k++) {
								r = sqrt(pow((x_eff - x_obs[j][k]), 2) + pow((y_eff - y_obs[j][k]), 2));
								if (r < sigma*pow(2, 1.0 / 6.0)) {
									theta = atan2((y_eff - y_obs[j][k]), (x_eff - x_obs[j][k]));
									F = 48 * eps*(pow(sigma, 12) / pow(r, 13) - pow(sigma, 6) / (2 * pow(r, 7)));
									Fx += F * cos(theta);
									Fy += F * sin(theta);

								}
							}
						}
					}
				}
				
				// Position Update
				dx = Fx * dt / zeta + sqrt(2 * kT*dt / zeta) * xi_x;
				dy = Fy * dt / zeta + sqrt(2 * kT*dt / zeta) * xi_y;

				x[i] = x[i] + dx;
				y[i] = y[i] + dy;

				// Periodic BC
				if (x[i] < 0) {
					x[i] = x[i] + L;
				}
				if (x[i] > L) {
					x[i] = x[i] - L;
				}
				if (y[i] < 0) // Particle hitting the bottom
				{
					if (flag[i] != 0) // Records time to compute MFPT
					{
						FPT += t * t; // Variable for sum of squared FPT
						flag[i] = 0; // Prevents duplicate for reflecting BC
						t_total += t;
						count--;
						printf("%d of %d particles reached the bottom.\n", N-count, N);
					}
					if (WALL == 1) // Reflecting BC at the bottom
					{
						y[i] = -y[i];
					}
				}
			}
			if (intsteps % 10 == 0) {
				fprintf(outfile, "a%i %f %f 0\n", i, x[i], y[i]);
			}
		}
		// Writes obstacle position to file
		if (intsteps % 10 == 0 && OBS == 1) {
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					fprintf(outfile, "a%i %f %f 0\n", 1000 + (j + 1) + i * 20, x_obs[i][j], y_obs[i][j]);
				}
			}
		}
		t += dt;
		intsteps++;
	}
	MFPT = t_total / N; // MFPT
	varFPT = FPT / N - MFPT * MFPT; // Variance of FPT
	sdFPT = sqrt(varFPT); // SD of FPT
	seFPT = sdFPT / sqrt(N - 1); // SE of FPT
	printf("MFPT is %f. Variance is %f. SD is %f. SE is %f.\n", MFPT, varFPT, sdFPT, seFPT);
	printf("Total time elapsed: %f seconds (%d integration steps).\n", t, intsteps);
	fclose(outfile);
}
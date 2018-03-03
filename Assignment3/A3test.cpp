// A1GEV.c for Assignment 1 Question 2 of MCSC6040
// Eric Ng (100446517)

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"

#define WALL 0 // Set to 1 if particles continue to move (elastic collision) when it hits the bottom. Set to 0 if particle stops moving altogether (much quicker computationally)
#define OBS 1 // Flag for obstacle. 1 to include obstacles, 0 to not include obstacles
#define N 100
#define n 20

int main() {

	int L = 50;
	double x[N], y[N];
	double x_obs[20][20], y_obs[20][20];	// positions of obstacles. max value of array is n = ceil(19.5/sigma)
	double offset = 1.25;
	double m = 1.0, zeta = 1.0, eps = 1.0, sigma = 1.0;
	double dt = 0.005, dx, dy;
	double t = 0.0;
	double r1, r2, xi_x, xi_y;
	double PI = 3.14159265359;
	double r, theta, phi, F, Fx, Fy;
	double flag[N];
	double t_total;
	int count; //Counts how many particles are in the domain
	int intsteps = 0;
	

	// File I/O
	FILE *outfile;
	if (OBS == 1) {
		outfile = fopen("A3obs.dat", "w");
	}
	else {
		outfile = fopen("A3no_obs.dat", "w");
	}

	for (double kT = 0.0; kT < 0.6; kT += 0.1)
	{
		printf("Working on kT = %f...\n", kT);
		if (OBS == 1)
			{
				// Initializing obstacle positions
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
				flag[i] = 1; // Re-initializes flag to detect when particles reach the bottom
			}

		t_total = 0;
		count = N;
		t = 0;
		intsteps = 0;
		while (count > 0)
			{
				// Position update with random motion (Brownian)
				for (int i = 0; i < N; i++) {
					if (y[i] > 0)
					{
						// Random Motion (Box-Muller transform)
						do {
							r1 = (double)rand() / RAND_MAX;
							r2 = (double)rand() / RAND_MAX;
						} while (r1 == 0);
						xi_x = sqrt(-2 * log(r1))*cos(2 * PI*r2);
						xi_y = sqrt(-2 * log(r1))*sin(2 * PI*r2);

						Fx = 0;
						Fy = -1;
						if (OBS == 1)
						{
							for (int j = 0; j < n; j++) {
								for (int k = 0; k < n; k++) {
									// Force (coordinate independent) from LJ/WCA Potential
									r = sqrt(pow((x[i] - x_obs[j][k]), 2) + pow((y[i] - y_obs[j][k]), 2));
									if (r < sigma*pow(2, 1.0 / 6.0)) {
										theta = atan2((y[i] - y_obs[j][k]), (x[i] - x_obs[j][k]));
										F = 48 * eps*(pow(sigma, 12) / pow(r, 13) - pow(sigma, 6) / (2 * pow(r, 7)));
										Fx += F * cos(theta);
										Fy += F * sin(theta);

									}
								}
							}
						}
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
								flag[i] = 0;
								t_total += t;
								count--;
								printf("%d of %d particles reached the bottom.\n", N - count, N);
							}
							if (WALL == 1)
							{
								y[i] = -y[i]; // Elastic collision at the bottom wall
							}
						}
					}
				}
				t += dt;
				intsteps++;
			}
		printf("MFPT is %f (%d integration steps).\n\n", t_total / N, intsteps);
		fprintf(outfile, "%f %f %d\n", kT, t_total/N, intsteps);

	}
	
	fclose(outfile);
}
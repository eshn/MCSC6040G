// A1GEV.c for Assignment 1 Question 2 of MCSC6040
// Eric Ng (100446517)

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#define N 50

int main() {
	int i, check_neighbor;
	double r, theta, phi, sigma = 1, eps = 1;
	double dt = 0.01;
	int intsteps = 0;
	int L = 10;
	double ax[N], ay[N], az[N], vx[N], vy[N], vz[N], x[N], y[N], z[N], vhx[N], vhy[N], vhz[N];
	double F = 0, t = 0, tmax = 1000;
	double m = 1.0, kT = 1.0;
	double aF;
	double Psim;
	double Pthe;
	double PI = 3.14159265359;
	double r1, r2, r3, r4;
	double Vp;
	int count_wall, count_particle; // counter for wall and particle collisions

	// Pointer for file I/O
	FILE *outt;
	outt = fopen("PvV.dat", "w");

	// Generates simulated data from L = 10 to 100 in increments of 10.
	for (int k = 0; k < 20; k++) {
		L = (k + 1) * 5;
		t = 0;
		aF = 0;
		intsteps = 0;
		count_wall = 0;
		count_particle = 0;
		// Initialization for Position, Velocity, Acceleration
		for (i = 0; i < N; i++) {
			check_neighbor = 0;
			// Initial Position
			if (i == 0) {
				x[i] = (double)rand() / RAND_MAX * L - L / 2.0;
				y[i] = (double)rand() / RAND_MAX * L - L / 2.0;
				z[i] = (double)rand() / RAND_MAX * L - L / 2.0;
			}
			// Checks for Overlap in Initialization
			else
			{
				while (check_neighbor == 0) {
					x[i] = (double)rand() / RAND_MAX * L - L / 2.0;
					y[i] = (double)rand() / RAND_MAX * L - L / 2.0;
					z[i] = (double)rand() / RAND_MAX * L - L / 2.0;
					for (int j = 0; j < i; j++) {
						r = sqrt(pow((x[i] - x[j]), 2) + pow((y[i] - y[j]), 2) + pow((z[i] - z[j]), 2)); // Calculates distance
						if (r < sigma) { // Re-generates if too close
							break;
						}
						else if (j == i - 1) {
							check_neighbor = 1; // Sets to 1 at end of list to stop while loop
						}
					}
				}
				check_neighbor = 0; // Resets flag
			}
			// Initial Velocity (Box Muller Transform)
			r1 = (double)rand() / RAND_MAX;
			r2 = (double)rand() / RAND_MAX;
			r3 = (double)rand() / RAND_MAX;
			r4 = (double)rand() / RAND_MAX;

			vx[i] = sqrt(kT / (m))*sqrt(-2 * log(r1))*cos(2 * PI*r2);
			vy[i] = sqrt(kT / (m))*sqrt(-2 * log(r1))*sin(2 * PI*r2);
			vz[i] = sqrt(kT / (m))*sqrt(-2 * log(r3))*cos(2 * PI*r4);

			// Initial Acceleration
			ax[i] = 0;
			ay[i] = 0;
			az[i] = 0;
		}

		while (t < tmax) {
			// Velocity Verlet
			// Position Update of every particle
			for (i = 0; i < N; i++) {
				vhx[i] = vx[i] + 0.5*ax[i] * dt;
				x[i] = x[i] + vhx[i] * dt;
				ax[i] = 0;

				vhy[i] = vy[i] + 0.5*ay[i] * dt;
				y[i] = y[i] + vhy[i] * dt;
				ay[i] = 0;

				vhz[i] = vz[i] + 0.5*az[i] * dt;
				z[i] = z[i] + vhz[i] * dt;
				az[i] = 0;
			}
			// Acceleration/second velocity update of every particle
			for (i = 0; i < N; i++) {
				for (int j = i + 1; j < N; j++) {
					// Force (coordinate independent) from LJ/WCA Potential
					r = sqrt(pow((x[i] - x[j]), 2) + pow((y[i] - y[j]), 2) + pow((z[i] - z[j]), 2));
					if (r < sigma*pow(2, 1.0 / 6.0)) {
						theta = atan2((y[i] - y[j]), (x[i] - x[j]));
						phi = acos((z[i] - z[j])/r);
						F = 48 * eps*(pow(sigma, 12) / pow(r, 13) - pow(sigma, 6) / (2 * pow(r, 7)));
						ax[i] += F * sin(phi)*cos(theta)/ m;
						ay[i] += F * sin(phi)*sin(theta) / m;
						az[i] += F * cos(phi) / m;
						ax[j] -= F * sin(phi)*cos(theta) / m;
						ay[j] -= F * sin(phi)*sin(theta) / m;
						az[j] -= F * cos(phi) / m;
						count_particle++;
					}
				}
				// Velocity update
				vx[i] = vhx[i] + 0.5*ax[i] * dt;
				vy[i] = vhy[i] + 0.5*ay[i] * dt;
				vz[i] = vhz[i] + 0.5*az[i] * dt;

				// Checks wall collision
				if (x[i]<-L / 2.0 || x[i]>L / 2.0) {
					vx[i] = -vx[i];
					ax[i] = -ax[i];
					aF += 2 * m*fabs(vx[i]) / dt;
					count_wall++;
				}
				if (y[i]<-L / 2.0 || y[i]>L / 2.0) {
					vy[i] = -vy[i];
					ay[i] = -ay[i];
					aF += 2 * m*fabs(vy[i]) / dt;
					count_wall++;
				}
				if (z[i]<-L / 2.0 || z[i]>L / 2.0) {
					vz[i] = -vz[i];
					az[i] = -az[i];
					aF += 2 * m*fabs(vz[i]) / dt;
					count_wall++;
				}
			}
			t += dt;
			intsteps++;
		}

		aF = aF / (double)intsteps;
		Psim = aF / (6 * L*L); // Surface area of cube

		Vp = (4.0 / 3.0)*PI*(pow(sigma / 2, 3)); // Excluded Volume
		Pthe = N * kT / (L*L*L - 4 * N*Vp); // VDW EOS
		printf("Particle Collision: %d \nCount: %d \nPsim: %f \nPthe: %f \n", count_particle, count_wall, Psim, Pthe);
		fprintf(outt, "%d %d %f %f\n", L*L*L, count_wall, Psim, Pthe);
	}
	fclose(outt);
}

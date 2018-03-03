/*
A2Langevin.cpp for Assignment 2 of MCSC6040G
Eric Ng (100446517)
February 12, 2018
*/
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"

#define N 1000000

int main() {
	// Uncomment for random seed
	//srand(time(NULL));

	////////////////////////// PART D, E, F PARAMETERS //////////////////////////
	double kT = 1;
	double zeta = 1;
	double m = 1;
	////////////////////////////////////////////////////////////////////////////

	int t, tmax = N;
	double dt = 0.01;
	double* x = new double[N];
	double* y = new double[N];
	double* z = new double[N];
	double vx = 1.0, vy = 1.0, vz = 1.0;
	double vhx, vhy, vhz;
	double ax, ay, az;
	double ran1, ran2, ran3, ran4;
	double gr1, gr2, gr3;
	double C1, C2, tau_relax;

	C1 = zeta / m;
	C2 = sqrt(2 * kT*C1 / (m*dt));
	tau_relax = 1 / C1;

	x[0] = 0;
	y[0] = 0;
	z[0] = 0;

	////////////////////////////////// PART A //////////////////////////////////
	/////////////////////////// 3D LANGEVIN DYNAMICS ///////////////////////////
	
	printf("--------------------PART A--------------------\n");
	// Initialization (Box Muller Transform)
	ran1 = (double)rand() / RAND_MAX;
	ran2 = (double)rand() / RAND_MAX;
	ran3 = (double)rand() / RAND_MAX;
	ran4 = (double)rand() / RAND_MAX;

	gr1 = (double)sqrt(-2 * log(ran1))*cos(2 * 3.14159265*ran2);
	gr2 = (double)sqrt(-2 * log(ran1))*sin(2 * 3.14159265*ran2);
	gr3 = (double)sqrt(-2 * log(ran3))*cos(2 * 3.14159265*ran4);

	ax = -C1 * vx + C2 * gr1;
	ay = -C1 * vy + C2 * gr2;
	az = -C1 * vz + C2 * gr3;

	for (t = 1; t < tmax; t++)
	{
		// Velocity Verlet
		vhx = vx + 0.5*ax*dt;
		vhy = vy + 0.5*ay*dt;
		vhz = vz + 0.5*az*dt;

		x[t] = x[t-1] + vhx * dt;
		y[t] = y[t-1] + vhy * dt;
		z[t] = z[t-1] + vhz * dt;

		ran1 = (double)rand() / RAND_MAX;
		ran2 = (double)rand() / RAND_MAX;
		ran3 = (double)rand() / RAND_MAX;
		ran4 = (double)rand() / RAND_MAX;

		// Generates new number if RNG gives 0 because this will give gr1, gr2 = +/- inf
		if (ran1 == 0)
		{
			ran1 = (double)rand() / RAND_MAX;
		}
		else if (ran2 == 0)
		{
			ran2 = (double)rand() / RAND_MAX;
		}
		else if (ran3 == 0)
		{
			ran3 = (double)rand() / RAND_MAX;
		}
		else if (ran4 == 0)
		{
			ran4 = (double)rand() / RAND_MAX;
		}

		gr1 = (double)sqrt(-2 * log(ran1))*cos(2 * 3.14159265*ran2);
		gr2 = (double)sqrt(-2 * log(ran1))*sin(2 * 3.14159265*ran2);
		gr3 = (double)sqrt(-2 * log(ran3))*cos(2 * 3.14159265*ran4);

		ax = -C1 * vx + C2 * gr1;
		ay = -C1 * vy + C2 * gr2;
		az = -C1 * vz + C2 * gr3;

		vx = vhx + 0.5*ax*dt;
		vy = vhy + 0.5*ay*dt;
		vz = vhz + 0.5*az*dt;

	}
	printf("Simulation Done.\n\n");
	
	////////////////////////////////// PART B //////////////////////////////////
	////////////////////////////// CALCULATING MSD /////////////////////////////
	printf("--------------------PART B--------------------\n");
	FILE *MSD;

	int DeltaT = 1;
	int count;
	int start, end; // index markers
	double xend, yend, zend;
	double MSDx, MSDy, MSDz; // MSD based on delta t
	double MSDx_total, MSDy_total, MSDz_total; // Overall MSD
	double Dx, Dy, Dz, D; // Diffusion Coefficient
	int thresh;

	thresh = tmax / 10;
	MSD = fopen("MSD.dat", "w");
	while (DeltaT < thresh)
	{
		count = 0; // counter for number of time steps as DeltaT varies
		MSDx = 0;
		MSDy = 0;
		MSDz = 0;
		for (int i = 0; i < tmax; i+=DeltaT)
		{
			if (i == 0) // First position is always marked
			{
				start = 0;
			}
			else if (i%DeltaT == 0) // Marks end of DeltaT
			{
				end = i;
				MSDx += pow(x[start] - x[end], 2.0);
				MSDy += pow(y[start] - y[end], 2.0);
				MSDz += pow(z[start] - z[end], 2.0);
				start = end; // Sets end of DeltaT to the start of a "new" Delta T
				count++;
			}
		}
		MSDx = MSDx / count;
		MSDy = MSDy / count;
		MSDz = MSDz / count;
		fprintf(MSD, "%d %f %f %f\n", DeltaT, MSDx, MSDy, MSDz);
		DeltaT++;
		if (DeltaT % 500 == 0) // This is just to keep track of progress
		{
			printf("Computing MSD... Progress: DeltaT: %d of %d \n", DeltaT, thresh);
		}
	}
	fclose(MSD);
	printf("MSD finished, moving on...\n\n");
}

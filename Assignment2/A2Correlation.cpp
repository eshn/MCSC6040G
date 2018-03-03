/*
A2Correlation.cpp for Assignment 2 of MCSC6040G
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

	///////////////////////////// PART E PARAMETERS /////////////////////////////
	double gamma = 1;
	double m = 3;
	////////////////////////////////////////////////////////////////////////////

	double zeta;
	double kT = 1;
	int t, tmax = N;
	double dt = 0.01;
	double* vx = new double[N];
	double* vy = new double[N];
	double* vz = new double[N];
	double x = 0, y = 0, z = 0;
	double vhx, vhy, vhz;
	double ax, ay, az;
	double ran1, ran2, ran3, ran4;
	double gr1, gr2, gr3;
	double C1, C2, tau_relax;

	zeta = gamma * m;

	C1 = gamma;
	C2 = sqrt(2 * kT*C1 / (m*dt));
	tau_relax = 1 / C1;

	vx[0] = 1.0;
	vy[0] = 1.0;
	vz[0] = 1.0;

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

	ax = -C1 * vx[0] + C2 * gr1;
	ay = -C1 * vy[0] + C2 * gr2;
	az = -C1 * vz[0] + C2 * gr3;

	for (t = 1; t < tmax; t++)
	{
		// Velocity Verlet
		vhx = vx[t-1] + 0.5*ax*dt;
		vhy = vy[t-1] + 0.5*ay*dt;
		vhz = vz[t-1] + 0.5*az*dt;

		x = x + vhx * dt;
		y = y + vhy * dt;
		z = z + vhz * dt;

		ran1 = (double)rand() / RAND_MAX;
		ran2 = (double)rand() / RAND_MAX;
		ran3 = (double)rand() / RAND_MAX;
		ran4 = (double)rand() / RAND_MAX;

		// Generates new number of RNG gives 0 because this will give gr1, gr2 = +/- inf
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

		ax = -C1 * vx[t-1] + C2 * gr1;
		ay = -C1 * vy[t - 1] + C2 * gr2;
		az = -C1 * vz[t - 1] + C2 * gr3;

		vx[t] = vhx + 0.5*ax*dt;
		vy[t] = vhy + 0.5*ay*dt;
		vz[t] = vhz + 0.5*az*dt;

	}
	printf("Simulation Done.\n\n");

	
	///////////////////////////////// PART B, C /////////////////////////////////
	/////////// VELOCITY AUTO CORRELATION FUNCTION, DIFFUSION CONSTANT /////////
	printf("--------------------PART B--------------------\n");
	FILE *VACF;

	int DeltaT = 1;
	int count;
	int start, end; // index markers
	double VACFx, VACFy, VACFz; // VACF based on delta t
	double Dx, Dy, Dz, D; // Diffusion Coefficient
	double VACFx_total = 0, VACFy_total = 0, VACFz_total = 0; // Total VACF for numerical integration
	double thresh;
	thresh = tmax / 500;
	VACF = fopen("VACF.dat", "w");
	while (DeltaT < thresh)
	{
		count = 0; // counter for number of time steps as DeltaT varies
		VACFx = 0;
		VACFy = 0;
		VACFz = 0;
		for (int i = 0; i < tmax; i++)
		{
			if (i == 0) // First position is always marked
			{
				start = 0;
			}
			else if ((i-start)%DeltaT == 0) // Marks end of DeltaT
			{
				end = i;
				VACFx += vx[start]*vx[end];
				VACFy += vy[start]*vy[end];
				VACFz += vz[start]*vz[end];
				start += 1; // Sets the start position to the next index
				count++;
			}
		}

		VACFx = VACFx / count;
		VACFy = VACFy / count;
		VACFz = VACFz / count;
		fprintf(VACF, "%d %f %f %f\n", DeltaT, VACFx, VACFy, VACFz);
		DeltaT++;
		VACFx_total += VACFx; // Total VACF for numerical integration with DeltaT = 1, dt = 0.01
		VACFy_total += VACFy;
		VACFz_total += VACFz;
		if (DeltaT % 500 == 0) // This is just to keep track of progress
		{
			printf("Computing Autocorrelation Function... Progress: DeltaT: %d of %d \n", DeltaT, tmax/500);
		}
	}
	fclose(VACF);
	printf("Autocorrelation finished, moving on...\n\n");
	
	// Computing diffusion coefficient
	Dx = VACFx_total * dt; // Total VACF for numerical integration with DeltaT = 1, dt = 0.01
	Dy = VACFy_total * dt;
	Dz = VACFz_total * dt;
	printf("The diffusion coefficients are:\n");
	printf("Dx: %f \n", Dx);
	printf("Dy: %f \n", Dy);
	printf("Dz: %f \n\n", Dz);

	////////////////////////////////// PART D //////////////////////////////////
	//////////////////////////////// SEE A2_2D.PY //////////////////////////////
}

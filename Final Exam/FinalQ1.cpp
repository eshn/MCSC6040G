/*
FinalQ1.cpp
Eric Ng (100446517)
April 20, 2018
MCSC6040G Final Exam
Code for Question 1 (MFPT)
Requires tools.cpp and tools.h to run
*/

#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "tools.h"

#define SEED 0 // Flag for seed. 0 - Off, 1 - On
#define N 1000
#define L 20

int main()
{
	double kT = 1.0, m = 1.0, z = 2.0;
	int t_total = 0, dt = 1;
	double x, t_avg;

	for (int i = 0; i < N; i++)
	{
		x = 0.1; // Initial position
		while (x < L) // Absorbing wall condition
		{
			x += sqrt(2 * kT * dt / z) * BMT(); // BMT() is an external procedure for the Box Muller Transform. See tools.cpp
			if (x < 0) // Reflecting BC at x = 0
			{
				x = -x;
			}
			t_total += dt; // Adds to total FPT
		}
	}
	t_avg = (double)t_total / N;
	printf("MFPT is %f time units.\n", t_avg);
	
}
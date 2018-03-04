#include "stdlib.h"
#include "stdio.h"
#include "math.h"

// RNG Procedure
double rnd()
{
	return (double)rand() / RAND_MAX;
}

// Probability of infection (by distance). IMPORTANT, THIS IS A PROBABILITY PER UNIT TIME.
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
	do {
		r1 = rnd();
	} while (r1 == 0);
	r2 = rnd();
	bm = (double)sqrt(-2 * v * log(r1))*cos(2 * PI*r2) + m;
	return bm;
}
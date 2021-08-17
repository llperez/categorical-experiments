#ifndef _MISC_H
#define _MISC_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>


/**
 * Returns the number of milliseconds of a gettimeofday call.
 */
long millis(struct timeval tx);

/**
 * Exp-max normalization for probabilities in log-space.
 * Allocates a new array with malloc() -- make sure to destroy with free().
 */
double *exp_max(double *p, unsigned int size);

/**
 * Generate an array of k unnormalized values in [0,1). 
 * Allocates a new array with malloc() -- make sure to destroy with free().
 */
double *gen_dist(unsigned int k, gsl_rng *r, bool log_p);

#endif // _MISC_H

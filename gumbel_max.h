#ifndef _GUMBEL_MAX_H
#define _GUMBEL_MAX_H

#include <gsl/gsl_rng.h>


/**
 * Structure for keeping precomputed for Gumbel samples.
 */
typedef struct gumbel_max_precomp_t {
  double *samples;
  size_t size;
} gumbel_max_precomp_t;


/**
 * Returns a sample using the Gumbel-max method.
 */
unsigned int gumbel_max(const double *p, unsigned int length, const gsl_rng *r);

/**
 * Returns a sample using the Gumbel-max method (with probabilities in log-scale).
 */
unsigned int gumbel_max_log_scale(const double *p, unsigned int length, const gsl_rng *r);

/**
 * Generates a precomputed Gumbel sample structure.
 */
gumbel_max_precomp_t *gumbel_max_precomp_create(size_t size, gsl_rng *r);

/**
 * Frees the precomputed Gumbel samples structure.
 */
void gumbel_max_precomp_free(gumbel_max_precomp_t *gmp);

/**
 * Returns a sample using the Gumbel-max method, with precomputed samples.
 */
unsigned int gumbel_max_precomp(const double *p, unsigned int length, const gsl_rng *r, const gumbel_max_precomp_t *gmp);

/**
 * Returns a sample using the Gumbel-max method, with precomputed
 * samples and probabilities in the log-scale.
 */
unsigned int gumbel_max_precomp_log_scale(const double *p, unsigned int length, const gsl_rng *r, const gumbel_max_precomp_t *gmp);

#endif // _GUMBEL_MAX_H

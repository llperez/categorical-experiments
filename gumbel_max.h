#ifndef _GUMBEL_MAX_H
#define _GUMBEL_MAX_H

#include <gsl/gsl_rng.h>

typedef struct gumbel_max_precomp_t {
  double *samples;
  size_t size;
} gumbel_max_precomp_t;

unsigned int gumbel_max(const double *p, unsigned int length, const gsl_rng *r);
unsigned int gumbel_max_log_scale(const double *p, unsigned int length, const gsl_rng *r);

gumbel_max_precomp_t *gumbel_max_precomp_create(size_t size, gsl_rng *r);
void gumbel_max_precomp_free(gumbel_max_precomp_t *gmp);

unsigned int gumbel_max_precomp(const double *p, unsigned int length, const gsl_rng *r, const gumbel_max_precomp_t *gmp);
unsigned int gumbel_max_precomp_log_scale(const double *p, unsigned int length, const gsl_rng *r, const gumbel_max_precomp_t *gmp);

#endif // _GUMBEL_MAX_H

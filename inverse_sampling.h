#ifndef _INVERSE_SAMPLING_H
#define _INVERSE_SAMPLING_H

#include <gsl/gsl_rng.h>

void inverse_method_build(const double *p, unsigned int size, double *out);
unsigned int inverse_method_draw(double *s, unsigned int size, const gsl_rng *r);
unsigned int inverse_method(const double *p, unsigned int size, const gsl_rng *r);
unsigned int inverse_method_s_preallocated(const double *p, double *s, unsigned int size, const gsl_rng *r);

#endif // _INVERSE_SAMPLING_H

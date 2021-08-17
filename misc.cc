#include "misc.h"

long millis(struct timeval tx) {
  return (long)tx.tv_sec * 1000 + (long)tx.tv_usec / 1000;
}

double *exp_max(double *p, unsigned int size) {

  // find max
  double max_p = p[0];
  for (unsigned int i=1;i<size;i++) {
    if (p[i] > max_p) {
      max_p = p[i];
    }
  }

  double *out_p = (double *)malloc(sizeof(double) * size);
  for (unsigned int i=0;i<size;i++) {
    out_p[i] = exp(p[i] - max_p);
  }

  return(out_p);
}

// generate an array of k unnormalized probabilities.
double *gen_dist(unsigned int k, gsl_rng *r, bool log_p) {
  double *out = (double *)malloc(sizeof(double) * k);

  for (unsigned int i=0;i<k;i++) {
    out[i] = gsl_rng_uniform(r);
    if (log_p) {
      out[i] = log(out[i]);
    }
  }

  return(out);
}

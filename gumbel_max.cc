#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include "gumbel_max.h"

unsigned int gumbel_max(const double *p, unsigned int length, const gsl_rng *r) {

  int which_max = 0;
  double max_val = log(p[0]) - log(-log(gsl_rng_uniform(r)));
  for (unsigned int i=1;i<length;i++) {

    // gumbel draw
    double c_val = log(p[i]) - log(-log(gsl_rng_uniform(r)));

    // update max
    if (c_val > max_val) {
      which_max = i;
      max_val = c_val;
    }
  }

  return(which_max);
}

unsigned int gumbel_max_log_scale(const double *p, unsigned int length, const gsl_rng *r) {

  int which_max = 0;
  double max_val = p[0] - log(-log(gsl_rng_uniform(r)));
  for (unsigned int i=1;i<length;i++) {

    // gumbel draw
    double c_val = p[i]  - log(-log(gsl_rng_uniform(r)));

    // update max
    if (c_val > max_val) {
      which_max = i;
      max_val = c_val;
    }
  }

  return(which_max);
}

gumbel_max_precomp_t *gumbel_max_precomp_create(size_t size, gsl_rng *r) {

  // allocate structure
  gumbel_max_precomp_t *out = (gumbel_max_precomp_t *)malloc(sizeof(gumbel_max_precomp_t));
  out->samples = (double *)malloc(sizeof(double) * size);
  out->size = size;


  // add the gumbel draws
  for (size_t i=0;i<size;i++) {
    out->samples[i] = log(-log(gsl_rng_uniform(r)));
  }

  return(out);
}

void gumbel_max_precomp_free(gumbel_max_precomp_t *gmp) {
  free(gmp->samples);
  free(gmp);
}

unsigned int gumbel_max_precomp(const double *p, unsigned int length, const gsl_rng *r, const gumbel_max_precomp_t *gmp) {


  // pick the starting position in the Gumbel-max structure at random
  int s_pos = gsl_rng_uniform_int(r, gmp->size);
  int which_max = 0;
  double max_val = log(p[0]) - gmp->samples[s_pos++];

  for (unsigned int i=1;i<length;i++) {

    // get next random Gumbel draw and move forward
    double gum_val = gmp->samples[(s_pos % gmp->size)];
    double c_val = log(p[i]) - gum_val;
    s_pos++;

    // update max
    if (c_val > max_val) {
      which_max = i;
      max_val = c_val;
    }
  }

  return(which_max);
}


unsigned int gumbel_max_precomp_log_scale(const double *p, unsigned int length, const gsl_rng *r, const gumbel_max_precomp_t *gmp) {

  // pick the starting position in the Gumbel-max structure at random
  int s_pos = gsl_rng_uniform_int(r, gmp->size);
  int which_max = 0;
  double max_val = p[0] - gmp->samples[s_pos++];

  for (unsigned int i=1;i<length;i++) {

    // get the next random Gumbel draw and move forward
    double gum_val = gmp->samples[(s_pos % gmp->size)];
    double c_val = p[i] - gum_val;
    s_pos++;


    // update max
    if (c_val > max_val) {
      which_max = i;
      max_val = c_val;
    }
  }

  return(which_max);
}



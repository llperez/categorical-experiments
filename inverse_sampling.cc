#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>

void inverse_method_build(const double *p, unsigned int size, double *out) {

  double csum = 0;
  for (unsigned int i=0;i<size;i++) {
    csum += p[i];
    out[i] = csum;
  }
}


unsigned int inverse_method_draw(double *s, unsigned int size, const gsl_rng *r) {

  double ux = gsl_rng_uniform(r) * s[size-1];
  int bot = 0;
  int top = size - 1;

  while (bot <= top) {
    int pos = bot + ((top - bot) / 2);
    
    if (s[pos] >= ux) {
      top = pos - 1;
    }
    else {
      bot = pos + 1;
    }
  }

  return(bot);
}

unsigned int inverse_method(const double *p, unsigned int size, const gsl_rng *r) {

  double *s = (double *)malloc(sizeof(double) * size); 
  inverse_method_build(p, size, s); 
  unsigned int draw = inverse_method_draw(s, size, r);
  free(s);
  return(draw);
}

unsigned int inverse_method_s_preallocated(const double *p, double *s, unsigned int size, const gsl_rng *r) {
  inverse_method_build(p, size, s); 
  unsigned int draw = inverse_method_draw(s, size, r);
  return(draw);
}


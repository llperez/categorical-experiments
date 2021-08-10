#include "gumbel_max.h"
#include "inverse_method.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// returns the number of milliseconds of a gettimeofday call. 
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


// Gumbel-max (linear scale)
void gum_lin(unsigned int k, int num_trials, int num_repeats, gsl_rng *r) {

  struct timeval time_x;
  struct timeval time_y;
  for (unsigned int i=0;i<num_trials;i++) {

    double *dist = gen_dist(k, r, false);
    
    long sum = 0;
    for (unsigned int j=0;j<num_repeats;j++) {
      gettimeofday(&time_x, NULL);

      gumbel_max(dist, k, r);

      gettimeofday(&time_y, NULL);      
      sum += millis(time_y) - millis(time_x);
    }
    
    printf("%f\n", (double)sum / num_repeats);
    free(dist);
    
  }
}

void gum_lin_precomp(unsigned int k, int num_trials, int num_repeats, gsl_rng *r) {

  struct timeval time_x;
  struct timeval time_y;

  long time_base = millis(time_y) - millis(time_x);
  
  for (unsigned int i=0;i<num_trials;i++) {

    double *dist = gen_dist(k, r, false);
    gumbel_max_precomp_t *gmp = gumbel_max_precomp_create(100000, r);
        
    long sum = 0;
    for (unsigned int j=0;j<num_repeats;j++) {
      gettimeofday(&time_x, NULL);

      gumbel_max_precomp(dist, k, r, gmp);

      gettimeofday(&time_y, NULL);      
      sum += millis(time_y) - millis(time_x);
    }
    
    printf("%f\n", (double)sum / num_repeats);
    free(dist);
    gumbel_max_precomp_free(gmp);
    
  }
}



// Gumbel-max (log scale)
void gum_log(unsigned int k, int num_trials, int num_repeats, gsl_rng *r) {

  struct timeval time_x;
  struct timeval time_y;
  for (unsigned int i=0;i<num_trials;i++) {

    double *dist = gen_dist(k, r, true);
    
    long sum = 0;
    for (unsigned int j=0;j<num_repeats;j++) {
      gettimeofday(&time_x, NULL);

      gumbel_max_log_scale(dist, k, r);

      gettimeofday(&time_y, NULL);      
      sum += millis(time_y) - millis(time_x);
    }
    
    printf("%f\n", (double)sum / num_repeats);
    free(dist);
    
  }
}

void gum_log_precomp(unsigned int k, int num_trials, int num_repeats, gsl_rng *r) {

  struct timeval time_x;
  struct timeval time_y;

  long time_base = millis(time_y) - millis(time_x);
  
  for (unsigned int i=0;i<num_trials;i++) {

    double *dist = gen_dist(k, r, true);
    gumbel_max_precomp_t *gmp = gumbel_max_precomp_create(100000, r);
        
    long sum = 0;
    for (unsigned int j=0;j<num_repeats;j++) {
      gettimeofday(&time_x, NULL);

      gumbel_max_precomp_log_scale(dist, k, r, gmp);

      gettimeofday(&time_y, NULL);      
      sum += millis(time_y) - millis(time_x);
    }
    
    printf("%f\n", (double)sum / num_repeats);
    free(dist);
    gumbel_max_precomp_free(gmp);
    
  }
}

// inverse-sampling (linear scale)
void inv_lin(unsigned int k, int num_trials, int num_repeats, gsl_rng *r) {

  double *s = (double *)malloc(sizeof(double) * k);
  struct timeval time_x;
  struct timeval time_y;
  for (unsigned int i=0;i<num_trials;i++) {

    double *dist = gen_dist(k, r, false);
    long sum = 0;
    for (unsigned int j=0;j<num_repeats;j++) {
      gettimeofday(&time_x, NULL);

      inverse_method(dist, k, r);

      gettimeofday(&time_y, NULL);
      sum += millis(time_y) - millis(time_x);
    }
    
    printf("%f\n", (double)sum / num_repeats);
    free(dist);    
  }

  free(s);
}

// inverse-sampling (log scale)
void inv_log(unsigned int k, int num_trials, int num_repeats, gsl_rng *r) {

  double *s = (double *)malloc(sizeof(double) * k);
  struct timeval time_x;
  struct timeval time_y;
  for (unsigned int i=0;i<num_trials;i++) {

    double *dist = gen_dist(k, r, true);
    long sum = 0;
    for (unsigned int j=0;j<num_repeats;j++) {
      gettimeofday(&time_x, NULL);

      double *dist_exp = exp_max(dist, k);
      inverse_method(dist_exp, k, r);
      free(dist_exp);

      gettimeofday(&time_y, NULL);
      sum += millis(time_y) - millis(time_x);
    }
    
    printf("%f\n", (double)sum / num_repeats);
    free(dist);    
  }

  free(s);
}


// gsl (linear scale)
void gsl_lin(unsigned int k, int num_trials, int num_repeats, gsl_rng *r) {
  struct timeval time_x;
  struct timeval time_y;
  for (unsigned int i=0;i<num_trials;i++) {

    double *dist = gen_dist(k, r, false);
    
    long sum = 0;
    for (unsigned int j=0;j<num_repeats;j++) {
      gettimeofday(&time_x, NULL);

      gsl_ran_discrete_t *grd = gsl_ran_discrete_preproc(k, dist);
      gsl_ran_discrete(r, grd);
      gsl_ran_discrete_free(grd);

      gettimeofday(&time_y, NULL);      
      sum += millis(time_y) - millis(time_x);
    }
    
    printf("%f\n", (double)sum / num_repeats);
    free(dist);    
  }  
}

void gsl_log(unsigned int k, int num_trials, int num_repeats, gsl_rng *r) {
  struct timeval time_x;
  struct timeval time_y;
  for (unsigned int i=0;i<num_trials;i++) {

    double *dist = gen_dist(k, r, false);
    
    long sum = 0;
    for (unsigned int j=0;j<num_repeats;j++) {
      gettimeofday(&time_x, NULL);

      double *dist_exp = exp_max(dist, k);
      gsl_ran_discrete_t *grd = gsl_ran_discrete_preproc(k, dist_exp);
      gsl_ran_discrete(r, grd);
      gsl_ran_discrete_free(grd);
      free(dist_exp);

      gettimeofday(&time_y, NULL);      
      sum += millis(time_y) - millis(time_x);
    }
    
    printf("%f\n", (double)sum / num_repeats);
    free(dist);    
  }  
}

int main(int argc, char **argv) {

  gsl_rng_env_setup();  
  const gsl_rng_type * T = gsl_rng_default;
  gsl_rng * r = gsl_rng_alloc (T);

  /**
  
  double p[20] = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.02, 0.02, 0.02, 0.02, 0.02};

  for (int i=0;i<20;i++) {
    p[i] = log(p[i]);
  }

  gumbel_max_precomp_t *gmp = gumbel_max_precomp_create(1000000, r);
  for (int i=0;i<1000000;i++) {
    int dx = gumbel_max_precomp(p, 20, r, gmp);
    printf("%d\n", dx);
  }

  gumbel_max_precomp_free(gmp);
  
  return(-1);
  **/
  
  if (argc < 5) {
    return(-1);
  }

  int k = atoi(argv[2]);    
  int num_trials = atoi(argv[3]);
  int num_repeats = atoi(argv[4]);
  
  
  if (!strcmp(argv[1], "gum_lin")) {    
    gum_lin(k, num_trials, num_repeats, r);
  }

  if (!strcmp(argv[1], "gum_log")) {    
    gum_log(k, num_trials, num_repeats, r);
  }
  

  if (!strcmp(argv[1], "inv_lin")) {    
    inv_lin(k, num_trials, num_repeats, r);
  }

  if (!strcmp(argv[1], "inv_log")) {    
    inv_log(k, num_trials, num_repeats, r);
  }

  if (!strcmp(argv[1], "gsl_lin")) {    
    gsl_lin(k, num_trials, num_repeats, r);
  }

  if (!strcmp(argv[1], "gsl_log")) {    
    gsl_log(k, num_trials, num_repeats, r);
  }

  if (!strcmp(argv[1], "gum_lin_precomp")) {    
    gum_lin_precomp(k, num_trials, num_repeats, r);
  }
  
  if (!strcmp(argv[1], "gum_log_precomp")) {    
    gum_log_precomp(k, num_trials, num_repeats, r);
  }
  
  return(-1);
}

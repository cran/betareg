# include <stdio.h>
# include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* A naive implementation of the 3_F_2 required for XBeta derivatives */
/* Author: Ioannis Kosmidis */
/* Date: 03/03/2017 */
void h3f2 (double* a, double* b, double* z, int* n, int* maxiter, double* eps, double* out)
{
  double factor;
  double outold;
  int i;
  int j;
  for (j=0; j < *n ; j++) {
    factor = 1;
    out[j] = 1;
    outold = 100000;
    i = 1;
    while ((i <= *maxiter)  && (fabs(outold - out[j]) > *eps)) {
      /* printf("%d\t\t%f\t\t%f\n", i, factor, out[j]); */
      outold = out[j];
      factor = factor * (i - b[j]) * (*z / i);
      out[j] += pow((a[j] / (a[j] + i)), 2) * factor;
      i+=1;
    }
  }
}


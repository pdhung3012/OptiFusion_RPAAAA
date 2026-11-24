/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* correlation.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "correlation.h"


/* Array initialization. */
static
void init_array (int m,
		 int n,
		 DATA_TYPE *float_n,
		 DATA_TYPE POLYBENCH_2D(data,N,M,n,m))
{
  int i, j;

  *float_n = (DATA_TYPE)N;

  for (i = 0; i < N; i++)
    for (j = 0; j < M; j++)
      data[i][j] = (DATA_TYPE)(i*j)/M + i;

}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int m,
		 DATA_TYPE POLYBENCH_2D(corr,M,M,m,m))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("corr");
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++) {
      if ((i * m + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, corr[i][j]);
    }
  POLYBENCH_DUMP_END("corr");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_correlation(int m, int n,
			DATA_TYPE float_n,
			DATA_TYPE POLYBENCH_2D(data,N,M,n,m),
			DATA_TYPE POLYBENCH_2D(corr,M,M,m,m),
			DATA_TYPE POLYBENCH_1D(mean,M,m),
			DATA_TYPE POLYBENCH_1D(stddev,M,m))
{
  int i, j, k;

  DATA_TYPE eps = SCALAR_VAL(0.1);


#pragma scop
#define CORR_TILE 32
  DATA_TYPE sqrt_fn = SQRT_FUN(float_n);

  /* Fuse mean, standard deviation, and column normalization into one
     loop over columns to improve locality.  Each column is fully
     independent, so the outer j loop is parallel. */
#pragma omp parallel for default(none)				\
  shared(m, n, data, mean, stddev, float_n, eps, sqrt_fn)	\
  private(i) schedule(static)
  for (j = 0; j < _PB_M; j++)
    {
      DATA_TYPE mean_j = SCALAR_VAL(0.0);

      for (i = 0; i < _PB_N; i++)
	mean_j += data[i][j];

      mean_j /= float_n;
      mean[j] = mean_j;

      DATA_TYPE stddev_j = SCALAR_VAL(0.0);
      for (i = 0; i < _PB_N; i++)
        {
          DATA_TYPE diff = data[i][j] - mean_j;
          stddev_j += diff * diff;
        }

      stddev_j /= float_n;
      stddev_j = SQRT_FUN(stddev_j);
      stddev_j = stddev_j <= eps ? SCALAR_VAL(1.0) : stddev_j;
      stddev[j] = stddev_j;

      DATA_TYPE norm = sqrt_fn * stddev_j;
      for (i = 0; i < _PB_N; i++)
	data[i][j] = (data[i][j] - mean_j) / norm;
    }

  /* Set the diagonal of the correlation matrix. */
#pragma omp parallel for default(none) shared(m, corr) schedule(static)
  for (i = 0; i < _PB_M; i++)
    corr[i][i] = SCALAR_VAL(1.0);

  /* Calculate the m * m correlation matrix (off-diagonal entries).
     Tile the i/j iteration space to improve reuse of the column
     vectors loaded from data[:,i] and data[:,j]. */
#pragma omp parallel for collapse(2) default(none)		\
  shared(m, n, data, corr) private(i, j, k) schedule(guided)
  for (ii = 0; ii < _PB_M - 1; ii += CORR_TILE)
    for (jj = 0; jj < _PB_M; jj += CORR_TILE)
      for (i = ii; i < ii + CORR_TILE && i < _PB_M - 1; i++)
	for (j = (jj > i + 1 ? jj : i + 1);
	     j < jj + CORR_TILE && j < _PB_M; j++)
	  {
	    DATA_TYPE acc = SCALAR_VAL(0.0);

#pragma omp simd reduction(+:acc)
	    for (k = 0; k < _PB_N; k++)
	      acc += data[k][i] * data[k][j];

	    corr[i][j] = acc;
	    corr[j][i] = acc;
	  }
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;
  int m = M;

  /* Variable declaration/allocation. */
  DATA_TYPE float_n;
  POLYBENCH_2D_ARRAY_DECL(data,DATA_TYPE,N,M,n,m);
  POLYBENCH_2D_ARRAY_DECL(corr,DATA_TYPE,M,M,m,m);
  POLYBENCH_1D_ARRAY_DECL(mean,DATA_TYPE,M,m);
  POLYBENCH_1D_ARRAY_DECL(stddev,DATA_TYPE,M,m);

  /* Initialize array(s). */
  init_array (m, n, &float_n, POLYBENCH_ARRAY(data));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_correlation (m, n, float_n,
		      POLYBENCH_ARRAY(data),
		      POLYBENCH_ARRAY(corr),
		      POLYBENCH_ARRAY(mean),
		      POLYBENCH_ARRAY(stddev));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, POLYBENCH_ARRAY(corr)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(data);
  POLYBENCH_FREE_ARRAY(corr);
  POLYBENCH_FREE_ARRAY(mean);
  POLYBENCH_FREE_ARRAY(stddev);

  return 0;
}
/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* covariance.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "covariance.h"


/* Array initialization. */
static
void init_array (int m, int n,
		 DATA_TYPE *float_n,
		 DATA_TYPE POLYBENCH_2D(data,N,M,n,m))
{
  int i, j;

  *float_n = (DATA_TYPE)n;

  for (i = 0; i < N; i++)
    for (j = 0; j < M; j++)
      data[i][j] = ((DATA_TYPE) i*j) / M;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int m,
		 DATA_TYPE POLYBENCH_2D(cov,M,M,m,m))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("cov");
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++) {
      if ((i * m + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, cov[i][j]);
    }
  POLYBENCH_DUMP_END("cov");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_covariance(int m, int n,
		       DATA_TYPE float_n,
		       DATA_TYPE POLYBENCH_2D(data,N,M,n,m),
		       DATA_TYPE POLYBENCH_2D(cov,M,M,m,m),
		       DATA_TYPE POLYBENCH_1D(mean,M,m))
{
  int i, j, k;

#pragma scop
  /* 1. Column means: parallel over independent columns. */
  #pragma omp parallel for default(none) shared(data, mean, float_n) private(i) schedule(static)
  for (j = 0; j < _PB_M; j++)
    {
      DATA_TYPE sum = SCALAR_VAL(0.0);
      #pragma omp simd reduction(+:sum)
      for (i = 0; i < _PB_N; i++)
        sum += data[i][j];
      mean[j] = sum / float_n;
    }

  /* 2. Center the data matrix: parallel over rows (contiguous access). */
  #pragma omp parallel for default(none) shared(data, mean) private(j) schedule(static)
  for (i = 0; i < _PB_N; i++)
    {
      #pragma omp simd
      for (j = 0; j < _PB_M; j++)
        data[i][j] -= mean[j];
    }

  /* 3. Covariance: interchange so k is outer (cov += data[k][i]*data[k][j]),
        tile the upper-triangle (i,j) loops, and fuse initialization,
        accumulation, normalization and symmetrization per row tile. */
  {
    const int TI = 32;
    const int TJ = 32;
    int ii, jj;

    #pragma omp parallel for default(none) shared(data, cov, float_n, TI, TJ) private(jj, k, i, j) schedule(guided)
    for (ii = 0; ii < _PB_M; ii += TI)
      {
        int i_end = ii + TI < _PB_M ? ii + TI : _PB_M;

        /* Initialize the upper-triangle rows owned by this tile. */
        for (i = ii; i < i_end; i++)
          for (j = i; j < _PB_M; j++)
            cov[i][j] = SCALAR_VAL(0.0);

        /* Tiled accumulation: k outer, (i,j) tile inner. */
        for (jj = ii; jj < _PB_M; jj += TJ)
          {
            int j_end = jj + TJ < _PB_M ? jj + TJ : _PB_M;
            for (k = 0; k < _PB_N; k++)
              {
                for (i = ii; i < i_end; i++)
                  {
                    DATA_TYPE dki = data[k][i];
                    int j_start = jj > i ? jj : i;
                    #pragma omp simd
                    for (j = j_start; j < j_end; j++)
                      cov[i][j] += dki * data[k][j];
                  }
              }
          }

        /* Normalize and symmetrize the rows owned by this tile. */
        for (i = ii; i < i_end; i++)
          {
            #pragma omp simd
            for (j = i; j < _PB_M; j++)
              {
                cov[i][j] /= (float_n - SCALAR_VAL(1.0));
                cov[j][i] = cov[i][j];
              }
          }
      }
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
  POLYBENCH_2D_ARRAY_DECL(cov,DATA_TYPE,M,M,m,m);
  POLYBENCH_1D_ARRAY_DECL(mean,DATA_TYPE,M,m);


  /* Initialize array(s). */
  init_array (m, n, &float_n, POLYBENCH_ARRAY(data));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_covariance (m, n, float_n,
		     POLYBENCH_ARRAY(data),
		     POLYBENCH_ARRAY(cov),
		     POLYBENCH_ARRAY(mean));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, POLYBENCH_ARRAY(cov)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(data);
  POLYBENCH_FREE_ARRAY(cov);
  POLYBENCH_FREE_ARRAY(mean);

  return 0;
}
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
{
  enum { TK = 32, TJ = 64 };

  #pragma omp parallel default(none) shared(mean, data, cov, float_n)
  {
    int ii, jj, kk;

    #pragma omp for schedule(static)
    for (jj = 0; jj < _PB_M; jj++)
      {
        int i;
        DATA_TYPE sum = SCALAR_VAL(0.0);
        for (i = 0; i < _PB_N; i++)
          sum += data[i][jj];
        mean[jj] = sum / float_n;
      }

    #pragma omp for schedule(static)
    for (ii = 0; ii < _PB_N; ii++)
      {
        int j;
        #pragma omp simd
        for (j = 0; j < _PB_M; j++)
          data[ii][j] -= mean[j];
      }

    #pragma omp for schedule(guided)
    for (ii = 0; ii < _PB_M; ii++)
      {
        int j, k;
        int kend, jend;

        for (jj = ii; jj < _PB_M; jj++)
          cov[ii][jj] = SCALAR_VAL(0.0);

        for (kk = 0; kk < _PB_N; kk += TK)
          {
            kend = kk + TK;
            if (kend > _PB_N) kend = _PB_N;
            for (jj = ii; jj < _PB_M; jj += TJ)
              {
                jend = jj + TJ;
                if (jend > _PB_M) jend = _PB_M;
                for (k = kk; k < kend; k++)
                  {
                    DATA_TYPE dki = data[k][ii];
                    #pragma omp simd
                    for (j = jj; j < jend; j++)
                      cov[ii][j] += dki * data[k][j];
                  }
              }
          }

        for (jj = ii; jj < _PB_M; jj++)
          {
            cov[ii][jj] /= (float_n - SCALAR_VAL(1.0));
            cov[jj][ii] = cov[ii][jj];
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
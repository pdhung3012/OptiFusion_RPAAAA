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
  {
    enum { BJ = 32 };
    const DATA_TYPE denom = SQRT_FUN(float_n);

#pragma omp parallel shared(mean, stddev, data, corr, float_n, eps, denom)
    {
      /* Compute mean[j] for each column, tiled over j. */
#pragma omp for schedule(static)
      for (int jj = 0; jj < _PB_M; jj += BJ)
	{
	  int jend = (jj + BJ < _PB_M) ? jj + BJ : _PB_M;
	  for (int j = jj; j < jend; j++)
	    mean[j] = SCALAR_VAL(0.0);
	  for (int i = 0; i < _PB_N; i++)
	    {
#pragma omp simd
	      for (int j = jj; j < jend; j++)
		mean[j] += data[i][j];
	    }
	  for (int j = jj; j < jend; j++)
	    mean[j] /= float_n;
	}

      /* Compute stddev[j] and center/reduce the column vectors,
	 fused and tiled over j. */
#pragma omp for schedule(static)
      for (int jj = 0; jj < _PB_M; jj += BJ)
	{
	  int jend = (jj + BJ < _PB_M) ? jj + BJ : _PB_M;
	  DATA_TYPE sdev[BJ];

	  for (int j = jj; j < jend; j++)
	    sdev[j - jj] = SCALAR_VAL(0.0);
	  for (int i = 0; i < _PB_N; i++)
	    {
#pragma omp simd
	      for (int j = jj; j < jend; j++)
		{
		  DATA_TYPE d = data[i][j] - mean[j];
		  sdev[j - jj] += d * d;
		}
	    }
	  for (int j = jj; j < jend; j++)
	    {
	      DATA_TYPE tmp = sdev[j - jj] / float_n;
	      tmp = SQRT_FUN(tmp);
	      tmp = (tmp <= eps) ? SCALAR_VAL(1.0) : tmp;
	      sdev[j - jj] = tmp;
	      stddev[j] = tmp;
	    }
	  for (int i = 0; i < _PB_N; i++)
	    {
#pragma omp simd
	      for (int j = jj; j < jend; j++)
		{
		  data[i][j] -= mean[j];
		  data[i][j] /= denom * sdev[j - jj];
		}
	    }
	}

      /* Calculate the m * m correlation matrix. */
#pragma omp for schedule(guided)
      for (int i = 0; i < _PB_M - 1; i++)
	{
	  corr[i][i] = SCALAR_VAL(1.0);

	  for (int jj = i + 1; jj < _PB_M; jj += BJ)
	    {
	      int jend = (jj + BJ < _PB_M) ? jj + BJ : _PB_M;
	      for (int j = jj; j < jend; j++)
		corr[i][j] = SCALAR_VAL(0.0);
	    }

	  for (int k = 0; k < _PB_N; k++)
	    {
	      DATA_TYPE x = data[k][i];
	      for (int jj = i + 1; jj < _PB_M; jj += BJ)
		{
		  int jend = (jj + BJ < _PB_M) ? jj + BJ : _PB_M;
#pragma omp simd
		  for (int j = jj; j < jend; j++)
		    corr[i][j] += x * data[k][j];
		}
	    }

	  for (int jj = i + 1; jj < _PB_M; jj += BJ)
	    {
	      int jend = (jj + BJ < _PB_M) ? jj + BJ : _PB_M;
#pragma omp simd
	      for (int j = jj; j < jend; j++)
		corr[j][i] = corr[i][j];
	    }
	}
    }

    corr[_PB_M-1][_PB_M-1] = SCALAR_VAL(1.0);
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
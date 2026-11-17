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
  DATA_TYPE sqrt_n = SQRT_FUN(float_n);
  const int TS = 32;
  DATA_TYPE denom[_PB_M];
  int ii, jj, kk, i_max, j_max, k_max;
  DATA_TYPE dik;

#pragma omp parallel private(i, j, k, ii, jj, kk, i_max, j_max, k_max, dik)
  {
    /* Compute mean and standard deviation for each column. */
#pragma omp for schedule(static)
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
	/* The following in an inelegant but usual way to handle
	   near-zero std. dev. values, which below would cause a zero-
	   divide. */
	stddev_j = stddev_j <= eps ? SCALAR_VAL(1.0) : stddev_j;
	stddev[j] = stddev_j;
	denom[j] = sqrt_n * stddev_j;
      }

    /* Center and reduce the column vectors. */
#pragma omp for schedule(static)
    for (i = 0; i < _PB_N; i++)
      {
#pragma omp simd
	for (j = 0; j < _PB_M; j++)
	  data[i][j] = (data[i][j] - mean[j]) / denom[j];
      }

    /* Calculate the m * m correlation matrix. */
#pragma omp for schedule(dynamic)
    for (ii = 0; ii < _PB_M; ii += TS)
      {
	i_max = ii + TS < _PB_M ? ii + TS : _PB_M;
	for (i = ii; i < i_max; i++)
	  {
	    corr[i][i] = SCALAR_VAL(1.0);
	    for (jj = i + 1; jj < _PB_M; jj += TS)
	      {
		j_max = jj + TS < _PB_M ? jj + TS : _PB_M;
#pragma omp simd
		for (j = jj; j < j_max; j++)
		  corr[i][j] = SCALAR_VAL(0.0);
		for (kk = 0; kk < _PB_N; kk += TS)
		  {
		    k_max = kk + TS < _PB_N ? kk + TS : _PB_N;
		    for (k = kk; k < k_max; k++)
		      {
			dik = data[k][i];
#pragma omp simd
			for (j = jj; j < j_max; j++)
			  corr[i][j] += dik * data[k][j];
		      }
		  }
		for (j = jj; j < j_max; j++)
		  corr[j][i] = corr[i][j];
	      }
	  }
      }
  }

  corr[_PB_M-1][_PB_M-1] = SCALAR_VAL(1.0);
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
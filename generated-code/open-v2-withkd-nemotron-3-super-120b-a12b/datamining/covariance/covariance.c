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
#define TILE_I 32
#define TILE_J 32
#define TILE_K 32

#pragma omp parallel for
for (int j = 0; j < _PB_M; j++) {
    mean[j] = SCALAR_VAL(0.0);
    for (int i = 0; i < _PB_N; i++) {
        mean[j] += data[i][j];
    }
    mean[j] /= float_n;
    for (int i = 0; i < _PB_N; i++) {
        data[i][j] -= mean[j];
    }
}

int num_tile_i = (_PB_M + TILE_I - 1) / TILE_I;
#pragma omp parallel for schedule(guided)
for (int ii_t = 0; ii_t < num_tile_i; ii_t++) {
    int ii = ii_t * TILE_I;
    int i_end = ii + TILE_I < _PB_M ? ii + TILE_I : _PB_M;
    for (int jj = ii; jj < _PB_M; jj += TILE_J) {
        int j_end = jj + TILE_J < _PB_M ? jj + TILE_J : _PB_M;
        DATA_TYPE tile_cov[TILE_I][TILE_J];
        for (int i = 0; i < TILE_I; i++) {
            for (int j = 0; j < TILE_J; j++) {
                tile_cov[i][j] = SCALAR_VAL(0.0);
            }
        }
        for (int kk = 0; kk < _PB_N; kk += TILE_K) {
            int k_end = kk + TILE_K < _PB_N ? kk + TILE_K : _PB_N;
            for (int i = ii; i < i_end; i++) {
                int j_start = (jj > i) ? jj : i;
                for (int j = j_start; j < j_end; j++) {
                    for (int k = kk; k < k_end; k++) {
                        tile_cov[i-ii][j-jj] += data[k][i] * data[k][j];
                    }
                }
            }
        }
        for (int i = ii; i < i_end; i++) {
            for (int j = jj; j < j_end; j++) {
                if (j < i) continue;
                DATA_TYPE val = tile_cov[i-ii][j-jj] / (float_n - SCALAR_VAL(1.0));
                cov[i][j] = val;
                cov[j][i] = val;
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

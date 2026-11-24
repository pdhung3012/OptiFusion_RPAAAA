/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* jacobi-2d.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "jacobi-2d.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		 DATA_TYPE POLYBENCH_2D(B,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
	A[i][j] = ((DATA_TYPE) i*(j+2) + 2) / n;
	B[i][j] = ((DATA_TYPE) i*(j+3) + 3) / n;
      }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("A");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      if ((i * n + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j]);
    }
  POLYBENCH_DUMP_END("A");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_jacobi_2d(int tsteps,
			    int n,
			    DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
			    DATA_TYPE POLYBENCH_2D(B,N,N,n,n))
{
  int t, i, j;

#pragma scop
#define TILE_SIZE 64
  for (t = 0; t < _PB_TSTEPS; t++)
    {
      // First stencil: update B from A
      #pragma omp parallel for collapse(2) private(ii, jj, i, j, i_min, i_max, j_min, j_max, i_start_load, i_end_load, j_start_load, j_end_load)
      for (ii = 1; ii < _PB_N - 1; ii += TILE_SIZE)
        for (jj = 1; jj < _PB_N - 1; jj += TILE_SIZE)
          {
            // Determine the tile boundaries for the computation (without halo)
            int i_min = ii;
            int i_max = (ii + TILE_SIZE - 1) < (_PB_N-2) ? (ii + TILE_SIZE - 1) : (_PB_N-2);
            int j_min = jj;
            int j_max = (jj + TILE_SIZE - 1) < (_PB_N-2) ? (jj + TILE_SIZE - 1) : (_PB_N-2);

            // Extend the tile by one in each direction for the stencil
            int i_start_load = (i_min > 1) ? i_min-1 : 1;
            int i_end_load = (i_max < _PB_N-2) ? i_max+1 : _PB_N-2;
            int j_start_load = (j_min > 1) ? j_min-1 : 1;
            int j_end_load = (j_max < _PB_N-2) ? j_max+1 : _PB_N-2;

            // Iterate over the extended tile
            for (i = i_start_load; i <= i_end_load; i++) {
              for (j = j_start_load; j <= j_end_load; j++) {
                // Only compute the inner points (the ones that are in the original tile without halo)
                if (i >= i_min && i <= i_max && j >= j_min && j <= j_max) {
                  B[i][j] = SCALAR_VAL(0.2) * (A[i][j] + A[i][j-1] + A[i][j+1] + A[i+1][j] + A[i-1][j]);
                }
              }
            }
          }

      // Second stencil: update A from B
      #pragma omp parallel for collapse(2) private(ii, jj, i, j, i_min, i_max, j_min, j_max, i_start_load, i_end_load, j_start_load, j_end_load)
      for (ii = 1; ii < _PB_N - 1; ii += TILE_SIZE)
        for (jj = 1; jj < _PB_N - 1; jj += TILE_SIZE)
          {
            int i_min = ii;
            int i_max = (ii + TILE_SIZE - 1) < (_PB_N-2) ? (ii + TILE_SIZE - 1) : (_PB_N-2);
            int j_min = jj;
            int j_max = (jj + TILE_SIZE - 1) < (_PB_N-2) ? (jj + TILE_SIZE - 1) : (_PB_N-2);

            int i_start_load = (i_min > 1) ? i_min-1 : 1;
            int i_end_load = (i_max < _PB_N-2) ? i_max+1 : _PB_N-2;
            int j_start_load = (j_min > 1) ? j_min-1 : 1;
            int j_end_load = (j_max < _PB_N-2) ? j_max+1 : _PB_N-2;

            for (i = i_start_load; i <= i_end_load; i++) {
              for (j = j_start_load; j <= j_end_load; j++) {
                if (i >= i_min && i <= i_max && j >= j_min && j <= j_max) {
                  A[i][j] = SCALAR_VAL(0.2) * (B[i][j] + B[i][j-1] + B[i][j+1] + B[i+1][j] + B[i-1][j]);
                }
              }
            }
          }
    }
#undef TILE_SIZE
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;
  int tsteps = TSTEPS;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(B, DATA_TYPE, N, N, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_jacobi_2d(tsteps, n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(A)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);

  return 0;
}

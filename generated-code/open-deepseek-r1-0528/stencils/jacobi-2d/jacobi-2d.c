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
#define TILE_I 64
#define TILE_J 64
  for (t = 0; t < _PB_TSTEPS; t++)
    {
      #pragma omp parallel for collapse(2) schedule(static) private(i, j)
      for (i = 1; i < _PB_N - 1; i += TILE_I)
        for (j = 1; j < _PB_N - 1; j += TILE_J)
          for (int ii = i; ii < i + TILE_I && ii < _PB_N - 1; ii++) {
            #pragma omp simd
            for (int jj = j; jj < j + TILE_J && jj < _PB_N - 1; jj++)
              B[ii][jj] = SCALAR_VAL(0.2) * (A[ii][jj] + A[ii][jj-1] + A[ii][jj+1] + A[ii+1][jj] + A[ii-1][jj]);
          }

      #pragma omp parallel for collapse(2) schedule(static) private(i, j)
      for (i = 1; i < _PB_N - 1; i += TILE_I)
        for (j = 1; j < _PB_N - 1; j += TILE_J)
          for (int ii = i; ii < i + TILE_I && ii < _PB_N - 1; ii++) {
            #pragma omp simd
            for (int jj = j; jj < j + TILE_J && jj < _PB_N - 1; jj++)
              A[ii][jj] = SCALAR_VAL(0.2) * (B[ii][jj] + B[ii][jj-1] + B[ii][jj+1] + B[ii+1][jj] + B[ii-1][jj]);
          }
    }
#undef TILE_I
#undef TILE_J
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
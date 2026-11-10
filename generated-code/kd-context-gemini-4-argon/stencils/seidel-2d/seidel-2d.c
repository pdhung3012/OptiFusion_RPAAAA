/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* seidel-2d.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "seidel-2d.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      A[i][j] = ((DATA_TYPE) i*(j+2) + 2) / n;
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
void kernel_seidel_2d(int tsteps,
		      int n,
		      DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int t, i, j;

#pragma scop
  if (_PB_TSTEPS >= 1 && _PB_N >= 3) {
    const int B = 32;
    int nt = _PB_TSTEPS;
    int ni = (_PB_N - 2 + B - 1) / B;
    int nj = (_PB_N - 2 + B - 1) / B;
    int max_w = (nt - 1) + (ni - 1) + (nj - 1);

    #pragma omp parallel default(none) shared(A, tsteps, n, nt, ni, nj, max_w, B) private(t, i, j)
    {
      for (int w = 0; w <= max_w; w++) {
        int tb_min = 0;
        if (w - (ni - 1) - (nj - 1) > tb_min) tb_min = w - (ni - 1) - (nj - 1);
        int tb_max = nt - 1;
        if (w < tb_max) tb_max = w;

        #pragma omp for schedule(static)
        for (int tb = tb_min; tb <= tb_max; tb++) {
          int rem = w - tb;
          int ib_min = (rem - (nj - 1) > 0) ? (rem - (nj - 1)) : 0;
          int ib_max = (rem < ni - 1) ? rem : (ni - 1);

          for (int ib = ib_min; ib <= ib_max; ib++) {
            int jb = rem - ib;
            int i_start = 1 + ib * B;
            int i_end = i_start + B - 1;
            if (i_end > _PB_N - 2) i_end = _PB_N - 2;

            int j_start = 1 + jb * B;
            int j_end = j_start + B - 1;
            if (j_end > _PB_N - 2) j_end = _PB_N - 2;

            t = tb;
            for (int ii = i_start; ii <= i_end; ii++) {
              i = ii + t;
              for (int jj = j_start; jj <= j_end; jj++) {
                j = jj + 2 * t + ii;
                if (i <= _PB_N - 2 + t && j <= _PB_N - 2 + 2 * t + ii) {
                  int orig_i = ii;
                  int orig_j = jj;
                  A[orig_i][orig_j] = (A[orig_i-1][orig_j-1] + A[orig_i-1][orig_j] + A[orig_i-1][orig_j+1]
                                     + A[orig_i][orig_j-1] + A[orig_i][orig_j] + A[orig_i][orig_j+1]
                                     + A[orig_i+1][orig_j-1] + A[orig_i+1][orig_j] + A[orig_i+1][orig_j+1]) / SCALAR_VAL(9.0);
                }
              }
            }
          }
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
  int tsteps = TSTEPS;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(A));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_seidel_2d (tsteps, n, POLYBENCH_ARRAY(A));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(A)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);

  return 0;
}
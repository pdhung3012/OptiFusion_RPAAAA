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
{
  const int B = 32;
  int L = _PB_N - 2;

  if (L > 0) {
    int nt_i = (L + B - 1) / B;
    int nt_j = (L + B - 1) / B;
    int s_max = 2 * (nt_i - 1) + (nt_j - 1);
    int s, ti, tj;
    int ti_min, ti_max;
    int i_start, i_end, j_start, j_end;

#pragma omp parallel private(t, s, ti, tj, i, j, ti_min, ti_max, i_start, i_end, j_start, j_end)
    {
      for (t = 0; t <= _PB_TSTEPS - 1; t++) {
        for (s = 0; s <= s_max; s++) {
          ti_min = s - (nt_j - 1);
          if (ti_min < 0)
            ti_min = 0;
          else
            ti_min = (ti_min + 1) / 2;

          ti_max = s / 2;
          if (ti_max > nt_i - 1)
            ti_max = nt_i - 1;

#pragma omp for schedule(dynamic, 1)
          for (ti = ti_min; ti <= ti_max; ti++) {
            tj = s - 2 * ti;

            i_start = 1 + ti * B;
            i_end = i_start + B - 1;
            if (i_end > _PB_N - 2)
              i_end = _PB_N - 2;

            j_start = 1 + tj * B;
            j_end = j_start + B - 1;
            if (j_end > _PB_N - 2)
              j_end = _PB_N - 2;

            for (i = i_start; i <= i_end; i++) {
              for (j = j_start; j <= j_end; j++) {
                A[i][j] = (A[i-1][j-1] + A[i-1][j] + A[i-1][j+1]
                           + A[i][j-1] + A[i][j] + A[i][j+1]
                           + A[i+1][j-1] + A[i+1][j] + A[i+1][j+1]) / SCALAR_VAL(9.0);
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
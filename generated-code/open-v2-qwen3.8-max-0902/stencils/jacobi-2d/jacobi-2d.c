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
  for (t = 0; t < _PB_TSTEPS; t++)
    {
      int ii, jj;
      #pragma omp parallel for collapse(2) schedule(static)
      for (ii = 1; ii < _PB_N - 1; ii += 32)
        for (jj = 1; jj < _PB_N - 1; jj += 32)
          {
            int i_max = (ii + 32 < _PB_N - 1) ? ii + 32 : _PB_N - 1;
            int j_max = (jj + 32 < _PB_N - 1) ? jj + 32 : _PB_N - 1;
            int ti, tj;
            for (ti = ii; ti < i_max; ti++)
              #pragma omp simd
              for (tj = jj; tj < j_max; tj++)
                B[ti][tj] = SCALAR_VAL(0.2) * (A[ti][tj] + A[ti][tj-1] + A[ti][1+tj] + A[1+ti][tj] + A[ti-1][tj]);
          }

      #pragma omp parallel for collapse(2) schedule(static)
      for (ii = 1; ii < _PB_N - 1; ii += 32)
        for (jj = 1; jj < _PB_N - 1; jj += 32)
          {
            int i_max = (ii + 32 < _PB_N - 1) ? ii + 32 : _PB_N - 1;
            int j_max = (jj + 32 < _PB_N - 1) ? jj + 32 : _PB_N - 1;
            int ti, tj;
            for (ti = ii; ti < i_max; ti++)
              #pragma omp simd
              for (tj = jj; tj < j_max; tj++)
                A[ti][tj] = SCALAR_VAL(0.2) * (B[ti][tj] + B[ti][tj-1] + B[ti][1+tj] + B[1+ti][tj] + B[ti-1][tj]);
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
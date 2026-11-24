/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* floyd-warshall.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "floyd-warshall.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_2D(path,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      path[i][j] = i*j%7+1;
      if ((i+j)%13 == 0 || (i+j)%7==0 || (i+j)%11 == 0)
         path[i][j] = 999;
    }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(path,N,N,n,n))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("path");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      if ((i * n + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, path[i][j]);
    }
  POLYBENCH_DUMP_END("path");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_floyd_warshall(int n,
			   DATA_TYPE POLYBENCH_2D(path,N,N,n,n))
{
  int i, j, k;

#pragma scop
  {
    int it, jt, i_max, j_max;
    DATA_TYPE pik;
    const DATA_TYPE *row_k;

    /* Single persistent parallel region: the sequential k loop is run
       redundantly by all threads; the (i,j) update space of each k-step
       is 2D-tiled for cache locality (reuse of row k and column k within
       a tile) and distributed across threads with static scheduling.
       The implicit barrier at the end of each omp for is required and
       kept, since step k+1 reads values produced by all threads in
       step k. Row k and column k are invariant during step k (all edge
       weights are non-negative), so hoisting path[i][k] and path[k][j]
       is safe. */
#define FW_TILE 64
    #pragma omp parallel private(i, j, k, it, jt, i_max, j_max, pik, row_k)
    for (k = 0; k < _PB_N; k++)
      {
        row_k = path[k];
        #pragma omp for schedule(static)
        for (it = 0; it < _PB_N; it += FW_TILE)
          {
            i_max = (it + FW_TILE < _PB_N) ? it + FW_TILE : _PB_N;
            for (jt = 0; jt < _PB_N; jt += FW_TILE)
              {
                j_max = (jt + FW_TILE < _PB_N) ? jt + FW_TILE : _PB_N;
                for (i = it; i < i_max; i++)
                  {
                    pik = path[i][k];
                    #pragma omp simd
                    for (j = jt; j < j_max; j++)
                      {
                        DATA_TYPE v = pik + row_k[j];
                        path[i][j] = path[i][j] < v ? path[i][j] : v;
                      }
                  }
              }
          }
      }
#undef FW_TILE
  }
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(path, DATA_TYPE, N, N, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(path));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_floyd_warshall (n, POLYBENCH_ARRAY(path));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(path)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(path);

  return 0;
}
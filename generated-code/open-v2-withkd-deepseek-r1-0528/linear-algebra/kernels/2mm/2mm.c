/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* 2mm.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "2mm.h"


/* Array initialization. */
static
void init_array(int ni, int nj, int nk, int nl,
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj),
		DATA_TYPE POLYBENCH_2D(C,NJ,NL,nj,nl),
		DATA_TYPE POLYBENCH_2D(D,NI,NL,ni,nl))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nk; j++)
      A[i][j] = (DATA_TYPE) ((i*j+1) % ni) / ni;
  for (i = 0; i < nk; i++)
    for (j = 0; j < nj; j++)
      B[i][j] = (DATA_TYPE) (i*(j+1) % nj) / nj;
  for (i = 0; i < nj; i++)
    for (j = 0; j < nl; j++)
      C[i][j] = (DATA_TYPE) ((i*(j+3)+1) % nl) / nl;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nl; j++)
      D[i][j] = (DATA_TYPE) (i*(j+2) % nk) / nk;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int ni, int nl,
		 DATA_TYPE POLYBENCH_2D(D,NI,NL,ni,nl))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("D");
  for (i = 0; i < ni; i++)
    for (j = 0; j < nl; j++) {
	if ((i * ni + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, D[i][j]);
    }
  POLYBENCH_DUMP_END("D");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_2mm(int ni, int nj, int nk, int nl,
		DATA_TYPE alpha,
		DATA_TYPE beta,
		DATA_TYPE POLYBENCH_2D(tmp,NI,NJ,ni,nj),
		DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj),
		DATA_TYPE POLYBENCH_2D(C,NJ,NL,nj,nl),
		DATA_TYPE POLYBENCH_2D(D,NI,NL,ni,nl))
{
  int i, j, k;

#pragma scop
#define TILE_I 32
#define TILE_J 32
#define TILE_K 32
  #pragma omp parallel default(none) shared(tmp, A, B, C, D, alpha, beta, ni, nj, nk, nl) firstprivate(TILE_I, TILE_J, TILE_K)
  {
    #pragma omp for collapse(2) schedule(static)
    for (int i0 = 0; i0 < ni; i0 += TILE_I)
      for (int j0 = 0; j0 < nj; j0 += TILE_J) {
        for (int i = i0; i < (i0 + TILE_I < ni ? i0 + TILE_I : ni); i++)
          for (int j = j0; j < (j0 + TILE_J < nj ? j0 + TILE_J : nj); j++)
            tmp[i][j] = SCALAR_VAL(0.0);
        
        for (int k0 = 0; k0 < nk; k0 += TILE_K) {
          for (int i = i0; i < (i0 + TILE_I < ni ? i0 + TILE_I : ni); i++) {
            for (int k = k0; k < (k0 + TILE_K < nk ? k0 + TILE_K : nk); k++) {
              DATA_TYPE temp = alpha * A[i][k];
              #pragma omp simd
              for (int j = j0; j < (j0 + TILE_J < nj ? j0 + TILE_J : nj); j++) {
                tmp[i][j] += temp * B[k][j];
              }
            }
          }
        }
      }
    
    #pragma omp for collapse(2) schedule(static)
    for (int i0 = 0; i0 < ni; i0 += TILE_I)
      for (int j0 = 0; j0 < nl; j0 += TILE_J) {
        for (int i = i0; i < (i0 + TILE_I < ni ? i0 + TILE_I : ni); i++)
          for (int j = j0; j < (j0 + TILE_J < nl ? j0 + TILE_J : nl); j++)
            D[i][j] *= beta;
        
        for (int k0 = 0; k0 < nj; k0 += TILE_K) {
          for (int i = i0; i < (i0 + TILE_I < ni ? i0 + TILE_I : ni); i++) {
            for (int k = k0; k < (k0 + TILE_K < nj ? k0 + TILE_K : nj); k++) {
              DATA_TYPE temp = tmp[i][k];
              #pragma omp simd
              for (int j = j0; j < (j0 + TILE_J < nl ? j0 + TILE_J : nl); j++) {
                D[i][j] += temp * C[k][j];
              }
            }
          }
        }
      }
  }
#undef TILE_I
#undef TILE_J
#undef TILE_K
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int ni = NI;
  int nj = NJ;
  int nk = NK;
  int nl = NL;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  DATA_TYPE beta;
  POLYBENCH_2D_ARRAY_DECL(tmp,DATA_TYPE,NI,NJ,ni,nj);
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,NI,NK,ni,nk);
  POLYBENCH_2D_ARRAY_DECL(B,DATA_TYPE,NK,NJ,nk,nj);
  POLYBENCH_2D_ARRAY_DECL(C,DATA_TYPE,NJ,NL,nj,nl);
  POLYBENCH_2D_ARRAY_DECL(D,DATA_TYPE,NI,NL,ni,nl);

  /* Initialize array(s). */
  init_array (ni, nj, nk, nl, &alpha, &beta,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B),
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(D));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_2mm (ni, nj, nk, nl,
	      alpha, beta,
	      POLYBENCH_ARRAY(tmp),
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B),
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(D));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(ni, nl,  POLYBENCH_ARRAY(D)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(tmp);
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(D);

  return 0;
}
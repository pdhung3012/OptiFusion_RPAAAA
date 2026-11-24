/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* 3mm.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "3mm.h"


/* Array initialization. */
static
void init_array(int ni, int nj, int nk, int nl, int nm,
		DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj),
		DATA_TYPE POLYBENCH_2D(C,NJ,NM,nj,nm),
		DATA_TYPE POLYBENCH_2D(D,NM,NL,nm,nl))
{
  int i, j;

  for (i = 0; i < ni; i++)
    for (j = 0; j < nk; j++)
      A[i][j] = (DATA_TYPE) ((i*j+1) % ni) / (5*ni);
  for (i = 0; i < nk; i++)
    for (j = 0; j < nj; j++)
      B[i][j] = (DATA_TYPE) ((i*(j+1)+2) % nj) / (5*nj);
  for (i = 0; i < nj; i++)
    for (j = 0; j < nm; j++)
      C[i][j] = (DATA_TYPE) (i*(j+3) % nl) / (5*nl);
  for (i = 0; i < nm; i++)
    for (j = 0; j < nl; j++)
      D[i][j] = (DATA_TYPE) ((i*(j+2)+2) % nk) / (5*nk);
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int ni, int nl,
		 DATA_TYPE POLYBENCH_2D(G,NI,NL,ni,nl))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("G");
  for (i = 0; i < ni; i++)
    for (j = 0; j < nl; j++) {
	if ((i * ni + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, G[i][j]);
    }
  POLYBENCH_DUMP_END("G");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_3mm(int ni, int nj, int nk, int nl, int nm,
		DATA_TYPE POLYBENCH_2D(E,NI,NJ,ni,nj),
		DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj),
		DATA_TYPE POLYBENCH_2D(F,NJ,NL,nj,nl),
		DATA_TYPE POLYBENCH_2D(C,NJ,NM,nj,nm),
		DATA_TYPE POLYBENCH_2D(D,NM,NL,nm,nl),
		DATA_TYPE POLYBENCH_2D(G,NI,NL,ni,nl))
{
  int i, j, k;

#pragma scop

   We must not change that. We can put the #define after the #pragma scop.

   However, the problem says: preserve the marker lines. So we leave the #pragma scop as is and then add the #define.

   But note: the #define must be inside the scop? Actually, the scop is a region of code. We can put the #define at the beginning of the scop.

   Alternatively, we can avoid #define and use:

        const int TI = 64;
        const int TJ = 64;

   and then use TI and TJ in the loops. This is acceptable.

   We'll do:

        const int TI = 64;
        const int TJ = 64;

   Now, we write the three tiled and parallelized matrix multiplications.

   We must be cautious: the inner k loop is not tiled. We are only tiling the i and j loops.

   We'll write the code accordingly.

   Let's code the first multiplication (E = A*B):

        const int TI = 64;
        const int TJ = 64;

        // E = A*B
        #pragma omp parallel for collapse(2) private(i, j, k) schedule(static)
        for (int ii = 0; ii < _PB_NI; ii += TI)
          for (int jj = 0; jj < _PB_NJ; jj += TJ)
          {
            int i_end = (ii + TI) < _PB_NI ? (ii + TI) : _PB_NI;
            int j_end = (jj + TJ) < _PB_NJ ? (jj + TJ) : _PB_NJ;
            for (int i = ii; i < i_end; i++)
              for (int j = jj; j < j_end; j++)
              {
                E[i][j] = SCALAR_VAL(0.0);
                for (k = 0; k < _PB_NK; k++)
                    E[i][j] += A[i][k] * B[k][j];
              }
          }

   Similarly for the others.

   But note: the variables i, j, k must be declared. We are using them in the inner loops. We can declare them in the inner scope.

   We'll declare k inside the inner loop? Actually, we can declare it in the inner loop:

        for (int k = 0; k < _PB_NK; k++)

   But the original code declares k at the beginning of the function. We are allowed to restructure. We'll declare k in the inner loop to avoid having to declare it as private in the OpenMP directive? Actually, the OpenMP directive already has 'private(i, j, k)'. But if we declare k inside, then it is private by default? 

   However, the OpenMP specification says that variables declared inside the loop are private. So we can remove k from the private clause? But we are using k in the inner loop and it is declared inside the parallel region? Actually, the inner loop is inside the parallel region. We can declare k inside the inner loop and then it is private to the thread.

   We'll do:

        for (int k = 0; k < _PB_NK; k++)

   and remove k from the private clause? But the OpenMP directive is on the outer loop. The inner loop is inside the parallel region. The variable k is declared in the inner loop, so it is private to the thread. We don't need to list it.

   However, the OpenMP specification says that variables declared inside a parallel region are private. So we can remove k from the private clause.

   But note: the original code declares k at the function level. We are changing the scope. This is allowed because we are restructuring.

   We'll change the inner loop to declare k.

   Similarly, we can declare i and j in the inner loops? Actually, we are using i and j in the inner two loops. We can declare them in the inner loops? 

   We'll do:

        for (int i = ii; i < i_end; i++)
          for (int j = jj; j < j_end; j++)

   Then we don't need to list i and j as private? Because they are declared inside the parallel region. 

   However, the OpenMP directive is on the outer two loops (ii and jj). The inner loops are inside the parallel region. The variables i and j are declared in the inner loops, so they are private to the thread.

   Therefore, we can remove the private clause entirely? 

   But note: the OpenMP specification says that the loop variables of the collapsed loops (ii and jj) are private by default. And the inner loops (i and j) are declared inside, so they are private. And k is declared inside the inner loop, so it is private.

   So we can remove the private clause.

   However, to be safe, we can leave it? But it's redundant. We'll remove it.

   Alternatively, we can keep it for clarity? The problem says: do not add unsupported transformations. We'll remove it.

   So the directive becomes:

        #pragma omp parallel for collapse(2) schedule(static)

   Now, we write the entire scop region.

   We must be cautious: the original code uses _PB_NI, _PB_NJ, etc. which are defined in the polybench header. We assume they are integers.

   We'll write:

        #pragma scop
        const int TI = 64;
        const int TJ = 64;

        // E = A*B
        #pragma omp parallel for collapse(2) schedule(static)
        for (int ii = 0; ii < _PB_NI; ii += TI)
          for (int jj = 0; jj < _PB_NJ; jj += TJ)
          {
            int i_end = (ii + TI) < _PB_NI ? (ii + TI) : _PB_NI;
            int j_end = (jj + TJ) < _PB_NJ ? (jj + TJ) : _PB_NJ;
            for (int i = ii; i < i_end; i++)
              for (int j = jj; j < j_end; j++)
              {
                E[i][j] = SCALAR_VAL(0.0);
                for (int k = 0; k < _PB_NK; k++)
                    E[i][j] += A[i][k] * B[k][j];
              }
          }

        // F = C*D
        #pragma omp parallel for collapse(2) schedule(static)
        for (int ii = 0; ii < _PB_NJ; ii += TI)
          for (int jj = 0; jj < _PB_NL; jj += TJ)
          {
            int i_end = (ii + TI) < _PB_NJ ? (ii + TI) : _PB_NJ;
            int j_end = (jj + TJ) < _PB_NL ? (jj + TJ) : _PB_NL;
            for (int i = ii; i < i_end; i++)
              for (int j = jj; j < j_end; j++)
              {
                F[i][j] = SCALAR_VAL(0.0);
                for (int k = 0; k < _PB_NM; k++)
                    F[i][j] += C[i][k] * D[k][j];
              }
          }

        // G = E*F
        #pragma omp parallel for collapse(2) schedule(static)
        for (int ii = 0; ii < _PB_NI; ii += TI)
          for (int jj = 0; jj < _PB_NL; jj += TJ)
          {
            int i_end = (ii + TI) < _PB_NI ? (ii + TI) : _PB_NI;
            int j_end = (jj + TJ) < _PB_NL ? (jj + TJ) : _PB_NL;
            for (int i = ii; i < i_end; i++)
              for (int j = jj; j < j_end; j++)
              {
                G[i][j] = SCALAR_VAL(0.0);
                for (int k = 0; k < _PB_NJ; k++)
                    G[i][j] += E[i][k] * F[k][j];
              }
          }
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int ni = NI;
  int nj = NJ;
  int nk = NK;
  int nl = NL;
  int nm = NM;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(E, DATA_TYPE, NI, NJ, ni, nj);
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, NI, NK, ni, nk);
  POLYBENCH_2D_ARRAY_DECL(B, DATA_TYPE, NK, NJ, nk, nj);
  POLYBENCH_2D_ARRAY_DECL(F, DATA_TYPE, NJ, NL, nj, nl);
  POLYBENCH_2D_ARRAY_DECL(C, DATA_TYPE, NJ, NM, nj, nm);
  POLYBENCH_2D_ARRAY_DECL(D, DATA_TYPE, NM, NL, nm, nl);
  POLYBENCH_2D_ARRAY_DECL(G, DATA_TYPE, NI, NL, ni, nl);

  /* Initialize array(s). */
  init_array (ni, nj, nk, nl, nm,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B),
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(D));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_3mm (ni, nj, nk, nl, nm,
	      POLYBENCH_ARRAY(E),
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B),
	      POLYBENCH_ARRAY(F),
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(D),
	      POLYBENCH_ARRAY(G));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(ni, nl,  POLYBENCH_ARRAY(G)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(E);
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);
  POLYBENCH_FREE_ARRAY(F);
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(D);
  POLYBENCH_FREE_ARRAY(G);

  return 0;
}

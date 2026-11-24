/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* gemm.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "gemm.h"


/* Array initialization. */
static
void init_array(int ni, int nj, int nk,
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(C,NI,NJ,ni,nj),
		DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nj; j++)
      C[i][j] = (DATA_TYPE) ((i*j+1) % ni) / ni;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nk; j++)
      A[i][j] = (DATA_TYPE) (i*(j+1) % nk) / nk;
  for (i = 0; i < nk; i++)
    for (j = 0; j < nj; j++)
      B[i][j] = (DATA_TYPE) (i*(j+2) % nj) / nj;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int ni, int nj,
		 DATA_TYPE POLYBENCH_2D(C,NI,NJ,ni,nj))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("C");
  for (i = 0; i < ni; i++)
    for (j = 0; j < nj; j++) {
	if ((i * ni + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, C[i][j]);
    }
  POLYBENCH_DUMP_END("C");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_gemm(int ni, int nj, int nk,
		 DATA_TYPE alpha,
		 DATA_TYPE beta,
		 DATA_TYPE POLYBENCH_2D(C,NI,NJ,ni,nj),
		 DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		 DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj))
{
  int i, j, k;

//BLAS PARAMS
//TRANSA = 'N'
//TRANSB = 'N'
// => Form C := alpha*A*B + beta*C,
//A is NIxNK
//B is NKxNJ
//C is NIxNJ
#pragma scop

    We can add the constants after that.

    However, the problem says: "preserve these marker lines unchanged". So we cannot add code between the marker and the loop? Actually, the marker is a pragma and then the code follows. We can add the constants after the pragma.

    We'll do:

        #pragma scop
        const int TB = 64;
        const int TX = 64;
        #pragma omp parallel for default(none) shared(C, A, B, alpha, beta) private(i) schedule(static)
        for (i = 0; i < _PB_NI; i++) {
            ... // rest of the code
        }

    But note: the constants TB and TX are used in the inner loops. They are shared? Actually, they are constants and are the same for all threads. We can list them as shared? Or we can define them inside the parallel region? But we defined them outside the parallel region. They are shared by default. We don't need to list them in the shared clause because they are not used in the parallel for construct? They are used in the inner loops. But the OpenMP standard requires that we list all variables from the outer scope. We must list them as shared? Or we can avoid by not using 'default(none)'? We are using 'default(none)', so we must list them.

    We can list them as shared:

        shared(C, A, B, alpha, beta, TB, TX)

    But note: TB and TX are constants defined in the scop. They are shared.

    Alternatively, we can define them inside the parallel region? But then we would have to define them for each thread. We can do:

        #pragma omp parallel
        {
            const int TB = 64, TX = 64;
            #pragma omp for schedule(static)
            for (i=0; i<_PB_NI; i++) {
                ...
            }
        }

    But that would require restructuring the parallel region. We are currently using a combined 'parallel for'. We can change to a 'parallel' and then a 'for'. But the problem does not require a specific form.

    However, the problem says: "Do not change observable behavior". The behavior is the same.

    We'll do:

        #pragma scop
        #pragma omp parallel default(none) shared(C, A, B, alpha, beta, _PB_NI, _PB_NJ, _PB_NK)
        {
            const int TB = 64, TX = 64;
            #pragma omp for schedule(static)
            for (i = 0; i < _PB_NI; i++) {
                ... // scaling and tiled multiplication
            }
        }

    But note: the variables _PB_NI, _PB_NJ, _PB_NK are used in the loop bounds. We must list them as shared? They are macros, so they are not variables. We don't need to list them.

    Actually, the preprocessor replaces them. So we don't need to list them.

    We can write:

        #pragma omp parallel default(none) shared(C, A, B, alpha, beta)
        {
            const int TB = 64, TX = 64;
            #pragma omp for schedule(static)
            for (i = 0; i < _PB_NI; i++) {
                ... 
            }
        }

    But note: the variables i, j, k, etc. are declared inside the for loop? We are declaring i in the for loop. The other variables (j, k, kk, jj, k_end, j_end) are declared inside the i loop. So they are private to the thread.

    We don't need to list them.

    However, the variable i is the loop variable of the 'omp for', so it is private by default.

    This is acceptable.

    But note: the original code had the entire scop as a single loop. We are changing the structure to a parallel region with a for inside. This is still within the scop.

    We'll do:

        #pragma scop
        #pragma omp parallel default(none) shared(C, A, B, alpha, beta)
        {
            const int TB = 64, TX = 64;
            #pragma omp for schedule(static)
            for (i = 0; i < _PB_NI; i++) {
                // Scale the row i of C by beta
                for (j = 0; j < _PB_NJ; j++)
                    C[i][j] *= beta;

                // Tiled multiplication
                for (kk = 0; kk < _PB_NK; kk += TB) {
                    int k_end = (kk + TB) < _PB_NK ? (kk+TB) : _PB_NK;
                    for (jj = 0; jj < _PB_NJ; jj += TX) {
                        int j_end = (jj + TX) < _PB_NJ ? (jj+TX) : _PB_NJ;
                        for (k = kk; k < k_end; k++) {
                            // Preload A[i][k] if possible? The compiler might do it.
                            // Vectorize the inner j loop
                            #pragma omp simd
                            for (j = jj; j < j_end; j++) {
                                C[i][j] += alpha * A[i][k] * B[k][j];
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
  int ni = NI;
  int nj = NJ;
  int nk = NK;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  DATA_TYPE beta;
  POLYBENCH_2D_ARRAY_DECL(C,DATA_TYPE,NI,NJ,ni,nj);
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,NI,NK,ni,nk);
  POLYBENCH_2D_ARRAY_DECL(B,DATA_TYPE,NK,NJ,nk,nj);

  /* Initialize array(s). */
  init_array (ni, nj, nk, &alpha, &beta,
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_gemm (ni, nj, nk,
	       alpha, beta,
	       POLYBENCH_ARRAY(C),
	       POLYBENCH_ARRAY(A),
	       POLYBENCH_ARRAY(B));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(ni, nj,  POLYBENCH_ARRAY(C)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);

  return 0;
}

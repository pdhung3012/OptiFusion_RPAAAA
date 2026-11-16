We are given a code for the syrk kernel. The scop region is:

#pragma scop
  for (i = 0; i < _PB_N; i++) {
    for (j = 0; j <= i; j++)
      C[i][j] *= beta;
    for (k = 0; k < _PB_M; k++) {
      for (j = 0; j <= i; j++)
        C[i][j] += alpha * A[i][k] * A[j][k];
    }
  }
#pragma endscop

We are to optimize this region.

First, we note that the code has two parts for each i:
1. A scaling loop: for j from 0 to i: C[i][j] *= beta.
2. A triple nested loop: for k from 0 to _PB_M-1, and for j from 0 to i: C[i][j] += alpha * A[i][k] * A[j][k].

We are to consider loop transformations in the order: tiling, fusion, skewing.

1. Loop Tiling:
   The inner two loops (j and k) in the second part are candidates for tiling. However, note that the j loop runs from 0 to i (which is a triangular loop). Tiling a triangular loop is more complex. Also, the k loop is over the entire range of M. We might consider tiling the k loop and the j loop together. But note that the j loop is triangular and the k loop is rectangular. We can tile the k loop and the j loop, but we must be cautious because the j loop depends on i.

   Alternatively, we can consider tiling the i loop and the k loop? But the i loop is the outermost and the k loop is inside. The j loop is the innermost and is triangular.

   However, note that the entire computation for a fixed i is independent of other i's? Actually, the computation for C[i][j] only depends on the j<=i and k, and the inner loop for j and k only uses A[i][k] and A[j][k] (with j<=i). So the outer i loop is parallel? But note: the inner j loop for a fixed i writes to C[i][j] for j<=i. Since each i iteration writes to a distinct set of C (the lower triangle), the i loop is parallel.

   Therefore, we can parallelize the i loop. But before that, we can consider tiling the inner loops to improve cache locality.

   The inner two loops (k and j) in the second part: 
      for (k=0; k<_PB_M; k++)
        for (j=0; j<=i; j++)
          C[i][j] += ... 

   We can interchange the k and j loops? But note: the j loop is triangular and the k loop is rectangular. Interchanging would require a rectangular j loop? Actually, we cannot interchange because the j loop bound depends on i and the k loop does not. However, we can tile the k loop and the j loop? But the j loop is triangular and the tile would have to be adjusted for each i.

   Alternatively, we can tile the k loop and then the j loop? But the j loop is the innermost and we are accessing A[i][k] and A[j][k]. Note that A[i][k] is fixed for a given i and k, and A[j][k] is accessed for j from 0 to i. So for a fixed k, we are accessing a column of A (the k-th column) from row 0 to i. This is a contiguous access in j? Actually, in C/C++ the array is row-major, so A[j][k] for fixed k and varying j is contiguous? Actually, no: because the array is row-major, the elements A[j][k] for fixed k and consecutive j are separated by the row length (which is M). So it's not contiguous. However, we are accessing the same column k for many j, so we are striding by M. This is bad for cache.

   We can try to interchange the j and k loops? Then we would have:
        for (j=0; j<=i; j++)
          for (k=0; k<_PB_M; k++)
            C[i][j] += alpha * A[i][k] * A[j][k];

   Now, for fixed j, we are accessing A[i][k] and A[j][k] for k from 0 to M-1. This is contiguous in k (since k is the last index). So this interchange would improve spatial locality for A.

   But note: the original loop order was k then j. The interchange is legal? There is no dependence between k and j? The inner loop writes to C[i][j] and the reads are from A. The writes to C[i][j] are accumulated and the order of accumulation does not matter (addition is commutative). So the interchange is legal.

   Therefore, we can interchange the k and j loops in the second part.

   After interchange, the inner two loops become:
        for (j=0; j<=i; j++)
          for (k=0; k<_PB_M; k++)
            C[i][j] += alpha * A[i][k] * A[j][k];

   Now, we can consider tiling the j and k loops? But note that the j loop is triangular and the k loop is rectangular. We can tile the k loop (which is now the innermost) and the j loop? But the j loop is the outer of the two. We can tile the j loop? But the j loop runs from 0 to i, which is not a constant. We can use a rectangular tile for j? But the j loop is triangular. Alternatively, we can tile the k loop only? But the k loop is the innermost and we are accessing contiguous memory in k. Tiling the k loop might not be necessary because the inner loop is already contiguous and we can rely on vectorization.

   However, note that the j loop is over a triangular region. We can break the j loop into blocks? But the block size would have to be adjusted for each i. This is complex.

   Alternatively, we can consider tiling the i loop? The i loop is the outermost and we are going to parallelize it. Tiling the i loop would break the parallelism? We don't want to do that.

   So, for now, we do the interchange and then consider tiling the k loop? But the k loop is the innermost and we are accessing contiguous memory. We can leave it as is and rely on vectorization.

   But note: the k loop is of size M, which might be large. We can tile the k loop to improve cache locality? We are accessing two rows of A: row i and row j. For a fixed j, we are accessing the entire row j and the entire row i? Actually, we are accessing the entire row i for every j? That is, for a fixed i, we are accessing row i for every j and k. And for each j, we are accessing row j. So the entire inner loop for a fixed i will access:
        row i: M times (once per k) for each j -> total M * (i+1) times? Actually, for each j, we access row i for k=0 to M-1 -> that's M accesses per j. And there are i+1 j's -> total M*(i+1) accesses to row i? And for row j, we access it for k=0 to M-1 -> M accesses per j. And we do that for each j from 0 to i -> total M*(i+1) accesses to the set of rows j (but each row j is accessed only once per j).

   The total data accessed for a fixed i is:
        row i: M*(i+1) accesses? Actually, no: for a fixed i and j, we access row i for every k -> that's M accesses to row i per j. And there are i+1 j's -> so row i is accessed M*(i+1) times? That is a lot. And the same row i is being accessed repeatedly. We are also accessing each row j (for j<=i) M times.

   The entire array A is of size n x m. We are accessing the entire row i and the entire rows j for j<=i. The total data accessed for a fixed i is: 
        row i: M elements, and rows j (j=0..i): each of M elements -> total M + (i+1)*M = M*(i+2) elements? But note: row i is accessed for every j and k -> that's M*(i+1) accesses? But the same element A[i][k] is accessed for every j? So we are accessing the same element multiple times.

   We can restructure to avoid reloading A[i][k] for every j? We can preload A[i][k] for a fixed k and then use it for all j? That is, we can interchange the j and k loops? We already did that. Now, for a fixed k, we are accessing A[i][k] and then for each j we access A[j][k]. So we can preload A[i][k] and then use it for all j? That would reduce the number of times we load A[i][k] from M*(i+1) to M (once per k). And for A[j][k], we are accessing it once per j and k.

   So the interchange already helps.

   Now, we can consider tiling the k loop? We are accessing A[i][k] and A[j][k] for k in a tile. We can load a tile of A[i][k] for a block of k, and then for each j, we load a tile of A[j][k] for the same block of k? This would improve cache locality because we are reusing the same block of k for multiple j.

   We can do:
        for (j=0; j<=i; j++) {
          for (k=0; k<_PB_M; k+=TILE_K) {
            for (kk=0; kk<TILE_K; kk++) {
              C[i][j] += alpha * A[i][k+kk] * A[j][k+kk];
            }
          }
        }

   But note: the inner loop is over kk, and we are accumulating into C[i][j]. This is a reduction. We can do that.

   However, we must choose a tile size TILE_K that is a divisor of M? Or we can handle the remainder. But the problem says to choose legal tile sizes. We can use a constant tile size, say 64 or 128, and then handle the remainder with a cleanup loop.

   But note: the original code does not have a remainder. We are allowed to restructure the loops. We can do:

        for (j=0; j<=i; j++) {
          int k;
          for (k=0; k<_PB_M-TILE_K+1; k+=TILE_K) {
            int kk;
            for (kk=0; kk<TILE_K; kk++) {
              C[i][j] += alpha * A[i][k+kk] * A[j][k+kk];
            }
          }
          // remainder
          for (; k<_PB_M; k++) {
            C[i][j] += alpha * A[i][k] * A[j][k];
          }
        }

   But this adds overhead. Alternatively, we can use a tile size that divides M? But M is a variable. We cannot assume that. We can use a fixed tile size and then the remainder loop.

   However, the problem says: "Choose legal tile sizes". We can choose a tile size that is a power of two and less than or equal to the cache line size? But we don't know the cache line size. We can choose a tile size of 64 or 128? But we must be cautious because the inner loop is now over kk and we are doing a reduction. The compiler might vectorize the inner kk loop? That would be good.

   But note: the inner kk loop is a reduction and the order of the kk loop does not matter. So we can do the tiling.

   However, we must also consider that the j loop is triangular. We are not tiling the j loop. We are only tiling the k loop.

   Alternatively, we can tile both j and k? But the j loop is triangular and the tile would have to be adjusted for each i. This is complex and might not be worth it.

   So, we decide to do:
      - Interchange the k and j loops in the second part (so that k becomes the innermost).
      - Then tile the k loop (which is now the innermost) with a fixed tile size.

   But note: the original code has two parts: the scaling and the triple loop. We can also consider fusing the scaling and the triple loop? But the scaling is a separate loop. We can try to fuse the scaling with the triple loop? Actually, the scaling is done for each i and j, and then the triple loop is done for the same i and j. We can combine them? We can do:

        for (i=0; i<_PB_N; i++) {
          for (j=0; j<=i; j++) {
            C[i][j] *= beta;
            for (k=0; k<_PB_M; k++) {
              C[i][j] += alpha * A[i][k] * A[j][k];
            }
          }
        }

   But note: the original code does the scaling for all j in the lower triangle for a fixed i, and then does the triple loop for the same i and all k and j. The above fused loop does the scaling and then immediately the triple loop for the same j. This is valid? Yes, because the scaling is done first and then the accumulation. And the accumulation does not depend on the scaling of other j's? It only depends on the current j and k.

   However, the original code does the scaling for all j first and then the triple loop. The fused version does the scaling and then the triple loop for each j. This changes the order: in the original, the scaling for j=0..i is done first, then the triple loop for j=0..i and k=0..M-1. In the fused version, for a fixed j, we do the scaling and then the triple loop for that j. This is acceptable? The scaling for a given j is independent of the triple loop for other j's? Yes. And the triple loop for a given j is independent of the scaling for other j's? Yes. So it is legal.

   But note: the triple loop in the original code is:

        for (k=0; k<_PB_M; k++) {
          for (j=0; j<=i; j++)
            C[i][j] += ... 

   This order is k then j. We are fusing and then interchanging? Actually, in the fused version we are doing j then k. We already decided to interchange to j then k. So the fused version with j then k is:

        for (i=0; i<_PB_N; i++) {
          for (j=0; j<=i; j++) {
            C[i][j] *= beta;
            for (k=0; k<_PB_M; k++) {
              C[i][j] += alpha * A[i][k] * A[j][k];
            }
          }
        }

   This is the same as the fused version we wrote above.

   Now, we can tile the k loop? We can do:

        for (i=0; i<_PB_N; i++) {
          for (j=0; j<=i; j++) {
            C[i][j] *= beta;
            int k;
            for (k=0; k<_PB_M; k+=TILE_K) {
              int kk_end = (k+TILE_K <= _PB_M) ? TILE_K : _PB_M - k;
              int kk;
              for (kk=0; kk<kk_end; kk++) {
                C[i][j] += alpha * A[i][k+kk] * A[j][k+kk];
              }
            }
          }
        }

   But note: the inner kk loop is a reduction and we are accumulating into C[i][j]. This is safe.

   However, we must choose a tile size. Let's choose TILE_K = 64. We assume that 64 is a reasonable tile size for cache lines (typically 64 bytes, and if DATA_TYPE is double, then 8 bytes per element -> 8 elements per cache line). But we don't know the exact cache line size. We can choose 64 as a starting point.

   But note: the problem says "choose legal tile sizes". We can choose any positive integer. We choose 64.

   Alternatively, we can avoid the remainder by using a loop that goes to _PB_M and then use a mask? But that is complex. We do the remainder loop.

   However, we can also use a fixed tile size and then the compiler might unroll the inner loop? We can rely on the compiler to vectorize the inner kk loop? That would be good.

   But note: the inner kk loop is a reduction and the compiler might vectorize it with vector reduction.

   So, we do:
      - Fuse the scaling and the triple loop, and then interchange the k and j loops (so that k is innermost) and then tile the k loop.

   But wait: in the fused version, we have j then k. We are tiling the k loop. We are not changing the j loop.

   Now, we must also consider parallelism. The outer i loop is parallel. We can add an OpenMP parallel for to the i loop.

   Also, we can consider the j loop? The j loop is triangular and the iterations for different j are independent? Yes, because each j writes to a distinct C[i][j]. So we can also parallelize the j loop? But note: the j loop is inside i. We can collapse the i and j loops? But the j loop bound depends on i. We cannot collapse because the j loop is triangular. We can only parallelize the i loop.

   So we do:

        #pragma omp parallel for private(j,k) // we will add more clauses later
        for (i=0; i<_PB_N; i++) {
          for (j=0; j<=i; j++) {
            ... 
          }
        }

   But note: the j loop is not parallel because the j loop is inside the i loop and we are parallelizing i. Each thread will get a set of i's and then do the j loop for that i.

   Now, we must consider data sharing. The arrays C and A are shared. The variables i, j, k are private. We must also consider alpha and beta: they are shared and read-only. We can use 'shared(alpha, beta)' or 'firstprivate'? But they are scalars and read-only, so shared is fine. But to be safe, we can use 'default(none)' and then specify.

   Also, we must consider the reduction in the inner loop? The inner kk loop is a reduction on C[i][j]. But note: the entire j loop for a fixed i is done by one thread? Because the i loop is parallel and each i is assigned to one thread. And for a fixed i, the j loop is done sequentially by that thread. So there is no reduction across threads for C[i][j]? Each C[i][j] is updated only by the thread that handles i. So we don't need a reduction clause for C.

   Therefore, we can do:

        #pragma omp parallel for private(j,k) default(none) shared(C, A, alpha, beta, _PB_N, _PB_M)

   But note: the tile size TILE_K is a constant. We can define it as a macro? Or we can use a constant. We'll define it as a constant.

   However, the problem says: do not change code outside the scop. We can define the tile size inside the scop? We can do:

        #define TILE_K 64

   But we must be cautious: if the scop is in a function, we can define it at the beginning of the scop? But the scop is inside a function. We can define it inside the scop region? Yes.

   Alternatively, we can use a constant without defining a macro? We can use:

        const int TILE_K = 64;

   But note: the scop region is inside a function. We can declare a constant inside the scop.

   However, the problem says: preserve the marker lines. We can add the constant inside the scop.

   But note: the original code uses _PB_N and _PB_M. We are using them. We must be cautious because _PB_N and _PB_M are macros? They are defined in the polybench header. We can use them.

   Now, let's write the transformed code:

        #pragma scop
        #define TILE_K 64
        #pragma omp parallel for private(j,k) default(none) shared(C, A, alpha, beta)
        for (i = 0; i < _PB_N; i++) {
          for (j = 0; j <= i; j++) {
            C[i][j] *= beta;
            int k;
            for (k = 0; k < _PB_M; k += TILE_K) {
              int kk_end = (k + TILE_K <= _PB_M) ? TILE_K : _PB_M - k;
              int kk;
              for (kk = 0; kk < kk_end; kk++) {
                C[i][j] += alpha * A[i][k+kk] * A[j][k+kk];
              }
            }
          }
        }
        #undef TILE_K
        #pragma endscop

   But note: the original code had two separate loops for the scaling and the triple loop. We fused them. This is a loop fusion transformation.

   We also interchanged the k and j loops (which is a loop interchange) and then tiled the k loop.

   However, we must check: is the fusion legal? We already argued that it is.

   But note: the original code did the scaling for all j first and then the triple loop. The fused version does the scaling and then the triple loop for each j. This changes the order of the scaling and the triple loop for different j's? In the original, the scaling for j=0..i is done first, then the triple loop for j=0..i and k=0..M-1. In the fused version, for a fixed j, we do the scaling and then the triple loop for that j. This is a different order. However, the scaling for a given j is independent of the triple loop for other j's? Yes. And the triple loop for a given j is independent of the scaling for other j's? Yes. So it is legal.

   But note: the triple loop in the original code is:

        for (k=0; k<_PB_M; k++) {
          for (j=0; j<=i; j++)
            C[i][j] += ... 

   This order is k then j. We changed it to j then k. We already argued that the interchange is legal because the accumulation is commutative.

   Therefore, the transformation is legal.

   Now, we must consider the OpenMP directives. We added a 'parallel for' for the i loop. We used 'private(j,k)' and 'default(none)' and shared the necessary arrays and scalars.

   We did not use 'collapse' because the j loop bound depends on i and is triangular.

   We did not use any scheduling clause. The default is static, which is fine.

   We did not use 'nowait' because there is no next loop to avoid synchronization.

   We did not use 'atomic' or 'critical' because there is no race condition: each thread writes to a distinct set of C (the rows i that it is assigned).

   We did not use reduction because the reduction is per thread and per C[i][j] and is done sequentially within the thread.

   We did not use SIMD directives. We can add a 'simd' directive to the inner kk loop? The inner kk loop is a reduction and the compiler might vectorize it. We can add:

        #pragma omp simd

   to the kk loop? But note: the kk loop is inside the k loop and the j loop. We can do:

        for (kk=0; kk<kk_end; kk++) 
          C[i][j] += ... 

   We can add:

        #pragma omp simd reduction(+:C[i][j])

   But note: the reduction is on a scalar that is accumulated. The compiler might do that automatically. We can add the directive to ensure vectorization.

   However, the problem says: we can use 'simd' and 'reduction'. We can do:

        #pragma omp simd reduction(+:C[i][j])
        for (kk=0; kk<kk_end; kk++) {
          C[i][j] += alpha * A[i][k+kk] * A[j][k+kk];
        }

   But note: the reduction is on a single element C[i][j]. The compiler will create a temporary vector for the reduction and then combine. This is safe.

   We can also consider aligning the arrays? We don't know the alignment. We can skip that.

   Also, we can consider the thread management? We can leave it to the runtime.

   Therefore, we add the simd directive to the inner kk loop.

   But note: the inner kk loop is inside the k loop and the j loop. We are already in a parallel region. The simd directive is for the innermost loop.

   So the final code for the scop region:

        #pragma scop
        #define TILE_K 64
        #pragma omp parallel for private(j,k) default(none) shared(C, A, alpha, beta, _PB_N, _PB_M)
        for (i = 0; i < _PB_N; i++) {
          for (j = 0; j <= i; j++) {
            C[i][j] *= beta;
            int k;
            for (k = 0; k < _PB_M; k += TILE_K) {
              int kk_end = (k + TILE_K <= _PB_M) ? TILE_K : _PB_M - k;
              int kk;
              #pragma omp simd reduction(+:C[i][j])
              for (kk = 0; kk < kk_end; kk++) {
                C[i][j] += alpha * A[i][k+kk] * A[j][k+kk];
              }
            }
          }
        }
        #undef TILE_K
        #pragma endscop

   We must be cautious: the variable 'kk_end' is computed and then used. It is private to the thread because it is inside the parallel region and the i loop is distributed.

   Also, note that the original code used _PB_N and _PB_M. We are using them. We assume they are defined.

   We did not change the function signature or the program I/O.

   We must check: is the tile size 64 legal? It is a positive integer. And we handle the remainder. So it is legal.

   Therefore, we output the entire code with the scop region replaced.

   But note: the original code had:

        for (i = 0; i < _PB_N; i++) {
          for (j = 0; j <= i; j++)
            C[i][j] *= beta;
          for (k = 0; k < _PB_M; k++) {
            for (j = 0; j <= i; j++)
              C[i][j] += alpha * A[i][k] * A[j][k];
          }
        }

   We replaced it with the fused and tiled version.

   We also added the OpenMP directives.

   We must preserve the rest of the code.

   Let's write the entire code with the scop region replaced.

   Note: we must not change anything outside the scop.

   We also note that the original code uses _PB_N and _PB_M. We are using them.

   We also note that the tile size 64 is chosen arbitrarily. We could choose a different size, but 64 is a common choice.

   We also note that the inner kk loop might be vectorized by the compiler even without the simd directive? But we added it to be explicit.

   We also note that the reduction in the simd directive is for the accumulation. The compiler will handle it.

   Therefore, we output the code.
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "syrk.h"


/* Array initialization. */
static
void init_array(int n, int m,
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(C,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(A,N,M,n,m))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;
  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      A[i][j] = (DATA_TYPE) ((i*j+1)%n) / n;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      C[i][j] = (DATA_TYPE) ((i*j+2)%m) / m;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(C,N,N,n,n))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("C");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
	if ((i * n + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, C[i][j]);
    }
  POLYBENCH_DUMP_END("C");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_syrk(int n, int m,
		 DATA_TYPE alpha,
		 DATA_TYPE beta,
		 DATA_TYPE POLYBENCH_2D(C,N,N,n,n),
		 DATA_TYPE POLYBENCH_2D(A,N,M,n,m))
{
  int i, j, k;

//BLAS PARAMS
//TRANS = 'N'
//UPLO  = 'L'
// =>  Form  C := alpha*A*A**T + beta*C.
//A is NxM
//C is NxN
#pragma scop
#define TILE_K 64
#pragma omp parallel for private(j,k) default(none) shared(C, A, alpha, beta, _PB_N, _PB_M)
  for (i = 0; i < _PB_N; i++) {
    for (j = 0; j <= i; j++) {
      C[i][j] *= beta;
      for (k = 0; k < _PB_M; k += TILE_K) {
        int kk_end = (k + TILE_K <= _PB_M) ? TILE_K : _PB_M - k;
        int kk;
        #pragma omp simd reduction(+:C[i][j])
        for (kk = 0; kk < kk_end; kk++) {
          C[i][j] += alpha * A[i][k+kk] * A[j][k+kk];
        }
      }
    }
  }
#undef TILE_K
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;
  int m = M;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  DATA_TYPE beta;
  POLYBENCH_2D_ARRAY_DECL(C,DATA_TYPE,N,N,n,n);
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,N,M,n,m);

  /* Initialize array(s). */
  init_array (n, m, &alpha, &beta, POLYBENCH_ARRAY(C), POLYBENCH_ARRAY(A));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_syrk (n, m, alpha, beta, POLYBENCH_ARRAY(C), POLYBENCH_ARRAY(A));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(C)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(A);

  return 0;
}
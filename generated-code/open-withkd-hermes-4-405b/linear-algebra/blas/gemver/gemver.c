/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* gemver.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "gemver.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE *alpha,
		 DATA_TYPE *beta,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		 DATA_TYPE POLYBENCH_1D(u1,N,n),
		 DATA_TYPE POLYBENCH_1D(v1,N,n),
		 DATA_TYPE POLYBENCH_1D(u2,N,n),
		 DATA_TYPE POLYBENCH_1D(v2,N,n),
		 DATA_TYPE POLYBENCH_1D(w,N,n),
		 DATA_TYPE POLYBENCH_1D(x,N,n),
		 DATA_TYPE POLYBENCH_1D(y,N,n),
		 DATA_TYPE POLYBENCH_1D(z,N,n))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;

  DATA_TYPE fn = (DATA_TYPE)n;

  for (i = 0; i < n; i++)
    {
      u1[i] = i;
      u2[i] = ((i+1)/fn)/2.0;
      v1[i] = ((i+1)/fn)/4.0;
      v2[i] = ((i+1)/fn)/6.0;
      y[i] = ((i+1)/fn)/8.0;
      z[i] = ((i+1)/fn)/9.0;
      x[i] = 0.0;
      w[i] = 0.0;
      for (j = 0; j < n; j++)
        A[i][j] = (DATA_TYPE) (i*j % n) / n;
    }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_1D(w,N,n))
{
  int i;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("w");
  for (i = 0; i < n; i++) {
    if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, w[i]);
  }
  POLYBENCH_DUMP_END("w");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_gemver(int n,
		   DATA_TYPE alpha,
		   DATA_TYPE beta,
		   DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		   DATA_TYPE POLYBENCH_1D(u1,N,n),
		   DATA_TYPE POLYBENCH_1D(v1,N,n),
		   DATA_TYPE POLYBENCH_1D(u2,N,n),
		   DATA_TYPE POLYBENCH_1D(v2,N,n),
		   DATA_TYPE POLYBENCH_1D(w,N,n),
		   DATA_TYPE POLYBENCH_1D(x,N,n),
		   DATA_TYPE POLYBENCH_1D(y,N,n),
		   DATA_TYPE POLYBENCH_1D(z,N,n))
{
  int i, j;

#pragma scop
We are given the gemver.c code and are to optimize the region between #pragma scop and #pragma endscop.
 The original code in the scop region has four loops:

 1. A double loop updating A[i][j] = A[i][j] + u1[i] * v1[j] + u2[i] * v2[j]
 2. A double loop updating x[i] = x[i] + beta * A[j][i] * y[j]
 3. A single loop updating x[i] = x[i] + z[i]
 4. A double loop updating w[i] = w[i] + alpha * A[i][j] * x[j]

 We are to apply loop transformations (tiling, fusion, skewing) first, then OpenMP directives.

 Step 1: Loop Transformations
   - Tiling: We note that the first and fourth loops are over A and might benefit from tiling for cache locality.
     However, the second loop accesses A[j][i] (transposed) and the fourth loop accesses A[i][j] (normal). 
     Also, the second loop writes to x and the fourth loop reads from x. There is a loop-carried dependence on x between the second and third loops, and then the fourth loop uses x.

   - Fusion: We have four separate loops. We cannot fuse the first loop with the second because the second uses the updated A? 
     Actually, the first loop updates A, and the second loop uses A. So they are in sequence and cannot be fused arbitrarily.

   However, note that the third loop is a single loop that can be fused with the second loop? 
     The second loop: for i, for j: x[i] += ... 
     The third loop: for i: x[i] += z[i]

     We can fuse the second and third loops? But note: the second loop is a double loop and the third is a single loop. 
     We can do:

        for (i=0; i<N; i++) {
            for (j=0; j<N; j++) {
                x[i] += beta * A[j][i] * y[j];
            }
            x[i] += z[i];
        }

     This is legal because the inner j loop for a fixed i only updates x[i] and then we add z[i] to the same x[i]. 
     This fusion reduces one pass over the x array.

   Similarly, the first loop and the fourth loop are both double loops over A. However, the first loop writes to A and the fourth loop reads from A. 
   They are separated by the second and third loops (which use A and then update x). So we cannot fuse the first and fourth.

   But note: the fourth loop uses x, which is updated by the second and third loops. So the first loop must come first, then the second and third, then the fourth.

   So we can only fuse the second and third loops.

   Step 2: OpenMP
     We can parallelize the outer loops of the double loops and the fused loop.

     For the first loop: 
        for (i=0; i<N; i++)
          for (j=0; j<N; j++)
            A[i][j] = ... 

        We can parallelize the i loop because each i is independent? 
        But note: the inner loop writes to A[i][j] for a fixed i and all j. So each i is independent.

     For the fused loop (second and third):
        for (i=0; i<N; i++) {
            for (j=0; j<N; j++) {
                x[i] += ...   // reduction on x[i]
            }
            x[i] += z[i];
        }

        We can parallelize the i loop. However, note that the inner j loop is a reduction on x[i]. 
        We can use a reduction clause for x? But note: the reduction is per i, and then we add z[i] to the same x[i]. 
        We can do:

          #pragma omp parallel for private(j) reduction(+:x[:N]) ... 

        However, the reduction on an array section might be supported in OpenMP 4.5 and later. Alternatively, we can use a private array for the reduction and then combine? 
        But note: the inner loop is a reduction on a single element x[i] for each i. So we can do:

          #pragma omp parallel for private(j)
          for (i=0; i<N; i++) {
              double tmp = 0.0;
              for (j=0; j<N; j++) {
                  tmp += beta * A[j][i] * y[j];
              }
              x[i] = tmp + z[i];   // but note: the original x[i] was already initialized? 
          }

        However, the original code had x[i] initially set to 0 in init_array, and then the second loop does x[i] += ... and then the third loop does x[i] += z[i]. 
        So we are effectively doing: x[i] = (initial 0 + the double loop sum) + z[i]. 

        But wait: the original code does not reset x[i] to 0 at the beginning of the second loop. The init_array sets x[i]=0.0. 
        So we can do:

          #pragma omp parallel for private(j)
          for (i=0; i<N; i++) {
              double tmp = 0.0;
              for (j=0; j<N; j++) {
                  tmp += beta * A[j][i] * y[j];
              }
              x[i] = tmp + z[i];
          }

        This is equivalent and avoids the reduction on the array.

     For the fourth loop:
        for (i=0; i<N; i++)
          for (j=0; j<N; j++)
            w[i] += alpha * A[i][j] * x[j];

        Similarly, we can do:

          #pragma omp parallel for private(j)
          for (i=0; i<N; i++) {
              double tmp = 0.0;
              for (j=0; j<N; j++) {
                  tmp += alpha * A[i][j] * x[j];
              }
              w[i] = tmp;   // but note: w[i] was initialized to 0 and then we are adding? 
          }

        However, the original w[i] is set to 0 in init_array and then we do w[i] += ... so we are effectively setting w[i] to the sum.

        But note: the original code does not reset w[i] to 0 at the beginning of the fourth loop. So we must accumulate? 
        Actually, the original code does:

          w[i] = w[i] + ... 

        and w[i] was set to 0 in init_array. So we can do:

          #pragma omp parallel for private(j)
          for (i=0; i<N; i++) {
              double tmp = 0.0;
              for (j=0; j<N; j++) {
                  tmp += alpha * A[i][j] * x[j];
              }
              w[i] = tmp;   // because the initial w[i] is 0, so we are effectively doing w[i] = 0 + tmp.
          }

        However, what if the scop region is called multiple times? The problem says we must preserve semantics. 
        The original code does not reset w[i] to 0 at the beginning of the scop region. But the init_array sets w[i]=0.0. 
        And the scop region is called only once in the kernel. So we are safe.

        But note: the original code does:

          w[i] = w[i] + ... 

        and we are replacing it with:

          w[i] = tmp   (which is the sum of the products)

        This is equivalent because the initial w[i] is 0.

     However, we must be cautious: the original code might be called multiple times? The problem says we must preserve the original semantics. 
     The kernel_gemver function is called once in main. And the init_array sets w[i]=0.0. So we are safe.

     But what if the kernel is called multiple times? The problem does not specify. We must preserve the original behavior: 
        The original code accumulates into w[i]. So if the kernel is called multiple times, then w[i] would accumulate the results of each call.

     However, the problem says: "Preserve function signatures, program inputs/outputs, and original semantics." 
     And the original code does not reset w[i] to 0 at the beginning of the kernel. So we must accumulate.

     Therefore, we cannot simply set w[i] = tmp. We must do:

          w[i] += tmp;

     Similarly, for the fused loop, we must do:

          x[i] = (original x[i] which is 0) + tmp + z[i] -> but wait, the original x[i] is 0 at the beginning of the kernel? 
          The init_array sets x[i]=0.0. And the kernel does not reset it. So we are safe to set x[i] = tmp + z[i] because the initial value is 0.

     However, note: the first loop updates A, then the second and third update x, then the fourth updates w. 
     The original x[i] is set to 0 in init_array and then updated in the second and third loops. 
     The original w[i] is set to 0 in init_array and then updated in the fourth loop.

     So we can do:

        For the fused loop: 
          x[i] = (initial 0) + (sum over j) + z[i] -> so we set x[i] = tmp + z[i]

        For the fourth loop:
          w[i] = (initial 0) + (sum over j) -> so we set w[i] = tmp

     But wait: the original code for the fourth loop does:

          w[i] = w[i] + ... 

        and w[i] is 0 initially. So we are effectively setting w[i] to the sum. 

     However, if the kernel is called multiple times, then the original code would accumulate. 
     But the problem says the kernel is called once. And the init_array is called only once. 

     But the problem says: "Preserve function signatures, program inputs/outputs, and original semantics." 
     We must not change the behavior if the kernel is called multiple times. 

     How is the kernel used? The main function calls it once. But the problem does not specify that it is called only once. 

     We must be safe: the original code accumulates. So we should do:

        For the fused loop: 
          We are replacing the original two loops that do:
             x[i] += (double loop) 
             x[i] += z[i]
          with:
             x[i] = (double loop result) + z[i]

          This is equivalent to the original because the initial x[i] is 0? 
          But what if the kernel is called again? Then the original x[i] would have the value from the previous call. 
          We are setting x[i] to a new value, not accumulating. 

        This is a problem.

     Therefore, we must do:

        For the fused loop: 
          We must accumulate the result of the double loop and then add z[i] to the existing x[i]. 

        We can do:

          #pragma omp parallel for private(j)
          for (i=0; i<N; i++) {
              double tmp = 0.0;
              for (j=0; j<N; j++) {
                  tmp += beta * A[j][i] * y[j];
              }
              x[i] += tmp;   // accumulate the double loop part
          }
          // Then the third loop: 
          for (i=0; i<N; i++) {
              x[i] += z[i];
          }

        But then we break the fusion? 

        Alternatively, we can do:

          #pragma omp parallel for private(j)
          for (i=0; i<N; i++) {
              double tmp = 0.0;
              for (j=0; j<N; j++) {
                  tmp += beta * A[j][i] * y[j];
              }
              x[i] += tmp + z[i];   // but note: this is not the same as (x[i] += tmp) and then (x[i] += z[i])? 
          }

        Actually, it is the same: x[i] = x[i] + tmp + z[i] is the same as (x[i] = x[i] + tmp) and then (x[i] = x[i] + z[i]).

        So we can do:

          #pragma omp parallel for private(j)
          for (i=0; i<N; i++) {
              double tmp = 0.0;
              for (j=0; j<N; j++) {
                  tmp += beta * A[j][i] * y[j];
              }
              x[i] += tmp + z[i];
          }

        This is safe for multiple calls.

     Similarly, for the fourth loop:

          #pragma omp parallel for private(j)
          for (i=0; i<N; i++) {
              double tmp = 0.0;
              for (j=0; j<N; j++) {
                  tmp += alpha * A[i][j] * x[j];
              }
              w[i] += tmp;   // accumulate
          }

     But note: the original code for the fourth loop does:

          w[i] = w[i] + ... 

        and we are doing the same.

     However, we must be cautious: the inner loop for the fourth loop uses x[j] which is the entire x array. 
     And we are updating w[i] for each i. There is no dependence between different i, so it is parallel.

     But note: the inner loop for the fourth loop: j from 0 to N-1, and we are reading x[j] for all j. 
     This is read-only in the inner loop. So it is safe.

     Now, what about the first loop? 

        for (i=0; i<N; i++)
          for (j=0; j<N; j++)
            A[i][j] = A[i][j] + u1[i] * v1[j] + u2[i] * v2[j];

        We can parallelize the i loop. Each i is independent.

        We can do:

          #pragma omp parallel for private(j)
          for (i=0; i<N; i++) {
              for (j=0; j<N; j++) {
                  A[i][j] = A[i][j] + u1[i] * v1[j] + u2[i] * v2[j];
              }
          }

        This is safe.

     However, we might consider tiling the first and fourth loops for cache locality? 

        The first loop: 
          A[i][j] is written and read. The access is row-wise: A[i][j] for fixed i and j from 0 to N-1 -> contiguous in j.
          u1[i] and u2[i] are constant for the entire inner loop for a fixed i.
          v1[j] and v2[j] are accessed in order j.

        This has good spatial locality for A and the v arrays. 

        The fourth loop:
          A[i][j] is read in row-wise order (contiguous in j) for fixed i.
          x[j] is read in order j (contiguous).
          w[i] is written for each j in the inner loop -> but we are accumulating in a temporary and then writing once? 
          Actually, in our transformed version we are writing w[i] once per i.

        So the fourth loop also has good locality.

        Therefore, we might not need tiling? 

        But note: the first loop writes to A, and then the second loop (which we fused) reads A in column-wise order (A[j][i] for fixed i and j). 
        This is a transpose access and will cause cache misses. 

        We might consider changing the second loop to row-wise? 

        The second loop originally:

          for (i=0; i<N; i++)
            for (j=0; j<N; j++)
              x[i] += beta * A[j][i] * y[j];

        This accesses A[j][i] for fixed i and j: this is the i-th column of A. 

        We can change the loop order? 

          for (j=0; j<N; j++)
            for (i=0; i<N; i++)
              x[i] += beta * A[j][i] * y[j];

        Now, for fixed j, we are accessing A[j][i] for i from 0 to N-1 -> contiguous row. 
        And we are updating x[i] for each i. 

        But note: this is a reduction on x[i] for each i. We can do:

          for (j=0; j<N; j++) {
              #pragma omp parallel for
              for (i=0; i<N; i++) {
                  x[i] += beta * A[j][i] * y[j];
              }
          }

        But then we have a parallel region inside a j loop -> too expensive.

        Alternatively, we can do:

          #pragma omp parallel for private(j)
          for (i=0; i<N; i++) {
              double tmp = 0.0;
              for (j=0; j<N; j++) {
                  tmp += beta * A[j][i] * y[j];
              }
              x[i] += tmp + z[i];
          }

        This is the same as before. The access to A is by columns. 

        We might consider tiling the second loop? 

        We can tile the i and j loops for the second loop? 

          for (ii=0; ii<N; ii+=TILE_I)
            for (jj=0; jj<N; jj+=TILE_J)
              for (i=ii; i<ii+TILE_I; i++)
                for (j=jj; j<jj+TILE_J; j++)
                  ... 

        But note: the inner two loops are over i and j, and we are accessing A[j][i] -> which is not contiguous in either loop. 
        Tiling might help by reusing a tile of A in cache? 

        However, the access pattern is not straightforward. 

        Given the complexity and the fact that the problem does not specify the size of N, we might leave it as is and rely on the fused loop and the temporary to avoid multiple writes to x[i]. 

        We'll focus on the OpenMP parallelization and the fusion we have.

     Summary of transformations:

        We will:
          - Fuse the second and third loops into one loop that does:
                x[i] += (sum over j of beta * A[j][i] * y[j]) + z[i]
          - Use a temporary variable to hold the inner sum to avoid multiple reductions on x[i] (and to avoid false sharing?).

        We will parallelize the outer i loops of the first loop, the fused loop, and the fourth loop.

     Step 3: OpenMP clauses

        We'll use:
          #pragma omp parallel for private(j)   [for the first and fourth loops]
          #pragma omp parallel for private(j)   [for the fused loop]

        We might consider schedule(static) for load balancing? The work per i is the same (N inner iterations). 

        We might also consider collapse? But the inner loop is not parallel (it's a reduction in the fused loop and the fourth loop, and a simple assignment in the first loop). 

        We'll use:

          #pragma omp parallel for schedule(static) private(j)

        We'll also use default(none) to be safe? But the problem says we can use data sharing clauses. 

        We'll list the shared variables: 
          For the first loop: A, u1, v1, u2, v2, N (or _PB_N) are shared.
          For the fused loop: A, beta, y, z, x, N are shared.
          For the fourth loop: A, alpha, x, w, N are shared.

        We'll use:

          #pragma omp parallel for default(none) shared(A, u1, v1, u2, v2, _PB_N) private(i, j)
          for (i=0; i<_PB_N; i++) {
              for (j=0; j<_PB_N; j++) {
                  A[i][j] = A[i][j] + u1[i] * v1[j] + u2[i] * v2[j];
              }
          }

        But note: the original code uses _PB_N. We must use that.

        Similarly for the others.

     However, note: the original code uses:

          for (i = 0; i < _PB_N; i++)
            for (j = 0; j < _PB_N; j++)

        We are keeping the same.

     We'll write the code accordingly.

     But note: the fused loop and the fourth loop use a temporary variable. We must declare it inside the parallel region? 
        We can declare it inside the i loop, so it is private to the thread.

     We'll do:

        #pragma omp parallel for default(none) shared(A, beta, y, z, x, _PB_N) private(i, j, tmp)
        for (i=0; i<_PB_N; i++) {
            double tmp = 0.0;
            for (j=0; j<_PB_N; j++) {
                tmp += beta * A[j][i] * y[j];
            }
            x[i] += tmp + z[i];
        }

        Similarly for the fourth loop.

     However, note: the variable 'tmp' is declared inside the loop, so it is private by default? 
        But we list it in the private clause to be safe.

     Alternatively, we can declare it inside the parallel region? 

        #pragma omp parallel
        {
            #pragma omp for
            for (i=0; i<_PB_N; i++) {
                double tmp = 0.0;
                ...
            }
        }

        But that would require a parallel region for each loop. We are using parallel for, so we can do:

          #pragma omp parallel for private(tmp, j) ... 

        and then declare tmp inside the loop? 

        Actually, we can declare tmp inside the loop and it will be private because it is declared inside the parallel region. 

        But the OpenMP standard says: variables declared inside the loop body are private. 

        So we can do:

          #pragma omp parallel for default(none) shared(A, beta, y, z, x, _PB_N) private(j)
          for (i=0; i<_PB_N; i++) {
              double tmp = 0.0;
              for (j=0; j<_PB_N; j++) {
                  tmp += beta * A[j][i] * y[j];
              }
              x[i] += tmp + z[i];
          }

        This is safe.

     Similarly for the fourth loop.

     We'll do the same for the first loop? It doesn't need a temporary, but we can use the same pattern.

     Step 4: SIMD

        We might add simd to the inner loops? 

        For the first loop: the inner j loop is a simple assignment and can be vectorized.

        For the fused loop: the inner j loop is a reduction (tmp += ...) and can be vectorized.

        For the fourth loop: the inner j loop is a reduction and can be vectorized.

        We can add:

          #pragma omp simd

        to the inner j loops? 

        But note: the inner j loop in the first loop: 
            A[i][j] = ... 
          This is a store and might not be efficiently vectorized? But it is a contiguous store.

        We'll try:

          for (i=0; i<_PB_N; i++) {
              #pragma omp simd
              for (j=0; j<_PB_N; j++) {
                  A[i][j] = A[i][j] + u1[i] * v1[j] + u2[i] * v2[j];
              }
          }

        Similarly for the others.

        However, the fused loop and the fourth loop have a reduction on a scalar. We can use:

          #pragma omp simd reduction(+:tmp)

        But note: the reduction is on a scalar that is local to the thread. 

        We'll do:

          for (i=0; i<_PB_N; i++) {
              double tmp = 0.0;
              #pragma omp simd reduction(+:tmp)
              for (j=0; j<_PB_N; j++) {
                  tmp += beta * A[j][i] * y[j];
              }
              x[i] += tmp + z[i];
          }

        And for the fourth loop:

          for (i=0; i<_PB_N; i++) {
              double tmp = 0.0;
              #pragma omp simd reduction(+:tmp)
              for (j=0; j<_PB_N; j++) {
                  tmp += alpha * A[i][j] * x[j];
              }
              w[i] += tmp;
          }

     We'll add these.

     Step 5: Final code structure

        We'll have:

          #pragma scop

          // First loop
          #pragma omp parallel for default(none) shared(A, u1, v1, u2, v2, _PB_N) private(i, j)
          for (i=0; i<_PB_N; i++) {
              #pragma omp simd
              for (j=0; j<_PB_N; j++) {
                  A[i][j] = A[i][j] + u1[i] * v1[j] + u2[i] * v2[j];
              }
          }

          // Fused second and third loops
          #pragma omp parallel for default(none) shared(A, beta, y, z, x, _PB_N) private(i, j, tmp)
          for (i=0; i<_PB_N; i++) {
              double tmp = 0.0;
              #pragma omp simd reduction(+:tmp)
              for (j=0; j<_PB_N; j++) {
                  tmp += beta * A[j][i] * y[j];
              }
              x[i] += tmp + z[i];
          }

          // Fourth loop
          #pragma omp parallel for default(none) shared(A, alpha, x, w, _PB_N) private(i, j, tmp)
          for (i=0; i<_PB_N; i++) {
              double tmp = 0.0;
              #pragma omp simd reduction(+:tmp)
              for (j=0; j<_PB_N; j++) {
                  tmp += alpha * A[i][j] * x[j];
              }
              w[i] += tmp;
          }

          #pragma endscop

     We must be cautious: the variable 'tmp' in the fused and fourth loops is declared inside the i loop, so it is private to the thread. 
     We list it in the private clause for clarity, but it is not necessary because it is declared inside the loop.

     However, the OpenMP standard requires that variables in the private clause are declared in the outer scope? 
        Actually, we can list it and it will be treated as private. But if we don't list it, it is still private because it is declared inside the parallel region.

     We'll list it to be explicit.

     But note: the first loop does not have a 'tmp'. We only have i and j.

     We'll write the code accordingly.

     However, note: the original code uses _PB_N. We must use that.

     Let's write the code.
#pragma scop

  #pragma omp parallel for default(none) shared(A, u1, v1, u2, v2, _PB_N) private(i, j)
  for (i = 0; i < _PB_N; i++) {
    #pragma omp simd
    for (j = 0; j < _PB_N; j++) {
      A[i][j] = A[i][j] + u1[i] * v1[j] + u2[i] * v2[j];
    }
  }

  #pragma omp parallel for default(none) shared(A, beta, y, z, x, _PB_N) private(i, j, tmp)
  for (i = 0; i < _PB_N; i++) {
    double tmp = 0.0;
    #pragma omp simd reduction(+:tmp)
    for (j = 0; j < _PB_N; j++) {
      tmp += beta * A[j][i] * y[j];
    }
    x[i] += tmp + z[i];
  }

  #pragma omp parallel for default(none) shared(A, alpha, x, w, _PB_N) private(i, j, tmp)
  for (i = 0; i < _PB_N; i++) {
    double tmp = 0.0;
    #pragma omp simd reduction(+:tmp)
    for (j = 0; j < _PB_N; j++) {
      tmp += alpha * A[i][j] * x[j];
    }
    w[i] += tmp;
  }

#pragma endscop
#pragma endscop
}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  DATA_TYPE beta;
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);
  POLYBENCH_1D_ARRAY_DECL(u1, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(v1, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(u2, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(v2, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(w, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(x, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(y, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(z, DATA_TYPE, N, n);


  /* Initialize array(s). */
  init_array (n, &alpha, &beta,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(u1),
	      POLYBENCH_ARRAY(v1),
	      POLYBENCH_ARRAY(u2),
	      POLYBENCH_ARRAY(v2),
	      POLYBENCH_ARRAY(w),
	      POLYBENCH_ARRAY(x),
	      POLYBENCH_ARRAY(y),
	      POLYBENCH_ARRAY(z));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_gemver (n, alpha, beta,
		 POLYBENCH_ARRAY(A),
		 POLYBENCH_ARRAY(u1),
		 POLYBENCH_ARRAY(v1),
		 POLYBENCH_ARRAY(u2),
		 POLYBENCH_ARRAY(v2),
		 POLYBENCH_ARRAY(w),
		 POLYBENCH_ARRAY(x),
		 POLYBENCH_ARRAY(y),
		 POLYBENCH_ARRAY(z));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(w)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(u1);
  POLYBENCH_FREE_ARRAY(v1);
  POLYBENCH_FREE_ARRAY(u2);
  POLYBENCH_FREE_ARRAY(v2);
  POLYBENCH_FREE_ARRAY(w);
  POLYBENCH_FREE_ARRAY(x);
  POLYBENCH_FREE_ARRAY(y);
  POLYBENCH_FREE_ARRAY(z);

  return 0;
}

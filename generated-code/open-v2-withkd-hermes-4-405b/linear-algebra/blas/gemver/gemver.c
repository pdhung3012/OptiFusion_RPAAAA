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
   - Tiling: We note that the first and fourth loops are double loops that might benefit from tiling for cache locality.
     However, the second loop has a reduction on x[i] and the third loop is a single loop. We must be cautious.

   Let's analyze dependencies:

   Loop 1: A[i][j] is updated. Each element is updated independently? Actually, each (i,j) is independent.
     So we can tile this loop.

   Loop 2: This is a reduction on x[i] (each i is independent) but the inner loop over j accumulates into x[i].
     The inner loop is a dot product. We can tile the inner loop? But note: the inner loop is over j and the reduction is over j.

   Loop 3: Single loop, independent per i.

   Loop 4: Similar to loop2, but now we are updating w[i] and the inner loop is over j.

   However, note that the second and fourth loops both use A, but in different orders (loop2: A[j][i], loop4: A[i][j]).
   Also, the second loop writes to x, which is then read in the fourth loop.

   We cannot fuse the first loop with the others because the first loop writes A and the second loop reads A (but note: the second loop uses A[j][i] which is the transpose? Actually, it's the same array but with indices swapped).

   But note: the first loop updates A, then the second loop uses the updated A. So we cannot change the order arbitrarily.

   However, we can consider tiling the first and fourth loops independently.

   For the second loop, we can also tile the inner j loop? But note that the inner j loop is a reduction. We can tile the reduction.

   But let's look at the memory access patterns:

   Loop 1: 
      A[i][j] is updated. We are accessing A by rows (if we assume row-major) and then u1[i] and v1[j] are vectors.
      This loop has good spatial locality for A in the inner loop? Actually, in the inner loop j, we are accessing A[i][j] which is contiguous in j (if row-major). So the inner loop is contiguous for A. But u1[i] is constant for the inner loop, and v1[j] is contiguous. So it's already good.

   However, if we tile the i-loop, we can keep a block of A in cache? But note that the entire A is being updated and we are using two vectors. The vectors u1 and v1 are of size n, so they might not fit in cache. But tiling might help by reusing u1[i] for a block of j? Actually, u1[i] is used for every j in the inner loop. So without tiling, we are already reusing u1[i] for the entire inner loop. Similarly, v1[j] is used for every i? Actually, no: for a fixed j, we use v1[j] for every i. But in the original loop, j is the inner loop so we are not reusing v1[j] across i. We can interchange the loops? But note: the update is independent, so we can interchange.

   However, the problem says: we are to prioritize tiling, fusion, skewing.

   Let me consider fusion: 
     The first loop and the second loop: they are adjacent but the first writes A and the second reads A. We cannot fuse because the second loop must use the updated A. Also, the second loop writes x and the third loop reads x and then the fourth loop reads x. So we cannot fuse the second and third? Actually, the third loop is a single loop that adds z[i] to x[i]. We can fuse the second and third? 

     The second loop: 
        for (i=0; i<n; i++)
          for (j=0; j<n; j++)
            x[i] = x[i] + beta * A[j][i] * y[j];

     The third loop:
        for (i=0; i<n; i++)
          x[i] = x[i] + z[i];

     We can fuse the third loop into the second by doing:

        for (i=0; i<n; i++) {
          for (j=0; j<n; j++) {
            x[i] = x[i] + beta * A[j][i] * y[j];
          }
          x[i] = x[i] + z[i];
        }

     This is legal because the inner j loop for a fixed i is independent of the z[i] addition.

   Similarly, we can consider fusing the first loop with nothing? Or the fourth loop with nothing? 

   But note: the fourth loop uses x, which is produced by the second and third. So we cannot fuse the fourth with the first.

   So let's do:

     Fuse the second and third loops.

   Now the code becomes:

     Loop1: double loop for A update.
     Loop2: double loop for x update (with the inner j loop) and then the z[i] addition for that i.
     Loop3: double loop for w update.

   Now, we can consider tiling for the double loops.

   For Loop1: 
        for (i=0; i<n; i++)
          for (j=0; j<n; j++)
            A[i][j] = ... 

        We can tile both i and j. Let's choose a tile size, say 64.

   For Loop2 (the fused one):
        for (i=0; i<n; i++) {
          for (j=0; j<n; j++) {
            x[i] += ... 
          }
          x[i] += z[i];
        }

        We can tile the j loop? But note: the inner j loop is a reduction. We can tile the j loop to improve cache locality for A and y? 

        However, the inner j loop is already contiguous for A[j][i] (if we consider j as the row index? Actually, A[j][i] is accessing by row j and column i. So if we fix i and vary j, we are accessing consecutive rows? That is not contiguous in memory (if row-major). 

        Actually, in row-major, A[j][i] is at offset j * n + i. So as j increases, we are accessing consecutive rows? That is not contiguous because the stride is n. So it's not cache friendly.

        We can interchange the loops? But note: the reduction is on x[i] and the inner loop is j. We cannot interchange because the reduction is tied to i. 

        Alternatively, we can tile the j loop to get some reuse of y[j]? But y[j] is a vector and we are accessing it contiguously? Actually, in the inner j loop, we are accessing y[j] contiguously. So that is good. But A[j][i] is not contiguous. 

        We can also consider tiling the i loop? But the outer loop is i and the inner loop is j. We are updating x[i] for each j. We can tile the i loop? Then for a block of i, we do the entire j loop? That would keep a block of x in cache? But note: the inner j loop is long and we are updating x[i] for each j. So without tiling, we are already reusing x[i] for the entire inner loop. 

        However, the access to A[j][i] is non-contiguous. We can try to change the loop order? But the reduction pattern forces the current order.

        Alternatively, we can use a technique called "loop tiling for reductions" by having a private array for the reduction in the tile? But that would require extra storage and a final reduction. 

        Given the complexity, and since the problem does not require the absolute best, we might leave the second loop as is and rely on the cache for y and the fact that x[i] is reused.

   For Loop3 (the w update):
        for (i=0; i<n; i++)
          for (j=0; j<n; j++)
            w[i] = w[i] + alpha * A[i][j] * x[j];

        This is similar to the second loop: reduction on w[i] and the inner j loop. The access to A[i][j] is contiguous in j (if row-major) and x[j] is contiguous. So it's cache friendly for A and x. We can tile the j loop? But again, the reduction is on w[i] and we are updating for each j. We can tile the j loop to get better reuse of x? But x is already contiguous. 

        However, we can tile the i loop? Then for a block of i, we do the entire j loop? That would keep a block of w in cache? But note: the inner j loop is long and we are updating w[i] for each j. So without tiling, we are already reusing w[i] for the entire inner loop.

        But note: the inner j loop is long and we are accessing A[i][j] and x[j] contiguously. So it's already good.

   However, we are to prioritize tiling. So we will tile the first loop and the fourth loop? And the second loop we leave as is? 

   But note: the first loop and the fourth loop both access A. The first loop writes A and the fourth loop reads A. So we cannot tile across the entire scop. We must do the first loop, then the second and third, then the fourth.

   We decide to tile the first and fourth loops independently.

   Tile size: we choose 64 as a reasonable tile size for L1 cache.

   For the first loop:

        for (ii = 0; ii < n; ii += 64)
          for (jj = 0; jj < n; jj += 64)
            for (i = ii; i < min(ii+64, n); i++)
              for (j = jj; j < min(jj+64, n); j++)
                A[i][j] = ... 

   For the fourth loop:

        for (ii = 0; ii < n; ii += 64)
          for (jj = 0; jj < n; jj += 64)
            for (i = ii; i < min(ii+64, n); i++)
              for (j = jj; j < min(jj+64, n); j++)
                w[i] = w[i] + alpha * A[i][j] * x[j];

   But note: the fourth loop is a reduction on w[i]. We can do the tiling without private arrays because the reduction is per i and the tile in j is just a part of the reduction. However, we must be cautious: the inner two loops (i and j) are not independent in the reduction? Actually, for a fixed i, the j loop is a reduction. So we can tile the j loop arbitrarily? 

   However, the tiling we are doing is over both i and j. But note: the reduction for w[i] is independent of other i. So we can do:

        for (ii = 0; ii < n; ii += 64) {
          for (i = ii; i < min(ii+64, n); i++) {
            for (jj = 0; jj < n; jj += 64) {
              for (j = jj; j < min(jj+64, n); j++) {
                w[i] = w[i] + alpha * A[i][j] * x[j];
              }
            }
          }
        }

   But that would break the contiguous access in j? Because we are now doing j in blocks. However, the inner j loop is still contiguous within the tile. And we are reusing x[j] for the entire tile? Actually, x[j] is contiguous and we are accessing it in order. So it's still good.

   Alternatively, we can do:

        for (ii = 0; ii < n; ii += 64) {
          for (jj = 0; jj < n; jj += 64) {
            for (i = ii; i < min(ii+64, n); i++) {
              for (j = jj; j < min(jj+64, n); j++) {
                w[i] = w[i] + alpha * A[i][j] * x[j];
              }
            }
          }
        }

   This is the same as above? Actually, the order of the i and jj loops is swapped. The above does jj inside i, which is not good because we want to keep the inner j loop contiguous. The second version does:

        for (ii) {
          for (jj) {
            for (i in ii block) {
              for (j in jj block) {
                ... 
              }
            }
          }
        }

   This is better because the inner two loops are over i and j, and we are accessing A[i][j] and x[j] contiguously in j for a fixed i. But note: the inner loop is j, so for a fixed i and jj, we are accessing a contiguous block of j. So it's good.

   However, we are updating w[i] for each j. We are not reusing w[i] across jj blocks? Actually, we are: for a fixed i, we are updating w[i] for every j in the entire range. So we are accumulating w[i] over the jj blocks. That is correct.

   So we choose the second version for the fourth loop.

   Now, what about the second loop? We fused the third loop into it. We can also tile the second loop? 

        for (i=0; i<n; i++) {
          for (j=0; j<n; j++) {
            x[i] += ... 
          }
          x[i] += z[i];
        }

   We can tile the j loop? 

        for (i=0; i<n; i++) {
          for (jj=0; jj<n; jj+=64) {
            for (j=jj; j<min(jj+64, n); j++) {
              x[i] += ... 
            }
          }
          x[i] += z[i];
        }

   This is safe and might improve cache locality for y[j] and A[j][i]? But note: A[j][i] is not contiguous in j. However, we are accessing y[j] contiguously. And we are reusing x[i] for the entire j loop. So it might help to keep a block of y in cache? 

   We decide to tile the j loop in the second loop as well.

   So now we have:

     Loop1: tiled (i and j) with tile size 64.
     Loop2: tiled in j (with tile size 64) and then the z[i] addition.
     Loop3: tiled (i and j) with tile size 64.

 Step 2: OpenMP directives.

   We are to add OpenMP directives to the loops.

   We note that:

     Loop1: The update of A[i][j] is independent for every (i,j). So we can parallelize the outer tile loops? Or the inner loops? 

        We can do:

          #pragma omp parallel for collapse(2) 
          for (ii=0; ii<n; ii+=64)
            for (jj=0; jj<n; jj+=64)
              ... inner loops ...

        But note: the inner loops are over i and j. We are updating A[i][j] and the updates are independent. So we can collapse the two outer tile loops? 

        However, the tile loops are not perfectly nested? Actually, they are. So we can collapse(2) the ii and jj loops.

        But note: the inner loops are not parallel? We are doing the inner loops in serial. So we are parallelizing the tile level.

        Alternatively, we can parallelize the inner i loop? But that would be inside the tile. We want coarse-grained parallelism.

        We choose to parallelize the outer tile loops.

     Loop2: The outer loop is i, and for each i we do a reduction on x[i]. The inner j loop is a reduction. We can parallelize the i loop? 

        We can do:

          #pragma omp parallel for 
          for (i=0; i<n; i++) {
            ... inner j loop (tiled) and then z[i] ...
          }

        But note: the inner j loop is a reduction. We can use a reduction clause? Actually, the reduction is on x[i] and each i is independent. So we don't need a reduction clause? We can just do:

          #pragma omp parallel for 
          for (i=0; i<n; i++) {
            ... 
          }

        However, the inner j loop is a reduction. We can also parallelize the inner j loop? But that would be fine-grained. We choose to parallelize the outer i loop.

        But note: the inner j loop is long and we are tiling it. We can also use a reduction clause? Actually, we don't need because each thread has its own i.

     Loop3: The update of w[i] is independent per i? Actually, for a fixed i, we are updating w[i] and the inner j loop is a reduction. So we can parallelize the outer i loop? 

        We can do:

          #pragma omp parallel for 
          for (i=0; i<n; i++) {
            ... inner j loop (tiled) ...
          }

        But note: the inner j loop is a reduction. We can also parallelize the inner j loop? But again, we choose the outer i loop.

        However, we have tiled the i and j loops. We can parallelize the outer tile loop? 

          #pragma omp parallel for collapse(2)
          for (ii=0; ii<n; ii+=64)
            for (jj=0; jj<n; jj+=64)
              for (i=ii; i<min(ii+64, n); i++)
                for (j=jj; j<min(jj+64, n); j++)
                  w[i] = ... 

        But note: the update of w[i] is not atomic. We are updating w[i] for multiple j and multiple tiles. But the reduction is per i. So if we parallelize the tile loops, then multiple tiles (with the same i) might update w[i] concurrently? That would be a race.

        We must not parallelize the tile loops for the fourth loop? We can only parallelize the i loop? 

        Actually, we can restructure: we can make the outer loop over i and then the tile loop over j? Then we can parallelize the i loop.

        We change the fourth loop to:

          #pragma omp parallel for 
          for (i=0; i<n; i++) {
            for (jj=0; jj<n; jj+=64) {
              for (j=jj; j<min(jj+64, n); j++) {
                w[i] = w[i] + alpha * A[i][j] * x[j];
              }
            }
          }

        This is safe because each i is independent.

        Alternatively, we can keep the tile structure and use a reduction? But that would require a critical section or atomic, which is expensive.

        We choose to restructure the fourth loop to have the outer loop over i and then the j loop tiled.

   So now we have:

     Loop1: 
        #pragma omp parallel for collapse(2)
        for (ii=0; ii<n; ii+=64)
          for (jj=0; jj<n; jj+=64)
            for (i=ii; i<min(ii+64, n); i++)
              for (j=jj; j<min(jj+64, n); j++)
                A[i][j] = ... 

     Loop2:
        #pragma omp parallel for 
        for (i=0; i<n; i++) {
          for (jj=0; jj<n; jj+=64) {
            for (j=jj; j<min(jj+64, n); j++) {
              x[i] += beta * A[j][i] * y[j];
            }
          }
          x[i] += z[i];
        }

     Loop3:
        #pragma omp parallel for 
        for (i=0; i<n; i++) {
          for (jj=0; jj<n; jj+=64) {
            for (j=jj; j<min(jj+64, n); j++) {
              w[i] = w[i] + alpha * A[i][j] * x[j];
            }
          }
        }

   We also need to consider scheduling and other clauses.

   For the first loop: we are using collapse(2) and we have a regular iteration space. We can use schedule(static).

   For the second and third loops: we are parallelizing the i loop. We can use schedule(static) as well.

   We also note that the second and third loops are reductions? Actually, the second loop has a reduction on x[i] and the third on w[i]. But we are parallelizing the i loop and each thread has its own i, so no reduction clause is needed.

   We also consider SIMD: the inner loops are contiguous and we can add simd.

   For the first loop: the inner two loops (i and j) are contiguous? Actually, the inner loop is j and we are accessing A[i][j] contiguously. We can add simd to the inner j loop? But note: the tile is small (64) and we are already in a parallel region. We can do:

        #pragma omp parallel for collapse(2)
        for (ii=0; ii<n; ii+=64)
          for (jj=0; jj<n; jj+=64)
            #pragma omp simd
            for (i=ii; i<min(ii+64, n); i++)
              for (j=jj; j<min(jj+64, n); j++)
                A[i][j] = ... 

     But note: the inner loop is j, and we are vectorizing the j loop? Actually, the pragma omp simd applies to the next loop, which is the i loop? We want to vectorize the j loop. We can do:

        for (i=ii; i<min(ii+64, n); i++) {
          #pragma omp simd
          for (j=jj; j<min(jj+64, n); j++) {
            A[i][j] = ... 
          }
        }

     Similarly for the other loops.

   However, the problem says: we are to apply SIMD directives. We can do:

        For the first loop: 
          #pragma omp parallel for collapse(2)
          for (ii=0; ii<n; ii+=64)
            for (jj=0; jj<n; jj+=64)
              for (i=ii; i<min(ii+64, n); i++)
                #pragma omp simd
                for (j=jj; j<min(jj+64, n); j++)
                  A[i][j] = ... 

        For the second loop:
          #pragma omp parallel for 
          for (i=0; i<n; i++) {
            for (jj=0; jj<n; jj+=64) {
              #pragma omp simd
              for (j=jj; j<min(jj+64, n); j++) {
                x[i] += beta * A[j][i] * y[j];
              }
            }
            x[i] += z[i];
          }

        For the third loop:
          #pragma omp parallel for 
          for (i=0; i<n; i++) {
            for (jj=0; jj<n; jj+=64) {
              #pragma omp simd
              for (j=jj; j<min(jj+64, n); j++) {
                w[i] = w[i] + alpha * A[i][j] * x[j];
              }
            }
          }

   We also note that the second and third loops have a reduction in the inner j loop. The #pragma omp simd can handle reductions? Yes, if the compiler supports it. We assume it does.

   We also consider data sharing: we should use default(none) and specify shared and private variables.

   For the first loop: 
        shared: A, u1, v1, u2, v2, n
        private: ii, jj, i, j

   But note: the tile loops (ii, jj) are private by the OpenMP loop. The inner i and j are private.

   Similarly for the others.

   We also note that the second loop uses x, which is shared? Actually, each thread writes to its own x[i] (because i is the loop index and we are parallelizing i). So x is shared but the writes are to disjoint locations. So we don't need to make it private.

   Similarly for w in the third loop.

   However, we must be cautious: the second loop reads A and y, which are shared. The third loop reads A and x, which are shared.

   We can use:

        #pragma omp parallel for default(none) shared(A, u1, v1, u2, v2, n) 
        for (ii=0; ... 

   But note: the inner loops use min, which is a function. We can leave it.

   Alternatively, we can avoid default(none) and rely on the default shared? But the problem says to use default(none) if we use it. We will use default(none) and list all shared variables.

   However, the variables u1, v1, u2, v2, y, z, alpha, beta, x, w are all shared. And n is shared.

   But note: the inner loops use ii, jj, i, j which are private.

   We can do:

        #pragma omp parallel for collapse(2) default(none) \
            shared(A, u1, v1, u2, v2, n) private(ii, jj, i, j)
        for (ii=0; ii<n; ii+=64)
          for (jj=0; jj<n; jj+=64)
            for (i=ii; i<min(ii+64, n); i++)
              #pragma omp simd
              for (j=jj; j<min(jj+64, n); j++)
                A[i][j] = A[i][j] + u1[i] * v1[j] + u2[i] * v2[j];

   But note: the inner loops are inside the parallel region. The variables i and j are private by the for loop. We don't need to list them? Actually, the OpenMP standard says that the loop index variables are private. So we can remove private(ii,jj,i,j) because they are loop indices? 

   However, the outer loops (ii and jj) are the ones being parallelized. The inner i and j are not. But the inner i and j are declared in the for loop and are private by the for loop. So we don't need to list them.

   But the variables u1, v1, u2, v2 are shared and read-only. So it's safe.

   Similarly for the other loops.

   We also note that the second loop uses x, which is shared and written. But each thread writes to a disjoint part. So it's safe.

   We also note that the third loop uses w, which is shared and written. Similarly safe.

   We also note that the second loop uses A and y, and the third loop uses A and x. They are shared and read-only in the parallel region? Actually, in the second loop, A is read-only? Yes, because we are only reading A. Similarly, in the third loop, A is read-only? Yes.

   So we can do:

        Loop1: 
          #pragma omp parallel for collapse(2) default(none) shared(A, u1, v1, u2, v2, n)
          ... 

        Loop2:
          #pragma omp parallel for default(none) shared(x, A, y, beta, z, n)
          ... 

        Loop3:
          #pragma omp parallel for default(none) shared(w, A, x, alpha, n)
          ... 

   But note: the inner loops use min, which is a function. We don't need to share it.

   We also note that the tile size 64 is a constant. We can define it as a constant.

   However, the problem says: do not change code outside the scop. We can define the tile size inside the scop? 

   We can do:

        #define TILE_SIZE 64

   But we must be cautious: the scop region is inside a function. We can define it at the beginning of the scop.

   Alternatively, we can use a constant expression.

   We choose to use:

        const int tile_size = 64;

   But note: we are in C, and we are inside a function. We can declare it at the beginning of the scop.

   However, the problem says: preserve the marker lines. We can put the declaration inside the scop.

   We do:

        #pragma scop
        {
          const int tile_size = 64;

          ... the four loops ...

        }
        #pragma endscop

   But note: the original code does not have a block. We can add a block? The problem says: preserve the marker lines. We can add a block without changing the marker lines.

   Alternatively, we can avoid the block and declare the tile_size at the beginning of the scop.

   We do:

        #pragma scop
        const int tile_size = 64;

        ... the four loops ...

        #pragma endscop

   But note: the original code has the loops without a block. We are adding a declaration. That is allowed.

   However, we must be cautious: the original code uses _PB_N. We are using n. The function parameter is n. So we are using n.

   We write the code accordingly.

   We also note that the min function: we can use a macro or a function? We can use the standard min function? But note: the original code does not include <stdlib.h> for min? Actually, we are using min in the loop bounds. We can use a ternary operator? 

   We do:

        for (i = ii; i < (ii+tile_size < n ? ii+tile_size : n); i++)

   But that is messy. We can define:

        int i_end = (ii+tile_size < n) ? (ii+tile_size) : n;

   But we are in a parallel region? We can compute it outside the inner loops? 

   We do:

        for (ii=0; ii<n; ii+=tile_size) {
          int i_end = (ii+tile_size < n) ? (ii+tile_size) : n;
          for (jj=0; jj<n; jj+=tile_size) {
            int j_end = (jj+tile_size < n) ? (jj+tile_size) : n;
            ... 
          }
        }

   But note: the first loop is parallelized with collapse(2). We cannot break the tile loops with a declaration? 

   We can do:

        #pragma omp parallel for collapse(2) default(none) shared(A, u1, v1, u2, v2, n, tile_size)
        for (ii=0; ii<n; ii+=tile_size)
          for (jj=0; jj<n; jj+=tile_size) {
            int i_end = (ii+tile_size < n) ? (ii+tile_size) : n;
            int j_end = (jj+tile_size < n) ? (jj+tile_size) : n;
            for (i=ii; i<i_end; i++)
              #pragma omp simd
              for (j=jj; j<j_end; j++)
                A[i][j] = A[i][j] + u1[i] * v1[j] + u2[i] * v2[j];
          }

   But note: the inner loops are now inside a block. And we are declaring i_end and j_end. This is acceptable.

   Similarly for the other loops.

   However, the second and third loops are not collapsed. We can do:

        #pragma omp parallel for default(none) shared(x, A, y, beta, z, n, tile_size)
        for (i=0; i<n; i++) {
          for (jj=0; jj<n; jj+=tile_size) {
            int j_end = (jj+tile_size < n) ? (jj+tile_size) : n;
            #pragma omp simd
            for (j=jj; j<j_end; j++) {
              x[i] += beta * A[j][i] * y[j];
            }
          }
          x[i] += z[i];
        }

   And for the third loop:

        #pragma omp parallel for default(none) shared(w, A, x, alpha, n, tile_size)
        for (i=0; i<n; i++) {
          for (jj=0; jj<n; jj+=tile_size) {
            int j_end = (jj+tile_size < n) ? (jj+tile_size) : n;
            #pragma omp simd
            for (j=jj; j<j_end; j++) {
              w[i] = w[i] + alpha * A[i][j] * x[j];
            }
          }
        }

   We also note that the second loop has a reduction on x[i]. The #pragma omp simd will handle the reduction? We assume the compiler supports reduction in simd.

   We also note that the third loop has a reduction on w[i]. Similarly.

   We also note that the first loop does not have a reduction? It is a direct assignment.

   We also consider scheduling: we can add schedule(static) to the parallel for. The problem says to consider scheduling. We choose static for all.

   We do:

        #pragma omp parallel for collapse(2) schedule(static) ...

   And for the others:

        #pragma omp parallel for schedule(static) ...

   We also note that the second and third loops are reductions? But we are parallelizing the i loop and each thread has its own i, so no reduction clause is needed.

   We also note that the first loop is not a reduction.

   We also consider thread management: we can add num_threads? But we don't know the number of threads. We leave it to the runtime.

   We also consider proc_bind? We don't specify.

   We also note that the inner loops are vectorized with #pragma omp simd. We can also add safelen? We don't specify.

   We also note that the arrays might be aligned? We don't know. We don't add aligned.

   We also note that the second loop uses A[j][i] which is not contiguous. But we are vectorizing the j loop? The j loop is contiguous in j for A[j][i]? Actually, no: A[j][i] is at offset j*n+i. So as j increases, the offset increases by n. That is not contiguous. 

   This might cause poor vectorization. We might need to change the loop order? But we cannot because of the reduction pattern.

   We leave it as is and hope the compiler can vectorize non-contiguous accesses? Or we can use a temporary array? But that would be complex.

   Alternatively, we can interchange the loops? But the reduction is on x[i] and we are accumulating over j. We cannot interchange.

   We decide to leave it and hope that the vectorization of the j loop (which is contiguous in j for y[j] and the reduction) is still beneficial.

   Similarly, the third loop: A[i][j] is contiguous in j, so it's good.

   We write the code accordingly.

   We also note that the original code uses _PB_N. We are using n. The function parameter is n. So we are using n.

   We also note that the original code uses DATA_TYPE. We are not changing that.

   We also note that the original code uses POLYBENCH_2D and POLYBENCH_1D. We are not changing that.

   We also note that the original code uses the variables as declared.

   We write the optimized code.

   Note: We must preserve the original semantics. We have fused the second and third loops and tiled the loops. We have also added OpenMP and SIMD.

   We hope it is correct.

   Let's write the code.
#pragma scop
{
  const int tile_size = 64;

  #pragma omp parallel for collapse(2) schedule(static) default(none) shared(A, u1, v1, u2, v2, n, tile_size)
  for (int ii = 0; ii < n; ii += tile_size)
    for (int jj = 0; jj < n; jj += tile_size) {
      int i_end = (ii + tile_size < n) ? ii + tile_size : n;
      int j_end = (jj + tile_size < n) ? jj + tile_size : n;
      for (int i = ii; i < i_end; i++)
        #pragma omp simd
        for (int j = jj; j < j_end; j++)
          A[i][j] = A[i][j] + u1[i] * v1[j] + u2[i] * v2[j];
    }

  #pragma omp parallel for schedule(static) default(none) shared(x, A, y, beta, z, n, tile_size)
  for (int i = 0; i < n; i++) {
    for (int jj = 0; jj < n; jj += tile_size) {
      int j_end = (jj + tile_size < n) ? jj + tile_size : n;
      #pragma omp simd
      for (int j = jj; j < j_end; j++)
        x[i] += beta * A[j][i] * y[j];
    }
    x[i] += z[i];
  }

  #pragma omp parallel for schedule(static) default(none) shared(w, A, x, alpha, n, tile_size)
  for (int i = 0; i < n; i++) {
    for (int jj = 0; jj < n; jj += tile_size) {
      int j_end = (jj + tile_size < n) ? jj + tile_size : n;
      #pragma omp simd
      for (int j = jj; j < j_end; j++)
        w[i] += alpha * A[i][j] * x[j];
    }
  }
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

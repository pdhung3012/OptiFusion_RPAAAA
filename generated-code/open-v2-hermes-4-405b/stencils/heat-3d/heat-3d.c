/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* heat-3d.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "heat-3d.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_3D(A,N,N,N,n,n,n),
		 DATA_TYPE POLYBENCH_3D(B,N,N,N,n,n,n))
{
  int i, j, k;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++)
        A[i][j][k] = B[i][j][k] = (DATA_TYPE) (i + j + (n-k))* 10 / (n);
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_3D(A,N,N,N,n,n,n))

{
  int i, j, k;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("A");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++) {
         if ((i * n * n + j * n + k) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
         fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j][k]);
      }
  POLYBENCH_DUMP_END("A");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_heat_3d(int tsteps,
		      int n,
		      DATA_TYPE POLYBENCH_3D(A,N,N,N,n,n,n),
		      DATA_TYPE POLYBENCH_3D(B,N,N,N,n,n,n))
{
  int t, i, j, k;

#pragma scop
We are given a code for a 3D heat equation solver. The main computational kernel is within the scop markers.
 The kernel has two nested loops over time steps (t) and then three spatial dimensions (i, j, k). 
 There are two main triple-nested loops: one updating B from A, and then one updating A from B.

 We are to apply loop transformations and OpenMP directives.

 Step 1: Loop Transformations (in order: tiling, fusion, skewing)

 1. Loop Tiling: 
    The inner three loops (i, j, k) are over a 3D grid. We can tile these loops to improve cache locality.
    However, note that the two triple-nested loops are adjacent and both are inside the time loop.
    We cannot fuse them because they are separated by the time loop and also because they update different arrays (B then A) and the second loop uses the updated B to update A.

    But note: the two triple-nested loops are independent in terms of the time step? Actually, they are in the same time step and the second loop uses the B that was just computed in the first loop. So we cannot fuse the two triple-nested loops because they are not adjacent in the sense of having the same outer loop without the time loop.

    However, we can tile each of the triple-nested loops individually.

    Let's consider the first triple-nested loop (updating B). We can tile the i, j, k loops. Similarly for the second.

    But note: the stencil operation in each loop accesses neighbors in all three dimensions. Tiling in 3D might be beneficial.

    We choose a tile size that fits in cache. Since the problem size is N (which is typically large, e.g., 256 or more) and the stencil is 7-point, we can choose a tile size that is a multiple of the cache line and also divides the inner dimensions.

    However, the problem does not specify N, so we choose a tile size that is a power of two and not too large. Let's choose 32 for the tile size in each dimension? But note: the inner loops run from 1 to N-2. We have to adjust the tile boundaries to avoid going out of bounds.

    Alternatively, we can tile the two inner loops (j and k) and leave i as the outer loop? Or tile all three? Since the stencil is symmetric, we can tile in all three dimensions.

    We'll do:

        for (t = 1; t <= TSTEPS; t++) {
            // Tile the first triple loop
            for (i = 1; i < _PB_N-1; i += TI) {
                for (j = 1; j < _PB_N-1; j += TJ) {
                    for (k = 1; k < _PB_N-1; k += TK) {
                        // Inner loops for the tile
                        for (ii = i; ii < min(i+TI, _PB_N-1); ii++) {
                            for (jj = j; jj < min(j+TJ, _PB_N-1); jj++) {
                                for (kk = k; kk < min(k+TK, _PB_N-1); kk++) {
                                    // original body for B
                                }
                            }
                        }
                    }
                }
            }

            // Similarly for the second triple loop updating A
            for (i = 1; i < _PB_N-1; i += TI) {
                for (j = 1; j < _PB_N-1; j += TJ) {
                    for (k = 1; k < _PB_N-1; k += TK) {
                        for (ii = i; ii < min(i+TI, _PB_N-1); ii++) {
                            for (jj = j; jj < min(j+TJ, _PB_N-1); jj++) {
                                for (kk = k; kk < min(k+TK, _PB_N-1); kk++) {
                                    // original body for A
                                }
                            }
                        }
                    }
                }
            }
        }

    But note: the tile size must be chosen so that the tile fits in cache. We don't know the exact cache size, but we can choose a conservative tile size. Let's set TI, TJ, TK to 32. However, we must be cautious because the inner tile loops will have bounds that are not known at compile time (due to min). This might prevent some optimizations.

    Alternatively, we can use a fixed tile size and adjust the outer loops to cover the entire range. We can also use OpenMP on the outer tile loops.

 2. Loop Fusion: 
    We cannot fuse the two triple-nested loops because they are separated by the time loop and also because the second loop uses the updated B. So fusion is not applicable.

 3. Loop Skewing: 
    The inner loops do not have carried dependencies that would prevent parallelism? Actually, the stencil is a 7-point stencil and the computation for each point is independent. So we don't need skewing.

    However, note that the two triple-nested loops are independent in the sense that the first writes B and the second reads B and writes A. But within each triple-nested loop, the iterations are independent. So we can parallelize each triple-nested loop.

 Step 2: OpenMP Directives

 We can parallelize the outer time loop? But note: the time steps are sequential. So we cannot parallelize the time loop.

 Instead, we can parallelize the spatial loops. We have two triple-nested loops per time step. We can parallelize each of them.

 However, note that the two triple-nested loops are independent? Actually, the first loop writes B and the second loop reads B and writes A. So we cannot run them in parallel. But we can run each triple-nested loop in parallel.

 We can use:

    #pragma omp parallel for collapse(2)  ... 

 But note: the inner loops are three levels. We can collapse two of them? Or we can collapse all three? However, the number of iterations might be too large and cause overhead. Alternatively, we can collapse two and leave the innermost as serial.

 However, the stencil computation is heavy enough that we want to parallelize the outer levels.

 But note: the tile loops we introduced are now the outer loops. We can parallelize the tile loops? That would be at a coarser granularity.

 Alternatively, we can avoid tiling and just use OpenMP on the original loops? But the problem says to prioritize loop transformations first.

 We have two options:

   Option 1: Without tiling, we can do:

        for (t = 1; t <= TSTEPS; t++) {
            #pragma omp parallel for collapse(2)  // for i and j, then k is inner
            for (i = 1; i < _PB_N-1; i++) {
                for (j = 1; j < _PB_N-1; j++) {
                    for (k = 1; k < _PB_N-1; k++) {
                        // update B
                    }
                }
            }

            #pragma omp parallel for collapse(2)
            for (i = 1; i < _PB_N-1; i++) {
                for (j = 1; j < _PB_N-1; j++) {
                    for (k = 1; k < _PB_N-1; k++) {
                        // update A
                    }
                }
            }
        }

   Option 2: With tiling, we can parallelize the outer tile loops:

        for (t = 1; t <= TSTEPS; t++) {
            #pragma omp parallel for collapse(2)  // for i and j tile loops
            for (i = 1; i < _PB_N-1; i += TI) {
                for (j = 1; j < _PB_N-1; j += TJ) {
                    for (k = 1; k < _PB_N-1; k += TK) {
                        // inner tile loops
                    }
                }
            }

            // similarly for the second triple loop
        }

   The advantage of tiling is cache locality. Without tiling, the original loops might have poor locality because the array accesses are to a 3D grid and the inner loop (k) is the fastest changing, which is contiguous in memory? Actually, in C, the last dimension is contiguous. So the inner loop (k) is contiguous. But the stencil accesses neighbors in i, j, k. So without tiling, the j and i loops might cause cache misses.

   Therefore, we choose to do tiling and then parallelize the outer tile loops.

   We set the tile sizes to 32 in each dimension. We define:

        #define TI 32
        #define TJ 32
        #define TK 32

   But note: the problem says to preserve the original code outside the scop. We cannot define macros inside the scop. We can use constants? Or we can use variables? But the tile size must be known at compile time for the min function? Actually, we can use:

        int ti_end = min(i+TI, _PB_N-1);

   But that would be a variable, and the inner loop bounds would be variable. This might prevent some optimizations. Alternatively, we can use:

        for (ii = i; ii < i+TI; ii++) {
            if (ii >= _PB_N-1) break;
            ... 
        }

   But that is not efficient. Instead, we can use:

        int ti_end = (i+TI < _PB_N-1) ? i+TI : _PB_N-1;

   and then:

        for (ii = i; ii < ti_end; ii++)

   This is acceptable.

   However, note that the original code uses _PB_N which is a macro. We can use that.

   We'll do:

        for (t = 1; t <= TSTEPS; t++) {
            // First triple loop: update B
            #pragma omp parallel for collapse(2)  // we collapse the i and j tile loops? But note: the k tile loop is the next. We can collapse two levels: i and j tile loops.
            for (i = 1; i < _PB_N-1; i += TI) {
                for (j = 1; j < _PB_N-1; j += TJ) {
                    for (k = 1; k < _PB_N-1; k += TK) {
                        int ti_end = (i+TI < _PB_N-1) ? i+TI : _PB_N-1;
                        int tj_end = (j+TJ < _PB_N-1) ? j+TJ : _PB_N-1;
                        int tk_end = (k+TK < _PB_N-1) ? k+TK : _PB_N-1;
                        for (ii = i; ii < ti_end; ii++) {
                            for (jj = j; jj < tj_end; jj++) {
                                for (kk = k; kk < tk_end; kk++) {
                                    // update B[ii][jj][kk]
                                }
                            }
                        }
                    }
                }
            }

            // Second triple loop: update A
            #pragma omp parallel for collapse(2)
            for (i = 1; i < _PB_N-1; i += TI) {
                for (j = 1; j < _PB_N-1; j += TJ) {
                    for (k = 1; k < _PB_N-1; k += TK) {
                        int ti_end = (i+TI < _PB_N-1) ? i+TI : _PB_N-1;
                        int tj_end = (j+TJ < _PB_N-1) ? j+TJ : _PB_N-1;
                        int tk_end = (k+TK < _PB_N-1) ? k+TK : _PB_N-1;
                        for (ii = i; ii < ti_end; ii++) {
                            for (jj = j; jj < tj_end; jj++) {
                                for (kk = k; kk < tk_end; kk++) {
                                    // update A[ii][jj][kk]
                                }
                            }
                        }
                    }
                }
            }
        }

   We collapse the two outer tile loops (i and j) to reduce the overhead of parallel regions.

   We also need to consider data sharing. The arrays A and B are shared. The loop variables are private by default in OpenMP for worksharing loops.

   We also note that the inner tile loops are now sequential. This is acceptable because the tile size is small.

   We must be cautious: the original code uses the same variable names (i, j, k) for the outer loops and then reuses them for the inner loops? Actually, we are using ii, jj, kk for the inner loops. So we are safe.

   We also note that the original code uses the same variable names for the two triple loops. We are reusing i, j, k for the tile loops in both the first and second triple loops. This is acceptable because they are in separate scopes.

 Step 3: Other OpenMP clauses

   - We can add 'default(none)' and then list the shared variables: A, B, and the loop bounds (t, _PB_N, TI, TJ, TK). But note: TI, TJ, TK are constants? We defined them as 32. We can use:

        #define TI 32
        #define TJ 32
        #define TK 32

     But we cannot define macros inside the scop. We can use enum? Or we can use const int? But the problem says to preserve the original code outside the scop. We are only allowed to change within the scop.

     We can define the tile sizes as constants inside the scop? But that would be inside the function. We can do:

        const int TI = 32;
        const int TJ = 32;
        const int TK = 32;

     However, note that the scop is inside a function. We can declare these at the beginning of the scop.

   - We can also use 'schedule(static)' which is the default, so we don't need to specify.

   - We don't need 'nowait' because there is an implicit barrier at the end of the parallel for. And we need the barrier because the second triple loop depends on the first.

   - We don't have reductions.

   - We can consider SIMD? The inner loops (the tile loops) are small and the stencil computation is not vectorized easily because of the neighbor accesses? But we can try to use SIMD on the innermost loop? However, the innermost loop (kk) is over a tile of size TK (32). We can use:

        #pragma omp simd

     on the innermost loop? But note: the innermost loop is the kk loop. We can do:

        for (kk = k; kk < tk_end; kk++) {
            // update
        }

     We can add:

        #pragma omp simd

     But we must ensure that the loop is safe for SIMD: no dependencies and the number of iterations is a multiple of the vector length? We don't know. We can use 'safelen' to specify the maximum safe vector length? But we don't know the architecture.

     Alternatively, we can rely on the compiler to vectorize. We can add:

        #pragma omp simd

     to the innermost loop of the tile.

   - We can also use 'aligned' if we know the alignment of A and B? But we don't. So we skip.

   - We can use 'proc_bind' to control thread binding? But we don't have specific requirements.

   So we decide:

        For the first triple loop (update B) and the second (update A), we do:

            #pragma omp parallel for collapse(2) default(none) \
                shared(A, B, _PB_N, TI, TJ, TK) \
                private(ii, jj, kk, ti_end, tj_end, tk_end) \
                // Note: the outer tile loops (i, j, k) are private by the worksharing construct? Actually, the worksharing loop variables are private. But we are using i, j, k as the tile loop indices. They are private by the worksharing. The inner tile loops (ii, jj, kk) and the end variables (ti_end, etc.) are private to each thread.

            But note: the tile loop indices (i, j, k) are defined in the worksharing loop and are private. The inner loops are inside the parallel region and are private to each thread.

        We also add SIMD to the innermost loop:

            for (kk = k; kk < tk_end; kk++) {
                #pragma omp simd
                ... the update ...

            }

        But wait: the update for one element is independent? Yes, because each element is updated independently. So we can vectorize the innermost loop.

        However, note that the innermost loop is over kk, and the stencil accesses in the k dimension: A[i][j][k+1] and A[i][j][k-1]. This might cause alignment issues? But the compiler can handle it.

        We can also try to vectorize the jj loop? But the jj loop is the next outer. We can collapse the two inner loops? But the stencil in j and k might be more complex.

        We decide to only vectorize the innermost loop (kk) because it is contiguous in memory.

 Implementation:

    We define the tile sizes as constants at the beginning of the scop.

    We restructure the two triple loops with tiling and then add OpenMP.

    We also note: the original code uses the same variable names (i, j, k) for the outer loops and then for the inner loops? We are changing the inner loops to use ii, jj, kk. So we are safe.

    We also note: the original code uses the same variable names for the two triple loops. We are reusing i, j, k for the tile loops in both the first and second triple loops. This is acceptable because they are in separate scopes.

    We must be cautious: the original code uses the same variable names for the time loop (t) and the spatial loops. We are not changing the time loop.

    We also note: the original code uses the same array names (A and B) and we are not changing that.

    We write the code accordingly.

    However, note: the original code uses the macro _PB_N. We use that.

    We also note: the original code uses the macro SCALAR_VAL and DATA_TYPE. We preserve.

    We also note: the original code uses the macro TSTEPS. We preserve.

    We write the optimized code within the scop.

    We also note: the original code has two triple loops. We are tiling both.

    We choose tile sizes: 32 for each dimension. We define:

        const int TI = 32;
        const int TJ = 32;
        const int TK = 32;

    at the beginning of the scop.

    We then restructure the loops.

    We also note: the original code has the time loop from t=1 to TSTEPS. We preserve.

    We also note: the original code has the spatial loops from 1 to _PB_N-1. We preserve.

    We also note: the stencil computation is the same.

    We write the code accordingly.

    We also note: the original code uses the same variable names for the two triple loops. We are reusing the same variable names for the tile loops (i, j, k) in both the first and second triple loops. This is acceptable because they are in separate scopes (each triple loop is a block).

    We also note: the inner tile loops use ii, jj, kk. We declare them inside the tile block.

    We also note: the end variables (ti_end, tj_end, tk_end) are declared inside the tile block.

    We also note: the OpenMP directives are placed immediately before the for loops they control.

    We also note: we use 'collapse(2)' to collapse the two outer tile loops (i and j). The k tile loop is the next level and we don't collapse it because we want to vectorize the innermost loop (which is the kk loop). But note: the k tile loop is the third level. We are not collapsing it. So the parallel for is over the i and j tile loops, and the k tile loop is run sequentially by each thread? Actually, no: the k tile loop is inside the j tile loop. We are collapsing the i and j tile loops, so the k tile loop is inside the parallel region and each thread will iterate over the k tile loop for the assigned i and j tiles.

    This is acceptable.

    We also note: the inner tile loops (ii, jj, kk) are run sequentially by each thread.

    We also note: we add '#pragma omp simd' to the innermost loop (kk) to suggest vectorization.

    We also note: we must declare the variables in the OpenMP clause. We use 'default(none)' and then list the shared and private variables.

    We list:

        shared: A, B, _PB_N, TI, TJ, TK
        private: i, j, k, ii, jj, kk, ti_end, tj_end, tk_end

    But note: the tile loop indices (i, j, k) are defined in the worksharing loop and are private by default. The inner loop indices (ii, jj, kk) and the end variables are declared inside the parallel region and are private to each thread.

    However, we must list them in the private clause? Actually, the OpenMP standard says that variables defined inside the parallel region are private. So we don't need to list them. But we list the ones that are declared outside? We don't have any.

    We only need to list the shared variables that are used and are not loop variables. The loop variables of the worksharing loop (i, j, k) are private by default. The inner loops are inside the parallel region and their variables are private.

    So we can do:

        #pragma omp parallel for collapse(2) default(none) shared(A, B, _PB_N, TI, TJ, TK)

    But note: the tile sizes (TI, TJ, TK) are constants and are shared. The arrays A and B are shared. The macro _PB_N is a constant and shared.

    We don't have any other variables.

    We write the code accordingly.

    We also note: the time loop is outside the parallel region. We are creating a new parallel region for each triple loop in each time step. This might have overhead. We can move the parallel region outside the time loop? But then we would have to synchronize the two triple loops. We cannot because the second triple loop depends on the first. So we leave it as is.

    Alternatively, we can use:

        #pragma omp parallel
        {
            for (t = 1; t <= TSTEPS; t++) {
                #pragma omp for collapse(2) ... // for the first triple loop
                ... 

                #pragma omp for collapse(2) ... // for the second triple loop
                ...
            }
        }

    This would create one parallel region for the entire time loop. This reduces the overhead of creating and destroying threads. We can do that.

    We change:

        #pragma omp parallel
        {
            for (t = 1; t <= TSTEPS; t++) {
                // First triple loop: update B
                #pragma omp for collapse(2) nowait  // we can use nowait? But the second triple loop depends on the first. So we cannot use nowait. We need a barrier.
                ... 

                // We need a barrier here? Actually, the 'omp for' has an implicit barrier at the end. So we don't need an explicit barrier.

                // Second triple loop: update A
                #pragma omp for collapse(2) nowait
                ... 
            }
        }

    But note: the second triple loop in the same time step depends on the first. So we cannot use nowait for the first triple loop. We must have a barrier after the first triple loop. The implicit barrier at the end of the 'omp for' is sufficient.

    We can do:

        #pragma omp parallel default(none) shared(A, B, t, TSTEPS, _PB_N, TI, TJ, TK)
        {
            for (t = 1; t <= TSTEPS; t++) {
                // First triple loop: update B
                #pragma omp for collapse(2) schedule(static)
                for (i = 1; i < _PB_N-1; i += TI) {
                    ... 
                }

                // Second triple loop: update A
                #pragma omp for collapse(2) schedule(static)
                for (i = 1; i < _PB_N-1; i += TI) {
                    ... 
                }
            }
        }

    This is more efficient.

    We also note: the time loop variable t is shared? Actually, it is the loop variable of the time loop and is incremented by the master thread. We must make it private? Or we can use a schedule that doesn't require t to be shared? Actually, the time loop is run by all threads? No, the time loop is outside the worksharing. The worksharing is only for the spatial loops. The time loop is run by all threads? Actually, the entire parallel region is executed by all threads. The time loop is inside the parallel region and is executed by all threads? That would be incorrect because each thread would run the entire time loop.

    We must change: the time loop should be run only once. We can do:

        #pragma omp parallel
        {
            #pragma omp single
            {
                for (t = 1; t <= TSTEPS; t++) {
                    // First triple loop: update B
                    #pragma omp for collapse(2) nowait
                    ... 

                    // We need a barrier after the first triple loop? Yes, because the second triple loop depends on the first. But the 'omp for' has an implicit barrier. So we don't need an explicit barrier.

                    // Second triple loop: update A
                    #pragma omp for collapse(2) nowait
                    ... 
                }
            }
        }

    But then the time loop is run by a single thread. The worksharing loops are distributed. This is acceptable.

    Alternatively, we can use:

        #pragma omp parallel
        {
            for (t = 1; t <= TSTEPS; t++) {
                #pragma omp single
                {
                    // This is not what we want.
                }
            }
        }

    We want the time loop to be run in order, and each iteration of the time loop is sequential. The worksharing is only for the spatial loops.

    We can do:

        #pragma omp parallel
        {
            for (t = 1; t <= TSTEPS; t++) {
                // First triple loop: update B
                #pragma omp for collapse(2) nowait
                ... 

                // We must have a barrier after the first triple loop? Yes, because the second triple loop depends on the first. But the 'omp for' has an implicit barrier. So we don't need an explicit barrier.

                // Second triple loop: update A
                #pragma omp for collapse(2) nowait
                ... 
            }
        }

    But note: the time loop is run by all threads. Each thread will run the entire time loop. That is incorrect.

    We must have the time loop run only once. We can do:

        #pragma omp parallel
        {
            #pragma omp single
            {
                for (t = 1; t <= TSTEPS; t++) {
                    // First triple loop: update B
                    #pragma omp for collapse(2) nowait
                    ... 

                    // We need a barrier after the first triple loop? Yes, but the 'omp for' has an implicit barrier. So we don't need an explicit barrier.

                    // Second triple loop: update A
                    #pragma omp for collapse(2) nowait
                    ... 
                }
            }
        }

    This way, the time loop is run by a single thread, and the worksharing loops are distributed to all threads.

    But note: the 'omp for' inside the single construct is allowed. The single construct doesn't prevent the worksharing.

    We choose this.

    We also note: the variables t, i, j, k, ii, jj, kk, ti_end, tj_end, tk_end are private? 

        - t: is the loop variable of the time loop and is declared in the outer scope. We must make it private? Actually, the single construct is executed by one thread, so t is only used by that thread. The worksharing loops are inside the single construct and are executed by all threads. The variables i, j, k are the tile loop indices and are private to the worksharing loop. The inner variables are private to each thread.

    We also note: the arrays A and B are shared.

    We write:

        const int TI = 32;
        const int TJ = 32;
        const int TK = 32;

        #pragma omp parallel default(none) shared(A, B, TSTEPS, _PB_N, TI, TJ, TK)
        {
            #pragma omp single
            {
                for (t = 1; t <= TSTEPS; t++) {
                    // First triple loop: update B
                    #pragma omp for collapse(2) schedule(static)
                    for (i = 1; i < _PB_N-1; i += TI) {
                        for (j = 1; j < _PB_N-1; j += TJ) {
                            for (k = 1; k < _PB_N-1; k += TK) {
                                int ti_end = (i+TI < _PB_N-1) ? i+TI : _PB_N-1;
                                int tj_end = (j+TJ < _PB_N-1) ? j+TJ : _PB_N-1;
                                int tk_end = (k+TK < _PB_N-1) ? k+TK : _PB_N-1;
                                for (ii = i; ii < ti_end; ii++) {
                                    for (jj = j; jj < tj_end; jj++) {
                                        #pragma omp simd
                                        for (kk = k; kk < tk_end; kk++) {
                                            B[ii][jj][kk] = ... ;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // Second triple loop: update A
                    #pragma omp for collapse(2) schedule(static)
                    for (i = 1; i < _PB_N-1; i += TI) {
                        for (j = 1; j < _PB_N-1; j += TJ) {
                            for (k = 1; k < _PB_N-1; k += TK) {
                                int ti_end = (i+TI < _PB_N-1) ? i+TI : _PB_N-1;
                                int tj_end = (j+TJ < _PB_N-1) ? j+TJ : _PB_N-1;
                                int tk_end = (k+TK < _PB_N-1) ? k+TK : _PB_N-1;
                                for (ii = i; ii < ti_end; ii++) {
                                    for (jj = j; jj < tj_end; jj++) {
                                        #pragma omp simd
                                        for (kk = k; kk < tk_end; kk++) {
                                            A[ii][jj][kk] = ... ;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

    We also note: the variable t is declared in the outer scope. We must declare it inside the scop? Actually, the original code declares t at the beginning of the function. We are not changing that. But inside the scop, we are using t. We must make sure it is private to the single thread? It is, because the single construct is executed by one thread.

    We also note: the original code has the time loop variable t declared as:

        int t, i, j, k;

    We are reusing t. We are not changing the declaration.

    We also note: the inner loops use ii, jj, kk. We declare them inside the tile block.

    We also note: the stencil computation is the same.

    We write the code accordingly.

    We also note: the original code uses the macro TSTEPS. We use that.

    We also note: the original code uses the macro _PB_N. We use that.

    We also note: the original code uses the macro SCALAR_VAL. We use that.

    We also note: the original code uses the macro DATA_TYPE. We use that.

    We write the optimized code.

    We also note: the original code has the scop markers. We preserve them.

    We output the entire code? No, only the scop region. But the problem says: "Return ONLY the complete optimized C code". And we are to preserve the scop markers.

    We output the entire function? Actually, the problem says: "Optimize only code in regions delimited by '#pragma scope' and '#pragma endscop'". So we only change the code between these markers.

    We output the entire kernel_heat_3d function? But the problem says: "Return ONLY the complete optimized C code". And we are to preserve the rest of the code. But the problem says: "Do not modify code outside the marked regions". So we only change the scop region.

    We output the entire file? The problem says: "Return ONLY the complete optimized C code". And we are to preserve the rest. But the problem says: "Preserve function signatures, program inputs/outputs, and original semantics". So we output the entire file with the scop region replaced.

    However, the problem says: "Return ONLY the complete optimized C code". And we are to return the entire code? But the example input is the entire file. So we return the entire file with the scop region replaced.

    But note: the problem says: "You MAY restructure loops and adjust associated C statements within the marked regions". So we only change the scop region.

    We write the entire file with the scop region replaced by our optimized version.

    However, the problem says: "Return ONLY the complete optimized C code". So we return the entire code.

    We write the entire code with the scop region replaced.

    We also note: the original code has the scop markers. We preserve them.

    We write:

        ... (the rest of the code)

        static
        void kernel_heat_3d(int tsteps,
		      int n,
		      DATA_TYPE POLYBENCH_3D(A,N,N,N,n,n,n),
		      DATA_TYPE POLYBENCH_3D(B,N,N,N,n,n,n))
        {
          int t, i, j, k;

        #pragma scop
            const int TI = 32;
            const int TJ = 32;
            const int TK = 32;

            #pragma omp parallel default(none) shared(A, B, tsteps, n, TI, TJ, TK)
            {
                #pragma omp single
                {
                    for (t = 1; t <= tsteps; t++) {
                        // First triple loop: update B
                        #pragma omp for collapse(2) schedule(static)
                        for (i = 1; i < n-1; i += TI) {
                            for (j = 1; j < n-1; j += TJ) {
                                for (k = 1; k < n-1; k += TK) {
                                    int ti_end = (i+TI < n-1) ? i+TI : n-1;
                                    int tj_end = (j+TJ < n-1) ? j+TJ : n-1;
                                    int tk_end = (k+TK < n-1) ? k+TK : n-1;
                                    for (ii = i; ii < ti_end; ii++) {
                                        for (jj = j; jj < tj_end; jj++) {
                                            #pragma omp simd
                                            for (kk = k; kk < tk_end; kk++) {
                                                B[ii][jj][kk] =   SCALAR_VAL(0.125) * (A[ii+1][jj][kk] - SCALAR_VAL(2.0) * A[ii][jj][kk] + A[ii-1][jj][kk])
                                                                 + SCALAR_VAL(0.125) * (A[ii][jj+1][kk] - SCALAR_VAL(2.0) * A[ii][jj][kk] + A[ii][jj-1][kk])
                                                                 + SCALAR_VAL(0.125) * (A[ii][jj][kk+1] - SCALAR_VAL(2.0) * A[ii][jj][kk] + A[ii][jj][kk-1])
                                                                 + A[ii][jj][kk];
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        // Second triple loop: update A
                        #pragma omp for collapse(2) schedule(static)
                        for (i = 1; i < n-1; i += TI) {
                            for (j = 1; j < n-1; j += TJ) {
                                for (k = 1; k < n-1; k += TK) {
                                    int ti_end = (i+TI < n-1) ? i+TI : n-1;
                                    int tj_end = (j+TJ < n-1) ? j+TJ : n-1;
                                    int tk_end = (k+TK < n-1) ? k+TK : n-1;
                                    for (ii = i; ii < ti_end; ii++) {
                                        for (jj = j; jj < tj_end; jj++) {
                                            #pragma omp simd
                                            for (kk = k; kk < tk_end; kk++) {
                                                A[ii][jj][kk] =   SCALAR_VAL(0.125) * (B[ii+1][jj][kk] - SCALAR_VAL(2.0) * B[ii][jj][kk] + B[ii-1][jj][kk])
                                                                + SCALAR_VAL(0.125) * (B[ii][jj+1][kk] - SCALAR_VAL(2.0) * B[ii][jj][kk] + B[ii][jj-1][kk])
                                                                + SCALAR_VAL(0.125) * (B[ii][jj][kk+1] - SCALAR_VAL(2.0) * B[ii][jj][kk] + B[ii][jj][kk-1])
                                                                + B[ii][jj][kk];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        #pragma endscop

        ... (the rest of the code)

    We also note: the original code uses the parameter n, not _PB_N. We use n.

    We also note: the original code uses the parameter tsteps, not TSTEPS. We use tsteps.

    We also note: the original code uses the arrays A and B. We use them.

    We also note: the original code uses the variable t, i, j, k. We use them for the outer loops. We use ii, jj, kk for the inner tile loops.

    We also note: we declare ii, jj, kk inside the tile block. We must declare them. We do:

        int ii, jj, kk;

    But we can declare them in the for loops? In C99, we can. But the original code might be C89? The original code declares the loop variables at the beginning of the function. We are inside a block. We can declare them at the beginning of the tile block? Or we can declare them in the inner for loops.

    We choose to declare them in the inner for loops:

        for (int ii = i; ii < ti_end; ii++) {
            for (int jj = j; jj < tj_end; jj++) {
                #pragma omp simd
                for (int kk = k; kk < tk_end; kk++) {
                    ...
                }
            }
        }

    This is C99. The original code uses C99? The original code uses:

        for (i = 0; i < n; i++)

    without declaring i inside the for. But the function declares i at the beginning. We are changing the inner loops to declare the loop variables inside the for. This is acceptable in C99.

    We assume the compiler supports C99.

    We write accordingly.

    We also note: the original code uses the macro POLYBENCH_3D. We preserve.

    We output the entire code.

    However, the problem says: "Return ONLY the complete optimized C code". So we return the entire code with the scop region replaced.

    We write the entire code.

    But note: the problem says: "Do not modify code outside the marked regions". So we only change the scop region.

    We output the entire file with the scop region replaced.

    We also note: the original code has the scop markers. We preserve them.

    We write the code accordingly.

    We also note: the original code has the function parameters: tsteps, n, A, B. We use them.

    We also note: the original code uses the macro SCALAR_VAL. We use that.

    We also note: the original code uses the macro DATA_TYPE. We use that.

    We write the code.

    We also note: the original code uses the macro POLYBENCH_3D. We use that.

    We write the code.

    We also note: the original code uses the macro POLYBENCH_DUMP_START, etc. We preserve.

    We write the entire code.

    However, the problem says: "Return ONLY the complete optimized C code". So we return the entire code.

    We output the entire code with the scop region replaced.

    We also note: the original code has the main function. We preserve.

    We write the entire code.

    But note: the problem says: "You are an expert HPC compiler optimizing C code for CPU execution with OpenMP". We have applied tiling and OpenMP.

    We hope this is correct.

    We output the code.
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "heat-3d.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_3D(A,N,N,N,n,n,n),
		 DATA_TYPE POLYBENCH_3D(B,N,N,N,n,n,n))
{
  int i, j, k;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++)
        A[i][j][k] = B[i][j][k] = (DATA_TYPE) (i + j + (n-k))* 10 / (n);
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_3D(A,N,N,N,n,n,n))

{
  int i, j, k;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("A");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++) {
         if ((i * n * n + j * n + k) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
         fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j][k]);
      }
  POLYBENCH_DUMP_END("A");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_heat_3d(int tsteps,
		      int n,
		      DATA_TYPE POLYBENCH_3D(A,N,N,N,n,n,n),
		      DATA_TYPE POLYBENCH_3D(B,N,N,N,n,n,n))
{
  int t, i, j, k;

#pragma scop
  const int TI = 32;
  const int TJ = 32;
  const int TK = 32;

  #pragma omp parallel default(none) shared(A, B, tsteps, n, TI, TJ, TK)
  {
    #pragma omp single
    {
      for (t = 1; t <= tsteps; t++) {
        // First triple loop: update B
        #pragma omp for collapse(2) schedule(static)
        for (i = 1; i < n-1; i += TI) {
          for (j = 1; j < n-1; j += TJ) {
            for (k = 1; k < n-1; k += TK) {
              int ti_end = (i+TI < n-1) ? i+TI : n-1;
              int tj_end = (j+TJ < n-1) ? j+TJ : n-1;
              int tk_end = (k+TK < n-1) ? k+TK : n-1;
              for (int ii = i; ii < ti_end; ii++) {
                for (int jj = j; jj < tj_end; jj++) {
                  #pragma omp simd
                  for (int kk = k; kk < tk_end; kk++) {
                    B[ii][jj][kk] =   SCALAR_VAL(0.125) * (A[ii+1][jj][kk] - SCALAR_VAL(2.0) * A[ii][jj][kk] + A[ii-1][jj][kk])
                                 + SCALAR_VAL(0.125) * (A[ii][jj+1][kk] - SCALAR_VAL(2.0) * A[ii][jj][kk] + A[ii][jj-1][kk])
                                 + SCALAR_VAL(0.125) * (A[ii][jj][kk+1] - SCALAR_VAL(2.0) * A[ii][jj][kk] + A[ii][jj][kk-1])
                                 + A[ii][jj][kk];
                  }
                }
              }
            }
          }
        }

        // Second triple loop: update A
        #pragma omp for collapse(2) schedule(static)
        for (i = 1; i < n-1; i += TI) {
          for (j = 1; j < n-1; j += TJ) {
            for (k = 1; k < n-1; k += TK) {
              int ti_end = (i+TI < n-1) ? i+TI : n-1;
              int tj_end = (j+TJ < n-1) ? j+TJ : n-1;
              int tk_end = (k+TK < n-1) ? k+TK : n-1;
              for (int ii = i; ii < ti_end; ii++) {
                for (int jj = j; jj < tj_end; jj++) {
                  #pragma omp simd
                  for (int kk = k; kk < tk_end; kk++) {
                    A[ii][jj][kk] =   SCALAR_VAL(0.125) * (B[ii+1][jj][kk] - SCALAR_VAL(2.0) * B[ii][jj][kk] + B[ii-1][jj][kk])
                                 + SCALAR_VAL(0.125) * (B[ii][jj+1][kk] - SCALAR_VAL(2.0) * B[ii][jj][kk] + B[ii][jj-1][kk])
                                 + SCALAR_VAL(0.125) * (B[ii][jj][kk+1] - SCALAR_VAL(2.0) * B[ii][jj][kk] + B[ii][jj][kk-1])
                                 + B[ii][jj][kk];
                  }
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
  POLYBENCH_3D_ARRAY_DECL(A, DATA_TYPE, N, N, N, n, n, n);
  POLYBENCH_3D_ARRAY_DECL(B, DATA_TYPE, N, N, N, n, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_heat_3d (tsteps, n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

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
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;
  int tsteps = TSTEPS;

  /* Variable declaration/allocation. */
  POLYBENCH_3D_ARRAY_DECL(A, DATA_TYPE, N, N, N, n, n, n);
  POLYBENCH_3D_ARRAY_DECL(B, DATA_TYPE, N, N, N, n, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_heat_3d (tsteps, n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

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

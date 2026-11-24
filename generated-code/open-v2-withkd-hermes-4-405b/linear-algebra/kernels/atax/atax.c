/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* atax.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "atax.h"


/* Array initialization. */
static
void init_array (int m, int n,
		 DATA_TYPE POLYBENCH_2D(A,M,N,m,n),
		 DATA_TYPE POLYBENCH_1D(x,N,n))
{
  int i, j;
  DATA_TYPE fn;
  fn = (DATA_TYPE)n;

  for (i = 0; i < n; i++)
      x[i] = 1 + (i / fn);
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      A[i][j] = (DATA_TYPE) ((i+j) % n) / (5*m);
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_1D(y,N,n))

{
  int i;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("y");
  for (i = 0; i < n; i++) {
    if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, y[i]);
  }
  POLYBENCH_DUMP_END("y");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_atax(int m, int n,
		 DATA_TYPE POLYBENCH_2D(A,M,N,m,n),
		 DATA_TYPE POLYBENCH_1D(x,N,n),
		 DATA_TYPE POLYBENCH_1D(y,N,n),
		 DATA_TYPE POLYBENCH_1D(tmp,M,m))
{
  int i, j;

#pragma scop
We are given the atax.c code and are to optimize the region between #pragma scop and #pragma endscop.
 The original code in the scop region:
   for (i = 0; i < _PB_N; i++)
     y[i] = 0;
   for (i = 0; i < _PB_M; i++)
     {
       tmp[i] = SCALAR_VAL(0.0);
       for (j = 0; j < _PB_N; j++)
         tmp[i] = tmp[i] + A[i][j] * x[j];
       for (j = 0; j < _PB_N; j++)
         y[j] = y[j] + A[i][j] * tmp[i];
     }

 We are to apply loop transformations (tiling, fusion, skewing) first, then OpenMP directives.

 Step 1: Loop Transformations
   - Tiling: We have two nested loops. The inner loops are over j. We can tile the j loops to improve cache locality for A and x (in the first inner loop) and for A and y (in the second inner loop). However, note that the two inner loops are adjacent and we might consider fusing them? But they are separated by the tmp[i] computation.

   Let's look at the structure:
     Loop1: i from 0 to M-1
        tmp[i] = 0
        Loop2: j from 0 to N-1: tmp[i] += A[i][j] * x[j]
        Loop3: j from 0 to N-1: y[j] += A[i][j] * tmp[i]

   We cannot fuse Loop2 and Loop3 because they are not adjacent and have a dependence (tmp[i] is computed in Loop2 and used in Loop3). However, we can tile the inner j loops.

   We'll tile the j loops in both inner loops. We choose a tile size, say 64 (a common cache line multiple). We'll use TI as the tile size.

   We can also consider tiling the outer i loop? But note that the outer i loop is over M, and the inner j loops are over N. The array A is of size MxN, so we might tile both i and j? However, the two inner loops are independent in j, so we can tile j without tiling i.

   Alternatively, we can tile the j dimension in both inner loops. This will help in reusing the same elements of A[i][j] and x[j] (for the first inner loop) and A[i][j] and y[j] (for the second inner loop) within a tile.

   We'll do:
     for (i = 0; i < _PB_M; i++) {
        tmp[i] = 0.0;
        for (jj = 0; jj < _PB_N; jj += TI) {
            for (j = jj; j < min(jj+TI, _PB_N); j++) {
                tmp[i] += A[i][j] * x[j];
            }
        }
        for (jj = 0; jj < _PB_N; jj += TI) {
            for (j = jj; j < min(jj+TI, _PB_N); j++) {
                y[j] += A[i][j] * tmp[i];
            }
        }
     }

   But note: the first inner loop (computing tmp[i]) is a reduction. We can also consider fusing the two inner loops? Actually, we cannot because the entire tmp[i] must be computed before the second inner loop. However, we can try to interchange the order? Not easily.

   Alternatively, we can tile the outer i loop? But the outer i loop is not the one causing cache misses for A. The inner j loops are the ones that access A by rows. Since C is row-major, accessing A[i][j] for fixed i and varying j is contiguous. So the main issue is that the same row of A is used in two consecutive inner loops. We are already reusing the row of A in the two inner loops. The tiling of j will help in keeping a block of A[i][j] and x[j] (or y[j]) in cache.

   However, note that the two inner loops are over the same j, so we can try to combine them? Actually, we can compute the entire tmp[i] and then update y. We cannot change the order.

   Another idea: we can tile the i loop? But the i loop is over M, and each i iteration is independent? Actually, the i iterations are independent. So we can also tile the i loop to improve locality for tmp and for the updates to y? But note that y is being updated by every i, so we want to update y in a way that is cache friendly. The current update pattern for y is: for each i, we update every j. This is a streaming update. Tiling i might not help y because y is updated in a scattered manner (each i updates the entire y). However, if we tile i, then we can update a block of y multiple times (for each i in the tile) and then move on. This might be beneficial.

   We can do:
     for (ii = 0; ii < _PB_M; ii += TI) {
        for (i = ii; i < min(ii+TI, _PB_M); i++) {
            tmp[i] = 0.0;
            for (j = 0; j < _PB_N; j++) {
                tmp[i] += A[i][j] * x[j];
            }
        }
        for (i = ii; i < min(ii+TI, _PB_M); i++) {
            for (j = 0; j < _PB_N; j++) {
                y[j] += A[i][j] * tmp[i];
            }
        }
     }

   But note: the update of y[j] in the second part of the tile: we are updating the entire y for each i in the tile. This means that for a fixed j, we update it for every i in the tile. This is good because we are reusing the same y[j] and the same A[i][j] for the entire tile? Actually, A[i][j] is different for each i. But we are accessing A by rows, so for a fixed j, we are accessing A[i][j] for consecutive i (in the tile) which is contiguous in memory (if stored by row). So this is good.

   However, we are updating y[j] for every i in the tile. This is a reduction over i for each j. We can do that.

   But note: the original code does:
        for i: 
            compute tmp[i] (using A[i][*] and x[*])
            then update y[*] using A[i][*] and tmp[i]

   In the tiled version, we are doing:
        for a tile of i:
            compute tmp[i] for each i in the tile
        then for the same tile of i:
            update y for each i in the tile

   This is correct because the update of y in the second part uses the tmp[i] computed in the first part for the same i.

   However, we are updating y for the entire j range for each i in the tile. This means that we are updating the entire y array multiple times (once per i in the tile). This is the same as the original, but now we are doing it in a block of i's.

   The advantage: we are reusing the same row of A (for a fixed i) in the two inner loops? Actually, in the tiled version, we are not reusing the row of A for the same i in the two inner loops because we break the inner j loops? We are not. We are doing the entire j loop for the first part (computing tmp) and then the entire j loop for the second part (updating y). So we are still accessing the row of A twice. But we are doing it for a block of i's. The hope is that the block of rows of A (for the tile of i) will stay in cache.

   We can also tile the j loops? We can do both: tile i and then within each i, tile j? That might be too complex.

   Given the constraints, we choose to tile the outer i loop. We choose a tile size, say 64, for i.

   Step 2: OpenMP
     We can parallelize the outer i loop? But note: the outer i loop is now tiled. We can parallelize the tile loop (over ii). Each tile is independent? Yes, because the computation of tmp[i] for i in a tile does not depend on other tiles, and the update of y is a reduction over i (so we can do it in any order). However, the update of y is a reduction: multiple tiles will update the same y[j]. So we need to use a reduction or critical section? But note: the update is y[j] += ... and it is done for every i. We can use an OpenMP reduction? But the reduction is over the entire i loop, and we are breaking it into tiles. We can do:

        #pragma omp parallel for private(i,j) reduction(+:y[0:_PB_N])
        for (ii = 0; ii < _PB_M; ii += TI) {
            for (i = ii; i < min(ii+TI, _PB_M); i++) {
                tmp[i] = 0.0;
                for (j = 0; j < _PB_N; j++) {
                    tmp[i] += A[i][j] * x[j];
                }
            }
            for (i = ii; i < min(ii+TI, _PB_M); i++) {
                for (j = 0; j < _PB_N; j++) {
                    y[j] += A[i][j] * tmp[i];
                }
            }
        }

     However, note that the reduction clause in OpenMP for arrays is not standard in all versions. Alternatively, we can use a critical section for updating y? That would be too expensive.

     Another idea: we can use a private array for y in each thread and then do a reduction at the end? But that would require extra memory and a final reduction step.

     Alternatively, we can note that the update of y is a reduction over i. We can change the order: instead of updating y for each i, we can accumulate the contributions per j. But that is what we are doing.

     Actually, the standard way is to use a reduction clause for the entire array. In OpenMP 3.1 and later, we can do:

        reduction(+:y[:_PB_N])

     But note: the original code initializes y to zero at the beginning. We are doing that outside the parallel region. Then we are updating y in the parallel region. We must be cautious: the reduction clause will create a private copy for each thread and then combine them. This is safe.

     However, the tile loop is over ii, and we are parallelizing the tile loop. Each thread will get a set of tiles. The reduction will be done at the end of the parallel region.

     But note: the tile loop is the outermost loop. We are parallelizing the tile loop. This is good because the tile size (TI) can be chosen to be large enough to amortize the overhead.

     We also need to consider the initialization of y. The initialization is done before the parallel region. We can leave it as is.

     We also need to consider the tmp array. The tmp array is written per i and then read per i. It is private to the i iteration. We can make it private to the thread? Actually, we are writing to tmp[i] for i in the tile. Since each tile is processed by one thread (because we are parallelizing the tile loop), we don't need to make tmp private. But note: the tmp array is shared. We are writing to different i's in different tiles, so it is safe. However, we must ensure that the same i is not processed by two threads. And it isn't because the tile loop is partitioned.

     But note: the tile loop is over ii, and we are processing a block of i's. The tmp array is indexed by i, and each i is unique. So it is safe.

     However, we must declare the variables appropriately. We'll use:

        #pragma omp parallel for private(i,j) reduction(+:y[0:_PB_N]) schedule(static)

     We choose static scheduling because the work per tile is uniform.

     We also note that the inner j loops can be vectorized. We can add a simd directive? But the compiler might do that automatically.

     We can also consider collapsing the tile loop? But the tile loop is the only loop we are parallelizing. We don't need to collapse.

     We also note that the initialization of y is outside the parallel region. We can leave it as is.

   Step 3: Code structure

     We'll do:

        // Initialize y
        for (i = 0; i < _PB_N; i++)
            y[i] = 0;

        // Tiled and parallelized computation
        #pragma omp parallel for private(i,j) reduction(+:y[0:_PB_N]) schedule(static)
        for (ii = 0; ii < _PB_M; ii += TI) {
            int ii_end = min(ii+TI, _PB_M);
            for (i = ii; i < ii_end; i++) {
                tmp[i] = SCALAR_VAL(0.0);
                for (j = 0; j < _PB_N; j++) {
                    tmp[i] += A[i][j] * x[j];
                }
            }
            for (i = ii; i < ii_end; i++) {
                for (j = 0; j < _PB_N; j++) {
                    y[j] += A[i][j] * tmp[i];
                }
            }
        }

   We choose TI = 64. We can make it a constant.

   However, note: the original code uses _PB_M and _PB_N. We must be cautious with the min function.

   We'll define:
        #define TI 64

   But note: the problem sizes might be small. We can use a variable tile size? But the problem says to choose legal tile sizes. We'll use 64.

   Also note: the inner j loops are not tiled. We could tile them as well? But we are already tiling the i loop. We'll leave the j loops as is for now.

   Alternatively, we can tile the j loops? We can do:

        for (i = ii; i < ii_end; i++) {
            tmp[i] = 0.0;
            for (jj = 0; jj < _PB_N; jj += TJ) {
                int jj_end = min(jj+TJ, _PB_N);
                for (j = jj; j < jj_end; j++) {
                    tmp[i] += A[i][j] * x[j];
                }
            }
        }

        for (i = ii; i < ii_end; i++) {
            for (jj = 0; jj < _PB_N; jj += TJ) {
                int jj_end = min(jj+TJ, _PB_N);
                for (j = jj; j < jj_end; j++) {
                    y[j] += A[i][j] * tmp[i];
                }
            }
        }

   We can choose TJ = 64 as well. This might help in reusing x and y in cache? But note: x is being read in the first inner loop and y is being written in the second. The array A is being read in both. Tiling j might help in reusing a block of A[i][j] and x[j] (or y[j]) in the cache.

   However, we are already tiling i. We can do both. We'll do:

        #define TI 64
        #define TJ 64

        for (ii = 0; ii < _PB_M; ii += TI) {
            int ii_end = min(ii+TI, _PB_M);
            for (i = ii; i < ii_end; i++) {
                tmp[i] = 0.0;
                for (jj = 0; jj < _PB_N; jj += TJ) {
                    int jj_end = min(jj+TJ, _PB_N);
                    for (j = jj; j < jj_end; j++) {
                        tmp[i] += A[i][j] * x[j];
                    }
                }
            }
            for (i = ii; i < ii_end; i++) {
                for (jj = 0; jj < _PB_N; jj += TJ) {
                    int jj_end = min(jj+TJ, _PB_N);
                    for (j = jj; j < jj_end; j++) {
                        y[j] += A[i][j] * tmp[i];
                    }
                }
            }
        }

   This is a 2-level tiling.

   But note: the problem says to prioritize tiling. We are doing two levels. We'll do this.

   Step 4: OpenMP again: we are parallelizing the outer tile loop (ii). We'll keep the same OpenMP directive.

   Step 5: We must be cautious about the min function. We can use:

        int ii_end = (ii + TI < _PB_M) ? (ii + TI) : _PB_M;

   But we can write a macro for min? Or use a function. We'll use:

        #define MIN(a,b) ((a) < (b) ? (a) : (b))

   But note: the original code uses _PB_M and _PB_N. We'll use the same.

   However, the problem says to preserve the original semantics. We are not changing the order of computation, just tiling and parallelizing.

   We'll write the code accordingly.

   Step 6: We must also consider the initialization of tmp. We are setting tmp[i] to 0.0 for each i. This is done in the first inner loop. We are not changing that.

   Step 7: We must declare the variables ii, jj, ii_end, jj_end. We'll declare them inside the scop region.

   Final code for the scop region:

        for (i = 0; i < _PB_N; i++)
            y[i] = 0;

        #define TI 64
        #define TJ 64
        #define MIN(a,b) ((a) < (b) ? (a) : (b))

        #pragma omp parallel for private(i,j,ii,jj,ii_end,jj_end) reduction(+:y[0:_PB_N]) schedule(static)
        for (ii = 0; ii < _PB_M; ii += TI) {
            int ii_end = MIN(ii+TI, _PB_M);
            for (i = ii; i < ii_end; i++) {
                tmp[i] = SCALAR_VAL(0.0);
                for (jj = 0; jj < _PB_N; jj += TJ) {
                    int jj_end = MIN(jj+TJ, _PB_N);
                    for (j = jj; j < jj_end; j++) {
                        tmp[i] += A[i][j] * x[j];
                    }
                }
            }
            for (i = ii; i < ii_end; i++) {
                for (jj = 0; jj < _PB_N; jj += TJ) {
                    int jj_end = MIN(jj+TJ, _PB_N);
                    for (j = jj; j < jj_end; j++) {
                        y[j] += A[i][j] * tmp[i];
                    }
                }
            }
        }

        #undef TI
        #undef TJ
        #undef MIN

   But note: the reduction clause for an array section: y[0:_PB_N] is supported in OpenMP 3.1 and later. We assume that.

   Alternatively, if we are concerned about portability, we could use a different method? But the problem does not specify the OpenMP version. We'll use this.

   Also note: the private clause: we list the variables that are private to each thread: i, j, ii, jj, ii_end, jj_end. The loop variables ii, i, jj, j are private by default in the worksharing loop? Actually, in OpenMP, the loop variable of the worksharing loop (ii) is private. The inner loops: i, j, jj are declared inside the parallel region? We are declaring them in the private clause to be safe.

   However, note: the inner loops are inside the parallel region. The variables i, j, jj, ii_end, jj_end are declared in the outer scope? Actually, we are declaring them in the for loops? We are not. We are using them without declaring. We must declare them. We'll declare them at the beginning of the scop region? But the original code declares i and j at the beginning of the function. We are allowed to restructure within the scop region.

   We'll declare:

        int i, j, ii, jj, ii_end, jj_end;

   inside the scop region? But note: the initialization of y uses i. We can declare i and j at the beginning of the scop region and then reuse them.

   We'll do:

        int i, j;
        for (i = 0; i < _PB_N; i++)
            y[i] = 0;

        #define TI 64
        #define TJ 64
        #define MIN(a,b) ((a) < (b) ? (a) : (b))

        int ii, jj, ii_end, jj_end;
        #pragma omp parallel for private(i,j,ii,jj,ii_end,jj_end) reduction(+:y[0:_PB_N]) schedule(static)
        for (ii = 0; ii < _PB_M; ii += TI) {
            ii_end = MIN(ii+TI, _PB_M);
            for (i = ii; i < ii_end; i++) {
                tmp[i] = SCALAR_VAL(0.0);
                for (jj = 0; jj < _PB_N; jj += TJ) {
                    jj_end = MIN(jj+TJ, _PB_N);
                    for (j = jj; j < jj_end; j++) {
                        tmp[i] += A[i][j] * x[j];
                    }
                }
            }
            for (i = ii; i < ii_end; i++) {
                for (jj = 0; jj < _PB_N; jj += TJ) {
                    jj_end = MIN(jj+TJ, _PB_N);
                    for (j = jj; j < jj_end; j++) {
                        y[j] += A[i][j] * tmp[i];
                    }
                }
            }
        }

        #undef TI
        #undef TJ
        #undef MIN

   But note: the variable i and j are declared at the beginning of the scop region. We are using them in the initialization and then in the parallel region. In the parallel region, we have a private clause for i and j, so each thread will have its own private copy. This is safe.

   However, the initialization of y uses the original i. That is fine because it is before the parallel region.

   We are done.

   But note: the problem says to preserve the original semantics. We are not changing the order of computation, just tiling and parallelizing. The reduction is done correctly.

   We'll output the entire code with the scop region replaced.

   However, note: the original code uses _PB_N and _PB_M. We are using the same.

   We'll write the code accordingly.

   Let's code accordingly.
#pragma scop
  for (i = 0; i < _PB_N; i++)
    y[i] = 0;
  #define TI 64
  #define TJ 64
  #define MIN(a,b) ((a) < (b) ? (a) : (b))
  int ii, jj, ii_end, jj_end;
  #pragma omp parallel for private(i,j,ii,jj,ii_end,jj_end) reduction(+:y[0:_PB_N]) schedule(static)
  for (ii = 0; ii < _PB_M; ii += TI) {
    ii_end = MIN(ii + TI, _PB_M);
    for (i = ii; i < ii_end; i++) {
      tmp[i] = SCALAR_VAL(0.0);
      for (jj = 0; jj < _PB_N; jj += TJ) {
        jj_end = MIN(jj + TJ, _PB_N);
        for (j = jj; j < jj_end; j++) {
          tmp[i] += A[i][j] * x[j];
        }
      }
    }
    for (i = ii; i < ii_end; i++) {
      for (jj = 0; jj < _PB_N; jj += TJ) {
        jj_end = MIN(jj + TJ, _PB_N);
        for (j = jj; j < jj_end; j++) {
          y[j] += A[i][j] * tmp[i];
        }
      }
    }
  }
  #undef TI
  #undef TJ
  #undef MIN
#pragma endscop
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int m = M;
  int n = N;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, M, N, m, n);
  POLYBENCH_1D_ARRAY_DECL(x, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(y, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(tmp, DATA_TYPE, M, m);

  /* Initialize array(s). */
  init_array (m, n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(x));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_atax (m, n,
	       POLYBENCH_ARRAY(A),
	       POLYBENCH_ARRAY(x),
	       POLYBENCH_ARRAY(y),
	       POLYBENCH_ARRAY(tmp));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(y)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(x);
  POLYBENCH_FREE_ARRAY(y);
  POLYBENCH_FREE_ARRAY(tmp);

  return 0;
}

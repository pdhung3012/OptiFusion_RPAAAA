We are given a code for the FDTD-2D kernel. The scop region contains a time loop (t) that iterates TMAX times. Inside the time loop, there are four nested loops:

1. The first loop: sets the first row of `ey` (i=0) to `_fict_[t]` for all j.
2. The second loop: updates `ey` for i from 1 to NX-1 and all j, using neighboring `hz` values.
3. The third loop: updates `ex` for all i and j from 1 to NY-1, using neighboring `hz` values.
4. The fourth loop: updates `hz` for i from 0 to NX-2 and j from 0 to NY-2, using neighboring `ex` and `ey` values.

We are to optimize the code within the scop region.

HIGHEST PRIORITY: Loop transformations (tiling, fusion, skewing)

1. Loop Tiling: We have four separate loops. Tiling might be beneficial for the inner loops to improve cache locality. However, note that the loops are over the same arrays and have dependencies. We must be cautious.

   - The first loop is a simple assignment and is not computationally heavy. Tiling it might not be necessary.
   - The second and third loops are similar and update `ey` and `ex` respectively. They both use `hz` from the previous time step. The fourth loop updates `hz` using the updated `ex` and `ey` (from the current time step). Therefore, the fourth loop must come after the second and third.

   We cannot tile across the time loop because of the sequential nature of the time steps. However, we can consider tiling the inner loops (i and j) for each of the four loops. But note that the four loops are independent in terms of the same time step? Actually, the second and third loops can be done in any order? Let's check:

   - The second loop updates `ey` (for i>=1) and the third loop updates `ex` (for j>=1). They do not depend on each other. So they can be fused? But note: the fourth loop uses both `ex` and `ey` that are updated in the second and third loops. So the fourth loop must come after both.

   However, the first loop sets `ey[0][j]` and then the second loop uses `ey[i][j]` (for i>=1) and also `hz[i][j]` and `hz[i-1][j]`. The first loop does not affect the second loop because the second loop starts at i=1. So the first loop can be done at any time before the second loop? Actually, the first loop sets `ey[0][j]` and the second loop does not use `ey[0][j]` (it starts at i=1). So the first loop is independent of the second and third? But note: the fourth loop uses `ey[i][j]` for i from 0 to NX-2. The first loop sets `ey[0][j]` and the second loop sets `ey[i][j]` for i>=1. So the fourth loop uses the updated `ey` (both the first row and the rest). Therefore, the first loop must be done before the fourth loop, but it can be done at the beginning of the time step.

   So the order of the four loops in the time step is fixed: first, then second and third (which can be in any order relative to each other), then fourth.

   Now, for tiling: we can tile each of the inner loops (i and j) for the second, third, and fourth loops. The first loop is a simple row assignment and might not benefit from tiling.

   However, note that the second and third loops are over the entire domain (with boundaries) and the fourth loop is over a slightly smaller domain. We can tile the second and third loops together? But they are independent and can be fused? Actually, they are two separate loops. We can consider fusing the second and third loops? But they have different loop bounds and update different arrays. However, they both iterate over the entire i and j (with the second having i from 1 to NX-1 and the third having j from 1 to NY-1). We cannot fuse them because the iteration spaces are not the same. But note: the second loop does not update the boundaries in j? Actually, the second loop does j from 0 to NY-1. The third loop does j from 1 to NY-1. So they are not the same.

   Alternatively, we can tile each loop individually. We'll choose a tile size that fits in cache. Since we are not given specific cache sizes, we can choose a reasonable tile size (e.g., 64x64) and let the user adjust.

   However, note that the fourth loop uses the results of the second and third. So we cannot tile the fourth loop until the second and third are done. Therefore, we will tile each of the second, third, and fourth loops independently.

   But wait: the second and third loops are independent and can be done in any order. We can also consider fusing the second and third? Actually, they are two separate loops. We can try to fuse them by creating a loop that does both updates? But note the iteration spaces are different. We can do:

   for (i=0; i < NX; i++) {
     for (j=0; j < NY; j++) {
        if (i>=1) { ... update ey ... }
        if (j>=1) { ... update ex ... }
     }
   }

   But note: the original second loop does i from 1 to NX-1 and j from 0 to NY-1. The third loop does i from 0 to NX-1 and j from 1 to NY-1. So if we do a fused loop for i from 0 to NX-1 and j from 0 to NY-1, then we can conditionally update. However, this might introduce conditionals that could be costly. Alternatively, we can split the fused loop into boundary and inner parts? But that might complicate.

   Given the complexity and the fact that the two loops are independent, we might leave them separate and tile each.

   However, the problem says: "loop fusion: investigate merging adjacent loops with compatible iteration spaces". The iteration spaces are not exactly the same, but they are adjacent and we can extend the iteration space to cover both and then conditionally execute. But note: the second loop does not update when i==0, and the third loop does not update when j==0. So we can do:

   for (i=0; i < NX; i++) {
     for (j=0; j < NY; j++) {
        if (i > 0) {
            ey[i][j] = ...;
        }
        if (j > 0) {
            ex[i][j] = ...;
        }
     }
   }

   This would replace the second and third loops. But note: the original second loop does i from 1 to NX-1 (so i>0 and i < NX) and j from 0 to NY-1. The condition i>0 is the same. Similarly, the third loop does j>0. So this fused loop would cover the same updates. However, we are now iterating over the entire domain (i=0 to NX-1, j=0 to NY-1) and conditionally updating. This might be acceptable and even beneficial because we are doing one pass over the data. Also, the conditionals are simple and might be predicted well.

   But note: the original second loop does not update the last row? Actually, the second loop does i from 1 to NX-1 (so i < NX) and j from 0 to NY-1. The fused loop does i from 0 to NX-1 and j from 0 to NY-1. So we are including i=0 and j=0, but we skip the update for i=0 in the ey part and j=0 in the ex part. So it's safe.

   Also, the fourth loop uses the updated ex and ey. So we can do the fused loop and then the fourth loop.

   Now, what about the first loop? The first loop sets ey[0][j] for all j. We can do that before the fused loop? Yes, because the fused loop for ey update starts at i=1. So the first loop is independent and can be done at the beginning.

   So the new structure for the time step would be:

   for(t=0; t<TMAX; t++) {
      // First loop: set ey[0][j] = _fict_[t] for all j.
      for (j=0; j<NY; j++)
          ey[0][j] = _fict_[t];

      // Fused loop for ey (i>=1) and ex (j>=1)
      for (i=0; i<NX; i++) {
          for (j=0; j<NY; j++) {
              if (i > 0) {
                  ey[i][j] = ey[i][j] - 0.5*(hz[i][j] - hz[i-1][j]);
              }
              if (j > 0) {
                  ex[i][j] = ex[i][j] - 0.5*(hz[i][j] - hz[i][j-1]);
              }
          }
      }

      // Fourth loop: update hz for i from 0 to NX-2, j from 0 to NY-2.
      for (i=0; i<NX-1; i++) {
          for (j=0; j<NY-1; j++) {
              hz[i][j] = hz[i][j] - 0.7*(ex[i][j+1] - ex[i][j] + ey[i+1][j] - ey[i][j]);
          }
      }
   }

   Now, we can tile the fused loop and the fourth loop.

   For the fused loop: we are iterating over i and j. We can tile both i and j. Similarly for the fourth loop.

   We choose tile sizes: let's say TILE_I and TILE_J. We can set them to 64, but we must be cautious because the inner loop might be small. We'll use:

   #define TILE_I 64
   #define TILE_J 64

   But note: the problem says "choose legal tile sizes". We are not given the actual NX and NY, so we choose a fixed tile size that is reasonable.

   We'll tile the fused loop:

   for (ii = 0; ii < NX; ii += TILE_I) {
     for (jj = 0; jj < NY; jj += TILE_J) {
        for (i = ii; i < min(ii+TILE_I, NX); i++) {
          for (j = jj; j < min(jj+TILE_J, NY); j++) {
              if (i > 0) {
                  ey[i][j] = ...;
              }
              if (j > 0) {
                  ex[i][j] = ...;
              }
          }
        }
     }
   }

   Similarly, for the fourth loop:

   for (ii = 0; ii < NX-1; ii += TILE_I) {
     for (jj = 0; jj < NY-1; jj += TILE_J) {
        for (i = ii; i < min(ii+TILE_I, NX-1); i++) {
          for (j = jj; j < min(jj+TILE_J, NY-1); j++) {
              hz[i][j] = ...;
          }
        }
     }
   }

   However, note that the fourth loop uses ex[i][j+1] and ey[i+1][j]. This means that we need to be cautious at the tile boundaries. But since we are tiling the entire domain and the stencil is small (only one element away), we can handle the boundaries by extending the tile by one? Actually, we are iterating over the entire domain and the conditionals in the inner loop take care of the boundaries. But note: the fourth loop does not update the last row and last column. So we are safe.

   But wait: the fourth loop uses ex[i][j+1] and ey[i+1][j]. This means that for a tile, we need the next row of ey and the next column of ex. Therefore, we must be cautious when tiling the fourth loop: we need to include the next row and next column in the tile? Actually, we are reading from ex and ey that are updated in the fused loop. The fused loop updates the entire domain (except the boundaries that are conditionally skipped). So the data is available. However, when we tile the fourth loop, we are only updating the inner part. The reads for ex[i][j+1] and ey[i+1][j] are from the same time step and are already computed by the fused loop. So we don't need to change the tile.

   But note: the fourth loop does not update the last row and last column, so we don't have to worry about writing beyond the tile.

   However, the reads for ex[i][j+1] when j is the last in the tile: we need j+1 to be in the next tile? Actually, we are iterating j from 0 to NY-2. So j+1 is at most NY-1, which is within the array. Similarly, i+1 is at most NX-1. So we are safe.

   Therefore, we can tile the fused loop and the fourth loop.

   But note: the first loop is a simple row assignment and is not heavy. We can leave it as is.

   So the plan is:

   - Fuse the second and third loops into one double nested loop with conditionals.
   - Tile the fused loop and the fourth loop.

2. Loop Skewing: We don't see any carried dependencies that would require skewing. The time loop is sequential, but the inner loops are parallel. So we skip skewing.

SECOND PRIORITY: OpenMP directives.

We are to parallelize the inner loops. We note:

- The time loop (t) is sequential because each time step depends on the previous.
- The first loop: sets the first row of ey. It is a simple loop over j. We can parallelize it? But note: it's a small loop (only one row) and might not be worth parallelizing. Also, it's inside a time loop that runs TMAX times. We can parallelize the j loop? But the overhead might be too high. We'll leave it sequential.

- The fused loop: we have two nested loops (i and j). We can parallelize the outer loop (i) or the inner loop (j). Since the loop is over a 2D domain, we can use collapse(2) to parallelize both. But note: the conditionals inside are independent. So we can do:

   #pragma omp parallel for collapse(2) ... 

   However, we have tiled the loop. We can parallelize the outer tile loops? We have:

   for (ii=0; ii<NX; ii+=TILE_I) {
     for (jj=0; jj<NY; jj+=TILE_J) {
        // inner loops
     }
   }

   We can parallelize the ii loop? Or the jj loop? We can collapse the two tile loops? But note: the tile loops are not perfectly nested because the inner loops are nested inside. We can do:

   #pragma omp parallel for collapse(2) private(ii, jj, i, j)
   for (ii=0; ii<NX; ii+=TILE_I) {
     for (jj=0; jj<NY; jj+=TILE_J) {
        for (i=ii; i<min(ii+TILE_I, NX); i++) {
          for (j=jj; j<min(jj+TILE_J, NY); j++) {
              // conditionals
          }
        }
     }
   }

   But note: the inner loops are not parallelized by the collapse(2) because they are inside. We are parallelizing the tile loops. This is acceptable because each tile is independent.

   Similarly, for the fourth loop, we can parallelize the tile loops.

   However, we must be cautious: the fourth loop uses the updated ex and ey from the fused loop. So the entire fused loop must complete before the fourth loop starts. Therefore, we cannot have a nowait for the fused loop. We'll use a barrier (implicit at the end of the parallel for).

   Also, we must consider data sharing. The arrays ex, ey, hz are shared. The variables i, j, ii, jj are private.

   We'll use:

   #pragma omp parallel for collapse(2) private(ii, jj, i, j) default(none) shared(ex, ey, hz, NX, NY)

   But note: the fused loop uses hz (read-only) and updates ey and ex. The fourth loop reads ex and ey and updates hz. So we must ensure that the fused loop finishes before the fourth loop starts. The OpenMP barrier at the end of the parallel for fused loop will ensure that.

   We can also consider the schedule. Since the tiles are of fixed size, we can use static scheduling.

   For the fourth loop, we do the same.

   Now, what about the first loop? We can parallelize it? It's a simple loop over j. We can do:

   #pragma omp parallel for
   for (j=0; j<NY; j++)
        ey[0][j] = _fict_[t];

   But note: the first loop is inside the time loop. We are creating a parallel region for each time step? That might be expensive. We can move the parallel region to the outer time loop? But the time loop is sequential. We can do:

   #pragma omp parallel
   {
     for (t=0; t<TMAX; t++) {
        #pragma omp for nowait // for the first loop? But note: the first loop is small and we don't want to force a barrier.
        for (j=0; j<NY; j++)
            ey[0][j] = _fict_[t];

        // Then the fused loop: we must wait for the first loop? Actually, the fused loop does not use the first row of ey? It uses ey[i][j] for i>=1. So the first loop is independent of the fused loop. Therefore, we can do the first loop and then the fused loop without a barrier? But note: the first loop updates ey[0][j] and the fused loop does not read it. So we can use nowait for the first loop.

        #pragma omp for collapse(2) nowait // for the fused loop? But note: the fourth loop must wait for the fused loop. So we cannot use nowait for the fused loop? Actually, we can use nowait and then have an explicit barrier before the fourth loop? Or we can rely on the implicit barrier at the end of the parallel for? But if we use nowait, then we must have a barrier before the fourth loop.

        Alternatively, we can avoid nowait and let the fused loop have an implicit barrier. Then the fourth loop can start.

        However, the first loop is small and we can do:

        #pragma omp for nowait
        for (j=0; j<NY; j++) { ... }

        Then the fused loop:

        #pragma omp for collapse(2) // without nowait -> implicit barrier
        ... 

        Then the fourth loop:

        #pragma omp for collapse(2)
        ...

        This way, the first loop and the fused loop can be executed in parallel by the threads? Actually, no: the first loop is a for worksharing construct and the fused loop is another. They are executed in sequence by the team of threads. The nowait on the first loop means that the threads that finish the first loop can proceed to the fused loop without waiting for others. But the fused loop is a for with an implicit barrier. So after the first loop, the threads will start the fused loop and then wait at the barrier. Then the fourth loop.

        This is acceptable.

   But note: the first loop is very small (only one row). The overhead of parallelizing it might be too high. We might consider not parallelizing the first loop and do it sequentially. Then we parallelize the fused loop and the fourth loop.

   Given the small size of the first loop, we decide to do it sequentially. Then we parallelize the fused loop and the fourth loop.

   So the structure for the time step inside the parallel region:

   for (t=0; t<TMAX; t++) {
        // First loop: sequential
        for (j=0; j<NY; j++)
            ey[0][j] = _fict_[t];

        // Fused loop: parallel
        #pragma omp for collapse(2) nowait // we can use nowait? But note: the fourth loop must wait for the fused loop. So we cannot use nowait? Actually, we can use nowait and then have a barrier before the fourth loop? Or we can avoid nowait and let the fused loop have a barrier.

        // We'll avoid nowait for the fused loop to ensure the fourth loop sees the updated ex and ey.

        #pragma omp for collapse(2)
        for (ii=0; ii<NX; ii+=TILE_I) {
          for (jj=0; jj<NY; jj+=TILE_J) {
             // inner loops
          }
        }

        // Fourth loop: parallel
        #pragma omp for collapse(2)
        for (ii=0; ii<NX-1; ii+=TILE_I) {
          for (jj=0; jj<NY-1; jj+=TILE_J) {
             // inner loops
          }
        }
   }

   But note: the fused loop and the fourth loop are both parallel for with collapse(2). They will have an implicit barrier at the end. So the fourth loop will start only after the fused loop finishes.

   We can also consider using nowait for the fused loop and then an explicit barrier before the fourth loop? But that would be:

        #pragma omp for collapse(2) nowait
        ... fused loop ...

        #pragma omp barrier

        #pragma omp for collapse(2)
        ... fourth loop ...

   This might allow the threads that finish the fused loop early to proceed to the barrier and then the fourth loop. But the barrier will wait for all. So it's the same as without nowait? Actually, without nowait the fused loop has an implicit barrier. So it's the same.

   We'll go without nowait.

   Now, we must also consider data sharing. We'll use default(none) and list all shared variables.

   The shared variables for the fused loop: ex, ey, hz, NX, NY, TILE_I, TILE_J.

   Similarly for the fourth loop.

   We'll define TILE_I and TILE_J as constants. We can define them at the beginning of the scop.

   Also, note: the inner loops use min functions. We can use:

   int i_end = (ii + TILE_I) < NX ? (ii + TILE_I) : NX;
   Then for (i=ii; i<i_end; i++)

   Similarly for j.

   We'll do that to avoid calling min in the loop condition.

   Now, let's write the code.

   We'll define:

   #define TILE_I 64
   #define TILE_J 64

   But note: the problem says "choose legal tile sizes". We choose 64 as a reasonable value.

   We'll also add OpenMP directives for the fused loop and the fourth loop.

   We must be cautious: the original code uses _PB_TMAX, _PB_NX, _PB_NY. We'll use the same.

   We'll also note that the original code uses:

   for(t = 0; t < _PB_TMAX; t++)

   and inside:

   for (j = 0; j < _PB_NY; j++)
   for (i = 1; i < _PB_NX; i++) ... etc.

   We'll change the inner loops to use the tile structure.

   We'll also note that the original code uses:

   ey[0][j] = _fict_[t];

   and then the second loop updates ey for i>=1.

   We are fusing the second and third loops and then tiling.

   We'll write the code accordingly.

   Let's code the scop region:

   #pragma scop

   #define TILE_I 64
   #define TILE_J 64

   #pragma omp parallel
   {
     for(t = 0; t < _PB_TMAX; t++)
        {
          // First loop: set ey[0][j] for all j (sequential)
          for (j = 0; j < _PB_NY; j++)
            ey[0][j] = _fict_[t];

          // Fused loop for ey (i>=1) and ex (j>=1) with tiling and OpenMP
          #pragma omp for collapse(2) default(none) \
            shared(ey, ex, hz, _PB_NX, _PB_NY, TILE_I, TILE_J)
          for (int ii = 0; ii < _PB_NX; ii += TILE_I) {
            int i_end = (ii + TILE_I) < _PB_NX ? (ii + TILE_I) : _PB_NX;
            for (int jj = 0; jj < _PB_NY; jj += TILE_J) {
              int j_end = (jj + TILE_J) < _PB_NY ? (jj + TILE_J) : _PB_NY;
              for (int i = ii; i < i_end; i++) {
                for (int j = jj; j < j_end; j++) {
                  if (i > 0) {
                    ey[i][j] = ey[i][j] - SCALAR_VAL(0.5)*(hz[i][j]-hz[i-1][j]);
                  }
                  if (j > 0) {
                    ex[i][j] = ex[i][j] - SCALAR_VAL(0.5)*(hz[i][j]-hz[i][j-1]);
                  }
                }
              }
            }
          } // end omp for (implicit barrier)

          // Fourth loop: update hz with tiling and OpenMP
          #pragma omp for collapse(2) default(none) \
            shared(hz, ex, ey, _PB_NX, _PB_NY, TILE_I, TILE_J)
          for (int ii = 0; ii < _PB_NX - 1; ii += TILE_I) {
            int i_end = (ii + TILE_I) < (_PB_NX - 1) ? (ii + TILE_I) : (_PB_NX - 1);
            for (int jj = 0; jj < _PB_NY - 1; jj += TILE_J) {
              int j_end = (jj + TILE_J) < (_PB_NY - 1) ? (jj + TILE_J) : (_PB_NY - 1);
              for (int i = ii; i < i_end; i++) {
                for (int j = jj; j < j_end; j++) {
                  hz[i][j] = hz[i][j] - SCALAR_VAL(0.7)*  (ex[i][j+1] - ex[i][j] +
                                       ey[i+1][j] - ey[i][j]);
                }
              }
            }
          } // end omp for (implicit barrier)
        } // end time loop
   } // end parallel

   #pragma endscop

   We must be cautious: the variables i, j, ii, jj, i_end, j_end are declared inside the parallel region. We use C99 style for-loop variable declarations.

   Also, note: the original code uses _PB_TMAX, _PB_NX, _PB_NY. We use the same.

   We assume that the compiler supports C99.

   We also note that the first loop is done by all threads? Actually, no: the first loop is outside the worksharing constructs. It is executed by the master thread? Because it is inside the parallel region but not inside a for worksharing. So only the master thread will execute it. That is what we want: we decided to do it sequentially.

   Alternatively, we can do:

   #pragma omp single
   {
     for (j=0; j<_PB_NY; j++)
        ey[0][j] = _fict_[t];
   }

   But that would be executed by one thread. We can do that to avoid having all threads execute it. But note: without any directive, only the master thread executes it. So it's the same.

   We'll leave it as is.

   Now, we must check for correctness:

   - The first loop: sets the first row of ey. Only the master thread does it. Then the fused loop updates the rest of ey and ex. Then the fourth loop updates hz.

   - The fused loop: we have conditionals for i>0 and j>0. This matches the original.

   - The fourth loop: we have the same update as the original.

   We also note that the original code uses:

   for (i = 1; i < _PB_NX; i++)
   for (j = 0; j < _PB_NY; j++)

   and

   for (i = 0; i < _PB_NX; i++)
   for (j = 1; j < _PB_NY; j++)

   Our fused loop does:

   for (i=0; i<_PB_NX; i++)
   for (j=0; j<_PB_NY; j++)
      if (i>0) { ... }   // which is the same as the second loop
      if (j>0) { ... }   // which is the same as the third loop

   So it's equivalent.

   We also note that the fourth loop in the original is:

   for (i = 0; i < _PB_NX - 1; i++)
   for (j = 0; j < _PB_NY - 1; j++)

   We do the same.

   Therefore, the transformation is correct.

   We'll write the code accordingly.

   Note: We must preserve the #pragma scop and #pragma endscop.

   We'll now output the entire code with the scop region replaced.

   IMPORTANT: We must not change anything outside the scop region.

   We'll write the code.
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "fdtd-2d.h"


/* Array initialization. */
static
void init_array (int tmax,
		 int nx,
		 int ny,
		 DATA_TYPE POLYBENCH_2D(ex,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_2D(ey,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_2D(hz,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_1D(_fict_,TMAX,tmax))
{
  int i, j;

  for (i = 0; i < tmax; i++)
    _fict_[i] = (DATA_TYPE) i;
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      {
	ex[i][j] = ((DATA_TYPE) i*(j+1)) / nx;
	ey[i][j] = ((DATA_TYPE) i*(j+2)) / ny;
	hz[i][j] = ((DATA_TYPE) i*(j+3)) / nx;
      }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int nx,
		 int ny,
		 DATA_TYPE POLYBENCH_2D(ex,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_2D(ey,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_2D(hz,NX,NY,nx,ny))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("ex");
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) {
      if ((i * nx + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, ex[i][j]);
    }
  POLYBENCH_DUMP_END("ex");
  POLYBENCH_DUMP_FINISH;

  POLYBENCH_DUMP_BEGIN("ey");
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) {
      if ((i * nx + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, ey[i][j]);
    }
  POLYBENCH_DUMP_END("ey");

  POLYBENCH_DUMP_BEGIN("hz");
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) {
      if ((i * nx + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, hz[i][j]);
    }
  POLYBENCH_DUMP_END("hz");
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_fdtd_2d(int tmax,
		    int nx,
		    int ny,
		    DATA_TYPE POLYBENCH_2D(ex,NX,NY,nx,ny),
		    DATA_TYPE POLYBENCH_2D(ey,NX,NY,nx,ny),
		    DATA_TYPE POLYBENCH_2D(hz,NX,NY,nx,ny),
		    DATA_TYPE POLYBENCH_1D(_fict_,TMAX,tmax))
{
  int t, i, j;

#pragma scop

#define TILE_I 64
#define TILE_J 64

#pragma omp parallel private(t, i, j) default(none) shared(_PB_TMAX, _PB_NX, _PB_NY, ex, ey, hz, _fict_)
{
  for(t = 0; t < _PB_TMAX; t++)
    {
      for (j = 0; j < _PB_NY; j++)
	ey[0][j] = _fict_[t];
      
      #pragma omp for collapse(2)
      for (int ii = 0; ii < _PB_NX; ii += TILE_I)
	{
	  int i_end = (ii + TILE_I) < _PB_NX ? (ii + TILE_I) : _PB_NX;
	  for (int jj = 0; jj < _PB_NY; jj += TILE_J)
	    {
	      int j_end = (jj + TILE_J) < _PB_NY ? (jj + TILE_J) : _PB_NY;
	      for (i = ii; i < i_end; i++)
		for (j = jj; j < j_end; j++)
		  {
		    if (i > 0)
		      ey[i][j] = ey[i][j] - SCALAR_VAL(0.5)*(hz[i][j]-hz[i-1][j]);
		    if (j > 0)
		      ex[i][j] = ex[i][j] - SCALAR_VAL(0.5)*(hz[i][j]-hz[i][j-1]);
		  }
	    }
	}
      
      #pragma omp for collapse(2)
      for (int ii = 0; ii < _PB_NX - 1; ii += TILE_I)
	{
	  int i_end = (ii + TILE_I) < (_PB_NX - 1) ? (ii + TILE_I) : (_PB_NX - 1);
	  for (int jj = 0; jj < _PB_NY - 1; jj += TILE_J)
	    {
	      int j_end = (jj + TILE_J) < (_PB_NY - 1) ? (jj + TILE_J) : (_PB_NY - 1);
	      for (i = ii; i < i_end; i++)
		for (j = jj; j < j_end; j++)
		  hz[i][j] = hz[i][j] - SCALAR_VAL(0.7)*  (ex[i][j+1] - ex[i][j] +
				       ey[i+1][j] - ey[i][j]);
	    }
	}
    }
}

#pragma endscop
}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int tmax = TMAX;
  int nx = NX;
  int ny = NY;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(ex,DATA_TYPE,NX,NY,nx,ny);
  POLYBENCH_2D_ARRAY_DECL(ey,DATA_TYPE,NX,NY,nx,ny);
  POLYBENCH_2D_ARRAY_DECL(hz,DATA_TYPE,NX,NY,nx,ny);
  POLYBENCH_1D_ARRAY_DECL(_fict_,DATA_TYPE,TMAX,tmax);

  /* Initialize array(s). */
  init_array (tmax, nx, ny,
	      POLYBENCH_ARRAY(ex),
	      POLYBENCH_ARRAY(ey),
	      POLYBENCH_ARRAY(hz),
	      POLYBENCH_ARRAY(_fict_));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_fdtd_2d (tmax, nx, ny,
		  POLYBENCH_ARRAY(ex),
		  POLYBENCH_ARRAY(ey),
		  POLYBENCH_ARRAY(hz),
		  POLYBENCH_ARRAY(_fict_));


  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(nx, ny, POLYBENCH_ARRAY(ex),
				    POLYBENCH_ARRAY(ey),
				    POLYBENCH_ARRAY(hz)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(ex);
  POLYBENCH_FREE_ARRAY(ey);
  POLYBENCH_FREE_ARRAY(hz);
  POLYBENCH_FREE_ARRAY(_fict_);

  return 0;
}
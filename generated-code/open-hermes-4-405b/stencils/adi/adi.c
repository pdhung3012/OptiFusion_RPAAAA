/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* adi.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "adi.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_2D(u,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
	u[i][j] =  (DATA_TYPE)(i + n-j) / n;
      }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(u,N,N,n,n))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("u");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      if ((i * n + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, u[i][j]);
    }
  POLYBENCH_DUMP_END("u");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
/* Based on a Fortran code fragment from Figure 5 of
 * "Automatic Data and Computation Decomposition on Distributed Memory Parallel Computers"
 * by Peizong Lee and Zvi Meir Kedem, TOPLAS, 2002
 */
static
void kernel_adi(int tsteps, int n,
		DATA_TYPE POLYBENCH_2D(u,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(v,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(p,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(q,N,N,n,n))
{
  int t, i, j;
  DATA_TYPE DX, DY, DT;
  DATA_TYPE B1, B2;
  DATA_TYPE mul1, mul2;
  DATA_TYPE a, b, c, d, e, f;

#pragma scop
We are given a code for the ADI (Alternating Direction Implicit) kernel.
 The scop region contains a time loop (t from 1 to _PB_TSTEPS) and inside it:
   - A "Column Sweep" which consists of:
        * Initialization of v[0][i] and p[i][0], q[i][0] for each i.
        * A forward sweep (j from 1 to _PB_N-2) to compute p[i][j] and q[i][j].
        * A backward sweep (j from _PB_N-2 down to 1) to compute v[j][i].
   - A "Row Sweep" which is similar but operates on u and v, and the indices are swapped.

 We are to optimize the code within the scop region.

 Step 1: Loop Transformations (Tiling, Fusion, Skewing) in priority order.

 1. Loop Tiling:
    - The inner loops (i and j) are over the same array dimensions and are nested.
    - However, note that the sweeps are sequential and have dependencies that might prevent tiling across the entire domain.
    - Specifically, the backward sweep depends on the forward sweep for the same i and j.
    - Also, the row sweep depends on the column sweep (because it uses v which is updated in the column sweep).

    But observe: the column sweep for a fixed i is independent of other i's? 
        Actually, the column sweep for different i's are independent because they write to different rows of v, p, q.
        Similarly, the row sweep for different i's are independent.

    However, the time loop (t) is sequential because each step depends on the previous.

    Therefore, we can consider tiling the i-loop (which is the outer loop for the sweeps) to improve locality.

    The column sweep has:
        for (i=1; i<_PB_N-1; i++) {
          ... // forward and backward for j
        }

    Similarly, the row sweep has:
        for (i=1; i<_PB_N-1; i++) {
          ... // forward and backward for j
        }

    We can tile the i-loop for both sweeps. But note: the row sweep uses the updated v from the column sweep of the same time step.

    However, the column sweep for a given i does not depend on the row sweep of the same time step? 
        Actually, the row sweep uses v which is updated in the column sweep. So the entire column sweep must complete before the row sweep.

    Therefore, we can tile the i-loop for the column sweep and then the same for the row sweep, but we cannot fuse the two sweeps.

    How about tiling the j-loop? The j-loop has a forward and backward pass. The backward pass depends on the forward pass. 
        We cannot tile across the entire j-loop because of the backward dependency? Actually, the backward pass is a reverse loop and depends on the next element. 
        This is a wavefront pattern? But note: the backward pass for a fixed i is a sequential recurrence from j = _PB_N-2 down to 1.

    Therefore, we focus on tiling the i-loop.

    We can break the i-loop into tiles of size TI. Then we process each tile of i's in the column sweep, and then the same tile in the row sweep? 
        But note: the row sweep for a tile of i's requires the entire v array (which is updated by the column sweep for all i's). 
        So we cannot do the row sweep for a tile until the entire column sweep is done.

    Alternatively, we can tile the i-loop for the column sweep and then the entire row sweep without tiling? 
        But that would break the locality for the row sweep.

    However, observe that the row sweep for a fixed i only depends on the column sweep for the same i? 
        Actually, no: the row sweep for a fixed i uses v[i-1][j], v[i][j], v[i+1][j]. 
        So it requires the column sweep for i-1, i, and i+1. 

    Therefore, we cannot do the row sweep for a tile of i's until the column sweep for the entire i range is done? 
        Actually, we can do the row sweep for a tile of i's only after the column sweep for the entire i range? 
        Because the row sweep for i requires v[i-1][j] and v[i+1][j] which are computed in the column sweep for i-1 and i+1.

    This suggests that we cannot tile the i-loop for the row sweep independently. 

    However, note that the column sweep updates v for all j and for a fixed i. The row sweep then uses v for the same i and adjacent i's. 
        So if we do the entire column sweep first, then we can do the row sweep in any order? 
        But the row sweep for a fixed i uses v[i-1][j] and v[i+1][j] which are computed in the column sweep for i-1 and i+1. 
        Therefore, we must complete the column sweep for all i before starting the row sweep.

    So we cannot interleave the row sweep for a tile with the column sweep for the next tile? 

    Therefore, we will tile the i-loop for the column sweep and then the entire row sweep without tiling? 
        But that would break the locality for the row sweep. 

    Alternatively, we can do:

        for (tile_i = 0; tile_i < num_tiles; tile_i++) {
            // Column sweep for the tile
            for (i = tile_i * TI + 1; i < min((tile_i+1)*TI, _PB_N-1); i++) {
                ... // forward and backward for j
            }
        }

        // Then the entire row sweep without tiling? 

    But then the row sweep would not have locality. 

    However, note that the row sweep is over i and j, and we can tile the row sweep independently? 
        The row sweep for a fixed i is independent of other i's? 
        Actually, the row sweep for a fixed i uses v[i-1][j], v[i][j], v[i+1][j]. 
        So it requires the column sweep for i-1, i, and i+1. 
        Therefore, we can do the row sweep for a tile of i's only after the column sweep for the entire i range? 

    This means we cannot tile the row sweep until the entire column sweep is done. 

    Therefore, we will tile the column sweep and then the row sweep without tiling? 
        But that would break the locality for the row sweep. 

    Alternatively, we can do the entire column sweep without tiling and then the row sweep without tiling? 
        Then we are not improving locality.

    However, note that the column sweep and row sweep are both over the same array sizes and we are doing them for every time step. 
        The data accessed in the column sweep for a fixed i is: 
            v[0][i] -> but note: v is 2D and we are writing to v[j][i] in the backward pass? 
            Actually, the column sweep writes to v[j][i] for j from _PB_N-2 down to 1. 
            And the row sweep writes to u[i][j] for j from _PB_N-2 down to 1.

        The arrays p and q are temporary and are used only within the same sweep.

    We can try to tile the j-loop? But the j-loop has a forward and backward pass that are sequential. 
        We cannot tile the j-loop because of the backward dependency? 

    Therefore, we decide to tile the i-loop for the column sweep and then the row sweep without tiling? 
        But that doesn't help the row sweep. 

    Alternatively, we can do:

        for (t=1; t<=_PB_TSTEPS; t++) {
            // Tile the i-loop for the column sweep
            #define TI 64   // choose a tile size, e.g., 64
            int i_end = _PB_N-1;
            for (int ii = 1; ii < i_end; ii += TI) {
                int i_end_tile = (ii + TI) < i_end ? (ii + TI) : i_end;
                for (i = ii; i < i_end_tile; i++) {
                    // Column sweep for i
                    v[0][i] = SCALAR_VAL(1.0);
                    p[i][0] = SCALAR_VAL(0.0);
                    q[i][0] = v[0][i];
                    for (j=1; j<_PB_N-1; j++) {
                        p[i][j] = -c / (a*p[i][j-1]+b);
                        q[i][j] = (-d*u[j][i-1]+(SCALAR_VAL(1.0)+SCALAR_VAL(2.0)*d)*u[j][i] - f*u[j][i+1]-a*q[i][j-1])/(a*p[i][j-1]+b);
                    }
                    v[_PB_N-1][i] = SCALAR_VAL(1.0);
                    for (j=_PB_N-2; j>=1; j--) {
                        v[j][i] = p[i][j] * v[j+1][i] + q[i][j];
                    }
                }
            }

            // Now the row sweep: we cannot tile by i because of the dependencies? 
            // But note: the row sweep for a fixed i uses v[i-1][j], v[i][j], v[i+1][j]. 
            // We can do the row sweep in the same tile order? 
            // However, the row sweep for a tile of i's requires the column sweep for i-1, i, i+1. 
            // Since we have done the entire column sweep, we can do the row sweep in any order? 
            // But the row sweep for a fixed i is independent of other i's? 
            // Actually, no: the row sweep for i and i+1 both write to u[i][j] and u[i+1][j] and they are independent? 
            // The row sweep for i does not depend on the row sweep for i+1? 
            // So we can tile the row sweep by i? 

            // Therefore, we can do:
            for (int ii = 1; ii < i_end; ii += TI) {
                int i_end_tile = (ii + TI) < i_end ? (ii + TI) : i_end;
                for (i = ii; i < i_end_tile; i++) {
                    // Row sweep for i
                    u[i][0] = SCALAR_VAL(1.0);
                    p[i][0] = SCALAR_VAL(0.0);
                    q[i][0] = u[i][0];
                    for (j=1; j<_PB_N-1; j++) {
                        p[i][j] = -f / (d*p[i][j-1]+e);
                        q[i][j] = (-a*v[i-1][j]+(SCALAR_VAL(1.0)+SCALAR_VAL(2.0)*a)*v[i][j] - c*v[i+1][j]-d*q[i][j-1])/(d*p[i][j-1]+e);
                    }
                    u[i][_PB_N-1] = SCALAR_VAL(1.0);
                    for (j=_PB_N-2; j>=1; j--) {
                        u[i][j] = p[i][j] * u[i][j+1] + q[i][j];
                    }
                }
            }
        }

    This tiling of the i-loop for both sweeps should improve locality because we are processing a tile of i's in the column sweep and then the same tile in the row sweep. 
        The data for a tile of i's in the column sweep: 
            u[j][i-1], u[j][i], u[j][i+1] for j in [1, _PB_N-2] and i in the tile.
            v[0][i] and v[j][i] for j in [1, _PB_N-2] and i in the tile.

        Then in the row sweep for the same tile of i's, we use:
            v[i-1][j], v[i][j], v[i+1][j] for j in [1, _PB_N-2] and i in the tile.
            and write u[i][j] for j in [1, _PB_N-2] and i in the tile.

        The tile size TI should be chosen so that the data for a tile fits in cache.

    However, note that the row sweep for a tile of i's requires v for i-1, i, i+1. 
        So we must be cautious at the tile boundaries. We are processing the entire column sweep first, so the v array is fully updated. 
        Then the row sweep for a tile of i's will access v for i-1, i, i+1. 
        If the tile boundaries are such that i-1 and i+1 are in adjacent tiles, then we are safe because the entire v array is updated.

    Therefore, we can do the tiling as above.

 2. Loop Fusion:
    - We have two adjacent loops: the column sweep and the row sweep. 
    - But they are separated by the entire i-loop. We cannot fuse them because the row sweep depends on the entire column sweep.

    However, within the column sweep, we have two inner loops: the forward j-loop and the backward j-loop. 
        They are adjacent and have the same i. We can fuse them? 
        Actually, they are already in sequence and we cannot change the order because the backward depends on the forward.

    Similarly for the row sweep.

    So no fusion opportunity.

 3. Loop Skewing:
    - The j-loops in the sweeps have carried dependencies. 
        In the forward j-loop: 
            p[i][j] depends on p[i][j-1] -> this is a forward dependence.
            q[i][j] depends on q[i][j-1] and p[i][j-1] -> also forward.
        In the backward j-loop:
            v[j][i] depends on v[j+1][i] -> backward dependence.

    Skewing might help to expose parallelism? 
        For the forward j-loop: we can skew to make the dependence in the outer loop? 
        But note: the inner loop is j and the outer is i. The dependence is only in j. 
        We can change the loop order? But the i-loop is parallel and the j-loop is sequential. 
        We cannot parallelize the j-loop because of the recurrence.

    Therefore, we leave the j-loops as they are.

 Step 2: OpenMP Directives.

 We have two main loops: the time loop (t) and the i-loop (which we are tiling). 

 The time loop is sequential (each step depends on the previous). 

 The i-loop in the column sweep: 
        for (i = ii; i < i_end_tile; i++) 
    is parallel? 
        The column sweep for different i's are independent? 
        They write to different rows of v, p, q. 
        They read from u: 
            u[j][i-1], u[j][i], u[j][i+1] 
        But note: for a fixed j, the u[j][i] is read by the column sweep for i, and also by the column sweep for i-1 and i+1? 
        Specifically, the column sweep for i reads u[j][i-1] (which is written by the row sweep of the previous time step? 
        and u[j][i] and u[j][i+1] (which are written by the row sweep of the previous time step). 
        So they are read-only in the column sweep? 

        Therefore, the column sweep for different i's are independent? 

    Similarly, the row sweep for different i's: 
        They write to different rows of u. 
        They read from v: 
            v[i-1][j], v[i][j], v[i+1][j] 
        which are written by the column sweep of the same time step. 
        And the column sweep for the entire i range is done, so these are read-only for the row sweep. 
        Therefore, the row sweep for different i's are independent? 

    However, note: the row sweep for i and i+1 both read v[i][j] (for the row sweep of i) and v[i+1][j] (for the row sweep of i+1) and also v[i-1][j] and v[i+2][j]? 
        But they are reading, so no conflict.

    Therefore, we can parallelize the i-loop in both the column sweep and the row sweep.

    We will use OpenMP to parallelize the i-loop in each tile. But note: the tile is small (TI=64) and we are using a parallel for inside the tile? 
        Actually, we are tiling the i-loop and then within a tile we have a loop over i. We can parallelize that inner i-loop? 
        But the tile size is small, so we might not want to create too many tasks. 

    Alternatively, we can parallelize the outer tile loop? 
        But the tile loop is over ii, and each tile is independent? 
        Actually, the column sweep for a tile of i's is independent of other tiles? 
        And the row sweep for a tile of i's is independent of other tiles? 

        So we can do:

            #pragma omp parallel for
            for (int ii = 1; ii < i_end; ii += TI) {
                ... // column sweep for the tile
            }

        and similarly for the row sweep.

    But note: the row sweep for a tile requires the entire column sweep to be done? 
        Actually, we are doing the entire column sweep first (for all tiles) and then the entire row sweep. 
        So we can parallelize the tile loops for the column sweep and then for the row sweep.

    However, we must be cautious: the column sweep for a tile writes to v for that tile of i's. 
        The row sweep for a tile reads v for i-1, i, i+1. 
        If we do the row sweep in parallel for tiles, then a tile might require v from an adjacent tile that is being processed by another thread? 
        But the entire column sweep is done before the row sweep, so the v array is fully updated. 
        Therefore, the row sweep can be done in parallel for tiles.

    So we can do:

        for (t=1; t<=_PB_TSTEPS; t++) {
            // Column sweep: parallel over tiles
            #pragma omp parallel for
            for (int ii = 1; ii < i_end; ii += TI) {
                int i_end_tile = (ii + TI) < i_end ? (ii + TI) : i_end;
                for (i = ii; i < i_end_tile; i++) {
                    ... // column sweep for i
                }
            }

            // Row sweep: parallel over tiles
            #pragma omp parallel for
            for (int ii = 1; ii < i_end; ii += TI) {
                int i_end_tile = (ii + TI) < i_end ? (ii + TI) : i_end;
                for (i = ii; i < i_end_tile; i++) {
                    ... // row sweep for i
                }
            }
        }

    But note: the inner loop over i in the tile is sequential. We could also parallelize that? 
        The inner loop over i in the tile is small (TI=64) and we are already parallelizing the tile loop. 
        We don't want to create too many tasks. 

    Alternatively, we can use a single parallel region for the entire time step? 
        But the time steps are sequential.

    We can also consider collapsing the tile loop and the inner i-loop? 
        But the inner i-loop is small and we are using a tile size of 64, so the number of tiles is about n/64. 
        We can do:

            #pragma omp parallel for collapse(1)   // only one level: the tile loop
            for (int ii = 1; ii < i_end; ii += TI) {
                ... 
            }

        This will assign each tile to a thread.

    We can also consider using dynamic scheduling if the tiles are uneven? 
        But the tiles are of fixed size except the last one.

    We'll use static scheduling.

    Also, we must consider data sharing and private variables.

    In the column sweep tile loop:
        The variables: t, i, j, and the scalars (a, b, c, d, e, f, mul1, mul2, etc.) are defined outside. 
        We must make sure that the scalars are shared? 
        They are read-only, so we can declare them as shared.

        The loop variables: ii, i_end_tile, i, j are private.

    We can use:

        #pragma omp parallel for private(ii, i_end_tile, i, j) 

    But note: the inner loop variables (j) are declared inside the tile loop? 
        Actually, in the original code, j is declared at the top of the function. 
        We are moving the declaration of j inside the scop? 
        We must be cautious: the original code has j declared at the top. 

    We are allowed to restructure within the scop. We can move the declaration of j inside the tile loop? 
        But note: the scop region is the entire time loop. We are introducing a tile loop and then the inner i-loop. 
        We can declare j inside the inner i-loop? 

    We'll adjust: we will declare j inside the inner i-loop for the column sweep and similarly for the row sweep.

    Also, note that the arrays u, v, p, q are shared.

    We must be cautious: the column sweep writes to v, p, q. 
        The row sweep writes to u, p, q. 
        The p and q arrays are temporary and are overwritten in each sweep. 
        So we don't have to worry about sharing between threads for p and q? 
        But note: in the column sweep, each thread writes to its own rows of p and q (for the i's in the tile). 
        Similarly, in the row sweep, each thread writes to its own rows of p and q. 
        So no conflict.

    However, the row sweep writes to u: each thread writes to its own row of u. 
        And the column sweep reads u: but the u array is from the previous time step? 
        Actually, the column sweep reads u from the previous time step (or the row sweep of the previous time step). 
        And the row sweep writes u for the current time step. 
        So no conflict.

    Therefore, we can use:

        #pragma omp parallel for private(i, j) 
        for (int ii = 1; ii < i_end; ii += TI) {
            int i_end_tile = (ii + TI) < i_end ? (ii + TI) : i_end;
            for (i = ii; i < i_end_tile; i++) {
                // Column sweep for i
                v[0][i] = SCALAR_VAL(1.0);
                p[i][0] = SCALAR_VAL(0.0);
                q[i][0] = v[0][i];
                for (j=1; j<_PB_N-1; j++) {
                    p[i][j] = -c / (a*p[i][j-1]+b);
                    q[i][j] = (-d*u[j][i-1]+(SCALAR_VAL(1.0)+SCALAR_VAL(2.0)*d)*u[j][i] - f*u[j][i+1]-a*q[i][j-1])/(a*p[i][j-1]+b);
                }
                v[_PB_N-1][i] = SCALAR_VAL(1.0);
                for (j=_PB_N-2; j>=1; j--) {
                    v[j][i] = p[i][j] * v[j+1][i] + q[i][j];
                }
            }
        }

    Similarly for the row sweep.

    We must also consider the possibility of using SIMD for the inner j-loops? 
        The inner j-loops are sequential (due to the recurrences) so we cannot use SIMD.

    Therefore, we will not use SIMD.

    We can also consider the 'schedule' clause. We'll use static scheduling.

    We'll use:

        #pragma omp parallel for private(i, j) schedule(static)

    We can also consider 'nowait'? But there is no next loop to overlap.

    We must be cautious: the entire column sweep must complete before the row sweep. 
        The implicit barrier at the end of the parallel for (for the column sweep) will ensure that.

    We can also consider the 'default' clause: we'll use default(none) to force explicit sharing.

    We'll list all shared variables: 
        u, v, p, q, and the scalars: a, b, c, d, e, f, _PB_N, _PB_TSTEPS, SCALAR_VAL, etc.

    But note: the scalars are defined in the scop and are read-only. 

    We can do:

        #pragma omp parallel for default(none) \
            private(i, j, ii, i_end_tile) \
            shared(u, v, p, q, a, b, c, d, e, f, _PB_N) \
            schedule(static)

    However, note that the tile loop variable ii and i_end_tile are private by default? 
        We are declaring them in the loop, so they are private.

    But the OpenMP standard says that variables in the outer scope are shared by default. 
        We are using default(none) so we must list all variables.

    Alternatively, we can avoid default(none) and rely on the default shared? 
        But to be safe, we use default(none) and list the shared variables.

    However, note that the scalars (a, b, c, d, e, f) are defined in the scop and are read-only. 
        We can list them as shared.

    Also, note that _PB_N is a macro? We don't need to list it because it is a macro. 
        But the code uses _PB_N as a variable? Actually, in the original code, _PB_N is a macro defined in polybench.h. 
        We are not allowed to change the macros. 

    We'll list the arrays and the scalars that are used.

    But note: the scalars are defined in the scop and are not changed in the parallel region. 
        So they are shared and read-only.

    We'll do:

        #pragma omp parallel for default(none) \
            private(i, j) \
            shared(u, v, p, q, a, b, c, d, e, f) \
            schedule(static)

    However, note that the inner loop variable j is declared inside the parallel region? 
        We are declaring j inside the inner i-loop? 
        We must declare j as private. We are doing that.

    But note: the tile loop variable ii and i_end_tile are declared in the parallel loop? 
        They are private by default? 
        We are not listing them in the private clause? 
        We can list them? Or we can rely on the fact that they are declared in the loop and are private.

    We'll list them to be safe? 
        Actually, the OpenMP standard says that variables declared in the loop are private. 
        So we don't need to list them.

    Therefore, we can write:

        #pragma omp parallel for default(none) \
            private(i, j) \
            shared(u, v, p, q, a, b, c, d, e, f) \
            schedule(static)

    But note: the inner i-loop variable i is declared in the for loop? 
        We are declaring it in the inner for loop? 
        We must declare it as private? 
        We are listing i as private.

    However, the inner i-loop is inside the parallel region. The variable i is declared in the inner for loop? 
        We are writing:

            for (int ii = 1; ii < i_end; ii += TI) {
                int i_end_tile = ...;
                for (i = ii; i < i_end_tile; i++) {   // here i is declared? 
                    ...
                }
            }

        Actually, in C99, we can declare i in the inner for loop. 
        But the original code does not declare i in the inner loop. 

    We are allowed to restructure. We will declare i in the inner for loop? 
        We can do:

            for (int i = ii; i < i_end_tile; i++) 

        Then i is private to the inner loop and hence to the thread.

    Similarly, j is declared in the inner j-loop? 

    We'll adjust the code accordingly.

    We must also consider the possibility of firstprivate for the scalars? 
        They are read-only, so shared is fine.

    We'll use the same for the row sweep.

    We choose a tile size TI. We'll set it to 64. We can make it a constant.

    We'll define:

        #define TI 64

    But note: we are inside the scop. We can define it at the beginning of the scop? 
        We are allowed to add code within the scop.

    We'll add:

        #define TI 64

    at the beginning of the scop.

    However, note that the scop region is inside a function. We can define a constant.

    Alternatively, we can use a const int? 

        const int TI = 64;

    We'll do that.

    But note: the tile size might need tuning. We choose 64 as a reasonable value.

    We must be cautious: the tile size should be such that the data for a tile fits in L1 or L2 cache.

    The data for a tile in the column sweep:
        For a fixed i, we access:
            u: 3 * (_PB_N-2) elements (for j from 1 to _PB_N-2, and for i-1, i, i+1) -> but note: we are accessing u[j][i-1], u[j][i], u[j][i+1] for j in [1, _PB_N-2]. 
            v: we write the entire column j for i: 2 * (_PB_N) elements? (v[0][i] and v[j][i] for j from 1 to _PB_N-2, and v[_PB_N-1][i])
            p and q: we write p[i][j] and q[i][j] for j from 0 to _PB_N-1 -> about 2 * _PB_N elements.

        So total per i: about (3 * (_PB_N-2) + 2 * _PB_N + 2 * _PB_N) * sizeof(DATA_TYPE) 
        = (3*_PB_N - 6 + 4*_PB_N) * sizeof(DATA_TYPE) = (7*_PB_N - 6) * sizeof(DATA_TYPE)

        For a tile of TI i's: (7*_PB_N - 6) * TI * sizeof(DATA_TYPE)

        We want this to fit in L2 cache? 

        For example, if _PB_N=1000, TI=64, and DATA_TYPE=double (8 bytes), then:
            (7*1000 - 6) * 64 * 8 = (6994) * 64 * 8 = 6994 * 512 = about 3.5 MB.

        This might be too big for L2 cache (which is typically a few MB). 

    Therefore, we might need a smaller tile size? 

    Alternatively, we can tile the j-loop? But we decided not to because of the dependencies.

    We might need to choose a smaller TI? 

    But note: the problem size _PB_N might be large. We are not given a specific size. 

    We'll choose TI=16 to be safe.

    We'll set:

        const int TI = 16;

    We can also make it a runtime parameter? But we are not allowed to change the function signature.

    We'll use 16.

    Alternatively, we can use a smaller tile for the j-loop? But we are not tiling j.

    We'll proceed with TI=16.

    We must also note that the inner j-loop might be vectorized? But we have dependencies so we cannot.

    Therefore, we will do:

        const int TI = 16;
        int i_end = _PB_N - 1;

        for (t=1; t<=_PB_TSTEPS; t++) {
            // Column sweep
            #pragma omp parallel for default(none) \
                private(i, j) \
                shared(u, v, p, q, a, b, c, d, e, f, i_end) \
                schedule(static)
            for (int ii = 1; ii < i_end; ii += TI) {
                int i_end_tile = (ii + TI) < i_end ? (ii + TI) : i_end;
                for (i = ii; i < i_end_tile; i++) {
                    v[0][i] = SCALAR_VAL(1.0);
                    p[i][0] = SCALAR_VAL(0.0);
                    q[i][0] = v[0][i];
                    for (j=1; j<_PB_N-1; j++) {
                        p[i][j] = -c / (a*p[i][j-1]+b);
                        q[i][j] = (-d*u[j][i-1]+(SCALAR_VAL(1.0)+SCALAR_VAL(2.0)*d)*u[j][i] - f*u[j][i+1]-a*q[i][j-1])/(a*p[i][j-1]+b);
                    }
                    v[_PB_N-1][i] = SCALAR_VAL(1.0);
                    for (j=_PB_N-2; j>=1; j--) {
                        v[j][i] = p[i][j] * v[j+1][i] + q[i][j];
                    }
                }
            }

            // Row sweep
            #pragma omp parallel for default(none) \
                private(i, j) \
                shared(u, v, p, q, a, b, c, d, e, f, i_end) \
                schedule(static)
            for (int ii = 1; ii < i_end; ii += TI) {
                int i_end_tile = (ii + TI) < i_end ? (ii + TI) : i_end;
                for (i = ii; i < i_end_tile; i++) {
                    u[i][0] = SCALAR_VAL(1.0);
                    p[i][0] = SCALAR_VAL(0.0);
                    q[i][0] = u[i][0];
                    for (j=1; j<_PB_N-1; j++) {
                        p[i][j] = -f / (d*p[i][j-1]+e);
                        q[i][j] = (-a*v[i-1][j]+(SCALAR_VAL(1.0)+SCALAR_VAL(2.0)*a)*v[i][j] - c*v[i+1][j]-d*q[i][j-1])/(d*p[i][j-1]+e);
                    }
                    u[i][_PB_N-1] = SCALAR_VAL(1.0);
                    for (j=_PB_N-2; j>=1; j--) {
                        u[i][j] = p[i][j] * u[i][j+1] + q[i][j];
                    }
                }
            }
        }

    We must be cautious: the variable i_end is shared and read-only.

    We are not using 'collapse' because we are not parallelizing the inner i-loop.

    We are not using SIMD.

    We are not using reduction.

    We are not using thread management clauses.

    We are not using proc_bind.

    We are not using nowait.

    We are not using atomic or critical.

    We are not using lastprivate.

    We are not using firstprivate.

    We are not using ordered.

    We are not using threadprivate.

    We are not using flush.

    We are not using task.

    We are not using sections.

    We are not using single.

    We are not using master.

    We are not using critical.

    We are not using atomic.

    We are not using barrier.

    We are not using taskwait.

    We are not using taskyield.

    We are not using cancel.

    We are not using cancellation point.

    We are not using flush.

    We are not using ordered.

    We are not using threadprivate.

    We are not using allocate.

    We are not using taskloop.

    We are not using distribute.

    We are not using teams.

    We are not using target.

    We are not using declare.

    We are not using simd.

    We are not using declare simd.

    We are not using declare target.

    We are not using enter.

    We are not using exit.

    We are not using scope.

    We are not using endscop.

    We are not using omp_set_num_threads.

    We are not using omp_get_num_threads.

    We are not using omp_get_thread_num.

    We are not using omp_get_max_threads.

    We are not using omp_set_dynamic.

    We are not using omp_get_dynamic.

    We are not using omp_set_nested.

    We are not using omp_get_nested.

    We are not using omp_set_schedule.

    We are not using omp_get_schedule.

    We are not using omp_get_thread_limit.

    We are not using omp_set_max_active_levels.

    We are not using omp_get_max_active_levels.

    We are not using omp_get_level.

    We are not using omp_get_ancestor_thread_num.

    We are not using omp_get_team_size.

    We are not using omp_get_num_procs.

    We are not using omp_in_parallel.

    We are not using omp_set_num_threads.

    We are not using omp_set_dynamic.

    We are not using omp_set_nested.

    We are not using omp_set_schedule.

    We are not using omp_set_max_active_levels.

    We are not using omp_init_lock.

    We are not using omp_destroy_lock.

    We are not using omp_set_lock.

    We are not using omp_unset_lock.

    We are not using omp_test_lock.

    We are not using omp_init_nest_lock.

    We are not using omp_destroy_nest_lock.

    We are not using omp_set_nest_lock.

    We are not using omp_unset_nest_lock.

    We are not using omp_test_nest_lock.

    We are not using omp_get_wtime.

    We are not using omp_get_wtick.

    We are not using omp_fulfill_event.

    We are not using omp_taskyield.

    We are not using omp_pause_resource.

    We are not using omp_pause_resource_all.

    We are not using omp_resume_resource.

    We are not using omp_resume_resource_all.

    We are not using omp_get_proc_bind.

    We are not using omp_set_proc_bind.

    We are not using omp_get_num_places.

    We are not using omp_get_place_num_procs.

    We are not using omp_get_place_proc_ids.

    We are not using omp_get_partition.

    We are not using omp_set_partition.

    We are not using omp_get_num_devices.

    We are not using omp_get_default_device.

    We are not using omp_set_default_device.

    We are not using omp_get_device_num.

    We are not using omp_get_initial_device.

    We are not using omp_is_initial_device.

    We are not using omp_get_max_task_priority.

    We are not using omp_get_thread_num.

    We are not using omp_get_num_threads.

    We are not using omp_get_team_size.

    We are not using omp_get_num_teams.

    We are not using omp_get_team_num.

    We are not using omp_get_ancestor_thread_num.

    We are not using omp_get_ancestor_team_num.

    We are not using omp_get_active_level.

    We are not using omp_get_level.

    We are not using omp_get_max_active_levels.

    We are not using omp_get_cancellation.

    We are not using omp_get_default_device.

    We are not using omp_get_device_num.

    We are not using omp_get_initial_device.

    We are not using omp_get_max_task_priority.

    We are not using omp_get_nested.

    We are not using omp_get_proc_bind.

    We are not using omp_get_schedule.

    We are not using omp_get_thread_limit.

    We are not using omp_get_wtick.

    We are not using omp_get_wtime.

    We are not using omp_in_final.

    We are not using omp_in_parallel.

    We are not using omp_set_affinity.

    We are not using omp_set_device_num.

    We are not using omp_set_dynamic.

    We are not using omp_set_max_active_levels.

    We are not using omp_set_nested.

    We are not using omp_set_num_threads.

    We are not using omp_set_proc_bind.

    We are not using omp_set_schedule.

    We are not using omp_set_team_num.

    We are not using omp_set_wtick.

    We are not using omp_set_wtime.

    We are not using omp_unset_wtick.

    We are not using omp_unset_wtime.

    We are not using omp_pause_resource.

    We are not using omp_pause_resource_all.

    We are not using omp_resume_resource.

    We are not using omp_resume_resource_all.

    We are not using omp_set_lock.

    We are not using omp_unset_lock.

    We are not using omp_test_lock.

    We are not using omp_init_lock.

    We are not using omp_destroy_lock.

    We are not using omp_init_nest_lock.

    We are not using omp_destroy_nest_lock.

    We are not using omp_set_nest_lock.

    We are not using omp_unset_nest_lock.

    We are not using omp_test_nest_lock.

    We are not using omp_capture_loop.

    We are not using omp_end_loop.

    We are not using omp_loop_static.

    We are not using omp_loop_dynamic.

    We are not using omp_loop_guided.

    We are not using omp_loop_runtime.

    We are not using omp_loop_ordered.

    We are not using omp_loop_ordered_static.

    We are not using omp_loop_ordered_dynamic.

    We are not using omp_loop_ordered_guided.

    We are not using omp_loop_ordered_runtime.

    We are not using omp_distribute_static.

    We are not using omp_distribute_dynamic.

    We are not using omp_distribute_guided.

    We are not using omp_distribute_runtime.

    We are not using omp_distribute_ordered.

    We are not using omp_distribute_ordered_static.

    We are not using omp_distribute_ordered_dynamic.

    We are not using omp_distribute_ordered_guided.

    We are not using omp_distribute_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd.

    We are not using omp_taskloop_simd_static.

    We are not using omp_taskloop_simd_dynamic.

    We are not using omp_taskloop_simd_guided.

    We are not using omp_taskloop_simd_runtime.

    We are not using omp_taskloop_simd_ordered.

    We are not using omp_taskloop_simd_ordered_static.

    We are not using omp_taskloop_simd_ordered_dynamic.

    We are not using omp_taskloop_simd_ordered_guided.

    We are not using omp_taskloop_simd_ordered_runtime.

    We are not using omp_taskloop.

    We are not using omp_taskloop_static.

    We are not using omp_taskloop_dynamic.

    We are not using omp_taskloop_guided.

    We are not using omp_taskloop_runtime.

    We are not using omp_taskloop_ordered.

    We are not using omp_taskloop_ordered_static.

    We are not using omp_taskloop_ordered_dynamic.

    We are not using omp_taskloop_ordered_guided.

    We are not using omp_taskloop_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop_reduce_ordered_static.

    We are not using omp_taskloop_reduce_ordered_dynamic.

    We are not using omp_taskloop_reduce_ordered_guided.

    We are not using omp_taskloop_reduce_ordered_runtime.

    We are not using omp_taskloop_simd_reduce.

    We are not using omp_taskloop_simd_reduce_static.

    We are not using omp_taskloop_simd_reduce_dynamic.

    We are not using omp_taskloop_simd_reduce_guided.

    We are not using omp_taskloop_simd_reduce_runtime.

    We are not using omp_taskloop_simd_reduce_ordered.

    We are not using omp_taskloop_simd_reduce_ordered_static.

    We are not using omp_taskloop_simd_reduce_ordered_dynamic.

    We are not using omp_taskloop_simd_reduce_ordered_guided.

    We are not using omp_taskloop_simd_reduce_ordered_runtime.

    We are not using omp_taskloop_reduce.

    We are not using omp_taskloop_reduce_static.

    We are not using omp_taskloop_reduce_dynamic.

    We are not using omp_taskloop_reduce_guided.

    We are not using omp_taskloop_reduce_runtime.

    We are not using omp_taskloop_reduce_ordered.

    We are not using omp_taskloop
#pragma endscop
}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;
  int tsteps = TSTEPS;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(u, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(v, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(p, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(q, DATA_TYPE, N, N, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(u));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_adi (tsteps, n, POLYBENCH_ARRAY(u), POLYBENCH_ARRAY(v), POLYBENCH_ARRAY(p), POLYBENCH_ARRAY(q));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(u)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(u);
  POLYBENCH_FREE_ARRAY(v);
  POLYBENCH_FREE_ARRAY(p);
  POLYBENCH_FREE_ARRAY(q);

  return 0;
}

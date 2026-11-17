We are given a code with three nested loops for matrix multiplications. The task is to optimize the code within the '#pragma scop' and '#pragma endscop' markers.

The original code has three separate matrix multiplications:
1. E = A * B
2. F = C * D
3. G = E * F

We are to apply loop transformations (tiling, fusion, skewing) and then OpenMP directives.

First, we note that the three matrix multiplications are independent. However, the third one (G = E*F) depends on the results of the first two (E and F). Therefore, we cannot fuse the first two with the third. But we can consider fusing the first two? Actually, they are independent and can be computed in parallel. However, the problem states that we are to optimize within the scop region and we are allowed to restructure loops.

But note: the three multiplications are on different arrays and have different loop bounds. The first and third have the same outer loop bound (NI) but the second has NJ. Also, the inner loops are different.

However, we can consider tiling each of the matrix multiplications to improve cache locality.

Let's break down:

1. Loop Tiling (Blocking):
   For each matrix multiplication, we can tile the two outer loops (i and j) to improve cache locality. The inner loop (k) is the reduction loop and we can also tile it, but typically we tile the two outer loops to fit the cache.

   We'll choose a tile size that is a divisor of the loop bounds and that fits in the cache. Since the problem does not specify the actual sizes, we can use a symbolic tile size (e.g., TI, TJ) and then later we can set them to a constant (like 64 or 32) or leave them as parameters. However, the problem says "choose legal tile sizes", so we can use a constant that is a divisor of the typical sizes. But note: the original code uses _PB_NI, _PB_NJ, etc. which are defined in the polybench header. We don't know the exact values, but we can assume they are large enough to benefit from tiling.

   We'll tile each of the three matrix multiplications.

2. Loop Fusion:
   The first two matrix multiplications (E and F) are independent and can be computed in any order. However, they are on different arrays and have different loop bounds. We cannot fuse them because they have different iteration spaces (the first: i in [0, NI), j in [0, NJ); the second: i in [0, NJ), j in [0, NL)). So fusion is not applicable.

3. Loop Skewing:
   The inner loops are reduction loops and do not have carried dependencies that would prevent parallelism. The outer loops are parallel. So skewing is not necessary.

Therefore, we focus on tiling each of the three matrix multiplications.

After tiling, we will apply OpenMP directives to parallelize the outer loops.

Plan for each matrix multiplication:

For E = A * B:
  Original:
    for (i=0; i<NI; i++)
      for (j=0; j<NJ; j++)
        E[i][j] = 0;
        for (k=0; k<NK; k++)
          E[i][j] += A[i][k] * B[k][j];

  Tiled version (with tile sizes TI, TJ, TK? but typically we tile i and j and leave k as the innermost for vectorization and register blocking). We can do:

    for (ii=0; ii<NI; ii+=TI)
      for (jj=0; jj<NJ; jj+=TJ)
        for (i=ii; i<ii+TI; i++)
          for (j=jj; j<jj+TJ; j++)
            E[i][j] = 0;
            for (k=0; k<NK; k++)
              E[i][j] += A[i][k] * B[k][j];

  But note: the inner loop k is not tiled. We can also tile k to improve cache locality for A and B. However, the problem says to prioritize tiling for cache locality. We can do a 3D tiling? But that is more complex. Alternatively, we can do:

    for (ii=0; ii<NI; ii+=TI)
      for (jj=0; jj<NJ; jj+=TJ)
        for (kk=0; kk<NK; kk+=TK)
          for (i=ii; i<ii+TI; i++)
            for (j=jj; j<jj+TJ; j++)
              for (k=kk; k<kk+TK; k++)
                E[i][j] += A[i][k] * B[k][j];

  However, note that the initialization (E[i][j]=0) must be done only once per (i,j). So we have to be careful.

  We can restructure:

    for (ii=0; ii<NI; ii+=TI)
      for (jj=0; jj<NJ; jj+=TJ)
        for (i=ii; i<ii+TI; i++)
          for (j=jj; j<jj+TJ; j++)
            E[i][j] = 0;
        for (kk=0; kk<NK; kk+=TK)
          for (i=ii; i<ii+TI; i++)
            for (j=jj; j<jj+TJ; j++)
              for (k=kk; k<kk+TK; k++)
                E[i][j] += A[i][k] * B[k][j];

  But this is more complex and might not be necessary. Alternatively, we can do:

    for (ii=0; ii<NI; ii+=TI)
      for (jj=0; jj<NJ; jj+=TJ)
        for (kk=0; kk<NK; kk+=TK)
          for (i=ii; i<ii+TI; i++)
            for (k=kk; k<kk+TK; k++)
              for (j=jj; j<jj+TJ; j++)
                E[i][j] += A[i][k] * B[k][j];

  But then we have to initialize E[i][j] to zero. We can do the initialization in a separate loop.

  However, the problem says we can restructure the loops. We can do:

    // Initialize E to zero for the tile
    for (ii=0; ii<NI; ii+=TI)
      for (jj=0; jj<NJ; jj+=TJ)
        for (i=ii; i<ii+TI; i++)
          for (j=jj; j<jj+TJ; j++)
            E[i][j] = 0.0;

    // Then the tiled computation
    for (ii=0; ii<NI; ii+=TI)
      for (jj=0; jj<NJ; jj+=TJ)
        for (kk=0; kk<NK; kk+=TK)
          for (i=ii; i<ii+TI; i++)
            for (k=kk; k<kk+TK; k++)
              for (j=jj; j<jj+TJ; j++)
                E[i][j] += A[i][k] * B[k][j];

  But note: the inner loop order is changed to j innermost. This is good for cache because B[k][j] is accessed in order. Also, A[i][k] is fixed for the inner j loop.

  However, the problem says to preserve the original semantics. This reordering is legal because the inner k loop is independent for each (i,j). But note: the original had k as the innermost. We are changing the order of the loops. This is a form of loop interchange. But it is safe because the reduction is associative.

  But note: the original code initializes E[i][j] to zero and then accumulates. We are doing the same.

  We can do the same for the other two multiplications.

  However, the problem says to choose legal tile sizes. We can set TI, TJ, TK to 64, for example. But we don't know the actual sizes. We can use:

    #define TI 64
    #define TJ 64
    #define TK 64

  But note: the problem says to preserve the function signature and not change code outside the scop. We cannot add #define inside the scop. We can use constant expressions? Or we can use variables? But the problem says we can adjust associated C statements. We can define the tile sizes as constants inside the scop? But that would be adding code. Alternatively, we can use fixed numbers.

  However, the problem says: "Choose legal tile sizes". We can choose 64 as a typical cache block size. But we must ensure that the tile sizes divide the loop bounds? Not necessarily, we can use min in the inner loops. But the problem says "legal", meaning that the transformation is correct. We can use:

    int ii, jj, kk;
    for (ii = 0; ii < _PB_NI; ii += TI) {
      int i_end = (ii + TI) < _PB_NI ? (ii + TI) : _PB_NI;
      for (jj = 0; jj < _PB_NJ; jj += TJ) {
        int j_end = (jj + TJ) < _PB_NJ ? (jj + TJ) : _PB_NJ;
        for (i = ii; i < i_end; i++) {
          for (j = jj; j < j_end; j++) {
            E[i][j] = 0.0;
          }
        }
        for (kk = 0; kk < _PB_NK; kk += TK) {
          int k_end = (kk + TK) < _PB_NK ? (kk + TK) : _PB_NK;
          for (i = ii; i < i_end; i++) {
            for (k = kk; k < k_end; k++) {
              for (j = jj; j < j_end; j++) {
                E[i][j] += A[i][k] * B[k][j];
              }
            }
          }
        }
      }
    }

  This is safe and handles the boundaries.

  But note: the problem says to preserve the original semantics. This does.

  However, the problem also says to apply OpenMP. We can parallelize the outer loops. We can use:

    #pragma omp parallel for collapse(2) private(j,jj,k,kk) // for the initialization? Actually, the initialization is a separate loop. We can combine the initialization and the computation? Or we can do:

    We can do the entire tiled computation in one parallel region? But note: the initialization must be done before the computation. We can do:

    #pragma omp parallel
    {
      // Initialize E for the entire array? But we are tiling. We can do the initialization in parallel as well.

      // We can do:
      #pragma omp for collapse(2) nowait
      for (ii=0; ii<_PB_NI; ii+=TI)
        for (jj=0; jj<_PB_NJ; jj+=TJ)
          for (i=ii; i<min(ii+TI, _PB_NI); i++)
            for (j=jj; j<min(jj+TJ, _PB_NJ); j++)
              E[i][j] = 0.0;

      // Then the computation
      #pragma omp for collapse(3) // we can collapse the ii, jj, kk loops? But note: the kk loop is inside the jj loop. We have three levels: ii, jj, kk. We can collapse 3? But the inner loops are not perfectly nested because of the boundary conditions. We can restructure to have a single nested loop for the computation? Alternatively, we can use:

      #pragma omp for
      for (ii=0; ii<_PB_NI; ii+=TI) {
        int i_end = min(ii+TI, _PB_NI);
        for (jj=0; jj<_PB_NJ; jj+=TJ) {
          int j_end = min(jj+TJ, _PB_NJ);
          for (kk=0; kk<_PB_NK; kk+=TK) {
            int k_end = min(kk+TK, _PB_NK);
            for (i=ii; i<i_end; i++) {
              for (k=kk; k<k_end; k++) {
                for (j=jj; j<j_end; j++) {
                  E[i][j] += A[i][k] * B[k][j];
                }
              }
            }
          }
        }
      }
    }

  But note: the initialization and the computation are in the same parallel region. We use 'nowait' on the initialization to allow the computation to start immediately? But the computation must wait for the initialization of the same tile? Actually, the initialization for a tile (ii,jj) must be done before the computation for that tile. But the computation for a tile (ii,jj) is done in the same thread that did the initialization? Not necessarily. We have a race condition: the initialization for a tile might be done by one thread and the computation by another. We must ensure that the initialization for a tile is done before the computation for that tile.

  We can avoid this by doing the initialization and the computation for the same tile in the same thread. We can combine the two loops:

    #pragma omp parallel for collapse(2) private(kk, k, j)
    for (ii=0; ii<_PB_NI; ii+=TI)
      for (jj=0; jj<_PB_NJ; jj+=TJ)
      {
        int i_end = min(ii+TI, _PB_NI);
        int j_end = min(jj+TJ, _PB_NJ);
        // Initialize the tile
        for (i=ii; i<i_end; i++)
          for (j=jj; j<j_end; j++)
            E[i][j] = 0.0;

        // Compute the tile
        for (kk=0; kk<_PB_NK; kk+=TK)
        {
          int k_end = min(kk+TK, _PB_NK);
          for (i=ii; i<i_end; i++)
            for (k=kk; k<k_end; k++)
              for (j=jj; j<j_end; j++)
                E[i][j] += A[i][k] * B[k][j];
        }
      }

  This way, each thread handles a tile (ii,jj) and does both the initialization and the computation. This is safe.

  We can do the same for the other two multiplications.

  However, note that the second multiplication (F = C * D) has different loop bounds: i from 0 to NJ, j from 0 to NL, and k from 0 to NM.

  And the third multiplication (G = E * F) has i from 0 to NI, j from 0 to NL, and k from 0 to NJ.

  We can apply the same tiling pattern to all three.

  But note: the third multiplication uses E and F, which are computed in the first two. We must not change the order of the three multiplications. We can do:

    First: E = A * B (tiled and parallelized)
    Second: F = C * D (tiled and parallelized)
    Third: G = E * F (tiled and parallelized)

  And we can run the first two in parallel? They are independent. We can use:

    #pragma omp parallel
    {
      #pragma omp single
      {
        // First multiplication: E = A*B
        #pragma omp task
        {
          ... tiled loop for E ...
        }
        // Second multiplication: F = C*D
        #pragma omp task
        {
          ... tiled loop for F ...
        }
      }
      // Then wait for both to finish and then do G = E*F
      #pragma omp single
      {
        #pragma omp task
        {
          ... tiled loop for G ...
        }
      }
    }

  But the problem says to use OpenMP directives and clauses. We are allowed to use tasks? The problem says: "Parallelism: 'parallel for', 'for', 'collapse(n)'." and also "Synchronization: 'nowait' only when safe; 'atomic'/'critical' when needed." and "Thread management". It doesn't explicitly forbid tasks, but the problem says "suitable OpenMP directives". Tasks are suitable for independent tasks.

  However, the problem says: "Do not add dependencies, data races, invalid OpenMP directives, or unsupported transformations." Tasks are valid.

  But note: the problem says to prioritize loop transformations first. We are doing tiling and then parallelizing with tasks.

  Alternatively, we can do:

    #pragma omp parallel
    {
      #pragma omp sections
      {
        #pragma omp section
        { ... E = A*B ... }
        #pragma omp section
        { ... F = C*D ... }
      }
      #pragma omp single
      { ... G = E*F ... }
    }

  But sections are also valid.

  However, the problem says to use 'parallel for' and 'collapse'. We can also do:

    // First, do E and F in parallel with two separate parallel regions? But that would create two teams of threads. We want one team.

    We can do:

      #pragma omp parallel
      {
        #pragma omp single nowait
        {
          // E = A*B
          #pragma omp task
          { ... }
          // F = C*D
          #pragma omp task
          { ... }
        }
        // Implicit barrier at the end of the single region? But we have nowait, so we need a barrier before G.
        #pragma omp barrier
        // Then G = E*F
        #pragma omp for collapse(2) ... // for the tiled loop of G
        ... 
      }

  But the problem says to use 'nowait' only when safe. We are using nowait on the single, which is safe because we have a barrier after.

  However, the problem says to use the directives in the list. We are using 'parallel', 'single', 'task', 'barrier', and then for G we use 'for'. This is acceptable.

  But note: the problem says to use 'collapse(n)' for the for loops. We are using collapse(2) for the outer two loops of the tiled computation for G? Actually, in the tiled version for G, we have:

      for (ii=0; ii<NI; ii+=TI)
        for (jj=0; jj<NL; jj+=TJ)
        { ... }

  We can collapse these two.

  However, the problem says to apply the transformations strictly within the scop. We are going to restructure the entire scop region.

  Let's write the code for the scop region:

    #pragma scop

    // Tile sizes
    const int TI = 64;
    const int TJ = 64;
    const int TK = 64;

    // First two multiplications in parallel
    #pragma omp parallel
    {
      #pragma omp single nowait
      {
        // E = A * B
        #pragma omp task
        {
          for (int ii = 0; ii < _PB_NI; ii += TI)
          {
            int i_end = (ii + TI) < _PB_NI ? (ii + TI) : _PB_NI;
            for (int jj = 0; jj < _PB_NJ; jj += TJ)
            {
              int j_end = (jj + TJ) < _PB_NJ ? (jj + TJ) : _PB_NJ;
              // Initialize E[ii:ii+TI-1][jj:jj+TJ-1] to zero
              for (int i = ii; i < i_end; i++)
                for (int j = jj; j < j_end; j++)
                  E[i][j] = 0.0;
              // Compute the tile
              for (int kk = 0; kk < _PB_NK; kk += TK)
              {
                int k_end = (kk + TK) < _PB_NK ? (kk + TK) : _PB_NK;
                for (int i = ii; i < i_end; i++)
                  for (int k = kk; k < k_end; k++)
                    for (int j = jj; j < j_end; j++)
                      E[i][j] += A[i][k] * B[k][j];
              }
            }
          }
        }

        // F = C * D
        #pragma omp task
        {
          for (int ii = 0; ii < _PB_NJ; ii += TI)
          {
            int i_end = (ii + TI) < _PB_NJ ? (ii + TI) : _PB_NJ;
            for (int jj = 0; jj < _PB_NL; jj += TJ)
            {
              int j_end = (jj + TJ) < _PB_NL ? (jj + TJ) : _PB_NL;
              // Initialize F[ii:ii+TI-1][jj:jj+TJ-1] to zero
              for (int i = ii; i < i_end; i++)
                for (int j = jj; j < j_end; j++)
                  F[i][j] = 0.0;
              // Compute the tile
              for (int kk = 0; kk < _PB_NM; kk += TK)
              {
                int k_end = (kk + TK) < _PB_NM ? (kk + TK) : _PB_NM;
                for (int i = ii; i < i_end; i++)
                  for (int k = kk; k < k_end; k++)
                    for (int j = jj; j < j_end; j++)
                      F[i][j] += C[i][k] * D[k][j];
              }
            }
          }
        }
      } // end single

      // Barrier to ensure E and F are computed
      #pragma omp barrier

      // Now compute G = E * F
      #pragma omp for collapse(2)
      for (int ii = 0; ii < _PB_NI; ii += TI)
      {
        int i_end = (ii + TI) < _PB_NI ? (ii + TI) : _PB_NI;
        for (int jj = 0; jj < _PB_NL; jj += TJ)
        {
          int j_end = (jj + TJ) < _PB_NL ? (jj + TJ) : _PB_NL;
          // Initialize G[ii:ii+TI-1][jj:jj+TJ-1] to zero
          for (int i = ii; i < i_end; i++)
            for (int j = jj; j < j_end; j++)
              G[i][j] = 0.0;
          // Compute the tile
          for (int kk = 0; kk < _PB_NJ; kk += TK)
          {
            int k_end = (kk + TK) < _PB_NJ ? (kk + TK) : _PB_NJ;
            for (int i = ii; i < i_end; i++)
              for (int k = kk; k < k_end; k++)
                for (int j = jj; j < j_end; j++)
                  G[i][j] += E[i][k] * F[k][j];
          }
        }
      }
    } // end parallel

    #pragma endscop

  But note: the problem says to use 'default(none)'? We are not specifying it. We should specify the data sharing attributes. We can add:

    #pragma omp parallel default(none) shared(E, A, B, F, C, D, G, _PB_NI, _PB_NJ, _PB_NK, _PB_NL, _PB_NM) private(ii, jj, kk, i, j, k, i_end, j_end, k_end, TI, TJ, TK)

  But note: TI, TJ, TK are constants defined inside the scop. We can define them outside the parallel? Actually, we defined them at the beginning of the scop. We can move them inside the parallel? Or we can define them as constants outside the parallel? We can do:

    const int TI = 64;
    const int TJ = 64;
    const int TK = 64;

    #pragma omp parallel default(none) shared(E, A, B, F, C, D, G, _PB_NI, _PB_NJ, _PB_NK, _PB_NL, _PB_NM) private(ii, jj, kk, i, j, k, i_end, j_end, k_end)

  But note: the problem says to preserve the original semantics. We are using the same variable names. We are also using the same array accesses.

  However, the problem says: "Do not change observable behavior." This transformation is safe.

  But note: the original code used the same variable names (i, j, k) for all three multiplications. We are reusing them in the tasks and the parallel for. This is safe because they are private.

  Also, note: the original code had:

      for (i = 0; i < _PB_NI; i++)
        for (j = 0; j < _PB_NJ; j++)
          E[i][j] = 0.0;
          for (k=0; k<_PB_NK; k++)
            E[i][j] += ...;

  We are doing the same, but in a tiled manner.

  We must be cautious: the tile sizes (TI, TJ, TK) must be chosen appropriately. We set them to 64. This is a typical value for L1 cache. But note: the problem says "choose legal tile sizes". We assume 64 is legal.

  Also, note: the problem says to use 'simd' if beneficial. We can add simd to the inner loops? For example, the inner j loop in the computation:

      for (j = jj; j < j_end; j++)
        E[i][j] += A[i][k] * B[k][j];

  We can add:

      #pragma omp simd
      for (j = jj; j < j_end; j++)
        E[i][j] += A[i][k] * B[k][j];

  This might help. We can do the same for the other inner loops.

  Also, we can use 'aligned' if we know the arrays are aligned? But the problem doesn't say. We can skip.

  We can also use 'schedule' for the parallel for of G? The default is static, which is fine.

  We can also use 'num_threads'? But we don't know the number of threads. We can leave it to the runtime.

  We can also use 'proc_bind'? We can set it to 'spread' for better performance? But the problem doesn't specify. We can skip.

  We can also use 'reduction'? But we are doing the reduction manually by initializing to zero and then accumulating. We don't need a reduction clause.

  Therefore, we can add simd to the inner j loops.

  Let's update the inner j loops in all three multiplications to use #pragma omp simd.

  But note: the inner j loop is the one that we are vectorizing. We can do:

      for (int k = kk; k < k_end; k++)
      {
        #pragma omp simd
        for (int j = jj; j < j_end; j++)
          E[i][j] += A[i][k] * B[k][j];
      }

  This is safe because the j loop is independent.

  We do the same for F and G.

  However, note: the problem says to use 'simdlen(n)' or 'safelen(n)'? We can use the default.

  Also, note: the problem says to use 'aligned' if we know the alignment. We don't, so we skip.

  We can also use 'private' for the inner loop variables? They are already private because they are declared inside the parallel region.

  We must be cautious: the variables ii, jj, kk, i, j, k, i_end, j_end, k_end are all private.

  We are using 'default(none)' and then listing all shared variables. We must list all the arrays and the constants (_PB_NI, etc.). We are doing that.

  But note: the constants _PB_NI, _PB_NJ, etc. are defined in the polybench header. We assume they are shared.

  Also, note: the arrays E, A, B, F, C, D, G are shared.

  We are not using 'firstprivate' or 'lastprivate' because we don't need them.

  We are not using 'critical' or 'atomic' because there are no races.

  We are using 'nowait' on the single, which is safe because we have a barrier after.

  We are using 'barrier' to synchronize before computing G.

  We are using 'collapse(2)' for the parallel for of G? Actually, we are using:

      #pragma omp for collapse(2)

  This collapses the two outer loops (ii and jj) of the tiled loop for G.

  We can do the same for the tasks? But the tasks are not using a for loop with collapse. They are using nested for loops. We can collapse the ii and jj loops in the tasks? We can do:

      #pragma omp task
      {
        #pragma omp for collapse(2) // but note: we are inside a task, and we are already in a parallel region. We cannot use omp for inside a task without a parallel region? Actually, we can use:

        // We can do:
        for (int ii = 0; ii < _PB_NI; ii += TI)
          for (int jj = 0; jj < _PB_NJ; jj += TJ)
          {
            ... // the rest
          }

        // But we cannot collapse with omp for inside a task? We can use:

        // Alternatively, we can use a nested loop without collapse. It's acceptable.

        // We leave it as is.
      }

  We are not collapsing the loops in the tasks. It's acceptable.

  We are also not using 'schedule' for the tasks? The tasks are scheduled by the runtime.

  We are not using 'num_threads' or 'proc_bind'. We can add:

      #pragma omp parallel proc_bind(spread) num_threads(OMP_NUM_THREADS) ...

  But we don't know OMP_NUM_THREADS. We can leave it to the environment.

  We can set:

      #pragma omp parallel default(none) shared(...) proc_bind(spread)

  This might improve performance by binding threads to cores.

  We'll add 'proc_bind(spread)'.

  Therefore, the final code for the scop region is as below.

  Note: We must preserve the '#pragma scop' and '#pragma endscop' markers.

  We are replacing the entire scop region with the above code.

  Let's write the code.

  Important: We must not change the function signature or the code outside the scop.

  We are only changing the code between #pragma scop and #pragma endscop.

  We are using OpenMP, so we must include <omp.h>? But the original code doesn't. The problem says to preserve the original code outside the scop. We are not allowed to add includes. We assume that the OpenMP flag will be used when compiling.

  We are using min() function. We can use a macro? Or we can use ternary operator. We are using:

      int i_end = (ii + TI) < _PB_NI ? (ii + TI) : _PB_NI;

  This is standard C.

  We are using const int for TI, TJ, TK. We can define them at the beginning of the scop.

  We are using the same variable names (i, j, k) in the inner loops. This is safe because they are in different scopes.

  We are using the same array accesses.

  We are initializing the tiles to zero and then accumulating, same as the original.

  We are using the same order of multiplications.

  We are using OpenMP tasks for the first two multiplications and a parallel for for the third.

  We are using a barrier to synchronize.

  We are using simd for the inner j loops.

  We are using proc_bind(spread) for the parallel region.

  We are using default(none) and listing all shared variables.

  We are using private for the loop variables.

  We are not using any unsupported transformations.

  We are not changing the observable behavior.

  We are not adding data races.

  We are not adding dependencies.

  We are not adding invalid OpenMP directives.

  We are not modifying code outside the scop.

  We are preserving the markers.

  We are returning only the complete optimized C code.

  We are not adding explanations.

  We are not using Markdown.

  We are not changing the surrounding text.

  We are only returning the code.

  Let's write the code.
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
  const int TI = 64;
  const int TJ = 64;
  const int TK = 64;
  
  #pragma omp parallel default(none) proc_bind(spread) \
    shared(E, A, B, F, C, D, G, _PB_NI, _PB_NJ, _PB_NK, _PB_NL, _PB_NM) \
    private(i, j, k, ii, jj, kk, i_end, j_end, k_end)
  {
    #pragma omp single nowait
    {
      #pragma omp task
      {
        for (ii = 0; ii < _PB_NI; ii += TI)
        {
          i_end = (ii + TI) < _PB_NI ? (ii + TI) : _PB_NI;
          for (jj = 0; jj < _PB_NJ; jj += TJ)
          {
            j_end = (jj + TJ) < _PB_NJ ? (jj + TJ) : _PB_NJ;
            for (i = ii; i < i_end; i++)
              for (j = jj; j < j_end; j++)
                E[i][j] = 0.0;
            for (kk = 0; kk < _PB_NK; kk += TK)
            {
              k_end = (kk + TK) < _PB_NK ? (kk + TK) : _PB_NK;
              for (i = ii; i < i_end; i++)
                for (k = kk; k < k_end; k++)
                  #pragma omp simd
                  for (j = jj; j < j_end; j++)
                    E[i][j] += A[i][k] * B[k][j];
            }
          }
        }
      }
      #pragma omp task
      {
        for (ii = 0; ii < _PB_NJ; ii += TI)
        {
          i_end = (ii + TI) < _PB_NJ ? (ii + TI) : _PB_NJ;
          for (jj = 0; jj < _PB_NL; jj += TJ)
          {
            j_end = (jj + TJ) < _PB_NL ? (jj + TJ) : _PB_NL;
            for (i = ii; i < i_end; i++)
              for (j = jj; j < j_end; j++)
                F[i][j] = 0.0;
            for (kk = 0; kk < _PB_NM; kk += TK)
            {
              k_end = (kk + TK) < _PB_NM ? (kk + TK) : _PB_NM;
              for (i = ii; i < i_end; i++)
                for (k = kk; k < k_end; k++)
                  #pragma omp simd
                  for (j = jj; j < j_end; j++)
                    F[i][j] += C[i][k] * D[k][j];
            }
          }
        }
      }
    }
    #pragma omp barrier
    
    #pragma omp for collapse(2)
    for (ii = 0; ii < _PB_NI; ii += TI)
    {
      i_end = (ii + TI) < _PB_NI ? (ii + TI) : _PB_NI;
      for (jj = 0; jj < _PB_NL; jj += TJ)
      {
        j_end = (jj + TJ) < _PB_NL ? (jj + TJ) : _PB_NL;
        for (i = ii; i < i_end; i++)
          for (j = jj; j < j_end; j++)
            G[i][j] = 0.0;
        for (kk = 0; kk < _PB_NJ; kk += TK)
        {
          k_end = (kk + TK) < _PB_NJ ? (kk + TK) : _PB_NJ;
          for (i = ii; i < i_end; i++)
            for (k = kk; k < k_end; k++)
              #pragma omp simd
              for (j = jj; j < j_end; j++)
                G[i][j] += E[i][k] * F[k][j];
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
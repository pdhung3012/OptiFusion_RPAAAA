We are given a code for the symm kernel. The scop region is the main computational kernel.
 The original code has three nested loops: i, j, k.
 The loop structure is:
   for (i=0; i<_PB_M; i++)
     for (j=0; j<_PB_N; j++)
        temp2 = 0;
        for (k=0; k<i; k++) {
           C[k][j] += alpha * B[i][j] * A[i][k];
           temp2 += B[k][j] * A[i][k];
        }
        C[i][j] = beta * C[i][j] + alpha * B[i][j] * A[i][i] + alpha * temp2;

 We are to optimize only within the scop region.

 Step 1: Loop Transformations (highest priority: tiling, fusion, skewing)

 1. Loop Tiling: 
    - The innermost loop (k) is over a triangular part (k from 0 to i-1). 
    - The outer loops (i and j) are rectangular. 
    - We note that the inner loop has a dependence on C[k][j] which is written and then read in the same j iteration but for different k and i.
    - Also, the inner loop accumulates into temp2 and then uses it in the update of C[i][j].
    - The access pattern for C in the inner loop: C[k][j] is updated for each k and for a fixed j and i. 
      This means that for a fixed j, we are updating a column of C (from row 0 to row i-1) and then updating C[i][j].
    - The access to A: A[i][k] is read for k from 0 to i-1. This is a row of A (but only the lower triangular part because k<=i-1 and i is the current row).
    - The access to B: 
          B[i][j] is read multiple times (once per k) and also in the final update.
          B[k][j] is read for each k (so a column of B).

    - We are updating C[k][j] for k in [0, i-1] and then C[i][j]. This is a column-wise update of C.

    - Tiling might be beneficial for the two outer loops (i and j) to improve locality. However, note that the inner loop (k) is triangular and depends on i.

    - We can consider tiling the i and j loops. But note that the inner loop (k) goes from 0 to i-1, so the tile in i must be handled carefully.

    - Alternatively, we can try to tile the j loop and the k loop? But k is triangular and depends on i.

    - Another idea: interchange the i and j loops? But note that the update of C[k][j] in the inner loop for a fixed i and j: 
          C[k][j] is updated for k from 0 to i-1. 
      Then the same C[k][j] might be updated again in a different i (if k < i' for a larger i'). 
      So the order of i matters: we are updating C[k][j] for increasing i (and k from 0 to i-1). 
      This means that for a fixed j, we update the entire column of C from row0 to row (m-1) as i increases.

    - We can try to tile the i loop and the j loop. We'll break the i loop into blocks of size TI and the j loop into blocks of size TJ.

    - However, note that the inner loop (k) goes from 0 to i-1. In a tile for i in [i0, i0+TI-1] and j in [j0, j0+TJ-1], 
      the inner loop for k would go from 0 to i-1, which is not aligned with the tile. 

    - We can restructure the loop to have a rectangular tile in i and j, but then the inner loop (k) would have to be adjusted.

    - Alternatively, we can try to tile the j loop and the k loop? But k is bounded by i, which is the outer loop.

    - Given the complexity, we might skip tiling for now and look at fusion.

 2. Loop Fusion:
    - The scop region has only one loop nest. There are no adjacent loops to fuse.

 3. Loop Skewing:
    - The inner loop (k) has a carried dependence? 
        The update: C[k][j] += ... 
        This is a write to C[k][j] and then in the same j iteration, for a fixed i, we update multiple k. 
        But note: for a fixed j and k, the same C[k][j] is updated for every i such that i > k. 
        So the order of i matters: we are accumulating into C[k][j] for increasing i.

        The dependence: 
          The write to C[k][j] in iteration (i, j, k) and then the same C[k][j] is read in the same j iteration for a larger i? 
          Actually, no: the update for C[k][j] is only in the inner loop for a fixed i and j. 
          But note: the same C[k][j] is updated again in a different i (if i' > k) and same j. 
          So the order of i is important: we must process i in increasing order.

        Also, the inner loop for k: we are updating C[k][j] for k from 0 to i-1. 
        Then in the same j iteration, we update C[i][j] at the end.

        The update of C[i][j] uses the old value of C[i][j] (from the previous state) and then overwrites it. 
        And note that for a fixed j, the update of C[i][j] happens only once (for the current i).

        The inner loop (k) does not have a carried dependence that would prevent parallelization? 
        Actually, the inner loop (k) is a reduction on temp2 and also updates C[k][j] for different k. 
        But note: the updates to C[k][j] are to different k, so they are independent. 
        However, the same C[k][j] is updated for different i (if k < i). 

        The loop-carried dependence is on the outer i loop: 
          For a fixed k and j, the value of C[k][j] is updated for every i such that i > k. 
          So the i loop is not parallel.

        We cannot parallelize the i loop because of the update to C[k][j] (which is a location that is updated in multiple i iterations).

        However, note that the update for C[k][j] is: 
            C[k][j] += alpha * B[i][j] * A[i][k]
        This is an accumulation. So the order of i does not matter? 
        Actually, it does: because floating point addition is not associative. 
        But if we are allowed to change the order (and the problem does not specify strict floating point associativity) then we might reorder.

        But the problem says: "Do not change observable behavior". So we must preserve the order.

        Therefore, we cannot parallelize the i loop.

        What about the j loop? 
          For a fixed i, the j loop iterates over columns. 
          The update for C[k][j] in the inner loop: 
             C[k][j] is updated for a fixed k and j. 
          And the update for C[i][j] is for a fixed i and j. 
          There is no dependence between different j. 

          So the j loop is parallel.

        Therefore, we can parallelize the j loop.

        But note: the inner loop (k) updates C[k][j] for k from 0 to i-1. 
        For a fixed i, the j loop is parallel. 

        However, the update of C[k][j] for a fixed k and j: 
          It is updated for every i such that i>k. 
          But in the j loop, we are iterating j independently. 
          And for a fixed i, the j loop is independent.

        So we can parallelize the j loop.

        But note: the inner loop (k) is a reduction on temp2? 
          temp2 is a scalar that is private to the (i,j) iteration. 
          So we can parallelize j without any reduction on temp2.

        However, the update to C[k][j] is a shared array. 
          We are updating C[k][j] for k in [0, i-1] and for a fixed j. 
          And we are doing this for multiple j in parallel. 
          But note: for a fixed k, the j loop updates different j. 
          So there is no conflict between j.

        Therefore, we can parallelize the j loop.

        But wait: the update of C[k][j] is done for a fixed j and k, and we are updating the same C[k][j] for multiple i? 
          However, in the j loop, we are iterating j independently. 
          And for a fixed j, the update of C[k][j] happens for multiple i (but that is in the i loop, which is sequential). 
          In the j loop, we are processing one j at a time? No, we are going to parallelize j so multiple j are processed concurrently.

          For a fixed k, the j loop updates C[k][j] for different j. 
          So no two threads will update the same C[k][j] because j is different.

        Therefore, the j loop is parallel.

        We can add: 
            #pragma omp parallel for 
            for (j=0; j<_PB_N; j++)

        But note: the inner loop (k) updates C[k][j] for k from 0 to i-1. 
        And the update of C[i][j] at the end. 
        Both are safe for parallel j.

        However, we must be cautious: the array C is shared. 
        We are updating different j, so no two threads write to the same element.

        So we can parallelize the j loop.

        But note: the original code has the i loop outer and j inner. 
        We are going to parallelize the j loop. 

        However, the j loop is the inner loop. 
        We can interchange the i and j loops? 
          The update for C[k][j] in the inner loop: 
             C[k][j] is updated for k from 0 to i-1. 
          Then the update for C[i][j] at the end.

          The update of C[k][j] for a fixed k and j: 
             It is updated for every i such that i>k. 
          So if we interchange i and j, then for a fixed j, we would update C[k][j] for k from 0 to i-1 and for increasing i. 
          But note: the update of C[k][j] for a fixed k and j must be done for all i>k. 
          And the update of C[i][j] for a fixed i and j must be done after the inner loop for that i and j.

          The interchange of i and j would change the order of updates to C[k][j] for a fixed k and j? 
          Currently, for a fixed k and j, the updates happen for i from k+1 to m-1. 
          If we interchange, then for a fixed j, we would do:
             for i from 0 to m-1:
                 for k from 0 to i-1: update C[k][j]   -> but note: k < i, so for a fixed k, we update when i>k.

          So the set of updates to C[k][j] is the same. 
          However, the order of updates for a fixed k and j: 
             Currently: i = k+1, k+2, ... , m-1.
             After interchange: same order? 
                 For a fixed j, we iterate i from 0 to m-1. 
                 For a fixed k, the update happens when i>k. 
                 So the updates for a fixed k and j occur in increasing i. 
                 So the order is the same.

          Therefore, we can interchange i and j.

        Why interchange? 
          Because we want to parallelize j. 
          If we leave j inner, then we can parallelize the inner loop? 
          But the inner loop is j, and we can do:
             #pragma omp parallel for 
             for (j=0; j<_PB_N; j++)
          But note: the inner loop (k) is inside j. 
          We can do:
             for (i=0; i<_PB_M; i++) {
                 #pragma omp parallel for 
                 for (j=0; j<_PB_N; j++) {
                     ... inner k loop ...
                 }
             }

          This is acceptable.

        However, we can also interchange to have j outer and then parallelize j. 
          Then the code becomes:
             for (j=0; j<_PB_N; j++) {
                 for (i=0; i<_PB_M; i++) {
                     ... inner k loop ...
                 }
             }

          And then we can parallelize the j loop.

        Which is better? 
          The original order: i outer, j inner. 
          The j loop is of size n, which might be large. 
          The i loop is of size m.

          The inner loop (k) is triangular: the total work is O(m^2 * n). 

          If we parallelize the inner j loop (without interchange) then we have a parallel region for each i. 
          That might be expensive if m is large.

          Alternatively, if we interchange and then parallelize the outer j loop, we have one parallel region for the entire j loop. 
          This is more efficient.

        So we decide to interchange the i and j loops.

        Steps for interchange:
          Original: 
            for i: 0..m-1
              for j: 0..n-1
          After interchange:
            for j: 0..n-1
              for i: 0..m-1

          And then we can parallelize the j loop.

        But note: the inner loop (k) is from 0 to i-1. 
          This is independent of j, so it's safe.

        Also, the update of C[i][j] at the end: 
          This is now for a fixed j and i. 
          And we are updating C[i][j] for each i in order. 
          But note: the update of C[i][j] uses the old value of C[i][j] (which is the value from the previous state, before the i loop started) and then overwrites it. 
          And the inner loop (k) updates C[k][j] for k from 0 to i-1. 
          For a fixed j, the update of C[k][j] for k<i happens in the same j iteration and for increasing i. 
          And the update of C[i][j] happens at the end of the i iteration for that j.

          The update of C[i][j] does not depend on the inner loop of the same i? 
          Actually, the inner loop updates C[k][j] for k<i, and then we update C[i][j]. 
          And the update of C[i][j] uses the old value of C[i][j] (which is the value from the initial state, or from a previous j? but j is fixed) and then overwrites it. 
          So it's safe.

        Therefore, we can interchange.

        After interchange, the code becomes:

          for (j = 0; j < _PB_N; j++)
            for (i = 0; i < _PB_M; i++) {
              temp2 = 0;
              for (k = 0; k < i; k++) {
                 C[k][j] += alpha * B[i][j] * A[i][k];
                 temp2 += B[k][j] * A[i][k];
              }
              C[i][j] = beta * C[i][j] + alpha * B[i][j] * A[i][i] + alpha * temp2;
            }

        Then we can parallelize the j loop.

        But note: the inner loop (k) updates C[k][j] for k from 0 to i-1. 
        For a fixed j, we are updating the same column of C (for rows 0 to i-1) and then row i. 
        And we are doing this for increasing i. 
        So the entire column j is updated in order of i. 

        And the j loop is parallel: each thread handles one j.

        This is safe.

 Step 2: OpenMP directives.

  We will add:
      #pragma omp parallel for 
      for (j=0; j<_PB_N; j++)

  But note: the inner loop (k) is a reduction on temp2? 
      temp2 is a scalar that is private to the (i,j) iteration. 
      So we don't need a reduction clause for temp2.

  However, we must consider the data sharing of the arrays.

  The arrays A, B, C are shared. 
  The variables i, k, temp2 are private by default in the parallel region? 
      But note: the j loop is the only loop in the parallel region. 
      The i and k loops are inside. 
      We must make sure that i, k, temp2 are private.

  We can use:
      #pragma omp parallel for private(i, k, temp2)

  But note: the variable alpha and beta are shared? 
      They are read-only, so shared is fine.

  Also, the array indices: we are using the original _PB_M and _PB_N.

  We can also consider scheduling. 
      The j loop: the work per j is O(m^2) (because the inner k loop goes up to i, and i goes up to m-1). 
      So the work is not constant per j. 
      We can use dynamic scheduling? 
      But note: the work per j is the same for every j? 
          For a fixed j, the work is: 
             for i from 0 to m-1: 
                 the inner k loop: k from 0 to i-1 -> about i/2 iterations on average? 
                 total for j: sum_{i=0}^{m-1} i = m*(m-1)/2.
          So the work per j is the same. 
          Therefore, static scheduling is fine.

  We can use: schedule(static)

  Also, we can consider SIMD for the inner k loop? 
      The inner k loop: 
          C[k][j] += alpha * B[i][j] * A[i][k];
          temp2 += B[k][j] * A[i][k];
      This is a short loop? The length is i, which varies. 
      We can try to vectorize the k loop.

      We can add: 
          #pragma omp simd
          for (k=0; k<i; k++)

      But note: the loop bound is i, which is not constant. 
      Also, the loop is not always vectorizable? 
          The update of C[k][j] is a scattered write? 
          Actually, k is the inner loop index, so the writes to C[k][j] are contiguous in k? 
          Because k is the row index and j is fixed. 
          So the array C is stored in row-major order? 
          The original code uses: 
              POLYBENCH_2D(C,M,N,m,n)
          The macro POLYBENCH_2D is defined in polybench.h. 
          Typically, it is row-major: C[i][j] is at offset i*N+j.

          Therefore, the access to C[k][j] for fixed j and varying k is contiguous? 
          Actually, no: because k is the row index, and j is fixed. 
          The array is stored by rows: 
             row0: C[0][0], C[0][1], ... 
             row1: C[1][0], C[1][1], ...
          So for fixed j, the elements C[k][j] for k=0,1,... are separated by N elements. 
          So it is a column access -> not contiguous.

      Therefore, the k loop is not contiguous in memory for C[k][j]. 
      But note: the inner loop also accesses A[i][k] and B[k][j]. 
          A[i][k]: for fixed i and varying k -> contiguous (row access).
          B[k][j]: for fixed j and varying k -> column access (not contiguous).

      So the k loop has two non-contiguous accesses (C and B) and one contiguous (A). 
      Vectorization might be limited.

      However, we can still try to use SIMD. 
      We can add: 
          #pragma omp simd
          for (k=0; k<i; k++)

      But note: the loop bound is i, which is a variable. 
      The compiler might not vectorize if the bound is not constant? 
      But OpenMP simd allows the bound to be a variable.

      We can try.

  Also, we can consider the entire j loop as a parallel for and then the inner k loop as simd.

  We can do:
      #pragma omp parallel for private(i, k, temp2)
      for (j=0; j<_PB_N; j++) {
          for (i=0; i<_PB_M; i++) {
              temp2 = 0;
              #pragma omp simd
              for (k=0; k<i; k++) {
                 C[k][j] += alpha * B[i][j] * A[i][k];
                 temp2 += B[k][j] * A[i][k];
              }
              C[i][j] = beta * C[i][j] + alpha * B[i][j] * A[i][i] + alpha * temp2;
          }
      }

  But note: the update of C[k][j] in the inner loop: 
        This is a reduction? 
        Actually, it is an accumulation: C[k][j] += ... 
        And the same C[k][j] is updated for multiple i (if k < i). 
        But in the inner k loop for a fixed i and j, we are updating C[k][j] for k from 0 to i-1. 
        And for a fixed j and k, the update happens only once per i (and i>k). 
        So in the inner k loop, each k is updated exactly once. 
        Therefore, the inner k loop does not have a reduction on C[k][j]? 
        It is just a sequence of independent updates to different k.

        However, the update of temp2 is a reduction. 
        But the inner k loop is a reduction on temp2? 
        Yes, and we are using a scalar temp2 that is private to the (i,j) iteration. 
        So the reduction on temp2 is local to the inner loop and does not require a reduction clause.

        But note: the inner k loop is being vectorized. 
        The reduction on temp2 can be handled by the vectorizer? 
        Yes, because it is a simple sum.

        However, the update of C[k][j] is a memory operation that is not a reduction in the vector sense? 
        It is a scattered write? 
        But note: the writes to C[k][j] are to different k, so they are independent. 
        The vectorizer can do them in any order? 
        But the order of floating point operations might change? 
        We are not allowed to change the order? 
        The problem says: "Do not change observable behavior". 
        So we must be cautious.

        However, the original order of the k loop is increasing k. 
        The vectorizer might change the order? 
        But the vectorizer typically preserves the order for the same iteration? 
        Actually, the vectorizer will do the operations in chunks, but the chunks are in increasing order. 
        So the order of the updates to C[k][j] for increasing k is preserved. 
        And the reduction on temp2 is also done in order? 
        The vectorizer will do the reduction in order? 
        Actually, the reduction on temp2 is a sum and the vectorizer will do it in a tree, which is not the same order. 
        But the problem does not specify strict floating point associativity? 
        We are allowed to change the order? 

        The problem says: "Do not change observable behavior". 
        So we must not change the result. 
        Therefore, we cannot use SIMD for the inner loop if it changes the order of the reduction? 

        However, the reduction on temp2 is a sum of positive terms? 
        But the problem does not specify. 

        We have two options:
          a) Do not use SIMD for the inner loop.
          b) Use SIMD and hope that the compiler does the reduction in a way that is acceptable? 
             But note: the OpenMP standard does not require strict associativity for reductions? 
             However, the problem says "Do not change observable behavior".

        We must be safe. 
        The original code does the reduction in increasing k order. 
        The vectorized version might do it in a different order. 
        So we cannot use SIMD for the inner loop.

        Therefore, we skip the SIMD for the inner loop.

  Alternatively, we can try to use SIMD only for the part that is safe? 
      The update of C[k][j] is independent and can be done in any order? 
      But note: the update of C[k][j] is an accumulation: 
          C[k][j] = C[k][j] + alpha * B[i][j] * A[i][k]
      This is a read and then write. 
      The order of these updates for different k does not matter? 
      Because they are independent. 
      So we can do them in any order. 
      But the reduction on temp2: 
          temp2 = temp2 + B[k][j] * A[i][k]
      This is a reduction that must be done in the original order? 
      The problem does not specify, but to be safe we avoid.

  We decide to not use SIMD for the inner loop.

  Therefore, we only parallelize the j loop.

  We can also consider the collapse clause? 
      The j loop is the only loop we are parallelizing. 
      We cannot collapse the i loop because the i loop is not parallel (due to the dependence on C[k][j] for increasing i).

  So we do:

      #pragma omp parallel for private(i, k, temp2)
      for (j=0; j<_PB_N; j++) {
          for (i=0; i<_PB_M; i++) {
              temp2 = 0;
              for (k=0; k<i; k++) {
                 C[k][j] += alpha * B[i][j] * A[i][k];
                 temp2 += B[k][j] * A[i][k];
              }
              C[i][j] = beta * C[i][j] + alpha * B[i][j] * A[i][i] + alpha * temp2;
          }
      }

  We must also consider data sharing: 
      The arrays A, B, C are shared. 
      The variables i, k, temp2 are private to each thread (each j iteration). 
      We have declared them as private in the pragma.

  Also, we can add default(none) to be safe? 
      But the problem does not require it. 
      We can do:

          #pragma omp parallel for default(none) private(i, k, temp2) shared(A, B, C, alpha, beta, _PB_M, _PB_N)

      But note: the original code uses _PB_M and _PB_N. 
      We must share them.

  However, the problem says: "Do not add dependencies, data races, invalid OpenMP directives, or unsupported transformations."

  We are not adding any data race because:
      - The update of C[k][j] for a fixed j and k: only one thread (the one handling j) updates it. 
        And for a fixed j and k, the update happens for multiple i (if k < i) but that is in the sequential i loop. 
        So no race.

      - The update of C[i][j] for a fixed i and j: only one thread (the one handling j) updates it. 
        And it is updated only once (at the end of the i iteration). 
        So no race.

  Therefore, the code is safe.

  But note: the original code had the i loop outer and j inner. 
  We interchanged to j outer and i inner. 
  This is a valid transformation.

  We output the code within the scop region.

  Note: the scop region markers must be preserved.

  The original scop region:

      #pragma scop
         for (i = 0; i < _PB_M; i++)
            for (j = 0; j < _PB_N; j++ )
            {
              temp2 = 0;
              for (k = 0; k < i; k++) {
                 C[k][j] += alpha*B[i][j] * A[i][k];
                 temp2 += B[k][j] * A[i][k];
              }
              C[i][j] = beta * C[i][j] + alpha*B[i][j] * A[i][i] + alpha * temp2;
           }
      #pragma endscop

  We replace it with:

      #pragma scop
         #pragma omp parallel for default(none) private(i, k, temp2) shared(A, B, C, alpha, beta, _PB_M, _PB_N)
         for (j = 0; j < _PB_N; j++) {
            for (i = 0; i < _PB_M; i++) {
               temp2 = 0;
               for (k = 0; k < i; k++) {
                  C[k][j] += alpha * B[i][j] * A[i][k];
                  temp2 += B[k][j] * A[i][k];
               }
               C[i][j] = beta * C[i][j] + alpha * B[i][j] * A[i][i] + alpha * temp2;
            }
         }
      #pragma endscop

  We added the OpenMP pragma and changed the loop order.

  But note: the original code used _PB_M and _PB_N. 
  We are using the same.

  Also, we added default(none) and explicitly listed the shared variables.

  We did not change the variable names.

  We must be cautious: the original code had the loops in the order i then j. 
  We changed to j then i. 
  This is a valid transformation as we argued.

  We did not do tiling because the inner loop is triangular and the interchange already gives us a good parallelization.

  We did not do fusion because there is only one loop.

  We did not do skewing because the inner loop is not parallel and we found a way to parallelize the j loop.

  Therefore, we output the code as below.

  Note: the entire function kernel_symm is:

      static
      void kernel_symm(int m, int n,
              DATA_TYPE alpha,
              DATA_TYPE beta,
              DATA_TYPE POLYBENCH_2D(C,M,N,m,n),
              DATA_TYPE POLYBENCH_2D(A,M,M,m,m),
              DATA_TYPE POLYBENCH_2D(B,M,N,m,n))
      {
        int i, j, k;
        DATA_TYPE temp2;

        #pragma scop
        #pragma omp parallel for default(none) private(i, k, temp2) shared(A, B, C, alpha, beta, _PB_M, _PB_N)
        for (j = 0; j < _PB_N; j++) {
            for (i = 0; i < _PB_M; i++) {
                temp2 = 0;
                for (k = 0; k < i; k++) {
                    C[k][j] += alpha * B[i][j] * A[i][k];
                    temp2 += B[k][j] * A[i][k];
                }
                C[i][j] = beta * C[i][j] + alpha * B[i][j] * A[i][i] + alpha * temp2;
            }
        }
        #pragma endscop

      }

  We removed the original loop variables i, j, k from the function body? 
      Actually, we are reusing them. 
      We declared them at the beginning: int i, j, k; 
      And then we use them in the loops. 
      This is safe.

  We also note that the original code had a variable temp2 declared. 
      We are using it inside the loop. 
      We made it private to the parallel region.

  We are done.

  However, note: the original code used _PB_M and _PB_N. 
      The function parameters are m and n. 
      But the scop region uses _PB_M and _PB_N. 
      This is because the PolyBench benchmarks use these macros to allow problem size changes. 
      We must use the same.

  We output the entire code? 
      The problem says: "Return ONLY the complete optimized C code."

  But note: the problem says: "Preserve the '#pragma scope' and '#pragma endscop' markers exactly as supplied."

  We have preserved them.

  We output the entire file? 
      The problem says: "Return ONLY the complete optimized C code."

  We must return the entire file? 
      But the problem says: "Optimize only code in regions delimited by '#pragma scope' and '#pragma endscop'."

  So we only change the scop region. 
  The rest of the file remains the same.

  We output the entire file with the scop region replaced.

  But note: the problem says: "Return ONLY the complete optimized C code. No explanations, Markdown code fences, or surrounding text."

  So we output the entire C code with the scop region optimized.

  We write the entire file.

  However, the problem says: "Do not modify code outside the marked regions."

  So we only change the scop region.

  We output the entire file with the scop region replaced by our optimized version.

  We do not change anything else.

  We write the code below.
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "symm.h"


/* Array initialization. */
static
void init_array(int m, int n,
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(C,M,N,m,n),
		DATA_TYPE POLYBENCH_2D(A,M,M,m,m),
		DATA_TYPE POLYBENCH_2D(B,M,N,m,n))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      C[i][j] = (DATA_TYPE) ((i+j) % 100) / m;
      B[i][j] = (DATA_TYPE) ((n+i-j) % 100) / m;
    }
  for (i = 0; i < m; i++) {
    for (j = 0; j <=i; j++)
      A[i][j] = (DATA_TYPE) ((i+j) % 100) / m;
    for (j = i+1; j < m; j++)
      A[i][j] = -999; //regions of arrays that should not be used
  }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int m, int n,
		 DATA_TYPE POLYBENCH_2D(C,M,N,m,n))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("C");
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
	if ((i * m + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, C[i][j]);
    }
  POLYBENCH_DUMP_END("C");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_symm(int m, int n,
		 DATA_TYPE alpha,
		 DATA_TYPE beta,
		 DATA_TYPE POLYBENCH_2D(C,M,N,m,n),
		 DATA_TYPE POLYBENCH_2D(A,M,M,m,m),
		 DATA_TYPE POLYBENCH_2D(B,M,N,m,n))
{
  int i, j, k;
  DATA_TYPE temp2;

//BLAS PARAMS
//SIDE = 'L'
//UPLO = 'L'
// =>  Form  C := alpha*A*B + beta*C
// A is MxM
// B is MxN
// C is MxN
//note that due to Fortran array layout, the code below more closely resembles upper triangular case in BLAS
#pragma scop
   #pragma omp parallel for default(none) private(i, k, temp2) shared(A, B, C, alpha, beta, _PB_M, _PB_N)
   for (j = 0; j < _PB_N; j++) {
      for (i = 0; i < _PB_M; i++) {
         temp2 = 0;
         for (k = 0; k < i; k++) {
            C[k][j] += alpha * B[i][j] * A[i][k];
            temp2 += B[k][j] * A[i][k];
         }
         C[i][j] = beta * C[i][j] + alpha * B[i][j] * A[i][i] + alpha * temp2;
      }
   }
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int m = M;
  int n = N;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  DATA_TYPE beta;
  POLYBENCH_2D_ARRAY_DECL(C,DATA_TYPE,M,N,m,n);
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,M,M,m,m);
  POLYBENCH_2D_ARRAY_DECL(B,DATA_TYPE,M,N,m,n);

  /* Initialize array(s). */
  init_array (m, n, &alpha, &beta,
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_symm (m, n,
	       alpha, beta,
	       POLYBENCH_ARRAY(C),
	       POLYBENCH_ARRAY(A),
	       POLYBENCH_ARRAY(B));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, n, POLYBENCH_ARRAY(C)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);

  return 0;
}
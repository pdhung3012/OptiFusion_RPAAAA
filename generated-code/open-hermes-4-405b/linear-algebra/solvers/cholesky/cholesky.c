We are given a code for Cholesky decomposition. The scop region is the main computational kernel.
 The original code has the following structure:

   for (i = 0; i < _PB_N; i++) {
      // j < i
      for (j = 0; j < i; j++) {
         for (k = 0; k < j; k++) {
            A[i][j] -= A[i][k] * A[j][k];
         }
         A[i][j] /= A[j][j];
      }
      // i==j case
      for (k = 0; k < i; k++) {
         A[i][i] -= A[i][k] * A[i][k];
      }
      A[i][i] = SQRT_FUN(A[i][i]);
   }

 We are to optimize the code within the scop region.

 Step 1: Loop Transformations (in order: tiling, fusion, skewing)

 1. Tiling: 
    - The inner loops are over j and k. The outer loop is i.
    - The inner two loops (j and k) are nested and have dependencies. 
    - However, note that the inner loop for k in the j-loop runs from 0 to j-1, and the j-loop runs from 0 to i-1.
    - The i-loop is sequential because each iteration depends on the previous ones (the lower triangle is updated and then used in subsequent iterations).
    - Therefore, we cannot tile the i-loop in a way that breaks the sequential order. But we can consider tiling the j and k loops? 
      However, the j-loop is not parallel and the k-loop is innermost and short. Tiling might not be beneficial for the inner loops because the inner loop is short and the data is small.

    - Alternatively, we can consider tiling the i-loop? But note that the i-loop is sequential. We can break the i-loop into tiles and then process each tile? 
      However, the entire i-loop is a reduction on the lower triangle. The problem is that the update for a given i depends on all previous j (which are less than i). 
      Therefore, we cannot parallelize the i-loop. 

    - Conclusion: Tiling does not seem applicable because the outer loop is sequential and the inner loops are short and have dependencies that prevent tiling.

 2. Loop Fusion:
    - We have two inner loops: one for j (with an inner k) and one for k (for the diagonal). 
    - The j-loop and the diagonal k-loop are adjacent and both are inside the i-loop. 
    - However, the j-loop updates A[i][j] for j from 0 to i-1, and the diagonal loop updates A[i][i]. 
    - There is no fusion opportunity because the two loops are independent and operate on different parts of the array (the j-loop on the off-diagonal and the diagonal loop on the diagonal). 
    - Also, the j-loop must complete before the diagonal loop because the diagonal loop uses the updated A[i][k] for k < i? Actually, the diagonal loop uses A[i][k] which is updated in the j-loop? 
      But note: the j-loop updates A[i][j] for j < i, and the diagonal loop uses A[i][k] for k < i. So the j-loop must complete before the diagonal loop. 
    - Therefore, we cannot fuse the j-loop and the diagonal loop because they are not adjacent in the sense of having the same iteration space and because the j-loop must run first.

 3. Loop Skewing:
    - The inner loops have dependencies. The j-loop has a dependence on k: the update of A[i][j] uses A[i][k] and A[j][k] for k<j. 
      This is a reduction on A[i][j] and the inner k-loop is a reduction. 
    - The diagonal loop is a reduction on A[i][i]. 
    - There are no carried dependencies that can be skewed to expose parallelism because the i-loop is sequential and the inner loops are reductions.

 Therefore, we do not apply any loop transformations.

 Step 2: OpenMP Directive-Level Optimizations

 We now consider OpenMP directives for the scop region.

 The outer i-loop is sequential. However, we can look for parallelism in the inner loops.

 1. The j-loop: 
      for (j = 0; j < i; j++) {
         for (k = 0; k < j; k++) {
            A[i][j] -= A[i][k] * A[j][k];
         }
         A[i][j] /= A[j][j];
      }

    - The j-loop iterates over j from 0 to i-1. Each j is independent? 
      Let's check dependencies: 
        The update for A[i][j] uses A[i][k] and A[j][k] for k<j. 
        Note that A[i][k] is read and written in the same i iteration? Actually, for a fixed i, the j-loop updates A[i][j] for each j. 
        But note: when updating A[i][j], we use A[i][k] for k<j. This means that for a fixed i, the j-loop is a reduction on A[i][j] and also uses A[i][k] that have been updated in previous j iterations? 
        Actually, the inner k-loop for a fixed j uses A[i][k] for k<j. But note that the j-loop runs from 0 to i-1. For a fixed i, the j-loop updates A[i][j] in increasing j. 
        However, when updating A[i][j], we use A[i][k] for k<j. But these A[i][k] are updated in the same i iteration for k<j? 
        This is a problem: the j-loop is not parallel because the update of A[i][j] depends on A[i][k] for k<j, which are updated in the same j-loop for j>k. 

    - Therefore, the j-loop is not parallel.

 2. The diagonal k-loop:
      for (k = 0; k < i; k++) {
         A[i][i] -= A[i][k] * A[i][k];
      }

    - This is a reduction on A[i][i]. The loop is parallelizable with a reduction.

    - We can use OpenMP to parallelize the k-loop for the diagonal.

 3. The inner k-loop in the j-loop:
      for (k = 0; k < j; k++) {
         A[i][j] -= A[i][k] * A[j][k];
      }

    - This is a reduction on A[i][j]. However, note that the j-loop is sequential. But for a fixed j, the k-loop is a reduction and can be parallelized? 
      However, the j-loop is sequential and the k-loop is innermost. We can parallelize the k-loop? But note: the k-loop is short (k from 0 to j-1) and j is at most i-1, which is less than n. 
      The overhead of parallelization might be too high. Also, we are already in a sequential j-loop. 

    - Alternatively, we can try to parallelize the j-loop? But we saw that the j-loop is not parallel because of the dependence on A[i][k] for k<j. 

    - However, note that the j-loop updates A[i][j] and the update for a particular j uses A[i][k] for k<j. This is a sequential dependence. 

    - Therefore, we cannot parallelize the j-loop. 

    - But we can parallelize the inner k-loop? The k-loop is a reduction and independent for each k. So we can use a reduction on the k-loop? 

    - However, the k-loop is inside the j-loop and the j-loop is sequential. We can try to collapse the j and k loops? But note the bounds: j from 0 to i-1, and k from 0 to j-1. 
      This is a triangular loop. We cannot collapse because the inner loop bound depends on the outer loop index.

    - Alternatively, we can use a reduction on the inner k-loop? But note: the reduction is on A[i][j] and the k-loop is the only loop that updates A[i][j]. 
      We can do:

          #pragma omp parallel for reduction(+:A[i][j]) 
          for (k=0; k<j; k++)

      However, we must be cautious: the variable A[i][j] is shared and we are reducing into it. But note that the j-loop is sequential and the k-loop is the inner loop. 
      We can do:

          for (j=0; j<i; j++) {
             #pragma omp parallel for reduction(+:A[i][j])
             for (k=0; k<j; k++) {
                 A[i][j] -= A[i][k] * A[j][k];
             }
             A[i][j] /= A[j][j];
          }

      But note: the division by A[j][j] is done after the reduction. This is safe.

    - However, the inner k-loop might be too short to benefit from parallelization. Also, the overhead of creating a parallel region for each j might be high.

    - We can try to move the parallel region outside the j-loop? But then we would have to make the j-loop parallel, which we cannot because of the dependence.

    - Alternatively, we can use a single parallel region for the entire j-loop and then use a reduction on the k-loop? But the j-loop is sequential.

    - We decide to parallelize the inner k-loop for the j-loop and the diagonal k-loop.

 4. The diagonal k-loop: 
        for (k=0; k<i; k++) {
            A[i][i] -= A[i][k] * A[i][k];
        }

      We can do:

          #pragma omp parallel for reduction(+:A[i][i])
          for (k=0; k<i; k++) {
              A[i][i] -= A[i][k] * A[i][k];
          }

 5. We must be cautious about data sharing. We are updating A[i][j] and A[i][i]. We must ensure that the updates are safe.

 6. We also note that the entire scop region is inside a sequential i-loop. We are only parallelizing the inner k-loops.

 7. We must use 'default(none)' and specify the shared variables. The array A is shared. The indices i, j, k are private by default in the inner parallel regions? 
    But note: in the j-loop, j is the outer loop index and we are creating a parallel region for the k-loop. The variable j is shared in the parallel region? 
    Actually, in the j-loop, j is fixed for the inner k-loop. So we can make j shared. Similarly, i is fixed for the entire j-loop and diagonal loop. 

    We can do:

        for (j=0; j<i; j++) {
            #pragma omp parallel for default(none) shared(A, i, j) reduction(+:A[i][j])
            for (k=0; k<j; k++) {
                A[i][j] -= A[i][k] * A[j][k];
            }
            A[i][j] /= A[j][j];
        }

        #pragma omp parallel for default(none) shared(A, i) reduction(+:A[i][i])
        for (k=0; k<i; k++) {
            A[i][i] -= A[i][k] * A[i][k];
        }

 8. However, note that the inner k-loop in the j-loop uses A[j][k]. The j-loop index j is fixed for the inner k-loop, so it is safe to share j. 
    But note: the inner k-loop also uses A[i][k] and A[j][k]. The array A is shared. 

 9. We must be cautious about false sharing and cache performance. But the inner k-loop is short and the reduction is on a single element. 
    The reduction will be implemented by creating a private copy for each thread and then combining. 

 10. We can also consider SIMD for the inner k-loops? The inner k-loops are reductions and the operations are independent. 
     We can add 'simd' to the inner k-loops? 

     For the j-loop inner k-loop:

         #pragma omp parallel for simd reduction(+:A[i][j])

     But note: the reduction with simd might be supported. However, the reduction variable is a single element. 
     The compiler might be able to vectorize the reduction. 

     Similarly for the diagonal k-loop.

 11. We can also consider scheduling. Since the inner k-loops are short, we can use 'schedule(static)' which is the default.

 12. We must be cautious about the number of threads. We are creating a parallel region for each j and for the diagonal. 
     This might create too many parallel regions. We can try to move the parallel region to the outer level? 

     However, the j-loop is sequential. We cannot parallelize the j-loop. 

     Alternatively, we can use a single parallel region for the entire i iteration? But the i-loop is sequential. 

     We decide to keep the parallel regions inside the j-loop and the diagonal loop. 

     But note: the overhead of creating a parallel region for each j might be high. We can use 'nowait' to avoid the barrier? 
     However, the j-loop is sequential and we must wait for the inner k-loop to finish before the next j. So we cannot use 'nowait'. 

     We can try to use a single parallel region for the entire j-loop? 

        #pragma omp parallel
        {
            for (j=0; j<i; j++) {
                #pragma omp for reduction(+:A[i][j])
                for (k=0; k<j; k++) {
                    A[i][j] -= A[i][k] * A[j][k];
                }
                // But note: the division must be done after the reduction and only by one thread? 
                // We cannot do the division inside the parallel region without synchronization.
            }
        }

     This is problematic because the division step must be done after the entire k-loop for that j. 

     We can do:

        #pragma omp parallel
        {
            for (j=0; j<i; j++) {
                #pragma omp single
                {
                    // We do the entire j iteration in a single thread? Then we lose parallelism.
                }
            }
        }

     This is not helpful.

 13. We decide to keep the parallel for inside the j-loop and the diagonal loop. 

 14. We can also consider the diagonal loop: it is a single loop and we can parallelize it without creating a parallel region for each i? 
     But note: the diagonal loop is inside the i-loop and the i-loop is sequential. 

 15. We must also consider the possibility of using SIMD without OpenMP parallel for? The compiler might vectorize the inner k-loops. 
     But we are already using OpenMP for parallel for. 

 16. We decide to use:

        for (i = 0; i < _PB_N; i++) {
           // j-loop
           for (j = 0; j < i; j++) {
              #pragma omp parallel for simd reduction(+:A[i][j])
              for (k = 0; k < j; k++) {
                 A[i][j] -= A[i][k] * A[j][k];
              }
              A[i][j] /= A[j][j];
           }
           // diagonal
           #pragma omp parallel for simd reduction(+:A[i][i])
           for (k = 0; k < i; k++) {
              A[i][i] -= A[i][k] * A[i][k];
           }
           A[i][i] = SQRT_FUN(A[i][i]);
        }

 17. We must specify the data sharing. We use 'default(none)' and list the shared variables. 

     For the j-loop inner k-loop: 
        shared: A, i, j
        reduction: A[i][j] (which is a single element)

     For the diagonal k-loop:
        shared: A, i
        reduction: A[i][i]

 18. We must be cautious: the array A is shared and we are updating different elements in the j-loop and the diagonal. 
     The j-loop updates A[i][j] for j from 0 to i-1, and the diagonal updates A[i][i]. 
     There is no overlap. 

 19. We also note that the inner k-loop in the j-loop uses A[j][k] for k from 0 to j-1. This part of the array is read-only in the j-loop? 
     Actually, for a fixed i, the j-loop runs for j from 0 to i-1. The element A[j][k] is updated in the j-loop when the outer i-loop index equals j? 
     But note: the outer i-loop runs from 0 to n-1. For a fixed j, the update of A[j][k] happens when i=j. 
     In the current i iteration, we are updating A[i][j] for j < i. The element A[j][k] for j < i is from a previous i iteration (when i=j) and is not updated in the current i iteration. 
     Therefore, it is safe to read A[j][k] in the current i iteration.

 20. We must also consider the possibility of using 'aligned' for the array A? The array A is allocated by the polybench framework. 
     We don't know the alignment. We can skip the aligned clause.

 21. We can also consider the schedule. The default schedule is static, which is fine for the inner k-loops.

 22. We can also consider the number of threads. We can let the runtime decide.

 23. We must be cautious about the reduction: the reduction on a single element might be inefficient. 
     The OpenMP reduction for a scalar is efficient. 

 24. We decide to use the above.

 However, note: the inner k-loop in the j-loop is very short (k from 0 to j-1, and j is at most i-1, which is less than n). 
 The overhead of creating a parallel region for each j might be too high. 

 Alternative approach: 

    We can try to parallelize the j-loop? But we saw that the j-loop is not parallel because of the dependence on A[i][k] for k<j. 

    However, note that the dependence is a reduction on A[i][j] and the inner k-loop is a reduction. 
    We can restructure the j-loop to use a reduction on the entire j-loop? 

    The update for A[i][j] is:

        A[i][j] = (A[i][j] - sum_{k=0}^{j-1} A[i][k] * A[j][k]) / A[j][j]

    We can precompute the entire lower triangle for a fixed i? 

    But note: the j-loop runs from 0 to i-1. The update for A[i][j] uses A[i][k] for k<j. 
    This means that the update for A[i][j] depends on the updated values of A[i][k] for k<j. 

    Therefore, we cannot change the order of j. 

    We can try to use a technique called "loop fission" to separate the reduction part? 

    We can do:

        for (j=0; j<i; j++) {
            // We cannot parallelize the j-loop because of the dependence on A[i][k] for k<j.
        }

    So we stick to parallelizing the inner k-loops.

 25. We decide to use the inner k-loop parallelization as described.

 26. We also note that the diagonal k-loop is a reduction and can be parallelized similarly.

 27. We write the code accordingly.

 28. We must preserve the original semantics and the scop markers.

 Final code for the scop region:

   for (i = 0; i < _PB_N; i++) {
      // j-loop
      for (j = 0; j < i; j++) {
         #pragma omp parallel for simd default(none) shared(A, i, j) reduction(+:A[i][j])
         for (k = 0; k < j; k++) {
            A[i][j] -= A[i][k] * A[j][k];
         }
         A[i][j] /= A[j][j];
      }
      // diagonal
      #pragma omp parallel for simd default(none) shared(A, i) reduction(+:A[i][i])
      for (k = 0; k < i; k++) {
         A[i][i] -= A[i][k] * A[i][k];
      }
      A[i][i] = SQRT_FUN(A[i][i]);
   }

 29. We must check the data types: the reduction is on a DATA_TYPE. The reduction operation is addition (for the subtraction, we are doing a reduction with subtraction? 
     Actually, the reduction is for the entire expression: 

        A[i][j] = A[i][j] - (A[i][k] * A[j][k])

     This is equivalent to:

        A[i][j] -= A[i][k] * A[j][k]

     So the reduction is a subtraction. But note: the reduction clause in OpenMP for subtraction is not directly supported. 
     We can use:

        reduction(-:A[i][j])

     However, the OpenMP standard (since OpenMP 3.1) supports user-defined reductions? But we are using a built-in type. 
     Alternatively, we can change the reduction to use addition by negating the product? 

        A[i][j] += (-A[i][k] * A[j][k])

     But that is the same as subtraction. 

     The OpenMP reduction clause for C/C++ supports the following operators: +, *, -, &, |, ^, &&, ||. 
     So we can use:

        reduction(-:A[i][j])

     However, note that the initial value for a subtraction reduction is 0? Actually, the initial value for a reduction with subtraction is 0, and then we do:

        A[i][j] = 0 - (sum of (A[i][k]*A[j][k]) for each thread) 

     But that is not what we want. We want:

        A[i][j] = original_A[i][j] - (sum of ...)

     How does OpenMP reduction work? 

        The reduction clause creates a private copy for each thread, initialized to the identity value for the operation. 
        For subtraction, the identity is 0? Then the private copy is 0, and then we do:

            private_A[i][j] -= A[i][k] * A[j][k]

        Then at the end, the original A[i][j] is updated by:

            A[i][j] = A[i][j] - private_A[i][j]   ??? 

        Actually, the reduction operation is applied in the same way as the operator. For subtraction, the reduction is:

            A[i][j] = A[i][j] - (private_A[i][j] of thread0 + private_A[i][j] of thread1 + ...)

        But note: the private copy is initialized to 0, and then we do:

            private_A[i][j] = private_A[i][j] - (A[i][k]*A[j][k])

        Then the reduction step is:

            A[i][j] = A[i][j] - private_A[i][j]

        This is not the same as the original. The original is:

            A[i][j] = A[i][j] - (A[i][k]*A[j][k])

        But in the parallel version, we are doing:

            A[i][j] = A[i][j] - (sum of (A[i][k]*A[j][k]) for all k)

        So it is the same.

        However, note that the reduction variable is the same for all threads? Actually, the reduction variable is shared and the private copies are combined by the reduction operator. 
        The reduction operator for subtraction is not commutative? But the order of subtraction does not matter because subtraction is not associative? 

        Actually, the reduction for subtraction is not associative. But note: the inner k-loop is a reduction and the order of the k-loop does not matter? 
        The operation is:

            A[i][j] = A[i][j] - (A[i][k1]*A[j][k1] + A[i][k2]*A[j][k2] + ...)

        So it is a linear operation and the order of the terms does not matter. 

        Therefore, we can use:

            reduction(-:A[i][j])

        But note: the OpenMP specification says that for C/C++, the reduction operator for subtraction is allowed. 

        However, we must be cautious: the initial value for the private copy is 0, and then we subtract the product. Then the reduction step subtracts the private copy from the original. 
        This is equivalent to:

            A[i][j] = A[i][j] - (sum of products)

        So it is correct.

 30. Alternatively, we can change the code to use a temporary variable for the reduction? 

        double tmp = 0;
        #pragma omp parallel for reduction(+:tmp)
        for (k=0; k<j; k++) {
            tmp += A[i][k] * A[j][k];
        }
        A[i][j] -= tmp;

     This avoids the issue with the reduction operator. And it is more efficient? 

     We can do:

        for (j=0; j<i; j++) {
            DATA_TYPE tmp = 0;
            #pragma omp parallel for simd reduction(+:tmp)
            for (k=0; k<j; k++) {
                tmp += A[i][k] * A[j][k];
            }
            A[i][j] = (A[i][j] - tmp) / A[j][j];
        }

     Similarly for the diagonal:

        DATA_TYPE tmp = 0;
        #pragma omp parallel for simd reduction(+:tmp)
        for (k=0; k<i; k++) {
            tmp += A[i][k] * A[i][k];
        }
        A[i][i] -= tmp;

     This is cleaner and avoids the issue with the reduction operator for subtraction. 

     We choose this method.

 31. We rewrite the scop region accordingly.

 32. We must declare the temporary variable. We can declare it inside the j-loop and the diagonal loop.

 33. We must be cautious about the data type: DATA_TYPE.

 34. The code becomes:

      for (i = 0; i < _PB_N; i++) {
         // j-loop
         for (j = 0; j < i; j++) {
            DATA_TYPE tmp = 0;
            #pragma omp parallel for simd default(none) shared(A, i, j) reduction(+:tmp)
            for (k = 0; k < j; k++) {
               tmp += A[i][k] * A[j][k];
            }
            A[i][j] = (A[i][j] - tmp) / A[j][j];
         }
         // diagonal
         DATA_TYPE tmp = 0;
         #pragma omp parallel for simd default(none) shared(A, i) reduction(+:tmp)
         for (k = 0; k < i; k++) {
            tmp += A[i][k] * A[i][k];
         }
         A[i][i] -= tmp;
         A[i][i] = SQRT_FUN(A[i][i]);
      }

 35. We must check the data sharing: 
        - In the j-loop: A, i, j are shared. tmp is reduced.
        - In the diagonal: A, i are shared. tmp is reduced.

 36. We must also note that the inner k-loop in the j-loop uses A[j][k]. The j-loop index j is fixed and shared. 

 37. We must be cautious about the performance: creating a parallel region for each j and for the diagonal might be expensive. 
     We can try to use a single parallel region for the entire i iteration? But the i-loop is sequential. 

     Alternatively, we can use orphaned directives and a single parallel region outside the i-loop? 

        #pragma omp parallel
        {
            for (i = 0; i < _PB_N; i++) {
               // j-loop: sequential
               for (j = 0; j < i; j++) {
                  DATA_TYPE tmp = 0;
                  #pragma omp for simd reduction(+:tmp) nowait
                  for (k = 0; k < j; k++) {
                     tmp += A[i][k] * A[j][k];
                  }
                  // But note: the j-loop is sequential and we must wait for the k-loop to finish? 
                  // We cannot use nowait because the next j depends on the previous j? 
                  // Actually, the j-loop is sequential and we must wait for the entire j iteration to finish? 
                  // We can use a barrier? But that would be expensive.

                  // Instead, we do the reduction and then update A[i][j] without a barrier? 
                  // But the update of A[i][j] is done by the master thread? 
                  // We can do:

                  #pragma omp single
                  {
                      A[i][j] = (A[i][j] - tmp) / A[j][j];
                  }
               }
               // diagonal
               DATA_TYPE tmp = 0;
               #pragma omp for simd reduction(+:tmp) nowait
               for (k = 0; k < i; k++) {
                  tmp += A[i][k] * A[i][k];
               }
               #pragma omp single
               {
                   A[i][i] -= tmp;
                   A[i][i] = SQRT_FUN(A[i][i]);
               }
            }
        }

     This avoids creating a parallel region for each j and the diagonal. But note: the j-loop is sequential and we are using a single thread for the j-loop? 
     The j-loop is inside the parallel region and we are using a 'single' construct for the update. The k-loop is distributed among threads. 

     However, the j-loop is sequential and we are not parallelizing the j-loop. The k-loop is parallelized. 

     But note: the j-loop runs for j from 0 to i-1. The k-loop for a fixed j is distributed. Then we do a single for the update. 
     Then we move to the next j. 

     This might be efficient. 

     But note: the entire i-loop is inside the parallel region. The i-loop is sequential and we are not parallelizing it. 

     We can do:

        #pragma omp parallel
        {
            for (i = 0; i < _PB_N; i++) {
               // j-loop: sequential, so we do it in a single thread? 
               for (j = 0; j < i; j++) {
                  DATA_TYPE tmp = 0;
                  #pragma omp for simd reduction(+:tmp) nowait
                  for (k = 0; k < j; k++) {
                     tmp += A[i][k] * A[j][k];
                  }
                  // We must wait for the k-loop to finish? The nowait means we don't wait, but the next j might start immediately? 
                  // But the j-loop is sequential and we are in a worksharing region? 
                  // Actually, the j-loop is not a worksharing region. We are in a parallel region and the j-loop is executed by all threads? 
                  // This is not what we want.

                  // We must make the j-loop sequential. We can use:

                  #pragma omp single
                  {
                      for (j = 0; j < i; j++) {
                         ... 
                      }
                  }

                  // But then we lose the parallelism for the k-loop? 

                  // Alternatively, we can use a barrier after the k-loop? 

                  // We decide to use a barrier after the k-loop? 

                  // But note: the j-loop is sequential and we must do the j iterations in order. 

                  // We can do:

                  #pragma omp single
                  for (j = 0; j < i; j++) {
                     DATA_TYPE tmp = 0;
                     #pragma omp for simd reduction(+:tmp) nowait
                     for (k = 0; k < j; k++) {
                         tmp += A[i][k] * A[j][k];
                     }
                     // We must wait for the k-loop to finish? We can use a barrier? 
                     #pragma omp barrier
                     A[i][j] = (A[i][j] - tmp) / A[j][j];
                  }

                  // But the barrier inside the j-loop would be expensive.

               }
            }
        }

     This is complex and might not be efficient.

 38. We decide to stick to the original plan of having a parallel region for each j and for the diagonal. 
     We hope that the inner k-loops are long enough to amortize the overhead. 

     But note: the inner k-loop for the j-loop has j iterations, and j is at most i-1, which is less than n. 
     The diagonal k-loop has i iterations. 

     We can set a threshold: only parallelize if the number of iterations is above a certain value? 
     But the problem size n might be large. 

     We decide to use the version with the temporary variable and without the single parallel region.

 39. We write the code as:

        for (i = 0; i < _PB_N; i++) {
           for (j = 0; j < i; j++) {
              DATA_TYPE tmp = 0;
              #pragma omp parallel for simd default(none) shared(A, i, j) reduction(+:tmp)
              for (k = 0; k < j; k++) {
                 tmp += A[i][k] * A[j][k];
              }
              A[i][j] = (A[i][j] - tmp) / A[j][j];
           }
           DATA_TYPE tmp = 0;
           #pragma omp parallel for simd default(none) shared(A, i) reduction(+:tmp)
           for (k = 0; k < i; k++) {
              tmp += A[i][k] * A[i][k];
           }
           A[i][i] -= tmp;
           A[i][i] = SQRT_FUN(A[i][i]);
        }

 40. We must be cautious about the data type of tmp: it is DATA_TYPE.

 41. We must also note that the original code uses SQRT_FUN, which is a macro for sqrt or sqrtf.

 42. We preserve the scop markers.

 43. We output the entire code with the scop region replaced.

 However, note: the original code uses _PB_N. We must not change that.

 Let's write the entire kernel_cholesky function with the optimized scop region.
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "cholesky.h"


/* Array initialization. */
static
void init_array(int n,
		DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    {
      for (j = 0; j <= i; j++)
	A[i][j] = (DATA_TYPE)(-j % n) / n + 1;
      for (j = i+1; j < n; j++) {
	A[i][j] = 0;
      }
      A[i][i] = 1;
    }

  /* Make the matrix positive semi-definite. */
  int r,s,t;
  POLYBENCH_2D_ARRAY_DECL(B, DATA_TYPE, N, N, n, n);
  for (r = 0; r < n; ++r)
    for (s = 0; s < n; ++s)
      (POLYBENCH_ARRAY(B))[r][s] = 0;
  for (t = 0; t < n; ++t)
    for (r = 0; r < n; ++r)
      for (s = 0; s < n; ++s)
	(POLYBENCH_ARRAY(B))[r][s] += A[r][t] * A[s][t];
    for (r = 0; r < n; ++r)
      for (s = 0; s < n; ++s)
	A[r][s] = (POLYBENCH_ARRAY(B))[r][s];
  POLYBENCH_FREE_ARRAY(B);

}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("A");
  for (i = 0; i < n; i++)
    for (j = 0; j <= i; j++) {
    if ((i * n + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j]);
  }
  POLYBENCH_DUMP_END("A");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_cholesky(int n,
		     DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int i, j, k;


#pragma scop
  for (i = 0; i < _PB_N; i++) {
     for (j = 0; j < i; j++) {
        DATA_TYPE tmp = 0;
        #pragma omp parallel for simd default(none) shared(A, i, j) reduction(+:tmp)
        for (k = 0; k < j; k++) {
           tmp += A[i][k] * A[j][k];
        }
        A[i][j] = (A[i][j] - tmp) / A[j][j];
     }
     DATA_TYPE tmp = 0;
     #pragma omp parallel for simd default(none) shared(A, i) reduction(+:tmp)
     for (k = 0; k < i; k++) {
        tmp += A[i][k] * A[i][k];
     }
     A[i][i] -= tmp;
     A[i][i] = SQRT_FUN(A[i][i]);
  }
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);

  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(A));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_cholesky (n, POLYBENCH_ARRAY(A));

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
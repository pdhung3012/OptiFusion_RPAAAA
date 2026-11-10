void kernel_doitgen(int nr, int nq, int np,
		    DATA_TYPE POLYBENCH_3D(A,NR,NQ,NP,nr,nq,np),
		    DATA_TYPE POLYBENCH_2D(C4,NP,NP,np,np),
		    DATA_TYPE POLYBENCH_1D(sum,NP,np))
{
  int r, q, p, s;

#pragma scop
#pragma omp parallel private(r, q, p, s)
  {
    DATA_TYPE sum_loc[NP];

#pragma omp for collapse(2) schedule(static)
    for (r = 0; r < _PB_NR; r++)
      for (q = 0; q < _PB_NQ; q++) {
        for (p = 0; p < _PB_NP; p++)
          sum_loc[p] = SCALAR_VAL(0.0);
        for (s = 0; s < _PB_NP; s++) {
          DATA_TYPE a = A[r][q][s];
#pragma omp simd
          for (p = 0; p < _PB_NP; p++)
            sum_loc[p] += a * C4[s][p];
        }
        for (p = 0; p < _PB_NP; p++)
          A[r][q][p] = sum_loc[p];
        if (r == _PB_NR - 1 && q == _PB_NQ - 1)
          for (p = 0; p < _PB_NP; p++)
            sum[p] = sum_loc[p];
      }
  }
#pragma endscop

}
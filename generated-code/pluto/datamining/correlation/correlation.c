#include <omp.h>
#include <math.h>
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* correlation.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "correlation.h"


/* Array initialization. */
static
void init_array (int m,
		 int n,
		 DATA_TYPE *float_n,
		 DATA_TYPE POLYBENCH_2D(data,N,M,n,m))
{
  int i, j;

  *float_n = (DATA_TYPE)N;

  for (i = 0; i < N; i++)
    for (j = 0; j < M; j++)
      data[i][j] = (DATA_TYPE)(i*j)/M + i;

}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int m,
		 DATA_TYPE POLYBENCH_2D(corr,M,M,m,m))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("corr");
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++) {
      if ((i * m + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, corr[i][j]);
    }
  POLYBENCH_DUMP_END("corr");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_correlation(int m, int n,
			DATA_TYPE float_n,
			DATA_TYPE POLYBENCH_2D(data,N,M,n,m),
			DATA_TYPE POLYBENCH_2D(corr,M,M,m,m),
			DATA_TYPE POLYBENCH_1D(mean,M,m),
			DATA_TYPE POLYBENCH_1D(stddev,M,m))
{
  int i, j, k;

  DATA_TYPE eps = SCALAR_VAL(0.1);


  int t1, t2, t3, t4, t5, t6, t7, t8;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
corr[m - 1][m - 1] = 1.;;
lbp=0;
ubp=floord(m-2,32);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8)
for (t2=lbp;t2<=ubp;t2++) {
  for (t3=t2;t3<=floord(m-1,32);t3++) {
    for (t4=32*t2;t4<=min(min(m-2,32*t2+31),32*t3+30);t4++) {
      lbv=max(32*t3,t4+1);
      ubv=min(m-1,32*t3+31);
#pragma ivdep
#pragma vector always
      for (t5=lbv;t5<=ubv;t5++) {
        corr[t4][t5] = 0.;;
      }
    }
  }
}
lbp=0;
ubp=floord(m-1,32);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8)
for (t2=lbp;t2<=ubp;t2++) {
  lbv=32*t2;
  ubv=min(m-2,32*t2+31);
#pragma ivdep
#pragma vector always
  for (t3=lbv;t3<=ubv;t3++) {
    corr[t3][t3] = 1.;;
    stddev[t3] = 0.;;
    mean[t3] = 0.;;
  }
  if (t2 >= ceild(m-32,32)) {
    stddev[(m-1)] = 0.;;
    mean[(m-1)] = 0.;;
  }
}
if (n >= 1) {
  lbp=0;
  ubp=floord(m-1,32);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8)
  for (t2=lbp;t2<=ubp;t2++) {
    for (t3=0;t3<=floord(n-1,32);t3++) {
      for (t4=32*t3;t4<=(min(n-1,32*t3+31))-7;t4+=8) {
        lbv=32*t2;
        ubv=min(m-1,32*t2+31);
#pragma ivdep
#pragma vector always
        for (t5=lbv;t5<=ubv;t5++) {
          mean[t5] += data[t4][t5];;
          mean[t5] += data[(t4+1)][t5];;
          mean[t5] += data[(t4+2)][t5];;
          mean[t5] += data[(t4+3)][t5];;
          mean[t5] += data[(t4+4)][t5];;
          mean[t5] += data[(t4+5)][t5];;
          mean[t5] += data[(t4+6)][t5];;
          mean[t5] += data[(t4+7)][t5];;
        }
      }
      for (;t4<=min(n-1,32*t3+31);t4++) {
        lbv=32*t2;
        ubv=min(m-1,32*t2+31);
#pragma ivdep
#pragma vector always
        for (t5=lbv;t5<=ubv;t5++) {
          mean[t5] += data[t4][t5];;
        }
      }
    }
  }
}
lbp=0;
ubp=floord(m-1,32);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8)
for (t2=lbp;t2<=ubp;t2++) {
  lbv=32*t2;
  ubv=min(m-1,32*t2+31);
#pragma ivdep
#pragma vector always
  for (t3=lbv;t3<=ubv;t3++) {
    mean[t3] /= float_n;;
  }
}
if (n >= 1) {
  lbp=0;
  ubp=floord(m-1,32);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8)
  for (t2=lbp;t2<=ubp;t2++) {
    for (t3=0;t3<=floord(n-1,32);t3++) {
      for (t4=32*t3;t4<=(min(n-1,32*t3+31))-7;t4+=8) {
        lbv=32*t2;
        ubv=min(m-1,32*t2+31);
#pragma ivdep
#pragma vector always
        for (t5=lbv;t5<=ubv;t5++) {
          stddev[t5] += ((data[t4][t5] - mean[t5]) * (data[t4][t5] - mean[t5]));;
          data[t4][t5] -= mean[t5];;
          stddev[t5] += ((data[(t4+1)][t5] - mean[t5]) * (data[(t4+1)][t5] - mean[t5]));;
          data[(t4+1)][t5] -= mean[t5];;
          stddev[t5] += ((data[(t4+2)][t5] - mean[t5]) * (data[(t4+2)][t5] - mean[t5]));;
          data[(t4+2)][t5] -= mean[t5];;
          stddev[t5] += ((data[(t4+3)][t5] - mean[t5]) * (data[(t4+3)][t5] - mean[t5]));;
          data[(t4+3)][t5] -= mean[t5];;
          stddev[t5] += ((data[(t4+4)][t5] - mean[t5]) * (data[(t4+4)][t5] - mean[t5]));;
          data[(t4+4)][t5] -= mean[t5];;
          stddev[t5] += ((data[(t4+5)][t5] - mean[t5]) * (data[(t4+5)][t5] - mean[t5]));;
          data[(t4+5)][t5] -= mean[t5];;
          stddev[t5] += ((data[(t4+6)][t5] - mean[t5]) * (data[(t4+6)][t5] - mean[t5]));;
          data[(t4+6)][t5] -= mean[t5];;
          stddev[t5] += ((data[(t4+7)][t5] - mean[t5]) * (data[(t4+7)][t5] - mean[t5]));;
          data[(t4+7)][t5] -= mean[t5];;
        }
      }
      for (;t4<=min(n-1,32*t3+31);t4++) {
        lbv=32*t2;
        ubv=min(m-1,32*t2+31);
#pragma ivdep
#pragma vector always
        for (t5=lbv;t5<=ubv;t5++) {
          stddev[t5] += ((data[t4][t5] - mean[t5]) * (data[t4][t5] - mean[t5]));;
          data[t4][t5] -= mean[t5];;
        }
      }
    }
  }
}
lbp=0;
ubp=floord(m-1,32);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8)
for (t2=lbp;t2<=ubp;t2++) {
  lbv=32*t2;
  ubv=min(m-1,32*t2+31);
#pragma ivdep
#pragma vector always
  for (t3=lbv;t3<=ubv;t3++) {
    stddev[t3] /= float_n;;
    stddev[t3] = sqrt(stddev[t3]);;
    stddev[t3] = ((stddev[t3] <= eps) ? 1. : stddev[t3]);;
  }
}
lbp=0;
ubp=floord(n-1,32);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8)
for (t2=lbp;t2<=ubp;t2++) {
  for (t3=0;t3<=floord(m-1,32);t3++) {
    for (t4=32*t2;t4<=(min(n-1,32*t2+31))-7;t4+=8) {
      lbv=32*t3;
      ubv=min(m-1,32*t3+31);
#pragma ivdep
#pragma vector always
      for (t5=lbv;t5<=ubv;t5++) {
        data[t4][t5] /= (sqrt(float_n) * stddev[t5]);;
        data[(t4+1)][t5] /= (sqrt(float_n) * stddev[t5]);;
        data[(t4+2)][t5] /= (sqrt(float_n) * stddev[t5]);;
        data[(t4+3)][t5] /= (sqrt(float_n) * stddev[t5]);;
        data[(t4+4)][t5] /= (sqrt(float_n) * stddev[t5]);;
        data[(t4+5)][t5] /= (sqrt(float_n) * stddev[t5]);;
        data[(t4+6)][t5] /= (sqrt(float_n) * stddev[t5]);;
        data[(t4+7)][t5] /= (sqrt(float_n) * stddev[t5]);;
      }
    }
    for (;t4<=min(n-1,32*t2+31);t4++) {
      lbv=32*t3;
      ubv=min(m-1,32*t3+31);
#pragma ivdep
#pragma vector always
      for (t5=lbv;t5<=ubv;t5++) {
        data[t4][t5] /= (sqrt(float_n) * stddev[t5]);;
      }
    }
  }
}
if (n >= 1) {
  lbp=0;
  ubp=floord(m-2,32);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8)
  for (t2=lbp;t2<=ubp;t2++) {
    for (t3=t2;t3<=floord(m-1,32);t3++) {
      for (t4=0;t4<=floord(n-1,32);t4++) {
        for (t5=32*t2;t5<=min(min(m-2,32*t2+31),32*t3+30);t5++) {
          for (t6=32*t4;t6<=(min(n-1,32*t4+31))-7;t6+=8) {
            lbv=max(32*t3,t5+1);
            ubv=min(m-1,32*t3+31);
#pragma ivdep
#pragma vector always
            for (t7=lbv;t7<=ubv;t7++) {
              corr[t5][t7] += (data[t6][t5] * data[t6][t7]);;
              corr[t5][t7] += (data[(t6+1)][t5] * data[(t6+1)][t7]);;
              corr[t5][t7] += (data[(t6+2)][t5] * data[(t6+2)][t7]);;
              corr[t5][t7] += (data[(t6+3)][t5] * data[(t6+3)][t7]);;
              corr[t5][t7] += (data[(t6+4)][t5] * data[(t6+4)][t7]);;
              corr[t5][t7] += (data[(t6+5)][t5] * data[(t6+5)][t7]);;
              corr[t5][t7] += (data[(t6+6)][t5] * data[(t6+6)][t7]);;
              corr[t5][t7] += (data[(t6+7)][t5] * data[(t6+7)][t7]);;
            }
          }
          for (;t6<=min(n-1,32*t4+31);t6++) {
            lbv=max(32*t3,t5+1);
            ubv=min(m-1,32*t3+31);
#pragma ivdep
#pragma vector always
            for (t7=lbv;t7<=ubv;t7++) {
              corr[t5][t7] += (data[t6][t5] * data[t6][t7]);;
            }
          }
        }
      }
    }
  }
}
lbp=0;
ubp=floord(m-2,32);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8)
for (t2=lbp;t2<=ubp;t2++) {
  for (t3=t2;t3<=floord(m-1,32);t3++) {
    for (t4=32*t2;t4<=min(min(m-2,32*t2+31),32*t3+30);t4++) {
      lbv=max(32*t3,t4+1);
      ubv=min(m-1,32*t3+31);
#pragma ivdep
#pragma vector always
      for (t5=lbv;t5<=ubv;t5++) {
        corr[t5][t4] = corr[t4][t5];;
      }
    }
  }
}

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;
  int m = M;

  /* Variable declaration/allocation. */
  DATA_TYPE float_n;
  POLYBENCH_2D_ARRAY_DECL(data,DATA_TYPE,N,M,n,m);
  POLYBENCH_2D_ARRAY_DECL(corr,DATA_TYPE,M,M,m,m);
  POLYBENCH_1D_ARRAY_DECL(mean,DATA_TYPE,M,m);
  POLYBENCH_1D_ARRAY_DECL(stddev,DATA_TYPE,M,m);

  /* Initialize array(s). */
  init_array (m, n, &float_n, POLYBENCH_ARRAY(data));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_correlation (m, n, float_n,
		      POLYBENCH_ARRAY(data),
		      POLYBENCH_ARRAY(corr),
		      POLYBENCH_ARRAY(mean),
		      POLYBENCH_ARRAY(stddev));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, POLYBENCH_ARRAY(corr)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(data);
  POLYBENCH_FREE_ARRAY(corr);
  POLYBENCH_FREE_ARRAY(mean);
  POLYBENCH_FREE_ARRAY(stddev);

  return 0;
}

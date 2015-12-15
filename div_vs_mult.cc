#include <omp.h>
#include <cstdio>

typedef double FLOAT;

const int n = 80000;
const int m = 10000;
FLOAT a[n];

int main()
{

  // multiplication
  double t_m = omp_get_wtime();
  for (int i = 0; i < m; ++i)
  {

    //#pragma omp simd
    //#pragma omp parallel for shared(a)
    for (int j = 0; j < n; ++j)
    {
      FLOAT x = (FLOAT)i;
      FLOAT y = (FLOAT)j;
      FLOAT z = (FLOAT)j;
      a[j] = (x * y + z * y + x * z) * (1.0 + z); // 7
    }
  }
  t_m = omp_get_wtime() - t_m;

  // division
  double t_d = omp_get_wtime();
  for (int i = 0; i < m; ++i)
  {
    //#pragma omp simd
    //#pragma omp parallel for shared(a)
    for (int j = 0; j < n; ++j)
    {
      FLOAT x = (FLOAT)i;
      FLOAT y = (FLOAT)j;
      FLOAT z = (FLOAT)j;
      a[j] = (x * y + z * y + x * z) / (1.0 + z); //
    }
  }
  t_d = omp_get_wtime() - t_d;

  double flops = (double)m * (double)n * 7.0;
  double gflops_s_m = flops/1.e9/t_m;
  double gflops_s_d = flops/1.e9/t_d;
  double ratio = gflops_s_m / gflops_s_d;
  double a = t_m / 7;
  double b = t_d - 6*t_m;

  printf("multiplication GLOPS/s: %f \n", gflops_s_m);
  printf("      division GLOPS/s: %f \n", gflops_s_d);
  printf("                 ratio: %f \n", ratio);
  printf("   time of mult or add: %f \n", a);
  printf("      time of division: %f \n", b);
  printf("                 ratio: %f \n", b/a);

}

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

int main()
{
  FILE *file;                       // 文件指针
  double x[1024];                   // 用于存储x数值的数组
  double y[1024];                   // 用于存储y数值的数组
  int count = 0;                    // 记录读入数据点的个数
  file = fopen("ode_out.dat", "r"); // 打开文件
  while (count < 1024 && fscanf(file, "%lf %lf", &x[count], &y[count]) == 2)// 读取文件中的浮点数
  {
    count++; // 成功读取一个浮点数后，计数器加1
  }
  fclose(file);
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, count);
  gsl_spline_init(spline, x, y, count);
  double x_spline, y_spline;
  for (x_spline = -4; x_spline <= 4; x_spline += 0.01)
  {
    y_spline = gsl_spline_eval(spline, x_spline, acc);
    printf("%f %f\n", x_spline, y_spline);
  }
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return 0;
}
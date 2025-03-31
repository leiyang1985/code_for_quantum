#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

double psi2(double x,void *par)
{
  void **par_arr=(void **)par;
  gsl_interp_accel *acc=(gsl_interp_accel *)(par_arr[0]);
  gsl_spline *spline=(gsl_spline *)(par_arr[1]);
  double psi=gsl_spline_eval(spline, x, acc);
  return psi*psi;
}
double psid2psi(double x,void *par)
{
  void **par_arr=(void **)par;
  gsl_interp_accel *acc=(gsl_interp_accel *)(par_arr[0]);
  gsl_spline *spline=(gsl_spline *)(par_arr[1]);
  double norm=*((double *)(par_arr[2]));
  double psi=gsl_spline_eval(spline, x, acc);
  double d2psi=gsl_spline_eval_deriv2(spline,x,acc);
  return -norm*norm*psi*d2psi;
}
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
  
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

  void *par[3];
  par[0]=(void *)acc;
  par[1]=(void *)spline;
  gsl_function F;
  F.function = &psi2;
  F.params = (void *)par;

  double norm,norm_err;
  gsl_integration_qags(&F, -4, 4, 1e-3, 1e-6, 1000, w, &norm, &norm_err);

  norm=1/sqrt(norm);
  par[2]=(void *)(&norm);
  F.function=&psid2psi;

  double T2,T2_err;
  gsl_integration_qags(&F, -4, 4, 1e-3, 1e-6, 1000, w, &T2, &T2_err);

  printf("<T^2>= %f\n",T2);

  gsl_integration_workspace_free(w);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return 0;
}
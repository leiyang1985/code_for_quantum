#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

int func (double x, const double y[], double f[],void *params)
{
  double E = *(double *)params;
  f[0] = y[1];
  f[1] = (x*x-E)*y[0];
  return GSL_SUCCESS;
}

int func_the (double x, const double y[], double f[],void *params)
{
  double E = *(double *)params;
  f[0] = y[1];
  f[1] = 2*tan(x)*y[1]+(1+tan(x)*tan(x))*(1+tan(x)*tan(x))*(tan(x)*tan(x)-E)*y[0];
  return GSL_SUCCESS;
}

int main ()
{
  double E=11;
    gsl_odeiv2_system sys = {func, NULL, 2, &E};
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-7, 1e-7, 0);

  double x0 = -10;
  double y[2] = { 0, 1e-14};
  while(x0<10)
  {
    gsl_odeiv2_driver_apply (d, &x0, x0+0.2, y);
    if(x0>-4.3&&x0<4.3)
    {
      printf("%f %f\n",x0,y[0]);
    }
  }
  gsl_odeiv2_driver_free (d);
  return 0;
}
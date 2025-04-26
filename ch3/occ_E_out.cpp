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

double y_end_out (double E)
{
    gsl_odeiv2_system sys = {func, NULL, 2, &E};
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-7, 1e-7, 0);

  double x0 = -5;
  double y[2] = { 0, 1};
  double x_end=5;
  int s= gsl_odeiv2_driver_apply (d, &x0, x_end, y);
  gsl_odeiv2_driver_free (d);
  return fabs(y[0]);
}

int main()
{
    for(double E=0;E<=10;E+=0.1)
    {
        printf("%f %f\n",E,y_end_out(E));
    }
}
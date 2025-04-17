#include <iostream>
using namespace std;

#include <stdio.h> //如需使用printf()，需引入该头文件
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

double sc(double x,int m)
{
    if(m%2==0)
    {
        return sin(M_PI/8.0*m*x);
    }
    else
    {
        return cos(M_PI/8.0*m*x);
    }
}
double dsc(double x,int m)
{
    if(m%2==0)
    {
        return M_PI/8.0*m*cos(M_PI/8.0*m*x);
    }
    else
    {
        return -M_PI/8.0*m*sin(M_PI/8.0*m*x);
    }
}
double ff(double x,void *par)
{
    void **par_arr=(void **)par;
    int m_l=*((int *)(par_arr[0]));
    int m_r=*((int *)(par_arr[1]));
    return sc(x,m_l)*sc(x,m_r);
}
double fxf(double x,void *par)
{
    void **par_arr=(void **)par;
    int m_l=*((int *)(par_arr[0]));
    int m_r=*((int *)(par_arr[1]));
    return sc(x,m_l)*x*sc(x,m_r);
}
double fx2f(double x,void *par)
{
    void **par_arr=(void **)par;
    int m_l=*((int *)(par_arr[0]));
    int m_r=*((int *)(par_arr[1]));
    return sc(x,m_l)*x*x*sc(x,m_r);
}
double fdf(double x,void *par)
{
    void **par_arr=(void **)par;
    int m_l=*((int *)(par_arr[0]));
    int m_r=*((int *)(par_arr[1]));
    return sc(x,m_l)*dsc(x,m_r);
}
double fd2f(double x,void *par)
{
    void **par_arr=(void **)par;
    int m_l=*((int *)(par_arr[0]));
    int m_r=*((int *)(par_arr[1]));
    return -M_PI/8.0*m_r*M_PI/8.0*m_r*sc(x,m_l)*sc(x,m_r);
}

void mat_out(double **x,double **x2,double **d,double **d2,double **ove,int L)
{
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    void *par[4];
    for(int i=0;i<L;i++)
    {
        int m_l=i+1;
        par[0]=(void *)(&m_l);
        for(int j=0;j<L;j++)
        {
            int m_r=j+1;
            par[1]=(void *)(&m_r);
            gsl_function F;
            F.function = &ff;
            F.params = (void *)par;
            double res,err;
            gsl_integration_qags(&F, -4, 4, 1e-3, 1e-6, 1000, w, &res, &err);
            ove[i][j]=res;
            F.function = &fxf;
            gsl_integration_qags(&F, -4, 4, 1e-3, 1e-6, 1000, w, &res, &err);
            x[i][j]=res;
            F.function = &fx2f;
            gsl_integration_qags(&F, -4, 4, 1e-3, 1e-6, 1000, w, &res, &err);
            x2[i][j]=res;
            F.function = &fdf;
            gsl_integration_qags(&F, -4, 4, 1e-3, 1e-6, 1000, w, &res, &err);
            d[i][j]=res;
            F.function = &fd2f;
            gsl_integration_qags(&F, -4, 4, 1e-3, 1e-6, 1000, w, &res, &err);
            d2[i][j]=res;
        }
    }
    gsl_integration_workspace_free(w);
}

double flu_pro(const gsl_vector *x, void *par)
{
    void **par_arr=(void **)par;
    int L=*((int *)(par_arr[0]));//三角函数个数
    double **x_mat=(double **)(par_arr[1]);//传递<i|x|j>矩阵元
    double **x2_mat=(double **)(par_arr[2]);//传递<i|x^2|j>矩阵元
    double **d_mat=(double **)(par_arr[3]);//传递<i|d/dx|j>矩阵元
    double **d2_mat=(double **)(par_arr[4]);//传递<i|d^2/dx^2|j>矩阵元
    double **ove=(double **)(par_arr[5]);//传递<i|j>重叠

    double x_exp=0;
    double x2_exp=0;
    double d_exp=0;
    double d2_exp=0;
    double norm=0;
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            x_exp+=x_mat[i][j]*gsl_vector_get(x,i)*gsl_vector_get(x,j);
            x2_exp+=x2_mat[i][j]*gsl_vector_get(x,i)*gsl_vector_get(x,j);
            d_exp+=d_mat[i][j]*gsl_vector_get(x,i)*gsl_vector_get(x,j);
            d2_exp+=d2_mat[i][j]*gsl_vector_get(x,i)*gsl_vector_get(x,j);
        }
        norm+=ove[i][i]*gsl_vector_get(x,i)*gsl_vector_get(x,i);
    }
    double res=(x2_exp/norm-x_exp/norm*x_exp/norm)*(d_exp/norm*d_exp/norm-d2_exp/norm)+(norm-1)*(norm-1)*100;//(norm-1)*(norm-1)*100主要是归一化
    return res;
}
int main(void)
{
    int L=20;
    cin>>L;
    double **x_mat=new double *[L];
    double **x2_mat=new double *[L];
    double **d_mat=new double *[L];
    double **d2_mat=new double *[L];
    double **ove=new double *[L];
    for(int i=0;i<L;i++)
    {
        x_mat[i]=new double [L];
        x2_mat[i]=new double [L];
        d_mat[i]=new double [L];
        d2_mat[i]=new double [L];
        ove[i]=new double [L];
    }
    mat_out(x_mat,x2_mat,d_mat,d2_mat,ove,L);
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            printf("%f ",ove[i][j]);
        }
        printf("\n");
    }
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            printf("%f ",x_mat[i][j]);
        }
        printf("\n");
    }
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            printf("%f ",x2_mat[i][j]);
        }
        printf("\n");
    }
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            printf("%f ",d_mat[i][j]);
        }
        printf("\n");
    }
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            printf("%f ",d2_mat[i][j]);
        }
        printf("\n");
    }


  const gsl_multimin_fminimizer_type *T =
    gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  int status;
  double size;

  /* Starting point */
  x = gsl_vector_alloc (L);
//   gsl_vector_set_all (x, 2.0);
  /* Set initial step sizes to 1 */
  ss = gsl_vector_alloc (L);
  
    void *par[6];
    par[0]=&L;
    par[1]=(void *)x_mat;
    par[2]=(void *)x2_mat;
    par[3]=(void *)d_mat;
    par[4]=(void *)d2_mat;
    par[5]=(void *)ove;
  /* Initialize method and iterate */
  minex_func.n = L;
  minex_func.f = flu_pro;
  minex_func.params = par;

  s = gsl_multimin_fminimizer_alloc (T, L);
  
  const gsl_rng_type *rngT=gsl_rng_ranlxs0;
  gsl_rng *r= gsl_rng_alloc(rngT);;
  gsl_rng_default_seed = ((unsigned long)(time(NULL)));
    double min=100;
    double vec[L];
    for(int repeat=0;repeat<1000;repeat++)
    {
        for(int i=0;i<L;i++)
        {
            gsl_vector_set (x, i,gsl_ran_gaussian(r,0.1));
            gsl_vector_set (ss,i,gsl_ran_gaussian(r,0.01));
        }
        gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
        size_t iter=0;
        do
            {
            iter++;
            status = gsl_multimin_fminimizer_iterate(s);

            if (status)
                break;

            size = gsl_multimin_fminimizer_size (s);
            status = gsl_multimin_test_size (size, 1e-8);

            // if (status == GSL_SUCCESS)
            //     {
            //     printf ("converged to minimum at\n");
            //     }

            // printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n",
            //         iter,
            //         gsl_vector_get (s->x, 0),
            //         gsl_vector_get (s->x, 1),
            //         s->fval, size);
            }
        while (status == GSL_CONTINUE && iter < 10000);
        if(min>s->fval)
        {
            min=s->fval;
            for(int i=0;i<L;i++)
            {
                vec[i]=gsl_vector_get(s->x,i);
            }
        }
    }
    cout<<L<<' '<<min<<endl;;
    
    // for(int i=0;i<L;i++)
    // {
    //     for(int k=0;k<L;k++)
    //     {
    //         cout<<vec[i][k]<<' ';
    //     }
    //     cout<<endl;
    // }
    for(double x=-4;x<4.1;x+=0.2)
    {
        cout<<x<<' ';
        double y=0;
        for(int k=0;k<L;k++)
        {
            y+=sc(x,k+1)*vec[k];
        }
        cout<<y<<endl;
    }
    
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
    return 0;
}
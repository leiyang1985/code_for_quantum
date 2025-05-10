#include <iostream>
using namespace std;

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_hermite.h>

double phi_n(int n,double x)
{
    double res=gsl_sf_hermite_phys(n,x)/sqrt(M_SQRTPI);
    for(int i=1;i<=n;i++)
    {
        res/=sqrt(2.0*i);
    }
    res*=exp(-x*x/2.0);
    return res;
}

double fket(double x,void *par)
{
    void **par_arr=(void **)par;
    int n=*((int *)par_arr[0]);
    double x0=*((double *)par_arr[1]);
    double sigma=*((double *)par_arr[2]);
    double res=gsl_sf_hermite_phys(n,x)/sqrt(M_SQRTPI);
    for(int i=1;i<=n;i++)
    {
        res/=sqrt(2.0*i);
    }
    res*=exp(-x*x/2.0);
    res*=1/sqrt(sqrt(2))/sqrt(M_SQRTPI)/sqrt(sigma)*exp(-(x-x0)*(x-x0)/4.0/sigma/sigma);
    return res;
}

void cn_out(double x0,double sigma,double *cn, int L)
{
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    void *par[3];
    par[1] = (void *)(&x0);
    par[2] = (void *)(&sigma);
    for(int i=0;i<L;i++)
    {
        par[0] = (void *)(&i);
        gsl_function F;
        F.function = &fket;
        F.params = (void *)par;
        double err;
        gsl_integration_qagi (&F, 1e-3, 1e-7, 1000, w, cn+i, &err);
    }
    gsl_integration_workspace_free(w);
}
double braket(double x,void *par)
{
    void **par_arr=(void **)par;
    int m=*((int *)par_arr[0]);
    int n=*((int *)par_arr[1]);
    double bra=gsl_sf_hermite_phys(m,x);
    for(int i=1;i<=m;i++)
    {
        bra/=sqrt(2.0*i);
    }
    double ket=gsl_sf_hermite_phys(n,x);
    for(int i=1;i<=n;i++)
    {
        ket/=sqrt(2.0*i);
    }
    double res=bra*ket/M_SQRTPI*exp(-x*x);
    return res;
}
double braxket(double x,void *par)
{
    void **par_arr=(void **)par;
    int m=*((int *)par_arr[0]);
    int n=*((int *)par_arr[1]);
    double bra=gsl_sf_hermite_phys(m,x);
    for(int i=1;i<=m;i++)
    {
        bra/=sqrt(2.0*i);
    }
    double ket=gsl_sf_hermite_phys(n,x);
    for(int i=1;i<=n;i++)
    {
        ket/=sqrt(2.0*i);
    }
    double res=bra*ket*x/M_SQRTPI*exp(-x*x);
    return res;
}
double brax2ket(double x,void *par)
{
    void **par_arr=(void **)par;
    int m=*((int *)par_arr[0]);
    int n=*((int *)par_arr[1]);
    double bra=gsl_sf_hermite_phys(m,x);
    for(int i=1;i<=m;i++)
    {
        bra/=sqrt(2.0*i);
    }
    double ket=gsl_sf_hermite_phys(n,x);
    for(int i=1;i<=n;i++)
    {
        ket/=sqrt(2.0*i);
    }
    double res=bra*ket*x*x/M_SQRTPI*exp(-x*x);
    return res;
}
void mat_out(double **norm,double **x,double **x2,int L)
{
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    void *par[2];
    for(int i=0;i<L;i++)
    {
        par[0] = (void *)(&i);
        for(int j=0;j<L;j++)
        {
            par[1]=(void *)(&j);
            gsl_function F;
            F.function = &braket;
            F.params = (void *)par;
            double err;
            gsl_integration_qagi (&F, 1e-3, 1e-7, 1000, w, &(norm[i][j]), &err);
            F.function = &braxket;
            gsl_integration_qagi (&F, 1e-3, 1e-7, 1000, w, &(x[i][j]), &err);
            F.function = &brax2ket;
            gsl_integration_qagi (&F, 1e-3, 1e-7, 1000, w, &(x2[i][j]), &err);
        }
    }
    gsl_integration_workspace_free(w);
}
/***************************|\phi(x,t)|^2 */
double phi(double *cn,double x,double t,int L)
{
    double res=0;
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            double delta=cn[i]*cn[j]*cos((i-j)*t)*phi_n(i,x)*phi_n(j,x);
            res+=delta;
        }
    }
    return res;
}

/*****************f(x)实函数期望：<f(x)> */
double mean(double *cn,double **f,double t,int L)
{
    double res=0;
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            res+=cn[i]*cn[j]*cos((i-j)*t)*f[i][j];
        }
    }
    return res;
}
int main()
{
    int L;
    cout<<"s.p. state number:";
    cin >> L;//输入单粒子维数
    double x0=1;
    cout<<"stating position:";
    cin >> x0;//输入单粒子维数
    double sigma=1/sqrt(2);
    // cout<<"sigma of Gaussion wave:";
    // cin >> sigma;//输入单粒子维数
    double cn[L];
    cn_out(x0,sigma,cn,L);
    cout<<"cn="<<endl;
    double cn2sum=0;
    for(int i=0;i<L;i++)
    {
        cout<<i<<' '<<cn[i]<<endl;
        cn2sum+=cn[i]*cn[i];
    }
    cout<<"cn2sum="<<cn2sum<<endl;
    for(int i=0;i<L;i++)
    {
        // cn[i]/=sqrt(cn2sum);
    }
    double **norm=new double *[L];
    double **x=new double *[L];
    double **x2=new double *[L];
    for(int i=0;i<L;i++)
    {
        norm[i]=new double [L];
        x[i]=new double [L];
        x2[i]=new double [L];
    }
    mat_out(norm,x,x2,L);
    cout<<"norm[][]="<<endl;
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            cout<<norm[i][j]<<' ';
        }
        cout<<endl;
    }
    cout<<"x[][]="<<endl;
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            cout<<x[i][j]<<' ';
        }
        cout<<endl;
    }
    cout<<"x2[][]="<<endl;
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            cout<<x2[i][j]<<' ';
        }
        cout<<endl;
    }
    cout<<"t <x>(t), Delta x(t)"<<endl;
    for(double t=0;t<=M_PI*2;t+=M_PI*2/20)
    {
        double phi_norm=mean(cn,norm,t,L);
        double x_mean=mean(cn,x,t,L);
        double x2_mean=mean(cn,x2,t,L);
        cout<<t<<' '<<phi_norm<<' '<<x_mean/phi_norm<<' '<<sqrt(x2_mean/phi_norm-x_mean/phi_norm*x_mean/phi_norm)<<endl;
    }
    cout<<"phi^2="<<endl;
    for(double x_cor=-4;x_cor<=4;x_cor+=0.2)
    {
        cout<<x_cor<<' ';
        for(double t=0;t<=M_PI*2;t+=M_PI*2/8)
        {
            cout<<phi(cn,x_cor,t,L)<<' ';
        }
        cout<<endl;
    }
    return 0;
}
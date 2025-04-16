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
//#include <gsl/gsl_randist.h>


inline void time2rand()
{
	srand((unsigned int)time(0));
}

//	"Polar" version without trigonometric calls
double average(double begin,double end)
{
	double res=((double)(rand()))/double(RAND_MAX)*(end-begin)+begin;
	return(res);
}

double gauss(double mu, double sigma)
{
	static bool deviateAvailable=false;	//	flag
	static float storedDeviate;			//	deviate from previous calculation
	//	If no deviate has been stored, the polar Box-Muller transformation is
	//	performed, producing two independent normally-distributed random
	//	deviates.  One is stored for the next round, and one is returned.
	if (!deviateAvailable)
	{
		double var1,var2,rsquared;

		//	choose pairs of uniformly distributed deviates, discarding those
		//	that don't fall within the unit circle
		do
		{
			var1=2.0*(double)(rand())/double(RAND_MAX) - 1.0;
			var2=2.0*(double)(rand())/double(RAND_MAX) - 1.0;
			rsquared=var1*var1+var2*var2;
		} while ( rsquared>=1.0);

		//	calculate polar tranformation for each deviate
		double polar=sqrt(-2.0*log(rsquared)/rsquared);

		//	store first deviate and set flag
		storedDeviate=var1*polar;
		deviateAvailable=true;

		//	return second deviate
		return var2*polar*sigma + mu;
	}

	//	If a deviate is available from a previous call to this function, it is
	//	returned, and the flag is set to false.
	else
	{
		deviateAvailable=false;
		return storedDeviate*sigma + mu;
	}
}

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
double ff(double x,void *par)
{
    void **par_arr=(void **)par;
    int m_l=*((int *)(par_arr[0]));
    int m_r=*((int *)(par_arr[1]));
    return sc(x,m_l)*sc(x,m_r);
}
double Hf(double x,void *par)
{
    void **par_arr=(void **)par;
    int m=*((int *)(par_arr[0]));
    return ((M_PI/8.0*m)*(M_PI/8.0*m)+x*x)*sc(x,m);
}
double dHf(double x,void *par)
{
    gsl_function F;
    F.function=Hf;
    F.params=par;
    double res,err;
    gsl_deriv_central(&F,x,1e-3,&res,&err);
    // cout<<"df("<<x<<")= "<<res<<endl;
    // printf("df(%15.16f)=%15.16f\n",x,res);
    return res;
}
double d2Hf(double x,void *par)
{
    gsl_function F;
    F.function=dHf;
    F.params=par;
    double res,err;
    gsl_deriv_central(&F,x,1e-3,&res,&err);
    // cout<<"d2f("<<x<<")= "<<res<<endl;
    return res;
}
double fHf(double x,void *par)
{
    void **par_arr=(void **)par;
    int m_l=*((int *)(par_arr[0]));
    int m_r=*((int *)(par_arr[1]));
    return sc(x,m_l)*((M_PI/8.0*m_r)*(M_PI/8.0*m_r)+x*x)*sc(x,m_r);
}
double fH2f(double x,void *par)
{
    void **par_arr=(void **)par;
    int m_l=*((int *)(par_arr[0]));
    int m_r=*((int *)(par_arr[1]));
    void *par_r[3];
    par_r[0]=(void *)(&m_r);
    double d2f_r=d2Hf(x,par_r);
    return sc(x,m_l)*(-d2f_r+x*x*Hf(x,par_r));
}
void mat_out(double **H,double **H2,double **ove,int L)
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
            F.function = &fHf;
            gsl_integration_qags(&F, -4, 4, 1e-3, 1e-6, 1000, w, &res, &err);
            H[i][j]=res;
            F.function = &fH2f;
            gsl_integration_qags(&F, -4, 4, 1e-3, 1e-6, 1000, w, &res, &err);
            H2[i][j]=res;
        }
    }
    gsl_integration_workspace_free(w);
}

double H_expect_oth(const gsl_vector *x, void *par)
{
    void **par_arr=(void **)par;
    int L=*((int *)(par_arr[0]));
    int index=*((int *)(par_arr[1]));
    double **vec_pre=(double **)(par_arr[2]);
    double **ove=(double **)(par_arr[3]);
    double **H=(double **)(par_arr[4]);

    double H_sum=0;
    double norm=0;
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            H_sum+=H[i][j]*gsl_vector_get(x,i)*gsl_vector_get(x,j);
        }
        norm+=ove[i][i]*gsl_vector_get(x,i)*gsl_vector_get(x,i);
    }
    double res=H_sum/norm+(norm-1)*(norm-1)*100;
    for(int i=0;i<index;i++)
    {
        double oth=0;
        for(int j=0;j<L;j++)
        {
            oth+=ove[j][j]*vec_pre[i][j]*gsl_vector_get(x,j);
        }
        res+=oth*oth*100;
    }
    return res;
}
double deltaH_oth(const gsl_vector *x, void *par)
{
    void **par_arr=(void **)par;
    int L=*((int *)(par_arr[0]));
    int index=*((int *)(par_arr[1]));
    double **vec_pre=(double **)(par_arr[2]);
    double **ove=(double **)(par_arr[3]);
    double **H=(double **)(par_arr[4]);
    double **H2=(double **)(par_arr[5]);

    double H_sum=0;
    double H2_sum=0;
    double norm=0;
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            H_sum+=H[i][j]*gsl_vector_get(x,i)*gsl_vector_get(x,j);
            H2_sum+=H2[i][j]*gsl_vector_get(x,i)*gsl_vector_get(x,j);
        }
        norm+=ove[i][i]*gsl_vector_get(x,i)*gsl_vector_get(x,i);
    }
    double res=H2_sum/norm-H_sum/norm*H_sum/norm;
    for(int i=0;i<index;i++)
    {
        double oth=0;
        for(int j=0;j<L;j++)
        {
            oth+=ove[j][j]*vec_pre[i][j]*gsl_vector_get(x,j);
        }
        res+=oth*oth*1000000;
    }
    return res;
}
double deltaH(const gsl_vector *c, void *par)
{
    void **par_arr=(void **)par;
    int L=*((int *)(par_arr[0]));
    int index=*((int *)(par_arr[1]));
    double **ove=(double **)(par_arr[3]);
    double **H=(double **)(par_arr[4]);
    double **H2=(double **)(par_arr[5]);

    double H_sum=0;
    double H2_sum=0;
    double norm=0;
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            H_sum+=H[i][j]*gsl_vector_get(c,i)*gsl_vector_get(c,j);
            H2_sum+=H2[i][j]*gsl_vector_get(c,i)*gsl_vector_get(c,j);
        }
        norm+=ove[i][i]*gsl_vector_get(c,i)*gsl_vector_get(c,i);
    }
    double res=H2_sum/norm-H_sum/norm*H_sum/norm;
    return res;
}
double H_expect(const gsl_vector *x, void *par)
{
    void **par_arr=(void **)par;
    int L=*((int *)(par_arr[0]));
    int index=*((int *)(par_arr[1]));
    double **vec_pre=(double **)(par_arr[2]);
    double **ove=(double **)(par_arr[3]);
    double **H=(double **)(par_arr[4]);

    double H_sum=0;
    double norm=0;
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            H_sum+=H[i][j]*gsl_vector_get(x,i)*gsl_vector_get(x,j);
        }
        norm+=ove[i][i]*gsl_vector_get(x,i)*gsl_vector_get(x,i);
    }
    double res=H_sum/norm;
    return res;
}
int main(void)
{
    time2rand();
    int L=20;
    cin>>L;
    double **H=new double *[L];
    double **H2=new double *[L];
    double **ove=new double *[L];
    double **vec=new double *[L];
    for(int i=0;i<L;i++)
    {
        H[i]=new double [L];
        H2[i]=new double [L];
        ove[i]=new double [L];
        vec[i]=new double [L];
    }
    mat_out(H,H2,ove,L);
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
            printf("%f ",H[i][j]);
        }
        printf("\n");
    }
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            printf("%f ",H2[i][j]);
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
    par[2]=(void *)vec;
    par[3]=(void *)ove;
    par[4]=(void *)H;
    par[5]=(void *)H2;
  /* Initialize method and iterate */
  minex_func.n = L;
  minex_func.f = H_expect_oth;
  minex_func.params = par;

  s = gsl_multimin_fminimizer_alloc (T, L);
  double E[L];
  double deltaE[L];
  int count=0;
  const gsl_rng_type *rngT=gsl_rng_ranlxs0;
  gsl_rng *r= gsl_rng_alloc(rngT);;
  gsl_rng_default_seed = ((unsigned long)(time(NULL)));
  int count_for_simplex=0;
    for(int index=0;index<L;index++)
    {
        par[1]=&index;
        double min=100;
        for(int repeat=0;repeat<1000;repeat++)
        {
            for(int i=0;i<L;i++)
            {
                gsl_vector_set (x, i,gauss(0,0.1));
                gsl_vector_set (ss,i,gauss(0,0.01));
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
                for(int i=0;i<L;i++)
                {
                    vec[index][i]=gsl_vector_get(s->x,i);
                }
                E[index]=H_expect(s->x,par);
            }
        }
        
        
    }
    cout<<L<<' ';
    for(int i=0;i<L;i++)
    {
        cout<<E[i]<<' ';
    }
    cout<<endl;
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
        for(int i=0;i<L;i++)
        {
            double y=0;
            for(int k=0;k<L;k++)
            {
                y+=sc(x,k+1)*vec[i][k];
            }
            cout<<y<<' ';
        }
        cout<<endl;
    }
    
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
    return 0;
}
#include <iostream>
using namespace std;

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_hermite.h>
double factor(int n)
{
    if(n==0)
    {
        return 1;
    }
    else
    {
        return n*factor(n-1);
    }
}
/*****************厄米多项式的递推 */
double Hp(double x,int n)
{
    if(n==0)
    {
        return 1;
    }
    else if(n==1)
    {
        return 2*x;
    }
    else
    {
        return 2*x*Hp(x,n-1)-2*(n-1)*Hp(x,n-2);
    }
}
/**********************厄米多项式的一阶导 */
double dHp(double x,int n)
{
    
    if(n==0)
    {
        return 0;
    }
    // else if(n==1)
    // {
    //     return 2;
    // }
    else
    {
        return 2*n*Hp(x,n-1);
        // return 2*Hp(x,n-1)+2*x*dHp(x,n-1)-2*(n-1)*dHp(x,n-2);
    }
}
double d2Hp(double x,int n)
{
    if(n==0)
    {
        return 0;
    }
    else if(n==1)
    {
        return 0;
    }
    else
    {
        return 4*n*(n-1)*Hp(x,n-2);
    }
}
/********未归一化的谐振子本征波函数 */
double occ(double x, int m)
{
    return exp(-0.5*x*x)*Hp(x,m);
}
/********未归一化的谐振子本征波函数的导数 */
double docc(double x,int m)
{
    return -x*exp(-0.5*x*x)*Hp(x,m)+exp(-0.5*x*x)*dHp(x,m);
}
double d2occ(double x,int m)
{
    return -exp(-0.5*x*x)*Hp(x,m)+x*x*exp(-0.5*x*x)*Hp(x,m)-x*exp(-0.5*x*x)*dHp(x,m)-x*exp(-0.5*x*x)*dHp(x,m)+exp(-0.5*x*x)*d2Hp(x,m);
}

double sc(double x, int m)
{
    if (m % 2 == 0)
    {
        return sin(M_PI / 8.0 * m * x);
    }
    else
    {
        return cos(M_PI / 8.0 * m * x);
    }
}
double d2sc(double x, int m)
{
    if (m % 2 == 0)
    {
        return -M_PI / 8.0 * m*M_PI / 8.0 * m*sin(M_PI / 8.0 * m * x);
    }
    else
    {
        return -M_PI / 8.0 * m*M_PI / 8.0 * m*cos(M_PI / 8.0 * m * x);
    }
}
double ff(double x, void *par)
{
    void **par_arr = (void **)par;
    int m_l = *((int *)(par_arr[0]));
    int m_r = *((int *)(par_arr[1]));
    return sc(x, m_l) * sc(x, m_r);
}

double fp2f(double x, void *par)
{
    void **par_arr = (void **)par;
    int m_l = *((int *)(par_arr[0]));
    int m_r = *((int *)(par_arr[1]));
    return -sc(x, m_l) * d2sc(x, m_r);
}
void mat_out(double **H,double **p, double **ove, int L)
{
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    void *par[4];
    for (int i = 0; i < L; i++)
    {
        int m_l = i+1;
        par[0] = (void *)(&m_l);
        for (int j = 0; j < L; j++)
        {
            int m_r = j+1;
            par[1] = (void *)(&m_r);
            gsl_function F;
            F.function = &ff;
            F.params = (void *)par;
            double res, err;
            gsl_integration_qags(&F, -4, 4, 1e-3, 1e-7, 1000, w, &res, &err);
            ove[i][j] = res;
            F.function = &fp2f;
            gsl_integration_qags(&F, -4, 4, 1e-3, 1e-7, 1000, w, &res, &err);
            H[i][j] = res;
        }
    }
    gsl_integration_workspace_free(w);
}

int main(void)
{
    int L = 20;
    cin >> L;//输入维数
    double **H = new double *[L];
    double **p=new double *[L];
    double **ove = new double *[L];
    double **vec = new double *[L];
    for (int i = 0; i < L; i++)
    {
        H[i] = new double[L];
        p[i]=new double [L];
        ove[i] = new double[L];
        vec[i] = new double[L];
    }
    mat_out(H, p,ove, L);
    cout<<"ove="<<endl;
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            printf("%f ", ove[i][j]);
        }
        printf("\n");
    }
    cout<<"normed ove="<<endl;
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            printf("%f ", ove[i][j]/sqrt(ove[i][i]*ove[j][j]));
        }
        printf("\n");
    }
    cout<<"p ele normed ="<<endl;
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            printf("%f ", p[i][j]/sqrt(ove[i][i]*ove[j][j]));
            // if(i==j-1)
            // {
            //     p[i][j]=-sqrt((i+1)*0.5);
            // }
            // else if(i==j+1)
            // {
            //     p[i][j]=sqrt((j+1)*0.5);
            // }
            // else
            // {
            //     p[i][j]=0;
            // }
        }
        printf("\n");
    }
    /********************H矩阵的归一化**************** */
    cout<<"H normed= "<<endl;
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            H[i][j] /= sqrt(ove[i][i] * ove[j][j]);
            // H[i][j]=0;
            // for(int k=0;k<L;k++)
            // {
            //     H[i][j]+=p[i][k]*p[k][j];
            // }
            // H[i][j]*=-1;
        }
    }
    cout<<"H="<<endl;
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            cout<<H[i][j]<<' ';
        }
        cout<<endl;
    }
    // for(int i=0;i<L;i++)
    // {
    //     for(int j=0;j<L;j++)
    //     {
    //         double H_ele=0;
    //         if(i==j+2)
    //         {
    //             H_ele=-sqrt((i+2)/2.0)*sqrt((i+1)/2.0);
    //         }
    //         else if(i==j)
    //         {
    //             H_ele=(2*i+1)/2.0;
    //         }
    //         else if(i==j-2)
    //         {
    //             H_ele=-sqrt((i-1)/2.0)*sqrt(i/2.0);
    //         }
    //         cout<<H[i][j]<<' ';
    //     }
    //     cout<<endl;
    // }
    /************************对角化的数据存储准备*/
    gsl_matrix *m=gsl_matrix_alloc(L,L);//对角化的对称矩阵
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            gsl_matrix_set(m,i,j,H[i][j]);//哈密顿矩阵元赋于对角化的对称矩阵
        }
    }
    gsl_vector *eval = gsl_vector_alloc(L);//本征值存储的数组
    gsl_matrix *evec = gsl_matrix_alloc(L, L);//存储本征态的二维数组，对也变换\Lambda矩阵

    gsl_eigen_symmv_workspace *w =gsl_eigen_symmv_alloc(L);
    gsl_eigen_symmv(m, eval, evec, w);
    gsl_eigen_symmv_free(w);
    gsl_eigen_symmv_sort(eval, evec,GSL_EIGEN_SORT_VAL_ASC);//GSL_EIGEN_SORT_VAL_ASC排本征值从小到大排序。

    for (int i = 0; i < L; i++)
    {
        double eval_i = gsl_vector_get(eval, i);
        gsl_vector_view evec_i = gsl_matrix_column(evec, i);

        printf("eigenvalue = %g c= %f\n", eval_i,sqrt(eval_i/.15421256880729318306));
        // printf("eigenvector = \n");
        // gsl_vector_fprintf(stdout,&evec_i.vector, "%g");
    }
    // gsl_sf_hermite_phys(10,4.0);
    // cout<<"in x,y plot:"<<endl;
    // for(double x=-4;x<4.1;x+=0.2)
    // {
    //     cout<<x<<' ';
    //     for (int i = 0; i < L; i++)
    //     {
    //         double eval_i = gsl_vector_get(eval, i);
    //         gsl_vector_view evec_i = gsl_matrix_column(evec, i);
    //         double y=0;
    //         for(int k=0;k<L;k++)
    //         {
    //             y+=sc(x,k+1)*gsl_vector_get(&(evec_i.vector),k);
    //         }
    //         cout<<y<<' ';
    //     }
    //     cout<<endl;
    // }
        

    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    return 0;
}
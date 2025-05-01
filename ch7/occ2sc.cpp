#include <iostream>
using namespace std;

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_hermite.h>

/*******************重叠的被积函数 */
double ff(double x, void *par)
{
    void **par_arr = (void **)par;
    int m_l = *((int *)(par_arr[0]));
    m_l=2*m_l+1;
    double zero_l = *((double *)(par_arr[1]));
    int m_r = *((int *)(par_arr[2]));
    m_r=2*m_r+1;
    double zero_r = *((double *)(par_arr[3]));
    return (gsl_sf_hermite_phys(m_l,x)*exp(-0.5*x*x)-x*zero_l)*(gsl_sf_hermite_phys(m_r,x)*exp(-0.5*x*x)-x*zero_r);
}
/************p^2算符积分函数 */
double fp2f(double x, void *par)
{
    void **par_arr = (void **)par;
    int m_l = *((int *)(par_arr[0]));
    m_l=2*m_l+1;
    double zero_l = *((double *)(par_arr[1]));
    int m_r = *((int *)(par_arr[2]));
    m_r=2*m_r+1;
    double zero_r = *((double *)(par_arr[3]));
    if(m_r==1)
    {
        return -(gsl_sf_hermite_phys(m_l,x)*exp(-0.5*x*x)-x*zero_l)*exp(-0.5*x*x)*2*x*(x*x-3);
    }
    else
    {
        return -(gsl_sf_hermite_phys(m_l,x)*exp(-0.5*x*x)-x*zero_l)*exp(-0.5*x*x)*((x*x-1)*gsl_sf_hermite_phys(m_r,x)-x*4*m_r*gsl_sf_hermite_phys(m_r-1,x)+4*m_r*(m_r-1)*gsl_sf_hermite_phys(m_r-2,x));
    }
}
void mat_out(double **H, double *her_zero,double **ove, int L)
{
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    void *par[4];
    for (int i = 0; i < L; i++)
    {
        int m_l = i;
        par[0] = (void *)(&m_l);
        par[1] = (void *)(her_zero+i);
        for (int j = 0; j < L; j++)
        {
            int m_r = j;
            par[2] = (void *)(&m_r);
            par[3] = (void *)(her_zero+j);
            gsl_function F;
            F.function = &ff;
            F.params = (void *)par;
            double res, err;
            gsl_integration_qags(&F, 0, 4, 1e-3, 1e-7, 1000, w, &res, &err);
            ove[i][j] = res;
            F.function = &fp2f;
            gsl_integration_qags(&F, 0, 4, 1e-3, 1e-7, 1000, w, &res, &err);
            H[i][j] = res;
        }
    }
    gsl_integration_workspace_free(w);
}

int main(void)
{
    int L = 20;
    cin >> L;//输入维数
    double *her_zero=new double [L];
    for(int i=0;i<L;i++)
    {
        her_zero[i]=gsl_sf_hermite_phys(2*i+1,4)*exp(-0.5*4*4)/4;
    }
    double **H = new double *[L];
    double **ove = new double *[L];
    double **vec = new double *[L];
    for (int i = 0; i < L; i++)
    {
        H[i] = new double[L];
        ove[i] = new double[L];
        vec[i] = new double[L];
    }
    mat_out(H, her_zero,ove, L);
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
    /********************归一化后H矩阵输出（注意非正交化）**************** */
    cout<<"normed H"<<endl;
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            printf("%f ", H[i][j]/sqrt(ove[i][i]*ove[j][j]));
        }
        printf("\n");
    }

    /************************对角化的数据存储准备*/
    gsl_matrix *A=gsl_matrix_alloc(L,L);//对角化的A矩阵，算符矩阵
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            gsl_matrix_set(A,i,j,H[i][j]);//哈密顿矩阵元赋于对角化的对称矩阵
        }
    }
    gsl_matrix *B=gsl_matrix_alloc(L,L);//对角化所用的B对称矩阵，重叠矩阵
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            gsl_matrix_set(B,i,j,ove[i][j]);//重叠矩阵元赋于ove对称矩阵
        }
    }
    gsl_vector *eval = gsl_vector_alloc(L);//本征值存储的数组
    gsl_matrix *evec = gsl_matrix_alloc(L, L);//存储本征态的二维数组，对也变换\Lambda矩阵
    gsl_eigen_gensymmv_workspace *w =gsl_eigen_gensymmv_alloc(L);
    gsl_eigen_gensymmv(A, B,eval, evec, w);
    gsl_eigen_gensymmv_free(w);
    gsl_eigen_gensymmv_sort(eval, evec,GSL_EIGEN_SORT_VAL_ASC);//GSL_EIGEN_SORT_VAL_ASC排本征值从小到大排序。

    for (int i = 0; i < L; i++)
    {
        double eval_i = gsl_vector_get(eval, i);
        gsl_vector_view evec_i = gsl_matrix_column(evec, i);

        printf("eigenvalue = %g c= %f\n", eval_i,sqrt(eval_i/.61685025402342250000));
        // printf("eigenvector = \n");
        // gsl_vector_fprintf(stdout,&evec_i.vector, "%g");
    }
    cout<<"in x,y plot:"<<endl;
    for(double x=0;x<4.1;x+=0.2)
    {
        cout<<x<<' ';
        for (int i = 0; i < L; i++)
        {
            double eval_i = gsl_vector_get(eval, i);
            gsl_vector_view evec_i = gsl_matrix_column(evec, i);
            double y=0;
            for(int k=0;k<L;k++)
            {
                y+=(gsl_sf_hermite_phys(2*k+1,x)*exp(-0.5*x*x)-x*her_zero[k])*gsl_vector_get(&(evec_i.vector),k);
            }
            cout<<y<<' ';
        }
        cout<<endl;
    }
    

    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    return 0;
}
#include <iostream>
using namespace std;

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

/**********哈密顿量矩阵元<ij|H|kl>的产生 */
void mat_out(double **H, int L)
{
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            for(int k=0;k<L;k++)
            {
                for(int l=0;l<L;l++)
                {
                    H[i*L+j][k*L+l]=0;
                    if(i==k&&j==l)
                    {
                        H[i*L+j][k*L+l]+=k+l+1;
                    }
                    if(i==k+2&&j==l)
                    {
                        H[i*L+j][k*L+l]+=0.25*sqrt((k+1)*(k+2));
                    }
                    if(i==k-2&&j==l)
                    {
                        H[i*L+j][k*L+l]+=0.25*sqrt(k*(k-1));
                    }
                    if(i==k&&j==l)
                    {
                        H[i*L+j][k*L+l]+=0.5*k;
                    }
                    H[i*L+j][k*L+l]+=0.25;
                    if(i==k&&j==l+2)
                    {
                        H[i*L+j][k*L+l]+=0.25*sqrt((l+1)*(l+2));
                    }
                    if(i==k&&j==l-2)
                    {
                        H[i*L+j][k*L+l]+=0.25*sqrt(l*(l-1));
                    }
                    if(i==k&&j==l)
                    {
                        H[i*L+j][k*L+l]+=0.5*l;
                    }
                    H[i*L+j][k*L+l]+=0.25;
                    if(i==k+1&&j==l+1)
                    {
                        H[i*L+j][k*L+l]-=0.5*sqrt((l+1)*(k+1));
                    }
                    if(i==k-1&&j==l-1)
                    {
                        H[i*L+j][k*L+l]-=0.5*sqrt(l*k);
                    }
                    if(i==k+1&&j==l-1)
                    {
                        H[i*L+j][k*L+l]-=0.5*sqrt((k+1)*l);
                    }
                    if(i==k-1&&j==l+1)
                    {
                        H[i*L+j][k*L+l]-=0.5*sqrt(k*(l+1));
                    }
                }
            }
        }
    }
}

int main(void)
{
    int L = 10;
    cout<<"s.p. state number:";
    cin >> L;//输入单粒子维数
    double **H = new double *[L*L];
    for (int i = 0; i < L*L; i++)
    {
        H[i] = new double[L*L];
    }
    mat_out(H, L);
    // for (int i = 0; i < L*L; i++)
    // {
    //     for (int j = 0; j < L*L; j++)
    //     {
    //         printf("%f ", H[i][j]);
    //     }
    //     printf("\n");
    // }
    /************************对角化的数据存储准备*/
    gsl_matrix *m=gsl_matrix_alloc(L*L,L*L);//对角化的对称矩阵
    for (int i = 0; i < L*L; i++)
    {
        for (int j = 0; j < L*L; j++)
        {
            gsl_matrix_set(m,i,j,H[i][j]);//哈密顿矩阵元赋于对角化的对称矩阵
        }
    }
    gsl_vector *eval = gsl_vector_alloc(L*L);//本征值存储的数组
    gsl_matrix *evec = gsl_matrix_alloc(L*L, L*L);//存储本征态的二维数组，对也变换\Lambda矩阵

    gsl_eigen_symmv_workspace *w =gsl_eigen_symmv_alloc(L*L);
    gsl_eigen_symmv(m, eval, evec, w);
    gsl_eigen_symmv_free(w);
    gsl_eigen_symmv_sort(eval, evec,GSL_EIGEN_SORT_VAL_ASC);//GSL_EIGEN_SORT_VAL_ASC排本征值从小到大排序。

    int print_max=10;
    if(L*L<print_max)
    {
        print_max=L*L;
    }
    for (int i = 0; i < print_max; i++)
    {
        double eval_i = gsl_vector_get(eval, i);
        gsl_vector_view evec_i = gsl_matrix_column(evec, i);

        /****************交换粒子而后计算重叠 */
        double rev_vec[L*L];
        for(int k=0;k<L;k++)
        {
            for(int l=0;l<L;l++)
            {
                rev_vec[l*L+k]=gsl_vector_get(&(evec_i.vector),k*L+l);
            }
        }
        double sum=0;
        for(int k=0;k<L;k++)
        {
            for(int l=0;l<L;l++)
            {
                sum+=gsl_vector_get(&(evec_i.vector),k*L+l)*rev_vec[k*L+l];
            }
        }

        printf("eigenvalue = %g rev_norm= %g \n", eval_i,sum);
        // printf("eigenvector = \n");
        // for(int k=0;k<L;k++)
        // {
        //     for(int l=0;l<L;l++)
        //     {
        //         cout<<gsl_vector_get(&(evec_i.vector),k*L+l)<<"|"<<k<<','<<l<<">"<<endl;
        //     }
        // }
    }        

    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    return 0;
}
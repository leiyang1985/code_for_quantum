#include <iostream>
#include <cmath>
using namespace std;

#include "gsl/gsl_sf_coupling.h"

double cg_from_3j(int j1,int j2,int m1,int m2,int j)
{
    double res=gsl_sf_coupling_3j(j1, j2,j,m1,m2,-m1-m2);
    res*=sqrt(j+1);
    if(abs(j1-j2+m1+m2)/2%2!=0)
    {
        res*=-1;
    }
    return res;
}

int main()
{
    /*********************用户明确j1 j2,注意是两倍整数****************************** */
    int j1, j2;
    cin>>j1;
    cin>>j2;
    /**************************CG系数存储三维数组*************************** */
    double **cg_store[(j1+j2-abs(j1-j2))/2+1];
    for(int j=abs(j1-j2);j<=j1+j2;j+=2)
    {
        cg_store[(j-abs(j1-j2))/2]=new double *[j+1];
        for(int m=-j;m<=j;m+=2)
        {
            cg_store[(j-abs(j1-j2))/2][(m+j)/2]=new double [j1+1];
            for(int m1=-j1;m1<=j1;m1+=2)
            {
                cg_store[(j-abs(j1-j2))/2][(m+j)/2][(m1+j1)/2]=314159;
            }
        }
    }
    for(int j=j1+j2;j>=abs(j1-j2);j-=2)
    {
        for(int m=j;m>=-j;m-=2)
        {
            if(j==j1+j2&&m==j1+j2)
            {
                int m1=j1;
                int m2=j2;
                cg_store[(j-abs(j1-j2))/2][(m+j)/2][j1]=1;
                cout<<"<"<<j1<<"/2,"<<j2<<"/2;"<<m1<<"/2,"<<m2<<"/2|"<<j<<"/2,"<<m<<"/2>="<<cg_store[(j-abs(j1-j2))/2][(m+j)/2][j1]<<"; <gsl>="<<cg_from_3j(j1,j2,m1,m2,j)<<";delta="<<cg_store[(j-abs(j1-j2))/2][(m+j)/2][j1]-cg_from_3j(j1,j2,m1,m2,j)<<endl;
            }
            else if(j!=j1+j2&&m==j)
            {
                /****确定m1范围 ******************* */
                int m1_max=j1;
                int m1_min=-j1;
                if(m1_min<m-j2)
                {
                    m1_min=m-j2;
                }
                if(m1_max>m+j2)
                {
                    m1_max=m+j2;
                }
                /*****************预设CG系数初始值************ */
                for(int m1=m1_min;m1<=m1_max;m1+=2)
                {
                    cg_store[(j-abs(j1-j2))/2][(m+j)/2][(m1+j1)/2]=m1+j1;
                }
                /***************根据Schmit正交法，实施正交归一化*********** */
                double ove_with_j_other[(j1+j2-j)/2];
                for(int j_other=j1+j2;j_other>j;j_other-=2)
                {
                    ove_with_j_other[(j1+j2-j_other)/2]=0;
                    for(int m1=m1_min;m1<=m1_max;m1+=2)
                    {
                        ove_with_j_other[(j1+j2-j_other)/2]+=cg_store[(j-abs(j1-j2))/2][(m+j)/2][(m1+j1)/2]*cg_store[(j_other-abs(j1-j2))/2][(m+j_other)/2][(m1+j1)/2];
                    }
                }
                for(int j_other=j1+j2;j_other>j;j_other-=2)
                {
                    for(int m1=m1_min;m1<=m1_max;m1+=2)
                    {
                        cg_store[(j-abs(j1-j2))/2][(m+j)/2][(m1+j1)/2]-=cg_store[(j_other-abs(j1-j2))/2][(m+j_other)/2][(m1+j1)/2]*ove_with_j_other[(j1+j2-j_other)/2];
                    }
                }
                double norm=0;
                for(int m1=m1_min;m1<=m1_max;m1+=2)
                {
                    norm+=cg_store[(j-abs(j1-j2))/2][(m+j)/2][(m1+j1)/2]*cg_store[(j-abs(j1-j2))/2][(m+j)/2][(m1+j1)/2];
                }
                norm=1/sqrt(norm);
                for(int m1=m1_min;m1<=m1_max;m1+=2)
                {
                    cg_store[(j-abs(j1-j2))/2][(m+j)/2][(m1+j1)/2]*=norm;
                }
                /***************************<j,m=j|J_{1z}|j+1,m=j>>=0相位保证 */
                double m1_exp=0;
                for(int m1=m1_min;m1<=m1_max;m1+=2)
                {
                    m1_exp+=m1*cg_store[(j-abs(j1-j2))/2][(m+j)/2][(m1+j1)/2]*cg_store[(j-abs(j1-j2))/2+1][(m+j)/2][(m1+j1)/2];
                }
                if(m1_exp<0)
                {
                    for(int m1=m1_min;m1<=m1_max;m1+=2)
                    {
                        int m2=m-m1;
                        cg_store[(j-abs(j1-j2))/2][(m+j)/2][(m1+j1)/2]*=-1;
                        cout<<"<"<<j1<<"/2,"<<j2<<"/2;"<<m1<<"/2,"<<m2<<"/2|"<<j<<"/2,"<<m<<"/2>="<<cg_store[(j-abs(j1-j2))/2][(m+j)/2][(m1+j1)/2]<<";<gsl>="<<cg_from_3j(j1,j2,m1,m2,j)<<";delta="<<cg_store[(j-abs(j1-j2))/2][(m+j)/2][(m1+j1)/2]-cg_from_3j(j1,j2,m1,m2,j)<<endl;
                    }
                }
                else
                {
                    for(int m1=m1_min;m1<=m1_max;m1+=2)
                    {
                        int m2=m-m1;
                        cout<<"<"<<j1<<"/2,"<<j2<<"/2;"<<m1<<"/2,"<<m2<<"/2|"<<j<<"/2,"<<m<<"/2>="<<cg_store[(j-abs(j1-j2))/2][(m+j)/2][(m1+j1)/2]<<";<gsl>="<<cg_from_3j(j1,j2,m1,m2,j)<<";delta="<<cg_store[(j-abs(j1-j2))/2][(m+j)/2][(m1+j1)/2]-cg_from_3j(j1,j2,m1,m2,j)<<endl;
                    }
                }
            }
            else
            {
                //J-|jm+1>=J1-|m1+1m2> +J2-|m1m2+1
                /*****************方程左边因子************* */
                double f_l=sqrt(j/2.0*(j/2.0+1)-m/2.0*(m/2.0+1));
                int m1_max=j1;
                int m1_min=-j1;
                if(m1_min<m-j2)
                {
                    m1_min=m-j2;
                }
                if(m1_max>m+j2)
                {
                    m1_max=m+j2;
                }
                for(int m1=m1_min;m1<=m1_max;m1+=2)
                {
                    int m2=m-m1;
                    cg_store[(j-abs(j1-j2))/2][(m+j)/2][(m1+j1)/2]=0;
                    /**********************源自J1-|m+1>=sum J1-|m1+1,m2>*********** */
                    if(m1+2<=j1&&m2+2>=-j1&&m2<=j2&&m2>=-j2)//m1+2也要在-j1~j1范围内
                    {
                        cg_store[(j-abs(j1-j2))/2][(m+j)/2][(m1+j1)/2]+=cg_store[(j-abs(j1-j2))/2][(m+j)/2+1][(m1+j1)/2+1]*sqrt(j1/2.0*(j1/2.0+1)-m1/2.0*(m1/2.0+1));
                    }
                    /**********************源自J2-|m+1>=sum J2-|m1,m2+2>*********** */
                    if(m2+2<=j2&&m2+2>=-j2&&m1<=j1&&m1>=-j1)//m2+2也要在-j2~j2范围内
                    {
                        cg_store[(j-abs(j1-j2))/2][(m+j)/2][(m1+j1)/2]+=cg_store[(j-abs(j1-j2))/2][(m+j)/2+1][(m1+j1)/2]*sqrt(j2/2.0*(j2/2.0+1)-m2/2.0*(m2/2.0+1));
                    }
                    cg_store[(j-abs(j1-j2))/2][(m+j)/2][(m1+j1)/2]/=f_l;
                    cout<<"<"<<j1<<"/2,"<<j2<<"/2;"<<m1<<"/2,"<<m2<<"/2|"<<j<<"/2,"<<m<<"/2>="<<cg_store[(j-abs(j1-j2))/2][(m+j)/2][(m1+j1)/2]<<";<gsl>="<<cg_from_3j(j1,j2,m1,m2,j)<<";delta="<<cg_store[(j-abs(j1-j2))/2][(m+j)/2][(m1+j1)/2]-cg_from_3j(j1,j2,m1,m2,j)<<endl;
                }
            }
        }
    }
    for(int j=abs(j1-j2);j<=j1+j2;j+=2)
    {
        for(int m=-j;m<=j;m+=2)
        {
            delete [] cg_store[(j-abs(j1-j2))/2][(m+j)/2];
        }
        delete [] cg_store[(j-abs(j1-j2))/2];
    }
    return 0;
}
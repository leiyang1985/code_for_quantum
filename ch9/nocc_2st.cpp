#include <iostream>
using namespace std;

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

double ove_f(int *a,int n);

double ove_f(int *c,int nc,int *a,int na,int *b,int nb,int *d,int nd)
{
    //<|>=1
    if(nc==0&&na==0&&nb==0&&nd==0)
    {
        return 1;
    }
    //<|c+=0
    if(nc!=0&&c[0]>0)
    {
        return 0;
    }
    //d|>=0
    if(nd!=0&&d[nd-1]<0)
    {
        return 0;
    }
    //<|c1~cn[a1,b1]d1~dn|>
    if(na==1&&nb==1)
    {
        if(a[0]+b[0]==0)
        {
            int e[nc+nd];
            for(int i=0;i<nc;i++)
            {
                e[i]=c[i];
            }
            for(int i=0;i<nd;i++)
            {
                e[i+nc]=d[i];
            }
            return ove_f(e,nc+nd);
        }
        else
        {
            return 0;
        }
    }
    //<|c1~cn[a1,b1~bn]d1~dn|>=<|c1~cn[a1,b1]b2~bnd1~dn|>-<|c1~cnb1[a1,b2~bn]d1~dn|>
    if(na==1&&nb!=1)
    {
        double res=0;
        if(a[0]+b[0]==0)
        {
            int e[nc+nb-1+nd];
            for(int i=0;i<nc;i++)
            {
                e[i]=c[i];
            }
            for(int i=0;i<nb-1;i++)
            {
                e[i+nc]=b[i+1];
            }
            for(int i=0;i<nd;i++)
            {
                e[i+nc+nb-1]=d[i];
            }
            res+=ove_f(e,nc+nb-1+nd);
        }
        int e[nc+1];
        for(int i=0;i<nc;i++)
        {
            e[i]=c[i];
        }
        e[nc]=b[0];
        res-=ove_f(e,nc+1,a,1,b+1,nb-1,d,nd);
        return res;
    }
    //<|c1~cn[a1~an,b1~bn]d1~dn|>=<|c1~cna1~an-1[an,b1~bn]d1~dn|>+theta<|c1~cn[a1~an-1,b1~bn]and1~dn|>
    if(na>1)
    {
        double res=0;
        int e[nc+na-1];
        for(int i=0;i<nc;i++)
        {
            e[i]=c[i];
        }
        for(int i=0;i<na-1;i++)
        {
            e[i+nc]=a[i];
        }
        res+=ove_f(e,nc+na-1,a+na-1,1,b,nb,d,nd);
        int f[nd+1];
        f[0]=a[na-1];
        for(int i=0;i<nd;i++)
        {
            f[i+1]=d[i];
        }
        if(nb%2==0)
        {
            res+=ove_f(c,nc,a,na-1,b,nb,f,nd+1);
        }
        else
        {
            res-=ove_f(c,nc,a,na-1,b,nb,f,nd+1);
        }
        return res;
    }
}
//<|a1~an|>=<|[a1~ai,ai+1~an]> for ai|>=0
double ove_f(int *a,int n)
{
    if(n==0)
    {
        return 1;
    }
    if(a[0]>0)
    {
        return 0;
    }
    if(a[n-1]<0)
    {
        return 0;
    }
    for(int i=n-1;i>=0;i--)
    {
        if(a[i]<0)
        {
            int b[n-i-1];
            for(int k=0;k<=n-i-2;k++)
            {
                b[k]=a[i+1+k];
            }
            return ove_f(NULL,0,a,i+1,b,n-i-1,NULL,0);
        }
    }
}

double ove_b(int *a,int n);

double ove_b(int *c,int nc,int *a,int na,int *b,int nb,int *d,int nd)
{
    //<|>=1
    if(nc==0&&na==0&&nb==0&&nd==0)
    {
        return 1;
    }
    //<|c+=0
    if(nc!=0&&c[0]>0)
    {
        return 0;
    }
    //d|>=0
    if(nd!=0&&d[nd-1]<0)
    {
        return 0;
    }
    //<|c1~cn[a1,b1]d1~dn|>
    if(na==1&&nb==1)
    {
        if(a[0]+b[0]==0)
        {
            int e[nc+nd];
            for(int i=0;i<nc;i++)
            {
                e[i]=c[i];
            }
            for(int i=0;i<nd;i++)
            {
                e[i+nc]=d[i];
            }
            if(a[0]<0)
            {
                return ove_b(e,nc+nd);
            }
            else
            {
                return -ove_b(e,nc+nd);
            }
        }
        else
        {
            return 0;
        }
    }
    //<|c1~cn[a1,b1~bn]d1~dn|>=<|c1~cn[a1,b1]b2~bnd1~dn|>+<|c1~cnb1[a1,b2~bn]d1~dn|>
    if(na==1&&nb!=1)
    {
        double res=0;
        if(a[0]+b[0]==0)
        {
            int e[nc+nb-1+nd];
            for(int i=0;i<nc;i++)
            {
                e[i]=c[i];
            }
            for(int i=0;i<nb-1;i++)
            {
                e[i+nc]=b[i+1];
            }
            for(int i=0;i<nd;i++)
            {
                e[i+nc+nb-1]=d[i];
            }
            if(a[0]<0)
            {
                res+=ove_b(e,nc+nb-1+nd);
            }
            else
            {
                res-=ove_b(e,nc+nb-1+nd);
            }
        }
        int e[nc+1];
        for(int i=0;i<nc;i++)
        {
            e[i]=c[i];
        }
        e[nc]=b[0];
        res+=ove_b(e,nc+1,a,1,b+1,nb-1,d,nd);
        return res;
    }
    //<|c1~cn[a1~an,b1~bn]d1~dn|>=<|c1~cna1~an-1[an,b1~bn]d1~dn|>+<|c1~cn[a1~an-1,b1~bn]and1~dn|>
    if(na>1)
    {
        double res=0;
        int e[nc+na-1];
        for(int i=0;i<nc;i++)
        {
            e[i]=c[i];
        }
        for(int i=0;i<na-1;i++)
        {
            e[i+nc]=a[i];
        }
        res+=ove_b(e,nc+na-1,a+na-1,1,b,nb,d,nd);
        int f[nd+1];
        f[0]=a[na-1];
        for(int i=0;i<nd;i++)
        {
            f[i+1]=d[i];
        }
        res+=ove_b(c,nc,a,na-1,b,nb,f,nd+1);
        return res;
    }
}
//<|a1~an|>=<|[a1~ai,ai+1~an]> for ai|>=0
double ove_b(int *a,int n)
{
    if(n==0)
    {
        return 1;
    }
    if(a[0]>0)
    {
        return 0;
    }
    if(a[n-1]<0)
    {
        return 0;
    }
    for(int i=n-1;i>=0;i--)
    {
        if(a[i]<0)
        {
            int b[n-i-1];
            for(int k=0;k<=n-i-2;k++)
            {
                b[k]=a[i+1+k];
            }
            return ove_b(NULL,0,a,i+1,b,n-i-1,NULL,0);
        }
    }
}

void mat_out(double **H,double **ove, int **basis_store,int dim,double *f_val,int **f_index,double *g_val,int **g_index,int g_num,int L,int N,bool isf)
{
    for(int i=0;i<dim;i++)
    {
        for(int j=0;j<dim;j++)
        {
            H[i][j]=0;
            int ket[N];
            int a_index=0;
            for(int k=0;k<L;k++)
            {
                for(int l=0;l<basis_store[j][k];l++)
                {
                    ket[a_index]=k+1;
                    a_index++;
                }
            }
            int bra[N];
            a_index=N-1;
            for(int k=0;k<L;k++)
            {
                for(int l=0;l<basis_store[i][k];l++)
                {
                    bra[a_index]=-(k+1);
                    a_index--;
                }
            }
            int a_cluster[2*N+4];
            for(int k=0;k<N;k++)
            {
                a_cluster[k]=bra[k];
            }
            for(int k=0;k<N;k++)
            {
                a_cluster[k+N]=ket[k];
            }
            if(isf)
            {
                ove[i][j]=ove_f(a_cluster,2*N);
            }
            else
            {
                ove[i][j]=ove_b(a_cluster,2*N);
            }
            H[i][j]=0;
            for(int k=0;k<N;k++)
            {
                a_cluster[k+N+2]=ket[k];
            }
            for(int k=0;k<L;k++)
            {
                a_cluster[N]=f_index[k][0]+1;
                a_cluster[N+1]=-(f_index[k][1]+1);
                if(isf)
                {
                    H[i][j]+=f_val[k]*ove_f(a_cluster,2*N+2);
                }
                else
                {
                    H[i][j]+=f_val[k]*ove_b(a_cluster,2*N+2);
                }
            }
            for(int k=0;k<N;k++)
            {
                a_cluster[k+N+4]=ket[k];
            }
            for(int k=0;k<g_num;k++)
            {
                a_cluster[N]=g_index[k][0]+1;
                a_cluster[N+1]=g_index[k][1]+1;
                a_cluster[N+2]=-(g_index[k][3]+1);
                a_cluster[N+3]=-(g_index[k][2]+1);
                if(isf)
                {
                    H[i][j]+=0.5*g_val[k]*ove_f(a_cluster,2*N+4);
                }
                else
                {
                    H[i][j]+=0.5*g_val[k]*ove_b(a_cluster,2*N+4);
                }
            }
        }
    }
}
void loop_one_state(int L,bool isf,int which_state,int N_left,int *basis,int **&basis_store,int &dim,bool isstore)
{
    int p_max=N_left;
    if(isf)
    {
        if(p_max>1)
        {
            p_max=1;
        }
    }
    if(which_state==L-1)
    {
        if(isf&&N_left>1)
        {
            return;
        }
        basis[which_state]=N_left;
        if(isstore)
        {
            for(int k=0;k<L;k++)
            {
                basis_store[dim][k]=basis[k];
            }
        }
        dim++;
    }
    else
    {
        for(basis[which_state]=0;basis[which_state]<=p_max;basis[which_state]++)
        {
            loop_one_state(L,isf,which_state+1,N_left-basis[which_state],basis,basis_store,dim,isstore);
        }
    }
    
}
void basis_out(int L,int N,bool isf,int **&basis_store,int &dim)
{
    int basis_one[L];
    if(L>=1)
    {
        dim=0;
        loop_one_state(L,isf,0,N,basis_one,basis_store,dim,false);
        cout<<"the dim is "<<dim<<endl;
        cout<<"do you want to go on?(Yes:1; NO:0)";
        int isgoon=0;
        cin>>isgoon;
        if(isgoon==1)
        {
            basis_store=new int *[dim];
            for(int i=0;i<dim;i++)
            {
                basis_store[i]=new int [L];
            }
            dim=0;
            loop_one_state(L,isf,0,N,basis_one,basis_store,dim,true);
        }
        else
        {
            exit(0);
        }
    }
}


void f_out(double *f_val,int **f_index,int L)
{
    for(int i=0;i<L;i++)
    {
        f_val[i]=i+0.5;
        f_index[i][0]=i;
        f_index[i][1]=i;
    }
}

void g_out(double *&g_val,int **&g_index,int &num,int L,bool isf)
{
    num=0;
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            if(i==j&&isf)
            {
                continue;
            }
            for(int k=0;k<L;k++)
            {
                for(int l=0;l<L;l++)
                {
                    if(k==l&&isf)
                    {
                        continue;
                    }
                    num++;
                }
            }
        }
    }
    g_val=new double [num];
    g_index=new int *[num];
    for(int i=0;i<num;i++)
    {
        g_index[i]=new int [4];
    }
    num=0;
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            if(i==j&&isf)
            {
                continue;
            }
            for(int k=0;k<L;k++)
            {
                for(int l=0;l<L;l++)
                {
                    if(k==l&&isf)
                    {
                        continue;
                    }
                    double H=0;
                    if(i==k+2&&j==l)
                    {
                        H+=0.25*sqrt((k+1)*(k+2));
                    }
                    if(i==k-2&&j==l)
                    {
                        H+=0.25*sqrt(k*(k-1));
                    }
                    if(i==k&&j==l)
                    {
                        H+=0.5*k;
                    }
                    H+=0.25;
                    if(i==k&&j==l+2)
                    {
                        H+=0.25*sqrt((l+1)*(l+2));
                    }
                    if(i==k&&j==l-2)
                    {
                        H+=0.25*sqrt(l*(l-1));
                    }
                    if(i==k&&j==l)
                    {
                        H+=0.5*l;
                    }
                    H+=0.25;
                    if(i==k+1&&j==l+1)
                    {
                        H-=0.5*sqrt((l+1)*(k+1));
                    }
                    if(i==k-1&&j==l-1)
                    {
                        H-=0.5*sqrt(l*k);
                    }
                    if(i==k+1&&j==l-1)
                    {
                        H-=0.5*sqrt((k+1)*l);
                    }
                    if(i==k-1&&j==l+1)
                    {
                        H-=0.5*sqrt(k*(l+1));
                    }
                    g_val[num]=H;
                    g_index[num][0]=i;
                    g_index[num][1]=j;
                    g_index[num][2]=k;
                    g_index[num][3]=l;
                    num++;
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
    int N=2;
    cout<<"particle num:";
    cin>>N;
    bool isf=true;
    cout<<"Fermion or Boson (1/0):";
    int fob=1;
    cin>>fob;
    if(fob==1)
    {
        isf=true;
    }
    else
    {
        isf=false;
    }
    int **basis;
    int dim;
    basis_out(L,N,isf,basis,dim);
    double *f_val=new double [L];
    int **f_index=new int *[L];
    for(int i=0;i<L;i++)
    {
        f_index[i]=new int [2];
    }
    f_out(f_val,f_index,L);
    double *g_val;
    int **g_index;
    int g_num=0;
    g_out(g_val,g_index,g_num,L,isf);
    // cout<<"g_num="<<g_num<<endl;
    double **H = new double *[dim];
    double **ove=new double *[dim];
    for (int i = 0; i < dim; i++)
    {
        H[i] = new double[dim];
        ove[i]=new double [dim];
    }
    mat_out(H,ove, basis,dim,f_val,f_index,g_val,g_index,g_num,L,N,isf);
    // cout<<"ove="<<endl;
    // for (int i = 0; i < dim; i++)
    // {
    //     for (int j = 0; j < dim; j++)
    //     {
    //         printf("%f ", ove[i][j]);
    //     }
    //     printf("\n");
    // }
    // cout<<"H="<<endl;
    // for (int i = 0; i < dim; i++)
    // {
    //     for (int j = 0; j < dim; j++)
    //     {
    //         printf("%f ", H[i][j]/sqrt(ove[i][i]*ove[j][j]));
    //     }
    //     printf("\n");
    // }
    /************************对角化的数据存储准备*/
    gsl_matrix *A=gsl_matrix_alloc(dim,dim);//对角化的A矩阵，算符矩阵
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            gsl_matrix_set(A,i,j,H[i][j]);//哈密顿矩阵元赋于对角化的对称矩阵
        }
    }
    gsl_matrix *B=gsl_matrix_alloc(dim,dim);//对角化所用的B对称矩阵，重叠矩阵
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            gsl_matrix_set(B,i,j,ove[i][j]);//重叠矩阵元赋于ove对称矩阵
        }
    }
    gsl_vector *eval = gsl_vector_alloc(dim);//本征值存储的数组
    gsl_matrix *evec = gsl_matrix_alloc(dim, dim);//存储本征态的二维数组，对也变换\Lambda矩阵
    gsl_eigen_gensymmv_workspace *w =gsl_eigen_gensymmv_alloc(dim);
    gsl_eigen_gensymmv(A, B,eval, evec, w);
    gsl_eigen_gensymmv_free(w);
    gsl_eigen_gensymmv_sort(eval, evec,GSL_EIGEN_SORT_VAL_ASC);//GSL_EIGEN_SORT_VAL_ASC排本征值从小到大排序。
    int print_max=10;
    if(dim<print_max)
    {
        print_max=dim;
    }
    for (int i = 0; i < print_max; i++)
    {
        double eval_i = gsl_vector_get(eval, i);
        gsl_vector_view evec_i = gsl_matrix_column(evec, i);

        printf("eigenvalue = %g \n", eval_i);
        // printf("eigenvector = \n");
        // for(int k=0;k<dim;k++)
        // {
        //     cout<<gsl_vector_get(&(evec_i.vector),k)<<"|";
        //     for(int l=0;l<L;l++)
        //     {
        //         cout<<basis[k][l]<<',';
        //     }
        //     cout<<">"<<endl;
        // }
    }

    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    return 0;
}
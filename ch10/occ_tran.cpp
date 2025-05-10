#include <iostream>
#include <cmath>
#include <cstring>
using namespace std;

#include <complex>
using Complex = std::complex<double>;


#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_hermite.h>
 
// 微分方程组函数指针类型
using DerivativeFunc = void(*)(double t, double w,const Complex y[], Complex dy_dt[], int n);
// RK4核心函数（使用原始数组实现）
void rk4_step(double t, double w,const Complex y[], Complex next_y[], double h, 
             DerivativeFunc f, int n) {
    // 分配临时存储空间
    Complex *k1 = new Complex[n];
    Complex *k2 = new Complex[n];
    Complex *k3 = new Complex[n];
    Complex *k4 = new Complex[n];
    Complex *y_temp = new Complex[n];
 
    // 计算k1
    f(t, w,y, k1, n);
    
    // 计算k2
    for(int i=0; i<n; ++i){
        y_temp[i] = y[i] + 0.5*h*k1[i];
    }
    f(t + 0.5*h, w,y_temp, k2, n);
    
    // 计算k3
    for(int i=0; i<n; ++i){
        y_temp[i] = y[i] + 0.5*h*k2[i];
    }
    f(t + 0.5*h, w,y_temp, k3, n);
    
    // 计算k4
    for(int i=0; i<n; ++i){
        y_temp[i] = y[i] + h*k3[i];
    }
    f(t + h, w,y_temp, k4, n);
    
    // 计算加权平均得到下一步解
    for(int i=0; i<n; ++i){
        Complex dy = (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0;
        next_y[i] = y[i] + h*dy;
    }
 
    // 释放临时内存
    delete[] k1;
    delete[] k2;
    delete[] k3;
    delete[] k4;
    delete[] y_temp;
}
/***********谐振子本征波函数 */
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

/*******************phi_m cos(x) phi_n积分函数 */
double mcosn(double x,void *par)
{
    void **par_arr=(void **)par;
    int m=*((int *)par_arr[0]);
    int n=*((int *)par_arr[1]);
    double w=*((double *)par_arr[2]);
    double t=*((double *)par_arr[3]);
    double bra=phi_n(m,x);
    double ket=phi_n(n,x);
    double res=bra*ket*cos(x);
    return res;
}


/******************<m|cos(x)cos(wt)|n> */
double mat_out(int m,int n,double omega,double t)
{
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    void *par[4];
    par[0]=(void *)(&m);
    par[1]=(void *)(&n);
    par[2]=(void *)(&omega);
    par[3]=(void *)(&t);
    gsl_function F;
    F.function = &mcosn;
    F.params = (void *)par;
    double res,err;
    gsl_integration_qagi (&F, 1e-3, 1e-7, 1000, w, &res, &err);
    gsl_integration_workspace_free(w);
    return res*cos(omega*t);
}


// 微分方程组：dy/dt = f(y,t)
void h_derive(double t, double w,const Complex y[], Complex dy_dt[],int n) {
    constexpr Complex I(0.0, 1.0);  // 虚数单位
    
    for(int i=0;i<n;i++)
    {
        dy_dt[i]=Complex(0,0);
        for(int j=0;j<n;j++)
        {
            dy_dt[i]+=exp(I*Complex(i-j,0)*Complex(t,0))*Complex(mat_out(i,j,w,t),0)*y[j];
            // dy_dt[i]+=exp(I*Complex(i-j,0)*Complex(t,0))*ei_out(i,j,w,t)*y[j];
        }
        dy_dt[i]*=-I;
    }
}
 
// 打印复数数组
void print_complex_array(const Complex arr[], int n) {
    std::cout << " ";
    double E_sum=0;
    for(int i=0; i<n; ++i){
        std::cout << abs(arr[i])*abs(arr[i]);
        if(i < n-1) std::cout << " ";
        E_sum+=abs(arr[i])*abs(arr[i])*(i+0.5);
    }
    std::cout <<' '<<E_sum<< " ";
}
 
int main()
{
    // 参数设置
    int n = 1;         // 方程个数
    cout<<"s.p. state number:";
    cin >> n;//输入单粒子维数
    double w;
    cout<<"hbar*omega:";
    cin >> w;//输入能级间隔
    const double x0 = 0.0;   // 初始x
    const double xf = 8*M_PI; // 终止x
    const double h = 0.1;    // 步长
 
    // 初始条件（复数数组）t=0时，处于基态
    Complex* y0 = new Complex[n];
    y0[0] = Complex(1.0, 0.0);  // 初始值 y0 = 1 + 0i
    for(int i=1;i<n;i++)
    {
        y0[i]=Complex(0,0);
    }
 
    // 存储结果的数组
    Complex* y = new Complex[n];
    Complex* next_y = new Complex[n];
 
    // 初始化当前解
    for(int i=0; i<n; ++i){
        y[i] = y0[i];
    }
 
    // 主循环
    std::cout << "t\t\t|C_k(t)|^2\t <E> \n";
    std::cout << "--------------------------------------\n";
    for(int step=0; ; ++step){
        double x = x0 + step*h;
        
        // 打印当前解
        std::cout.precision(4);
        std::cout << x << "\t";
        print_complex_array(y, n);
        std::cout << "\n";
        
        // 终止条件
        if(x >= xf) break;
        
        // 执行RK4步进
        rk4_step(x, w,y, next_y, h, h_derive, n);
        
        // 更新解向量（使用memcpy进行高效内存复制）
        std::memcpy(y, next_y, n*sizeof(Complex));
    }
 
    // 释放内存
    delete[] y0;
    delete[] y;
    delete[] next_y;
 
    return 0;
}
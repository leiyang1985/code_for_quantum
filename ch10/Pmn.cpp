#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>
using namespace std;
 
using Complex = std::complex<double>;

double Pmn(double wnm,double w,double t)
{
    Complex phase_p=Complex(0,1)*Complex((wnm-w)*t,0);
    Complex first=(exp(phase_p)-Complex(1,0))*Complex(w+wnm,0);
    Complex phase_n=Complex(0,1)*Complex((wnm+w)*t,0);
    Complex second=(exp(phase_n)-Complex(1,0))*Complex(w-wnm,0);
    double res=abs(first+second);
    res*=res;
    res/=(wnm*wnm-w*w)*(wnm*wnm-w*w);
    return res;
}

int main()
{
    Pmn(1,0.99,1000);
    double wnm=1;
    for(double w=-2;w<=2;w+=0.1)
    {
        cout<<w<<' ';
        for(double t=10;t<=10000;t*=10)
        {
            cout<<Pmn(wnm,w,t)<<' ';
        }
        cout<<endl;
    }
    for(double deltaw=0.01;deltaw>=1e-7;deltaw*=0.1)
    {
        double w=1-deltaw;
        cout<<w<<' ';
        for(double t=10;t<=10000;t*=10)
        {
            cout<<Pmn(wnm,w,t)<<' ';
        }
        cout<<endl;
        w=1+deltaw;
        cout<<w<<' ';
        for(double t=10;t<=10000;t*=10)
        {
            cout<<Pmn(wnm,w,t)<<' ';
        }
        cout<<endl;
        w=-1-deltaw;
        cout<<w<<' ';
        for(double t=10;t<=10000;t*=10)
        {
            cout<<Pmn(wnm,w,t)<<' ';
        }
        cout<<endl;
        w=-1+deltaw;
        cout<<w<<' ';
        for(double t=10;t<=10000;t*=10)
        {
            cout<<Pmn(wnm,w,t)<<' ';
        }
        cout<<endl;
    }
    // double w=1.00001;
    // for(double t=50000;t<=150000;t+=10000)
    // {
    //     cout<<t<<' '<<Pmn(wnm,w,t)<<endl;
    // }
    return 0;
}
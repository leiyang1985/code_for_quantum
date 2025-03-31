#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

double psi2(double x,void *par)//求norm的积分函数
{
  void **par_arr=(void **)par;
  gsl_interp_accel *acc=(gsl_interp_accel *)(par_arr[0]);
  gsl_spline *spline=(gsl_spline *)(par_arr[1]);
  double psi=gsl_spline_eval(spline, x, acc);
  return psi*psi;
}
double phipsi(double x,void *par)//求<phi|psi>内积
{
  void **par_arr=(void **)par;
  gsl_interp_accel *acc_l=(gsl_interp_accel *)(par_arr[0]);
  gsl_spline *spline_l=(gsl_spline *)(par_arr[1]);
  double norm_l=*((double *)par_arr[2]);
  gsl_interp_accel *acc_r=(gsl_interp_accel *)(par_arr[3]);
  gsl_spline *spline_r=(gsl_spline *)(par_arr[4]);
  double norm_r=*((double *)par_arr[5]);
  double phi=gsl_spline_eval(spline_l, x, acc_l)*norm_l;
  double psi=gsl_spline_eval(spline_r, x, acc_r)*norm_r;
  return phi*psi;
}

double phidpsi(double x,void *par)//求<phi|\hat p|psi>矩阵元，注意虚数单元被略去。
{
  void **par_arr=(void **)par;
  gsl_interp_accel *acc_l=(gsl_interp_accel *)(par_arr[0]);
  gsl_spline *spline_l=(gsl_spline *)(par_arr[1]);
  double norm_l=*((double *)par_arr[2]);
  gsl_interp_accel *acc_r=(gsl_interp_accel *)(par_arr[3]);
  gsl_spline *spline_r=(gsl_spline *)(par_arr[4]);
  double norm_r=*((double *)par_arr[5]);
  double phi=gsl_spline_eval(spline_l, x, acc_l)*norm_l;
  double dpsi=gsl_spline_eval_deriv(spline_r, x, acc_r)*norm_r;
  return phi*dpsi;
}

int main()
{
  FILE *file;                       // 文件指针
  double x[64];                   // 用于存储x数值的数组
  double y[6][64];                   // 用于存储y数值的数组
  char temp[64];
  int count = 0;                    // 记录读入数据点的个数
  file = fopen("ode.out.csv", "r"); // 打开文件
  fscanf(file,"%s %s %s %s %s %s %s",temp,temp,temp,temp,temp,temp,temp);//舍去文件头x y1 y2...标注
  printf("%s\n",temp);
  while (count < 64 && fscanf(file, "%lf %lf %lf %lf %lf %lf %lf", &x[count], &(y[0][count]),&(y[1][count]),&(y[2][count]),&(y[3][count]),&(y[4][count]),&(y[5][count])) == 7)// 读取文件中的浮点数
  {
    printf("%lf %lf %lf %lf %lf %lf %lf\n",x[count], (y[0][count]),(y[1][count]),(y[2][count]),(y[3][count]),(y[4][count]),(y[5][count]));
    count++; // 成功读取一个x位置后，计数器加1
  }
  fclose(file);
  //*****************为不同能量波函数建立不同的插值内存空间**************************//
  gsl_interp_accel *acc[6];
  gsl_spline *spline[6]; 
  for(int i=0;i<6;i++)
  {
    acc[i]= gsl_interp_accel_alloc();
    spline[i]= gsl_spline_alloc(gsl_interp_cspline, count);
    gsl_spline_init(spline[i], x, &(y[i][0]), count);
  }
  //*****************为不同能量波函数建立不同的插值内存空间**************************//

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);//建立求积分的内存空间
  double norm[6];
  double err;
  void *par[6];
  gsl_function F;
  F.function = &psi2;
  F.params = (void *)par;
  for(int i=0;i<6;i++)
  {
    par[0]=(void *)(acc[i]);
    par[1]=(void *)(spline[i]);
    
    gsl_integration_qags(&F, -4, 4, 1e-3, 1e-6, 1000, w, &(norm[i]), &err);//积分求norm
    norm[i]=1/sqrt(norm[i]);
  }

  //*********************正交归一性验证****************************/
  F.function = &phipsi;
  double product;
  for(int i=0;i<6;i++)
  {
    par[0]=(void *)(acc[i]);
    par[1]=(void *)(spline[i]);
    par[2]=(void *)(&(norm[i]));
    for(int j=0;j<6;j++)
    {
      par[3]=(void *)(acc[j]);
      par[4]=(void *)(spline[j]);
      par[5]=(void *)(&(norm[j]));
      gsl_integration_qags(&F, -4, 4, 1e-3, 1e-6, 1000, w, &product, &err);
      printf("%f ",product);
    }
    printf("\n");
  }
  //*********************正交归一性验证****************************/

  double p[6][6];
  //*********************<p>矩阵元计算****************************/
  F.function = &phidpsi;
  for(int i=0;i<6;i++)
  {
    par[0]=(void *)(acc[i]);
    par[1]=(void *)(spline[i]);
    par[2]=(void *)(&(norm[i]));
    for(int j=0;j<6;j++)
    {
      par[3]=(void *)(acc[j]);
      par[4]=(void *)(spline[j]);
      par[5]=(void *)(&(norm[j]));
      gsl_integration_qags(&F, -4, 4, 1e-3, 1e-6, 1000, w, &product, &err);
      printf("i*%f ",product);
      p[j][i]=product;
    }
    printf("\n");
  }
  //*********************<p>矩阵元计算****************************/

  //**************p的指数作用****/
  double eta=2;//输入eta
  double yk[count];
  double ysum[count];
  for(int i=0;i<count;i++)
  {
    yk[i]=y[0][i]*norm[0];//初始化c系数为|n=0>项
    ysum[i]=y[0][i]*norm[0];//定义最终系数数组，已加入|n=0>
  }
  int k=0;
  double sum2;
  gsl_interp_accel *acck= gsl_interp_accel_alloc();
  gsl_spline *splinek=gsl_spline_alloc(gsl_interp_cspline, count);
  do
  {
    gsl_spline_init(splinek, x, yk, count);
    sum2=0;
    for(int i=0;i<count;i++)
    {
      yk[i]=eta/double (k+1)*gsl_spline_eval_deriv(splinek, x[i], acck);//对应于式（\ref{eq:k2k1}）
      ysum[i]+=yk[i];//将第k+1项计入最终波函数中
      sum2+=yk[i]*yk[i];//sum2用于新增波函数大小，作为收敛依据
    }
    k++;//k计数加1
  } while (sum2>count*1e-14);
  for(int i=0;i<count;i++)
  {
    printf("%f %f\n",x[i],ysum[i]);//输入出最终波函数
  }
  gsl_spline_free (splinek);
  gsl_interp_accel_free (acck);
  //**************p的指数作用****/


  //*************使用矩阵图像来做****//
  double ck[6]={1,0,0,0,0,0};//初始化c系数为|n=0>项
  double sumc[6]={1,0,0,0,0,0};//定义最终系数数组，已加入|n=0>
  double ck1[6];
  k=0;
  do
  {
    sum2=0;
    for(int m=0;m<6;m++)
    {
      ck1[m]=0;
      for(int l=0;l<6;l++)
      {
        ck1[m]+=ck[l]*p[l][m]*eta/double(k+1);//对应于式（\ref{eq:k2k1}）
      }
      sumc[m]+=ck1[m];//将第k+1项计入最终波函数中
      sum2+=ck1[m]*ck1[m];//sum2用于衡量系数的大小，作为收敛依据
    }
    k++;//k计数加1
  } while (sum2>1e-14);
  
  for(int i=0;i<count;i++)
  {
    double ysum=0;
    for(int k=0;k<6;k++)
    {
      ysum+=sumc[k]*y[k][i]*norm[k];
    }
    // printf("%f %f\n",x[i],ysum);
  }
  //*************使用矩阵图像来做****//


  gsl_integration_workspace_free(w);
  for(int i=0;i<6;i++)
  {
    gsl_spline_free (spline[i]);
    gsl_interp_accel_free (acc[i]);
  }
  
  return 0;
}
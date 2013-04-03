/**************************************************************
*
* 文 件 名： function.c
* 
* 文件描述：
*     该文件包含了解题函数的函数体
*
* 创建人 ： 谌阳平
*
* 创建时间：2012/11/2
* 
* Character Set：UTF-8
*
* 版  本： 1.0  
*
***************************************************************/
#ifndef __FUNCTION__
#include "function.h"
#endif

#include <stdio.h>
#include <math.h>


//宏定义区
#define tmax 2000   //模拟总时长
#define g 9.8       //重力加速度
#define den 1       //文件数据写入密度


//全局变量定义区


double theta1,theta2,dtheta1[4],dtheta2[4],ddtheta1[4],ddtheta2[4],Energy,Energy_bef,dEnergy,Energy0;
double m1=1,m2=1,l1,l2,m3;
double K1,K2,u;
int step=0;
double dt0=0.0033;
double dt;
double Di1=1e-12;
double Di2;
double Di3=1e-8;
int cnt1;
int cnt2;
int cnt3;

/**************************************************************
*
* 函数名称：dt_etest
*
* 函数输入变量：void
*
* 返回值：double
*
* 函数功能：角速度误差计算
*
* 版本：1.0
*
* 创建人：谌阳平 
*
* 创建时间：2012/12/11
*
***************************************************************/
double dt_etest(void);

/**************************************************************
*
* 函数名称：t_etest
*
* 函数输入变量：void
*
* 返回值：double
*
* 函数功能：角度误差计算
*
* 版本：1.0
*
* 创建人：谌阳平 
*
* 创建时间：2012/12/11
*
***************************************************************/
double t_etest(void);

/**************************************************************
*
* 函数名称：E_etest
*
* 函数输入变量：void
*
* 返回值：double
*
* 函数功能：能量误差计算
*
* 版本：1.0
*
* 创建人：谌阳平 
*
* 创建时间：2012/12/11
*
***************************************************************/
double E_etest(void);

void Initialize(void){
    double cost1,sint1,cost2,sint2,cost2mt1,sint2mt1;

    m1 = 1;
    m2 = 1;
    l1 = 1;
    l2 = 1;
    m3 = 1;
    u=0;
    theta1 = 5;
    theta2 = 5;
    dtheta1[0] = 0;
    dtheta2[0] = 0;


    theta1 = theta1 * M_PI / 180;
    theta2 = theta2 * M_PI / 180;
    dtheta1[0] = dtheta1[0] * M_PI / 180;
    dtheta2[0] = dtheta2[0] * M_PI / 180;
    ddtheta1[0]=0;
    ddtheta2[0]=0;

    cost1 = cos(theta1);
    sint1 = sin(theta1);
    cost2 = cos(theta2);
    sint2 = sin(theta2);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;

    Energy=1/6.*(-3*g*((m2+2*m3)*l2*cost2+(m1+2*(m2+m3))*l1*cost1)+(m1+3*(m2+m3))*l1*l1*dtheta1[0]*dtheta1[0]+(m2+3*m3)*l2*l2*dtheta2[0]*dtheta2[0]+3*(m2+2*m3)*l1*l2*dtheta1[0]*dtheta2[0]*cost2mt1);
    Energy_bef=Energy;
    Energy0=Energy;
    dEnergy=Energy-Energy_bef;
    Di2=Di1*1e-1;
    // Di3=Di1*1e5;

    dt=dt0;
}

void move(void){
    double cost1,sint1,cost2,sint2,cost2mt1,sint2mt1;
    double A11,A12,A21,A22;
    double b1,b2;
    double ddtheta1_bef,ddtheta2_bef,dtheta1_tmp,dtheta2_tmp;
    double theta1_tmp,theta2_tmp;
    int i;

    for (i = 3; i > 0 ; --i)
    {
        dtheta1[i]=dtheta1[i-1];
        ddtheta1[i]=ddtheta1[i-1];
        dtheta2[i]=dtheta2[i-1];
        ddtheta2[i]=ddtheta2[i-1];
    }
    theta1_tmp=theta1;
    theta2_tmp=theta2;

    cost1 = cos(theta1);
    sint1 = sin(theta1);
    cost2 = cos(theta2);
    sint2 = sin(theta2);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;

    A11=(m1+3*m2+3*m3)*l1*l1/3;
    A12=(m2+2*m3)*l1*l2*cost2mt1/2;
    A21=(m2+2*m3)*l1*l2*cost2mt1/2;
    A22=(m2+3*m3)*l2*l2/3;

    b1=-l1*(m1/2+m2+m3)*g*sint1+(m2/2+m3)*l1*l2*dtheta2[0]*dtheta2[0]*sint2mt1;
    b2=-l2*(m2/2+m3)*(g*sint2+l1*dtheta1[0]*dtheta1[0]*sint2mt1);


    ddtheta1_bef=(b1*A22-b2*A12)/(A11*A22-A21*A12);
    ddtheta2_bef=(b2*A11-b1*A21)/(A11*A22-A21*A12);


    theta1=theta1+dtheta1[0]*dt+ddtheta1_bef*dt*dt/2;
    theta2=theta2+dtheta2[0]*dt+ddtheta2_bef*dt*dt/2;
    dtheta1_tmp=dtheta1[0];
    dtheta2_tmp=dtheta2[0];
    dtheta1[0]=dtheta1[0]+ddtheta1_bef*dt;
    dtheta2[0]=dtheta2[0]+ddtheta2_bef*dt;
    // printf("%e %e %e %e\n", theta1,theta2,dtheta1[0],dtheta2[0]);


    cost1 = cos(theta1);
    sint1 = sin(theta1);
    cost2 = cos(theta2);
    sint2 = sin(theta2);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;

    A11=(m1+3*m2+3*m3)*l1*l1/3;
    A12=(m2+2*m3)*l1*l2*cost2mt1/2;
    A21=(m2+2*m3)*l1*l2*cost2mt1/2;
    A22=(m2+3*m3)*l2*l2/3;
 
    b1=-l1*(m1/2+m2+m3)*g*sint1+(m2/2+m3)*l1*l2*dtheta2[0]*dtheta2[0]*sint2mt1;
    b2=-l2*(m2/2+m3)*(g*sint2+l1*dtheta1[0]*dtheta1[0]*sint2mt1);

    // printf("%8.2e %8.2e %8.2e %8.2e %8.2e %8.2e\n",Fr,FMr,Tt,T_r,Ft,FMt);
    // getchar();


    ddtheta1[0]=(b1*A22-b2*A12)/(A11*A22-A21*A12);
    ddtheta2[0]=(b2*A11-b1*A21)/(A11*A22-A21*A12);
    // printf("%e %e\n", ddtheta1[0],ddtheta1_bef);
    dtheta1[0]=dtheta1_tmp+(ddtheta1[0]+ddtheta1_bef)*dt/2;
    dtheta2[0]=dtheta2_tmp+(ddtheta2[0]+ddtheta2_bef)*dt/2;
    ddtheta1[0]=ddtheta1_bef;
    ddtheta2[0]=ddtheta2_bef;

    // theta1=theta1_tmp+(dtheta1[0]+dtheta1[1])*dt/2+(ddtheta1[1]+ddtheta1[0])*dt*dt/4;
    // theta2=theta2_tmp+(dtheta2[0]+dtheta2[1])*dt/2+(ddtheta2[1]+ddtheta2[0])*dt*dt/4;


    Energy=1/6.*(-3*g*((m2+2*m3)*l2*cost2+(m1+2*(m2+m3))*l1*cost1)+(m1+3*(m2+m3))*l1*l1*dtheta1[0]*dtheta1[0]+(m2+3*m3)*l2*l2*dtheta2[0]*dtheta2[0]+3*(m2+2*m3)*l1*l2*dtheta1[0]*dtheta2[0]*cost2mt1);
    dEnergy=Energy-Energy_bef;
    Energy_bef=Energy;
}

void move_RK(double h){
    double cost1,sint1,cost2,sint2,cost2mt1,sint2mt1;
    double k[5][5];
    double A11,A12,A21,A22;
    double b1,b2;
    double theta1_tmp,theta2_tmp;
    int i;

    for (i = 3; i > 0 ; --i)
    {
        dtheta1[i]=dtheta1[i-1];
        ddtheta1[i]=ddtheta1[i-1];
        dtheta2[i]=dtheta2[i-1];
        ddtheta2[i]=ddtheta2[i-1];
    }

    theta1_tmp=theta1;
    theta2_tmp=theta2;

    cost1 = cos(theta1);
    sint1 = sin(theta1);
    cost2 = cos(theta2);
    sint2 = sin(theta2);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;

    A11=(m1+3*m2+3*m3)*l1*l1/3;
    A12=(m2+2*m3)*l1*l2*cost2mt1/2;
    A21=(m2+2*m3)*l1*l2*cost2mt1/2;
    A22=(m2+3*m3)*l2*l2/3;

    b1=-l1*(m1/2+m2+m3)*g*sint1+(m2/2+m3)*l1*l2*dtheta2[0]*dtheta2[0]*sint2mt1;
    b2=-l2*(m2/2+m3)*(g*sint2+l1*dtheta1[0]*dtheta1[0]*sint2mt1);

    k[1][1]=(b1*A22-b2*A12)/(A11*A22-A21*A12);
    k[2][1]=(b2*A11-b1*A21)/(A11*A22-A21*A12);
    k[3][1]=dtheta1[0];
    k[4][1]=dtheta2[0];


    dtheta1[0]=dtheta1[1]+h*k[1][1]/2;
    dtheta2[0]=dtheta2[1]+h*k[2][1]/2;
    theta1=theta1_tmp+k[3][1]*h/2;
    theta2=theta2_tmp+k[4][1]*h/2;
    // printf("%e %e %e %e\n", theta1,theta2,dtheta1[0],dtheta2[0]);

    cost1 = cos(theta1);
    sint1 = sin(theta1);
    cost2 = cos(theta2);
    sint2 = sin(theta2);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;

    A11=(m1+3*m2+3*m3)*l1*l1/3;
    A12=(m2+2*m3)*l1*l2*cost2mt1/2;
    A21=(m2+2*m3)*l1*l2*cost2mt1/2;
    A22=(m2+3*m3)*l2*l2/3;

    b1=-l1*(m1/2+m2+m3)*g*sint1+(m2/2+m3)*l1*l2*dtheta2[0]*dtheta2[0]*sint2mt1;
    b2=-l2*(m2/2+m3)*(g*sint2+l1*dtheta1[0]*dtheta1[0]*sint2mt1);

    k[1][2]=(b1*A22-b2*A12)/(A11*A22-A21*A12);
    k[2][2]=(b2*A11-b1*A21)/(A11*A22-A21*A12);
    k[3][2]=dtheta1[0];
    k[4][2]=dtheta2[0];

    dtheta1[0]=dtheta1[1]+h*k[1][2]/2;
    dtheta2[0]=dtheta2[1]+h*k[2][2]/2;
    theta1=theta1_tmp+k[3][2]*h/2;
    theta2=theta2_tmp+k[4][2]*h/2;

    cost1 = cos(theta1);
    sint1 = sin(theta1);
    cost2 = cos(theta2);
    sint2 = sin(theta2);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;

    A11=(m1+3*m2+3*m3)*l1*l1/3;
    A12=(m2+2*m3)*l1*l2*cost2mt1/2;
    A21=(m2+2*m3)*l1*l2*cost2mt1/2;
    A22=(m2+3*m3)*l2*l2/3;

    b1=-l1*(m1/2+m2+m3)*g*sint1+(m2/2+m3)*l1*l2*dtheta2[0]*dtheta2[0]*sint2mt1;
    b2=-l2*(m2/2+m3)*(g*sint2+l1*dtheta1[0]*dtheta1[0]*sint2mt1);

    k[1][3]=(b1*A22-b2*A12)/(A11*A22-A21*A12);
    k[2][3]=(b2*A11-b1*A21)/(A11*A22-A21*A12);
    k[3][3]=dtheta1[0];
    k[4][3]=dtheta2[0];

    dtheta1[0]=dtheta1[1]+h*k[1][3];
    dtheta2[0]=dtheta2[1]+h*k[2][3];
    theta1=theta1_tmp+k[3][3]*h;
    theta2=theta2_tmp+k[4][3]*h;

    cost1 = cos(theta1);
    sint1 = sin(theta1);
    cost2 = cos(theta2);
    sint2 = sin(theta2);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;

    A11=(m1+3*m2+3*m3)*l1*l1/3;
    A12=(m2+2*m3)*l1*l2*cost2mt1/2;
    A21=(m2+2*m3)*l1*l2*cost2mt1/2;
    A22=(m2+3*m3)*l2*l2/3;

    b1=-l1*(m1/2+m2+m3)*g*sint1+(m2/2+m3)*l1*l2*dtheta2[0]*dtheta2[0]*sint2mt1;
    b2=-l2*(m2/2+m3)*(g*sint2+l1*dtheta1[0]*dtheta1[0]*sint2mt1);

    k[1][4]=(b1*A22-b2*A12)/(A11*A22-A21*A12);
    k[2][4]=(b2*A11-b1*A21)/(A11*A22-A21*A12);
    k[3][4]=dtheta1[0];
    k[4][4]=dtheta2[0];


    dtheta1[0]=dtheta1[1]+h*(k[1][1]+2*k[1][2]+2*k[1][3]+k[1][4])/6;
    dtheta2[0]=dtheta2[1]+h*(k[2][1]+2*k[2][2]+2*k[2][3]+k[2][4])/6;
    theta1=theta1_tmp+h*(k[3][1]+2*k[3][2]+2*k[3][3]+k[3][4])/6;
    theta2=theta2_tmp+h*(k[4][1]+2*k[4][2]+2*k[4][3]+k[4][4])/6;

    // printf("%e %e %e %e\n", k[1][1],k[1][2],k[1][3],k[1][4]);
    // printf("%e\n\n", (k[1][1]+2*k[1][2]+2*k[1][3]+k[1][4])/6); 


    Energy=1/6.*(-3*g*((m2+2*m3)*l2*cost2+(m1+2*(m2+m3))*l1*cost1)+(m1+3*(m2+m3))*l1*l1*dtheta1[0]*dtheta1[0]+(m2+3*m3)*l2*l2*dtheta2[0]*dtheta2[0]+3*(m2+2*m3)*l1*l2*dtheta1[0]*dtheta2[0]*cost2mt1);
    dEnergy=Energy-Energy0;
    Energy_bef=Energy;
}

double t_calc(void){
    double Dc3;
    double dt_update;
    Dc3=E_etest();

    if(Dc3/Di3<1e-15) return 3*dt;
    dt_update=dt*pow(Di3/Dc3, 0.2)*0.8;


    return dt_update;
}

double dt_etest(void){
    double vt2;
    double vt3;
    double theta1_t,theta2_t;
    double dt1[4],dt2[4],ddt1[4],ddt2[4];
    double E_t1,E_t2;
    int i;
    for (i = 3; i >= 0 ; --i)
    {
        dt1[i]=dtheta1[i];
        ddt1[i]=ddtheta1[i];
        dt2[i]=dtheta2[i];
        ddt2[i]=ddtheta2[i];
    }
    theta1_t=theta1;
    theta2_t=theta2;
    E_t1=Energy;
    E_t2=Energy_bef;

    move_RK(dt);
    vt2=dtheta1[0];

    for (i = 0; i < 4 ; ++i)
    {
        dtheta1[i]=dt1[i];
        ddtheta1[i]=ddt1[i];
        dtheta2[i]=dt2[i];
        ddtheta2[i]=ddt2[i];
    }
    theta1=theta1_t;
    theta2=theta2_t;

    move_RK(dt/2);
    move_RK(dt/2);
    vt3=dtheta1[0];
    // printf("%e\n", vt3-vt2);
    // getchar();

    for (i = 0; i < 4 ; ++i)
    {
        dtheta1[i]=dt1[i];
        ddtheta1[i]=ddt1[i];
        dtheta2[i]=dt2[i];
        ddtheta2[i]=ddt2[i];
    }
    theta1=theta1_t;
    theta2=theta2_t;
    Energy=E_t1;
    Energy_bef=E_t2;
    return fabs(vt3-vt2);
}

double t_etest(void){
    double t2;
    double t3;
    double theta1_t,theta2_t;
    double dt1[4],dt2[4],ddt1[4],ddt2[4];
    double E_t1,E_t2;
    int i;
    for (i = 3; i >= 0 ; --i)
    {
        dt1[i]=dtheta1[i];
        ddt1[i]=ddtheta1[i];
        dt2[i]=dtheta2[i];
        ddt2[i]=ddtheta2[i];
    }
    theta1_t=theta1;
    theta2_t=theta2;
    E_t1=Energy;
    E_t2=Energy_bef;
    move_RK(dt);
    t2=theta1;

    for (i = 0; i < 4 ; ++i)
    {
        dtheta1[i]=dt1[i];
        ddtheta1[i]=ddt1[i];
        dtheta2[i]=dt2[i];
        ddtheta2[i]=ddt2[i];
    }
    theta1=theta1_t;
    theta2=theta2_t;

    move_RK(dt/2);
    move_RK(dt/2);
    t3=theta1;
    // printf("%e\n", vt3-vt2);
    // getchar();

    for (i = 0; i < 4 ; ++i)
    {
        dtheta1[i]=dt1[i];
        ddtheta1[i]=ddt1[i];
        dtheta2[i]=dt2[i];
        ddtheta2[i]=ddt2[i];
    }
    theta1=theta1_t;
    theta2=theta2_t;
    Energy=E_t1;
    Energy_bef=E_t2;
    return fabs(t3-t2);
}

double E_etest(void){
    double e2;
    double e3;
    double theta1_t,theta2_t;
    double dt1[4],dt2[4],ddt1[4],ddt2[4];
    double E_t1,E_t2;
    int i;
    for (i = 3; i >= 0 ; --i)
    {
        dt1[i]=dtheta1[i];
        ddt1[i]=ddtheta1[i];
        dt2[i]=dtheta2[i];
        ddt2[i]=ddtheta2[i];
    }
    theta1_t=theta1;
    theta2_t=theta2;
    E_t1=Energy;
    E_t2=Energy_bef;
    move_RK(dt);
    e2=Energy;

    for (i = 0; i < 4 ; ++i)
    {
        dtheta1[i]=dt1[i];
        ddtheta1[i]=ddt1[i];
        dtheta2[i]=dt2[i];
        ddtheta2[i]=ddt2[i];
    }
    theta1=theta1_t;
    theta2=theta2_t;
    Energy=E_t1;
    Energy_bef=E_t2;


    move_RK(dt/2);
    move_RK(dt/2);
    e3=Energy;
    // printf("%e\n", vt3-vt2);
    // getchar();

    for (i = 0; i < 4 ; ++i)
    {
        dtheta1[i]=dt1[i];
        ddtheta1[i]=ddt1[i];
        dtheta2[i]=dt2[i];
        ddtheta2[i]=ddt2[i];
    }
    theta1=theta1_t;
    theta2=theta2_t;
    Energy=E_t1;
    Energy_bef=E_t2;
    return fabs(e3-e2);
}
/**************************************************************
*
* 文 件 名： main.c
* 
* 文件描述：
*        计算机模拟物理第三次作业主程序
*
* 创建人 ： 谌阳平
*
* 创建时间：2012/11/2
*
* 版  本： 1.0   
* 
* Character Set：UTF-8
*
* 编译平台：Mac OS X 10.8.2 with Intel Core i5 and 8GB RAM
*
* 编译器版本：gcc version 4.2.1 (Based on Apple Inc. build 5658) (LLVM build 2336.11.00)
*
***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "var.h"


struct sphere
{
	double mass;
    double l;
	double theta;
	double d_theta;
	double dd_theta;
}s1,s2;



void Initialize(void);
void move(void);

int main(int argc, char const *argv[])
{
    int step=0;
	Initialize();
    FILE *fp = fopen("/Users/Fermi/Desktop/data_plot/data.dat", "w");
    // fprintf(fp, "#%-7s %-12s %-12s %-12s %-12s\n","t","s1.theta","s2.theta","s1.d_theta","s2.d_theta");
    while(step * dt < tmax){
        if(step%1==0) fprintf(fp, "%-8.3f %-12.5e %-12.5e %-12.5e %-12.5e\n",step*dt,s1.theta*180/M_PI,s2.theta*180/M_PI,s1.d_theta*180/M_PI,s2.d_theta*180/M_PI);
        // getchar();
        move();
        step++;
    }
    fclose(fp);

    printf("\nProject is Finished\n");
	
}

void Initialize(){
    /*printf("Please input the mass of sphere 1: ");
    scanf("%lf",&s1.mass);
    printf("Please input the length of line 1: ");
    scanf("%lf",&s1.l);
	printf("Please input initial theta of sphere 1: ");
    scanf("%lf",&s1.theta);
    printf("Please input initial d_theta of sphere 1: ");
    scanf("%lf",&s1.d_theta);
    
    printf("\n");
    printf("Please input the mass of sphere 2: ");
    scanf("%lf",&s2.mass);
    printf("Please input the length of line 2: ");
    scanf("%lf",&s2.l);
	printf("Please input initial theta of sphere 2: ");
    scanf("%lf",&s2.theta);
    printf("Please input initial d_theta of sphere 2: ");
    scanf("%lf",&s2.d_theta);*/
    s1.mass = 1;
    s2.mass = 0.2;
    s1.l = 1;
    s2.l = 0.3;
    s1.theta = 180;
    s2.theta = 90;
    s1.d_theta = 0;
    s2.d_theta = 0;
    
    s1.theta = s1.theta * M_PI / 180;
    s2.theta = s2.theta * M_PI / 180;
    s1.dd_theta=0;
    s2.dd_theta=0;
}

void move(void){
    double T2,cost1,sint1,cost2,sint2,cost2mt1,sint2mt1;
    double s1dd_theta,s2dd_theta,s1d_theta,s2d_theta;
    s1d_theta=s1.d_theta;
    s2d_theta=s2.d_theta;

    cost1 = cos(s1.theta);
    sint1 = sin(s1.theta);
    cost2 = cos(s2.theta);
    sint2 = sin(s2.theta);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;

    T2 = (s2.mass * g * sint1 * sint2mt1 +  s2.mass * g * cost2 + s2.mass * s1.l * s1.d_theta * s1.d_theta * cost2mt1 + s2.mass * s2.l * s2.d_theta * s2.d_theta) / (s1.mass + s2.mass * sint2mt1 * sint2mt1);
    // T1 = s1.mass * g * cost1 + s1.mass * T2 * cost2mt1 + s1.mass * s1.l * s1.d_theta *s1.d_theta;
    s1dd_theta = T2 * sint2mt1 / s1.l - g * sint1 / s1.l;
    s2dd_theta = -T2 * sint2mt1 / s2.l + g * sint1 / s2.l - g * sint2 / s2.l - s1.l * s1.d_theta * s1.d_theta * sint2mt1 / s2.l;

    s1.theta = s1.theta + s1.d_theta * dt + s1dd_theta * dt * dt / 2;
    s2.theta = s2.theta + s2.d_theta * dt + s2dd_theta * dt * dt / 2;

    s1.d_theta = s1.d_theta + s1dd_theta * dt;
    s2.d_theta = s2.d_theta + s2dd_theta * dt;

    cost1 = cos(s1.theta);
    sint1 = sin(s1.theta);
    cost2 = cos(s2.theta);
    sint2 = sin(s2.theta);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;

    T2 = (s2.mass * g * sint1 * sint2mt1 + s2.mass * g * cost2 + s2.mass * s1.l * s1.d_theta * s1.d_theta * cost2mt1 + s2.mass * s2.l * s2.d_theta * s2.d_theta) / (s1.mass + s2.mass * sint2mt1 * sint2mt1);
    // T1 = s1.mass * g * cost1 + s1.mass * T2 * cost2mt1 + s1.mass * s1.l * s1.d_theta *s1.d_theta;

    s1.dd_theta = T2 * sint2mt1 / s1.l - g * sint1 / s1.l;
    s2.dd_theta = -T2 * sint2mt1 / s2.l + g * sint1 / s2.l - g * sint2 / s2.l - s1.l * s1.d_theta * s1.d_theta * sint2mt1 / s2.l;

    s1.d_theta = s1d_theta + (s1.dd_theta + s1dd_theta) * dt / 2;
    s2.d_theta = s2d_theta + (s2.dd_theta + s2dd_theta) * dt / 2;
}
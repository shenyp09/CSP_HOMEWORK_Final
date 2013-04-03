/**************************************************************
*
* 文 件 名： function.h
* 
* 文件描述：
*     该文件包含了解题函数的函数声明
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
#define __FUNCTION__ 

#include "function.c"


/**************************************************************
*
* 函数名称：Initialize
*
* 函数输入变量：void
*
* 返回值：void
*
* 函数功能：数据初始化
*
* 版本：1.0
*
* 创建人：谌阳平 
*
* 创建时间：2012/12/11
*
***************************************************************/
void Initialize(void);

/**************************************************************
*
* 函数名称：move
*
* 函数输入变量：void
*
* 返回值：void
*
* 函数功能：Verlet算法
*
* 版本：1.0
*
* 创建人：谌阳平 
*
* 创建时间：2012/12/11
*
***************************************************************/
void move(void);

/**************************************************************
*
* 函数名称：move_RK
*
* 函数输入变量：double h
*
* 返回值：void
*
* 函数功能：经典四阶Runge-Kutta算法
*
* 版本：1.0
*
* 创建人：谌阳平 
*
* 创建时间：2012/12/11
*
***************************************************************/
void move_RK(double h);

/**************************************************************
*
* 函数名称：t_calc
*
* 函数输入变量：void
*
* 返回值：double
*
* 函数功能：更新时间步长
*
* 版本：1.0
*
* 创建人：谌阳平 
*
* 创建时间：2012/12/11
*
***************************************************************/
double t_calc(void);


#endif


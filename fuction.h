#ifndef FUCTION_H_INCLUDED
#define FUCTION_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265

double solid_rate(double sx,double sy,double sl,double cx,double cy,double r);
double s_tri(double dot[][2]);
double s_fan(double r,double dot[][2]);
void line_cir(double cx,double cy,double r,double dot1[],double dot2[]);
double s_tra(double sl,double dot1[],double dot2[],double dot3[],double dot4[]);
int compare_x(const void *p,const void *q);
int compare_y(const void *p,const void *q);
void cmp_two(double a[],double b[]);



#endif // FUCTION_H_INCLUDED

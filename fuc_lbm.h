#ifndef FUC_LBM_H_INCLUDED
#define FUC_LBM_H_INCLUDED

#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>

void init();
double feq(int k,double rho, double u[2]);
void evolution();
void output(int m);
void Error();
void output_ex1(int m);
double solid_rate(double sx,double sy,double sl,double cx,double cy,double r);


#endif // FUC_LBM_H_INCLUDED

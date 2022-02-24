#ifndef FUC_LBM_H_INCLUDED
#define FUC_LBM_H_INCLUDED

#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>

#define PI 3.14159265

using namespace std;
const int Q = 9;            //D2Q9ģ��
const int NX = 200;         //x����
const int NY = 800;         //y����
const int N = 2;            //������Ŀ
double cir_r = 12.5;  //�뾶

int e[Q][2] = { {0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1} };
int oppk[9] = { 0,3,4,1,2,7,8,5,6 };                                         //Ϊ�˱�ʾf[-i](x,t)������Ӧ�ķֲ��������෴�ķ���
double w[Q] = { 4.0 / 9,1.0 / 9,1.0 / 9,1.0 / 9,1.0 / 9,1.0 / 36,1.0 / 36,1.0 / 36,1.0 / 36 };
/*������ĩβ��0��ʾ��Ӧ����һʱ�䲽����Ϣ*/
/*          �����ܶ�                �����ٶ�        ��һʱ�䲽�������ٶ�      �ֲ�����           ��ײ��ķֲ�����*/
double rho[NX + 1][NY + 1], u[NX + 1][NY + 1][2], u0[NX + 1][NY + 1][2], f[NX + 1][NY + 1][Q], F[NX + 1][NY + 1][Q];
/*          ���������ٶ�    ��һʱ�䲽�������ٶ�     λ��          ���ٶ�             �Ƕ�             ��������ɢ��Ԫ���ٶ�*/
double velocity_central[N][2], velocity0_central[N][2], cir_x[N], cir_y[N], angular[N], angular0[N], angle[N], angle0[N], vel_particle[NX + 1][NY + 1][2];
/*          �������          rate��Ӧϵ��            ������ײ��                 ���������м���*/
double rate[NX + 1][NY + 1], rate0[N][NX + 1][NY + 1], B[NX + 1][NY + 1], B0[N][NX + 1][NY + 1], om[NX + 1][NY + 1][Q], OM[NX + 1][NY + 1][2] = { 0 };
/*          ��������    ת��*/
double Ff[N][2], Ff0[N][2], Tf[N], Tf0[N];               //������һ�ڵ����Ϣ
int i, j, k, n, a;
double aver_rho, Cl, Cd;
/*  �����ٶ�               �����ʼ�ܶ� �����ܶ�    ճ��  ��� */
double c, Re, dx, dy, Lx, Ly, dt, rho0, rhos, tau_f, niu, error;
/*           �������� ת������ ��������ϵ��*/
double g, K, W_, Up, mass, I, coefficient;
/*  */
double lx, Ch, density, Cden, viscosity, Ct, Cf, f_;
double x[N][2];
int pi, pj;
double  Fp[N][N][2], Fp_sum[N][2], Fw[N][2], Fcol[N][2];
double len, threshold;
double ep, ep_, ew, ew_, d1, d2, d3, d4;
double leftw, rightw, upw, downw;

void init();
double feq(int k, double rho, double u[2]);
void evolution();
void output(int m);
void Error();
double solid_rate(double i, double j, double cx, double cy);
void partical_collision(int i, int j);
double partical_wall(double d);

#endif // FUC_LBM_H_INCLUDED

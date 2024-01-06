#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

int Testnode(double x[2], double cx, double cy);
double lengthratio(double x1[2], double x2[2], double cx, double cy);

extern double cir_r;

double solid_rate(double i, double j, double cx, double cy)
{
    int k, flag[4], temp, m, n;
    double nodenear[4][2], q1, q2, l, epsilon;

    nodenear[0][0] = i - 0.5; nodenear[0][1] = j - 0.5;//左下
    nodenear[1][0] = i + 0.5; nodenear[1][1] = j - 0.5;//右下
    nodenear[2][0] = i + 0.5; nodenear[2][1] = j + 0.5;//右上
    nodenear[3][0] = i - 0.5; nodenear[3][1] = j + 0.5;//左上
    for (k = 0; k < 4; k++)
        flag[k] = Testnode(nodenear[k], cx, cy);
    temp = flag[3] + flag[2] + flag[1] + flag[0];   //在圆外的角点个数
    switch (temp)
    {
    case 0:
        epsilon = 1;
        break;
    case 1:
    {
        m = 3 * flag[3] + 2 * flag[2] + 1 * flag[1] + 0 * flag[0];
        n = 3 * flag[3] + 2 * flag[2] + 1 * flag[1] + 0 * flag[0];  //在圆外的点序号，与m相同
        q1 = lengthratio(nodenear[m], nodenear[(m - 1 + 4) % 4], cx, cy);    //该点与顺时针邻点成的线
        q2 = lengthratio(nodenear[n], nodenear[(n + 1 + 4) % 4], cx, cy);    //该点与逆时针邻点成的线
        l = sqrt(q1 * q1 + q2 * q2);
        epsilon = (1 - 0.5 * q1 * q2) + (cir_r * cir_r * asin(0.5 * l / cir_r) - 0.25 * l * sqrt(4 * cir_r * cir_r - l * l));
    }
    break;
    case 2:
    {
        m = (((0 * flag[0] + 4 * flag[0] * flag[3]) + 1 * flag[1] + 2 * flag[2] + 3 * flag[3] - 1) / 2) % 4;
        n = (((0 * flag[0] + 4 * flag[0] * flag[3]) + 1 * flag[1] + 2 * flag[2] + 3 * flag[3] + 1) / 2) % 4;
        q1 = lengthratio(nodenear[m], nodenear[(m - 1 + 4) % 4], cx, cy);
        q2 = lengthratio(nodenear[n], nodenear[(n + 1 + 4) % 4], cx, cy);
        l = sqrt(1 + (q2 - q1) * (q2 - q1));
        epsilon = 0.5 * (2 - q1 - q2) + (cir_r * cir_r * asin(0.5 * l / cir_r) - 0.25 * l * sqrt(4 * cir_r * cir_r - l * l));
    }
    break;
    case 3:
    {
        m = (((0 * flag[0] + 4 * flag[0] * flag[3] * flag[2] + 4 * flag[1] * flag[0] * flag[3]) + (1 * flag[1] + 4 * flag[1] * flag[0] * flag[3]) + 2 * flag[2] + 3 * flag[3] - 3) / 3) % 4;
        n = (((0 * flag[0] + 4 * flag[0] * flag[3] * flag[2] + 4 * flag[1] * flag[0] * flag[3]) + (1 * flag[1] + 4 * flag[1] * flag[0] * flag[3]) + 2 * flag[2] + 3 * flag[3] + 3) / 3) % 4;
        q1 = lengthratio(nodenear[m], nodenear[(m - 1 + 4) % 4], cx, cy);
        q2 = lengthratio(nodenear[n], nodenear[(n + 1 + 4) % 4], cx, cy);
        l = sqrt((1 - q1) * (1 - q1) + (1 - q2) * (1 - q2));
        epsilon = 0.5 * (1 - q1) * (1 - q2) + (cir_r * cir_r * asin(0.5 * l / cir_r) - 0.25 * l * sqrt(4 * cir_r * cir_r - l * l));
    }
    break;
    case 4:
        epsilon = 0;
        break;
    }
    return epsilon;
}
int Testnode(double x[2], double cx, double cy)
{
    int flag;
    if (sqrt((x[0] - cx) * (x[0] - cx) + (x[1] - cy) * (x[1] - cy)) - cir_r <= 0)
        flag = 0;
    else
        flag = 1;
    return flag;
}
double lengthratio(double x1[2], double x2[2], double cx, double cy)
{
    double q1, q2, a0, b0, c0;
    a0 = (x2[0] - x1[0]) * (x2[0] - x1[0]) + (x2[1] - x1[1]) * (x2[1] - x1[1]);
    b0 = 2 * (x2[0] - x1[0]) * (x1[0] - cx) + 2 * (x2[1] - x1[1]) * (x1[1] - cy);
    c0 = (x1[0] - cx) * (x1[0] - cx) + (x1[1] - cy) * (x1[1] - cy) - cir_r * cir_r;
    q1 = (-b0 + sqrt(b0 * b0 - 4 * a0 * c0)) / 2.0 / a0;
    q2 = (-b0 - sqrt(b0 * b0 - 4 * a0 * c0)) / 2.0 / a0;
    if (q1 > 0 && q1 <= 1)
        return q1;
    else
        return q2;
}
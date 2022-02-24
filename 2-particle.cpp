#include "2-particle.h"

//cd F:\0two_particles
//icl/fast -o p 2-particle.cpp solidrate.cpp

double sum,mm0,sum0;

void init()                           //初始化
{
	x[0][0]=99.9;x[0][1]=700;
	x[1][0]=100.1;x[1][1]=650;
    dx = 1;
    dy = 1;
    Lx = dx * double(NX);
    Ly = dy * double(NY);
    dt = 1;                     //时间步长
    coefficient = dx * dx / dt; //计算力和扭矩的系数
    c = dx / dt;                //格子速度
    /*for (i = 0; i < 9; i++)     //D2Q9速度模型*格子速度
        for (j = 1; j < 2; j++)
            e[i][j] *= c;*/
    rho0 = 1.0;                 //格子密度
    tau_f = 0.6;                //松弛时间
    //tau_f = tau_f / dt;       //无量纲松弛时间
    lx = 0.02;             //实际宽度
    Ch = lx / Lx;           //长度系数
    density = 1000;         //实际密度
    Cden = density / rho0;  //密度系数
    viscosity = 1e-5;       //实际粘度
    Ct = (tau_f - 0.5) * Ch * Ch / (3 * viscosity);     //时间系数
    Cf = Cden * Ch / (Ct * Ct);                         //中间系数
    f_ = density * 9.81 / Cf;
    g = f_ / rho0;                                      //重力加速度
    niu = (tau_f - 0.5) / 3;                            //格子粘度

    W_ = Lx / (2 * cir_r);      //通道宽度与颗粒直径比      
    //经验解
    K = 1 / (log(W_) - 0.9175 + 1.7244 * pow(W_, -2) - 1.7302 * pow(W_, -4) + 2.4056 * pow(W_, -6) - 4.5913 * pow(W_, -8));
    Up = 4 * cir_r * cir_r * (rho0 - rhos) * g / (16 * K * niu);        //颗粒稳定沉降的经验速度
    Re = rho0 * 2 * fabs(Up) * cir_r / niu;      //雷诺数
    std::cout << "雷诺数：Re= " << Re << endl;
    for (a = 0; a < N; a++) {
        velocity_central[a][0] = 0;        //颗粒初速度
        velocity_central[a][1] = 0;
        angular[a] = 0;      //初角速度
        angle[a] = 0;        //初始角度
    }
	
    mass = PI * cir_r * cir_r * rhos;   //颗粒质量
    I = 0.5 * mass * cir_r * cir_r;     //颗粒转动惯量

    threshold = dx;
    len = 2 * cir_r + threshold;
    ep = ep_ = 0.01;
    ew = ew_ = ep / 2;
    leftw = 0 - cir_r; rightw = NX + cir_r; upw = NY + cir_r; downw = 0 - cir_r;

    for (i = 0; i <= NX; i++)
        for (j = 0; j <= NY; j++)
        {
            u[i][j][0] = 0;       //流体速度初始化
            u[i][j][1] = 0;
            rho[i][j] = rho0;
            vel_particle[i][j][0] = 0;       //颗粒合成速度（平移速度+旋转速度）初始化
            vel_particle[i][j][1] = 0;
            for (k = 0; k < Q; k++)
                f[i][j][k] = feq(k, rho[i][j], u[i][j]);
            rate[i][j] = 0;
            for (a = 0; a < N; a++) {
                rate0[a][i][j] = solid_rate(i, j, x[a][0], x[a][1]);          //初始固体比率
                rate[i][j] += rate0[a][i][j];
                B0[a][i][j] = rate0[a][i][j] * (tau_f - 0.5) / (1 - rate0[a][i][j] + tau_f - 0.5);
            }
            B[i][j] = rate[i][j] * (tau_f - 0.5) / (1 - rate[i][j] + tau_f - 0.5);
        }
    mm0 = (-1 + rho0 / rhos) * g * mass;
}

int main()
{
    using namespace std;  
    std::cout << "请指定颗粒密度：rhos= ";
    std::cin >> rhos;               //颗粒密度
    init();
    /*将所要监测的数据输出到一个文件*/
    ostringstream outdata;
    outdata << "data of rhos=" << rhos<<" ep=" <<ep<< ".txt";
    ofstream out(outdata.str().c_str());
    out << "Title=\"other data of particle\"\n"
        << "rhos=" << rhos << ", Up= " << Up << ", Re= " << Re << ", Ct=" << Ct << "\n"
        //<<"cir_x="<<cir_x<<", cir_y="<<cir_y<<", r="<<cir_r<<"\n"
        << "VARIABLES=\"Time step\",\"颗粒1横坐标\",\"颗粒1纵坐标\",\"2横坐标\",\"2纵坐标\",\"1水平速度\",\"1竖直速度\",\" 2水平速度\",\" 2竖直速度\"\n" << endl;

    for ( n = 0;; n++)
    {
        evolution();
        sum = 0;
        sum0 = 0;
        for(i=0;i<=NX;i++)
            for (j = 0; j <= NY; j++) {
                sum0 += rate0[0][i][j];
                sum += rate[i][j];

            }
        sum0 = Ff[1][1] + Fp_sum[1][1] + Fw[1][1] + mm0;
		if(n%50==0){
			out << setw(6) << n << "  " << setprecision(6) << x[0][0] << "  " << x[0][1] << "  " << x[1][0] << "  " << x[1][1] << "  "
				<<velocity_central[0][0]<< "  "<<velocity_central[0][1]<<"  "<<velocity_central[1][0]<<" "<< velocity_central[1][1]<< endl;
		}
        if (n % 200 == 0)
        {
            Error();
            cout << "The" << n << "th computation result:" << endl;
            cout << "The u,v of point(NX/2,NY/2) is: " << setprecision(6)
                << u[NX / 2][NY / 2][0] << ", " << u[NX / 2][NY / 2][1] << endl;
            cout << "The max relative error of uv is: "
                << setiosflags(ios::scientific) << error << endl;
            output(n);
            if (error > 1) break;
            if (n > 100000) break;        //测试的时间步
			
        }
		if (x[0][1]<30) break;
    }
    out.close();
    return 0;
}


void evolution()	                //演化：碰撞，迁移，宏观量计算，边界处理
{
    //碰撞
    for (i = 0; i <= NX; i++) {
        for (j = 0; j <= NY; j++) {
            for (k = 0; k < Q; k++) {
                om[i][j][k]=f[i][j][oppk[k]]-f[i][j][k]+feq(k,rho[i][j],vel_particle[i][j])-feq(oppk[k],rho[i][j],u[i][j]);         
                F[i][j][k] = f[i][j][k] + (feq(k, rho[i][j], u[i][j]) - f[i][j][k]) * (1 - B[i][j]) / tau_f + B[i][j] * om[i][j][k];//+3*e[k][0]*P0*w[k]  
            }
        }
    }

    //迁移
    for (i = 1; i < NX; i++)
        for (j = 1; j < NY; j++) {
            f[i][j][0] = F[i][j][0];
            f[i][j][1] = F[i - 1][j][1];
            f[i][j][2] = F[i][j - 1][2];
            f[i][j][3] = F[i + 1][j][3];
            f[i][j][4] = F[i][j + 1][4];
            f[i][j][5] = F[i - 1][j - 1][5];
            f[i][j][6] = F[i + 1][j - 1][6];
            f[i][j][7] = F[i + 1][j + 1][7];
            f[i][j][8] = F[i - 1][j + 1][8];
        }

    //宏观量计算
    for (a = 0; a < N; a++) {
        Ff0[a][0] = Ff[a][0]; Ff0[a][1] = Ff[a][1]; Tf0[a] = Tf[a];           //保存上一时间步的颗粒对流体的阻力和转矩
        Ff[a][0] = 0; Ff[a][1] = 0; Tf[a] = 0;                        //阻力和扭矩累加前置零
    }

    aver_rho = 0;                                       //平均密度累加前置零

    for (i = 1; i < NX; i++)
        for (j = 1; j < NY; j++)
        {
            u0[i][j][0] = u[i][j][0]; //保存误差计算中上一时间步速度的数据
            u0[i][j][1] = u[i][j][1];

            rho[i][j] = 0;            //密度和速度累加前置零
            u[i][j][0] = 0;
            u[i][j][1] = 0;
            for (k = 0; k < Q; k++)        //计算密度和速度
            {
                rho[i][j] += f[i][j][k];
                u[i][j][0] += e[k][0] * f[i][j][k];       
                u[i][j][1] += e[k][1] * f[i][j][k];
            }
            u[i][j][0] /= rho[i][j];
            u[i][j][1] /= rho[i][j];
            aver_rho += rho[i][j];

            //计算阻力和扭矩
            OM[i][j][0] = 0;
            OM[i][j][1] = 0;
            if (rate[i][j] != 0) {
                for (k = 0; k < Q; k++) {
                    OM[i][j][0] += om[i][j][k] * e[k][0];
                    OM[i][j][1] += om[i][j][k] * e[k][1];
                }
                for (a = 0; a < N; a++) {
                    Ff[a][0] += B0[a][i][j] * OM[i][j][0] * coefficient;        //这里是颗粒对流体的力还是流体对颗粒的力？
                    Ff[a][1] += B0[a][i][j] * OM[i][j][1] * coefficient;        //co=dx*dx/dt

                    Tf[a] += coefficient * ((i - x[a][0]) * dx * B0[a][i][j] * OM[i][j][1] - (j - x[a][1]) * dx * B0[a][i][j] * OM[i][j][0]);
                }
            }
        }
    /*aver_rho /= (double)((NX - 1) * (NY - 1));         //平均密度
    Cl = 0;//Ff[1]/(aver_rho*U*U*cir_r); //升力系数,暂且不算
    Cd = 0;//Ff[0]/(aver_rho*U*U*cir_r); //阻力系数
    */

    //颗粒间碰撞
    for (pi = 0; pi < N; pi++) {
        Fp[pi][pi][0] = 0; Fp[pi][pi][1] = 0;
        for (pj = pi + 1; pj < N; pj++) {
            if ((fabs(x[pi][0] - x[pj][0]) > len) || (fabs(x[pi][1] - x[pj][1]) > len)) {
                Fp[pi][pj][0] = 0; Fp[pj][pi][0] = 0;
                Fp[pi][pj][1] = 0; Fp[pj][pi][1] = 0;
            }
            else
                partical_collision(pi, pj);
        }
    }

    for (pi = 0; pi < N; pi++) {
        Fp_sum[pi][0] = 0;
        Fp_sum[pi][1] = 0;
        for (pj = 0; pj < N; pj++) {
            Fp_sum[pi][0] += Fp[pi][pj][0];
            Fp_sum[pi][1] += Fp[pi][pj][1];
        }
    }

    //颗粒与墙壁的碰撞
    for (pi = 0; pi < N; pi++) {
        /*d1 = x[pi][0] - leftw; d2 = rightw - x[pi][0];
        d3 = x[pi][1] - downw; d4 = upw - x[pi][1];*/
        d1 = 2*x[pi][0]; d2 = 2*(NX - x[pi][0]);
        d3 = 2*x[pi][1]; d4 = 2*(NY - x[pi][1]);
        Fw[pi][0] = d1 * partical_wall(d1) - d2 * partical_wall(d2);
        Fw[pi][1] = d3 * partical_wall(d3) - d4 * partical_wall(d4);
    }

    //颗粒间+墙壁碰撞力
    for (pi = 0; pi < N; pi++) {
        Fcol[pi][0] = Fp_sum[pi][0] + Fw[pi][0];
        Fcol[pi][1] = Fp_sum[pi][1] + Fw[pi][1];
    }

    for (a = 0; a < N; a++) {
        Ff[a][0] = -Ff[a][0];
        Ff[a][1] = -Ff[a][1];
        Tf[a] = -Tf[a];
        //更新颗粒运动参数
        velocity0_central[a][0] = velocity_central[a][0]; velocity0_central[a][1] = velocity_central[a][1];  angular0[a] = angular[a]; angle0[a] = angle[a];  //储存上一时间步的颗粒速度信息

        velocity_central[a][0] = velocity0_central[a][0] + (Ff0[a][0] + Ff[a][0]) / (2 * mass) * dt+Fcol[a][0]/mass;                            //颗粒水平速度
        velocity_central[a][1] = velocity0_central[a][1] + ((Ff0[a][1] + Ff[a][1]) / (2 * mass) + (-1 + rho0 / rhos) * g+ Fcol[a][1] / mass) * dt; //颗粒竖直速度  密度不更新
        x[a][0] = x[a][0] + (velocity_central[a][0] + velocity0_central[a][0]) / 2 * dt;                                 //颗粒中心坐标
        x[a][1] = x[a][1] + (velocity_central[a][1] + velocity0_central[a][1]) / 2 * dt;
    }

    for (i = 1; i < NX; i++)
        for (j = 1; j < NY; j++) {
            rate[i][j] = 0;
            for (a = 0; a < N; a++) {
                rate0[a][i][j] = solid_rate(i, j,  x[a][0], x[a][1]);          //更新固体比率
                rate[i][j] += rate0[a][i][j];
                B0[a][i][j] = rate0[a][i][j] * (tau_f - 0.5) / (1 - rate0[a][i][j] + tau_f - 0.5);
            }
            B[i][j] = rate[i][j] * (tau_f - 0.5) / (1 - rate[i][j] + tau_f - 0.5);
        }

    for (a = 0; a < N; a++) {
        angular[a] = angular0[a] + (Tf[a] + Tf0[a]) / (2 * I) * dt;     //角速度
        angle[a] = angle0[a] + (angular[a] + angular0[a]) / 2 * dt;     //角度
    }


    for (i = 1; i < NX; i++)                         //更新颗粒的合成速度（平移速度+旋转速度）
        for (j = 1; j < NY; j++) {
            vel_particle[i][j][0] = 0;                 //固体比率为0，颗粒合成速度设为0
            vel_particle[i][j][1] = 0;
            for (a = 0; a < N; a++) {
                if (B0[a][i][j] == 0) {
                    ;
                }
                else {
                    vel_particle[i][j][0] = velocity_central[a][0] - (j - x[a][1]) * dx * angular[a];
                    vel_particle[i][j][1] = velocity_central[a][1] + (i - x[a][0]) * dx * angular[a];
                }
            }
        }



    //边界处理，均采用非平衡外推格式
    for (i = 1; i < NX; i++)	           //上下边界
        for (k = 0; k < Q; k++)
        {
            rho[i][0] = rho[i][1];
            u[i][0][0] = u[i][1][0];
            u[i][0][1] = u[i][1][1];
            f[i][0][k] = feq(k, rho[i][0], u[i][0]) + f[i][1][k] - feq(k, rho[i][1], u[i][1]);

            rho[i][NY] = rho[i][NY - 1];
            u[i][NY][0] = 0;
            u[i][NY][1] = 0;
            f[i][NY][k] = feq(k, rho[i][NY], u[i][NY]) + f[i][NY - 1][k] - feq(k, rho[i][NY - 1], u[i][NY - 1]);
        }

    for (j = 0; j <= NY; j++)	           //左右边界
        for (k = 0; k < Q; k++)
        {
            rho[NX][j] = rho[NX - 1][j];
            f[NX][j][k] = feq(k, rho[NX][j], u[NX][j]) + f[NX - 1][j][k] - feq(k, rho[NX - 1][j], u[NX - 1][j]);

            rho[0][j] = rho[1][j];
            f[0][j][k] = feq(k, rho[0][j], u[0][j]) + f[1][j][k] - feq(k, rho[1][j], u[1][j]);
        }
}

double feq(int k, double rho, double u[2])   //平衡分布函数
{
    double eu, uv, feq;
    eu = (e[k][0] * u[0] + e[k][1] * u[1]);
    uv = (u[0] * u[0] + u[1] * u[1]);
    feq = w[k] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * uv);
    return feq;
}

void Error()                               //误差计算
{
    double temp1, temp2;
    temp1 = 0;
    temp2 = 0;
    for (i = 1; i < NX; i++)
        for (j = 1; j < NY; j++)
        {
            temp1 += (u[i][j][0] - u0[i][j][0]) * (u[i][j][0] - u0[i][j][0]) + (u[i][j][1] - u0[i][j][1]) * (u[i][j][1] - u0[i][j][1]);
            temp2 += u[i][j][0] * u[i][j][0] + u[i][j][1] * u[i][j][1];
        }
    temp1 = sqrt(temp1);
    temp2 = sqrt(temp2);
    error = temp1 / (temp2 + 1e-30);
}

void output(int m)                          //输出数据
{
    ostringstream name;
    name << "rhos=" << rhos << " ep" << ep<<"_"<< m << ".dat";
    ofstream out(name.str().c_str());
    out << "Title=\"sedimentation of a circular particle\"\n"
        << "VARIABLES=\"X\",\"Y\",\"U\",\"V\"\n"
        << "ZONE T= \"BOX\", I= "
        << NX + 1 << ", J=" << NY + 1 << ", F=POINT" << endl;
    for (j = 0; j <= NY; j++)
        for (i = 0; i <= NX; i++)
        {
            out << i * dx << " " << j * dy << " "
                << u[i][j][0] << " " << u[i][j][1] << endl;
        }
    out.close();
}

void partical_collision(int i, int j)
{
    double para, d;

    d = sqrt((x[i][0] - x[j][0]) * (x[i][0] - x[j][0]) + (x[i][1] - x[j][1]) * (x[i][1] - x[j][1]));

    if (d > len)
        para = 0;
    else if (d <= 2 * cir_r)
        para = (2 * cir_r - d) / (d*ep_);
    else
        para = (len - d) * (len - d) / (d*ep);	//threshold为1，否则要除以平方

    Fp[i][j][0] = para * (x[i][0] - x[j][0]);
    Fp[i][j][1] = para * (x[i][1] - x[j][1]);

    Fp[j][i][0] = -Fp[i][j][0];
    Fp[j][i][1] = -Fp[i][j][1];

}
double partical_wall(double d)
{
    double para;

    if (d > len) para = 0;
    else if (d <= 2 * cir_r) para = (2 * cir_r - d) / (d*ew_);
    else para = (len - d) * (len - d) / (d*ew);

    return para;
}
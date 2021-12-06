#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include <algorithm>

void init();	//初始化函数
double feq(int k,double rho, double u[2]);		//平衡态分布函数
void evolution();			//演化函数
void output(int m);			//输出
void output_extra();
void Error();				//误差控制
void output0(int m);

using namespace std;
const int Q=9;          //D2Q9模型
const double lx=1;	    //y实际长度
const double ly=0.1;	//y实际长度
const int NY=20;        //y方向网格数
const int NX=NY*10;     //x方向网格数
double U;				//泊肃叶最大流速

double e[Q][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};	//离散速度模型(D2Q9模型)
double w[Q]={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};		//平衡态分布函数的权重
double rho[NX+1][NY+1],p[NX+1][NY+1],u[NX+1][NY+1][2],u0[NX+1][NY+1][2],f[NX+1][NY+1][Q];
double F[NX+1][NY+1][Q];
int i,j,k,n;
double c,Re,dx,dy,Lx,Ly,dt,rho0,P0,tau_f,niu,error,hh;
double Ch,density,Cden,viscosity,Ct,Cf,f_,g,gradp;
double rho_in,rho_out;

int main()
{
    using namespace std;
    init();				//初始化
    for(n=0;;n++)
    {
        evolution();	//分布函数的演化(碰撞、迁移、宏观量计算和边界处理)
        if(n%100==0)
        {
            Error();	//计算相邻时间步的流场速度的迭代误差
			//控制台屏幕显示
            cout<<"The"<<n<<"th computation result:"<<endl;
            cout<<"The u,v of point(NX/2,NY/2) is: "<<setprecision(6)
                <<u[NX/2][NY/2][0]<<", "<<u[NX/2][NY/2][1]<<endl;		
            cout<<"The max relative error of uv is: "
                <<setiosflags(ios::scientific)<<error<<endl;
				
			if(error>1) break;
            if(n>=1000)
            {
                if(n%1000==0) output(n);
                if(error<1.0e-8){		//收敛误差
                    output(n);			//输出流场的速度
					output0(n);			//输出部分流场的数据，可在matlab进行可视化(可以不用)
					output_extra();		//输出截面速度分布的数值解和精确解
                    break;
                }
            }
            if(n>100000) break;			//最大迭代步
        }
    }
    return 0;
}

void init()
{
    dx=lx/NX; //空间步长
    dy=ly/NY;
	dt=dx;	//时间步长
	c=dx/dt;//格子速度
	Re=0.1;
	tau_f=1; //无量纲松弛时间
	niu=(tau_f-0.5)/3*dt;	//格子粘度
	U=Re*niu/NY;			//特征速度
	rho0=2.7;				//流场的初始化密度
	rho_out=2.6996;			//出口边界的密度
	rho_in=2.7004;			//入口边界的密度  	压力梯度通过入出口的密度差实现

	//格子单位到物理单位的转换系数
    Ch = ly / NY;           //长度系数
    density = 1000;         //实际密度
    Cden = density / rho0;  //密度系数
    viscosity = 1e-6;       //实际粘度
    Ct = (tau_f - 0.5) * Ch * Ch / (3 * viscosity);     //时间系数
    Cf = Cden * Ch / (Ct * Ct);                         //中间系数
    f_ = density * 9.81 / Cf;
    g = f_ / rho0;                                     //格子重力加速度

    std::cout<<"tau_f= "<<tau_f<<endl;	//控制台屏幕显示

    for(i=0;i<=NX;i++)   //初始化
        for(j=0;j<=NY;j++)
        {
            u[i][j][0]=0;
            u[i][j][1]=0;	//流场速度
			
            rho[i][j]=rho0;	//流场密度
			p[i][j]=rho[i][j]*c*c/3;	//流场压力
        }
}

void evolution()	//计算平衡态分布函数
{
    //碰撞
    for(i=0;i<=NX;i++){
        for(j=0;j<=NY;j++){
            for(k=0;k<Q;k++){
                F[i][j][k]=f[i][j][k]+(feq(k,rho[i][j],u[i][j])-f[i][j][k])/tau_f;		
            }
        }
    }
    //迁移
    for(i=0;i<=NX;i++)
        for(j=0;j<=NY;j++)
            f[i][j][0]=F[i][j][0];
    for(j=0;j<=NY;j++){
        for(i=1;i<=NX;i++)
            f[i][j][1]=F[i-1][j][1];//右
        for(i=0;i<NX;i++)
            f[i][j][3]=F[i+1][j][3];//左
    }
    for(i=0;i<=NX;i++){
        for(j=1;j<=NY;j++)
            f[i][j][2]=F[i][j-1][2];//上
        for(j=0;j<NY;j++)
            f[i][j][4]=F[i][j+1][4];//下
    }
    for(j=1;j<=NY;j++){
        for(i=1;i<=NX;i++)
            f[i][j][5]=F[i-1][j-1][5];//右上
        for(i=0;i<NX;i++)
            f[i][j][6]=F[i+1][j-1][6];//左上
    }
    for(j=0;j<NY;j++){
        for(i=1;i<=NX;i++)
            f[i][j][8]=F[i-1][j+1][8];//右下
        for(i=0;i<NX;i++)
            f[i][j][7]=F[i+1][j+1][7];//左下
    }
	
	//宏观量计算
    for(i=1;i<NX;i++)	
        for(j=1;j<NY;j++)
        {
            u0[i][j][0]=u[i][j][0];			//累加前置零
            u0[i][j][1]=u[i][j][1];
			rho[i][j]=0;
            u[i][j][0]=0;
            u[i][j][1]=0;
            for(k=0;k<Q;k++)
            {
				rho[i][j]+=f[i][j][k];		//流场密度
                u[i][j][0]+=e[k][0]*f[i][j][k];
                u[i][j][1]+=e[k][1]*f[i][j][k];
            }
			p[i][j]=rho[i][j]*c*c/3;		//流场压力 
            u[i][j][0]/=rho[i][j];
            u[i][j][1]/=rho[i][j];			//流场速度
        }

	#if 1//边界处理
	for(j=1;j<NY;j++)		//左右边界  非平衡反弹格式(Zou-He格式)		
		{
			rho[NX][j]=rho_out;
			p[NX][j]=rho[NX][j]/3;
			u[NX][j][0]=(2*(f[NX][j][1]+f[NX][j][5]+f[NX][j][8])+f[NX][j][2]+f[NX][j][4]+f[NX][j][0])/rho[NX][j]-1;
			u[NX][j][1]=0;
			f[NX][j][3]=f[NX][j][1]-2/3*rho[NX][j]*u[NX][j][0];
			f[NX][j][7]=f[NX][j][5]+1/2*(f[NX][j][2]-f[NX][j][4])-1/6*rho[NX][j]*u[NX][j][0];
			f[NX][j][6]=f[NX][j][8]-1/2*(f[NX][j][2]-f[NX][j][4])-1/6*rho[NX][j]*u[NX][j][0];


			rho[0][j]=rho_in;
			p[0][j]=rho[0][j]/3;
			u[0][j][0]=1-(f[0][j][2]+f[0][j][4]+f[0][j][0]+2*(f[0][j][3]+f[0][j][6]+f[0][j][7]))/rho[0][j];
			u[0][j][1]=0;
			f[0][j][1]=f[0][j][3]+2/3*rho[0][j]*u[0][j][0];
			f[0][j][5]=f[0][j][7]-1/2*(f[0][j][2]-f[0][j][4])+1/6*rho[0][j]*u[0][j][0];
			f[0][j][8]=f[0][j][6]+1/2*(f[0][j][2]-f[0][j][4])+1/6*rho[0][j]*u[0][j][0];
		}
	#endif	
	#if 0
        for(j=1; j<NY; j++)   //左右边界备选 非平衡外推格式
            for(k=0; k<Q; k++)   
            {   
                u[NX][j][0]=u[NX-1][j][0];				
                u[NX][j][1]=u[NX-1][j][1];   
                rho[NX][j]=rho_out;
				p[NX][j]=rho[NX][j]*c*c/3;
                f[NX][j][k]=feq(k,rho[NX][j],u[NX][j])+f[NX-1][j][k]-feq(k,rho[NX-1][j],u[NX-1][j]);   
				
                u[0][j][0]=u[1][j][0]; 
                u[0][j][1]=u[1][j][1];   
                rho[0][j]=rho_in; 
				p[0][j]=rho[0][j]*c*c/3;				
                f[0][j][k]=feq(k,rho[0][j],u[0][j])+f[1][j][k]-feq(k,rho[1][j],u[1][j]);
				
            }  
	#endif
	for(i=0;i<=NX;i++)	//上下边界	非平衡外推格式
		for(k=0;k<Q;k++)
		{
			rho[i][0]=rho[i][1];
			f[i][0][k]=feq(k,rho[i][0],u[i][0])+f[i][1][k]-feq(k,rho[i][1],u[i][1]);

			rho[i][NY]=rho[i][NY-1];
			f[i][NY][k]=feq(k,rho[i][NY],u[i][NY])+f[i][NY-1][k]-feq(k,rho[i][NY-1],u[i][NY-1]);
		}
}

double feq(int k,double rho, double u[2])		//平衡态分布函数(这里采用的是Qian提出的格式)
{
    double eu,uv,feq;
    eu=(e[k][0]*u[0]+e[k][1]*u[1]);
    uv=(u[0]*u[0]+u[1]*u[1]);
    feq=w[k]*rho*(1.0+3.0*eu+4.5*eu*eu-1.5*uv);
	
    return feq;
}

void output(int m)  //输出
{
    ostringstream name;
    name<<"pos_"<<m<<".dat";
    ofstream out(name.str().c_str());
    out<<"Title=\"LBM Lid Driven Flow\"\n"
        <<"VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"p\"\n"
        <<"ZONE T= \"BOX\", I= "
        <<NX-1<<", J="<<NY-1<<", F=POINT"<<endl;
    for(j=1; j<NY; j++)
        for(i=1; i<NX; i++)
        {
            out<<double(i)/NX*NX/NY<<" "<<double(j)/NY<<" "
                <<u[i][j][0]<<" "<<u[i][j][1]<<" "<<p[i][j]<<endl;	//输出流场的速度与压力
        }
}
void output0(int m)  //输出局部流场参数(非必要)
{
    ostringstream name;
    name<<"matlab"<<m<<".dat";
	ofstream out(name.str().c_str());
    for(j=0; j<=NY; j++)
        for(i=150; i<170; i++)
        {
            out<<double(i)/NX*NX/NY<<" "<<double(j)/NY<<" "		//输出局部的流场参数
                <<u[i][j][0]<<endl;
        }
}
void output_extra()		//输出截面速度分布的数值解和精确解(方便进行数据处理，非必要)
{
	double u_real[NY+1],Err[NY+1],maxErr,umax;
	int test=3.0/4*NX;  //水平格子坐标(截面的水平位置)
	
	umax=(rho_in-rho_out)/3/lx*ly*ly/(8*rho0*niu);	//理论最大流速
	for(i=0;i<=NY;i++){
		u_real[i]=umax*4*(double (i)/NY)*(1-double (i)/NY);	//理论截面流速分布
		Err[i]=sqrt((u_real[i]-u[test][i][0])*(u_real[i]-u[test][i][0])+(0-u[test][i][1])*(0-u[test][i][1]))/umax;		//相对速度误差
	}
	maxErr = *max_element(Err,Err+NY+1);
	Re=umax*ly/niu;
    ostringstream outdata;
    outdata << "tau_f=" << tau_f << ",NY=" << NY << ".txt";
    ofstream out(outdata.str().c_str());
    out << "Title=\"泊肃叶\"\n"
        << "理论最大速度Umax=" << umax << ",最大相对速度误差maxerr= " << maxErr << ", 松弛时间tau_f= " << tau_f  << ",分子步长dx=" << dx << "\n"
        <<"Re="<< Re <<", 格子重力加速度g="<< g <<", 压力梯度="<<(rho_in-rho_out)/3/lx<<",粘度niu="<<niu<<"\n"
        << "VARIABLES=\"格子纵坐标\",\"x/H=150/40\",\" 理论解\",\"  \"\n" << endl;
	for(i=0;i<=NY;i++){
        out << setw(4) << i << "  " << setprecision(8) << u[test][i][0] << "  " << u_real[i] << endl;	//输出截面流速的数值解和理论解
	}
	out.close();
    ostringstream outdata_ex;
    outdata_ex << "tau_f NY and maxErr"<< ".txt";
    ofstream app(outdata_ex.str().c_str(),ofstream::app);
	app <<tau_f<<"  "<<NY<<"  "<<maxErr<<"  "<<n<<endl;		//输出最大相对误差
	app.close();
}
void Error()  //误差控制
{
    double temp1,temp2;
    temp1=0;
    temp2=0;
    for(i=1; i<NX; i++)
        for(j=1;j<NY;j++)
        {
            temp1 += (u[i][j][0]-u0[i][j][0])*(u[i][j][0]-u0[i][j][0])+(u[i][j][1]-u0[i][j][1])*(u[i][j][1]-u0[i][j][1]);
            temp2 += u[i][j][0]*u[i][j][0]+u[i][j][1]*u[i][j][1];
        }
        temp1=sqrt(temp1);
        temp2=sqrt(temp2);
        error=temp1/(temp2+1e-30);
}
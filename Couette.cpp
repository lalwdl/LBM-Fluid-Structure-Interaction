
#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include <algorithm>

void init();
double feq(int k,double rho, double u[2]);
void evolution();
void output(int m);
void output_extra();
void Error();

using namespace std;
const int Q=9;          //D2Q9模型
const double lx=1;	//y实际长度
const double ly=0.1;	//y实际长度
const int NY=20;       //y方向
const int NX=NY*10;       //x方向
double U;		//最大流速

double e[Q][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};
int k_[9]={0,1,4,3,2,7,8,5,6};
double w[Q]={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};
double rho[NX+1][NY+1],p[NX+1][NY+1],u[NX+1][NY+1][2],u0[NX+1][NY+1][2],f[NX+1][NY+1][Q];
double F[NX+1][NY+1][Q];
int i,j,k,n;
double c,Re,dx,dy,Lx,Ly,dt,rho0,P0,tau_f,niu,error,hh;
double Ch,density,Cden,viscosity,Ct,Cf,f_,g,gradp;
double rho_in,rho_out,umax,P;

//cd/d E:\0 poiseuille\success
//icl/fast -o p main.cpp
int main()
{
    using namespace std;
    init();
    for(n=0;;n++)
    {
        evolution();
        if(n%100==0)
        {
            Error();
            cout<<"The"<<n<<"th computation result:"<<endl;
            cout<<"The u,v of point(NX/2,NY/2) is: "<<setprecision(6)
                <<u[NX/2][NY/2][0]<<", "<<u[NX/2][NY/2][1]<<endl;
            cout<<"The max relative error of uv is: "
                <<setiosflags(ios::scientific)<<error<<endl;
			if(error>1) break;
            if(n>=1000)
            {
                if(n%1000==0) output(n);
				if(error>1) break;
                if(error<1.0e-8){
                    output(n);
					output_extra();
                    break;
                }
            }
            if(n>50000) break;
        }
    }
    return 0;
}


void init()
{
    dx=lx/NX;
    dy=ly/NY;
	dt=dx;
	c=dx/dt;
	Re=0.1;
	tau_f=1;
	niu=(tau_f-0.5)/3*dt;
	U=Re*niu/NY;
	/*gradp=-U*8*niu/(ly*ly);*/
	rho0=2.7;
	rho_out=2.6996;
	rho_in=2.7004;
	umax=(rho_in-rho_out)/3/lx*ly*ly/(8*rho0*niu);
	P=2;//压力梯度的无量纲数
	U=4*umax/P;
	//U=0.001;
	//umax=U;
	Re=U*lx/niu;
	
	/*tau_f=5;
	niu=(tau_f-0.5)/3;
	Re=0.1*Lx/niu;*/
	
    /*Re=100;
    niu=0.5*U*Ly/Re;
    tau_f=3.0*niu+0.5;*/
	
    Ch = ly / NY;           //长度系数
    density = 1000;         //实际密度
    Cden = density / rho0;  //密度系数
    viscosity = 1e-6;       //实际粘度
    Ct = (tau_f - 0.5) * Ch * Ch / (3 * viscosity);     //时间系数
    Cf = Cden * Ch / (Ct * Ct);                         //中间系数
    f_ = density * 9.81 / Cf;
    g = f_ / rho0;                                     //重力加速度
   // P0=0.00001;
	
	/*rho0=1.0;   
	gradp=16*U*niu*rho0/((NY-1)*(NY-1));
	rho_out=1.0;
	rho_in=rho_out+3*gradp*NX;*/
	/*rho_in=1.05;
	rho_out=1;
	rho0=(rho_in+rho_out)/2;
	niu=0.1;
	
	tau_f=1.09;//3*niu+0.5;
	
	dt=dx*dx*(tau_f-0.5)/(3*niu);
    c=dx/dt;
	U=(rho_in-rho_out)*c*c/3*NY*NY/(16*niu*rho0*NX);
	Re=U*NY/(2*niu);
	
	for(i=0;i<9;i++)
		for(j=0;j<2;j++){
			e[i][j]*=c;
		}*/
	//rho0=1.0;

	
    std::cout<<"tau_f= "<<tau_f<<endl;

    for(i=0;i<=NX;i++)   //初始化
        for(j=0;j<=NY;j++)
        {
            u[i][j][0]=0;
            u[i][j][1]=0;
			
			u[i][NY][0]=U;
			
            rho[i][j]=rho0;
			p[i][j]=rho[i][j]*c*c/3;
			for(k=0; k<Q;k++)
                f[i][j][k]=feq(k,rho[i][j],u[i][j]);
        }
	/*for(j=1;j<NY;j++){
		rho[0][j]=rho_in;rho[NX][j]=rho_out;
		//u[0][j][0]=4*U*double(j)/NY*(1-double(j)/NY);
	}*/
}


void evolution()	//计算平衡态分布函数
{
    //碰撞
    for(i=0;i<=NX;i++){
        for(j=0;j<=NY;j++){
            for(k=0;k<Q;k++){
                F[i][j][k]=f[i][j][k]+(feq(k,rho[i][j],u[i][j])-f[i][j][k])/tau_f;//-3*e[k][0]*gradp*0.0164*w[k];
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

    for(i=1;i<NX;i++)	//计算宏观量速度
        for(j=1;j<NY;j++)
        {
            u0[i][j][0]=u[i][j][0];
            u0[i][j][1]=u[i][j][1];
			rho[i][j]=0;
            u[i][j][0]=0;
            u[i][j][1]=0;
            for(k=0;k<Q;k++)
            {
				rho[i][j]+=f[i][j][k];
                u[i][j][0]+=e[k][0]*f[i][j][k];
                u[i][j][1]+=e[k][1]*f[i][j][k];
            }
			p[i][j]=rho[i][j]*c*c/3;		//未更新边界出口的压力
            u[i][j][0]/=rho[i][j];
            u[i][j][1]/=rho[i][j];
        }

	/*for(j=1;j<NY;j++)	//左右边界
		for(k=0;k<Q;k++)
		{
			u[NX][j][0]=u[NX-1][j][0];
			rho[NX][j]=rho[NX-1][j];
			f[NX][j][k]=feq(k,rho[NX][j],u[NX][j])+f[NX-1][j][k]-feq(k,rho[NX-1][j],u[NX-1][j]);

			u[0][j][0]=U*(1-(j-20)*(j-20)/400);
			rho[0][j]=rho[1][j];
			f[0][j][k]=feq(k,rho[0][j],u[0][j])+f[1][j][k]-feq(k,rho[1][j],u[1][j]);
		}*/
	#if 1//边界处理，采用非平衡态外推格式
	for(j=1;j<NY;j++)	//左右边界		
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
        for(j=1; j<NY; j++)   
            for(k=0; k<Q; k++)   
            {   
                u[NX][j][0]=u[NX-1][j][0];				
                u[NX][j][1]=u[NX-1][j][1];   
                rho[NX][j]=rho_out;//rho[NX][j]=rho[NX-1][j];
				p[NX][j]=rho[NX][j]*c*c/3;
                f[NX][j][k]=feq(k,rho[NX][j],u[NX][j])+f[NX-1][j][k]-feq(k,rho[NX-1][j],u[NX-1][j]);   
                u[0][j][0]=u[1][j][0]; 
				//u[0][j][0]=4*U*double(j)/NY*(1-double(j)/NY);
                u[0][j][1]=u[1][j][1];   
                rho[0][j]=rho_in;//rho[0][j]=rho[1][j]; 
				p[0][j]=rho[0][j]*c*c/3;				
                f[0][j][k]=feq(k,rho[0][j],u[0][j])+f[1][j][k]-feq(k,rho[1][j],u[1][j]);
				
				/*f[0][j][1]=f[NX][j][1];
				f[0][j][5]=f[NX][j][5];
				f[0][j][8]=f[NX][j][8];
				
				f[NX][j][3]=f[0][j][3];
				f[NX][j][6]=f[0][j][6];
				f[NX][j][7]=f[0][j][7];*/
            }  
	#endif
	for(i=0;i<=NX;i++)	//上下边界
		for(k=0;k<Q;k++)
		{
			u[i][NY][0]=U;
			rho[i][0]=rho[i][1];
			f[i][0][k]=feq(k,rho[i][0],u[i][0])+f[i][1][k]-feq(k,rho[i][1],u[i][1]);

			rho[i][NY]=rho[i][NY-1];
			f[i][NY][k]=feq(k,rho[i][NY],u[i][NY])+f[i][NY-1][k]-feq(k,rho[i][NY-1],u[i][NY-1]);
		}

	
}

double feq(int k,double rho, double u[2])
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
                <<u[i][j][0]<<" "<<u[i][j][1]<<" "<<p[i][j]<<endl;
        }
}
void output_extra()
{
	double u_real[NY+1],Err[NY+1],maxErr;
	int test=3.0/4*NX;  //水平格子坐标
	
	for(i=0;i<=NY;i++){
		u_real[i]=U*(P*(double (i)/NY)*(1-double (i)/NY)+ double (i)/NY);
		Err[i]=sqrt((u_real[i]-u[test][i][0])*(u_real[i]-u[test][i][0])+(0-u[test][i][1])*(0-u[test][i][1]))/umax;
	}
	maxErr = *max_element(Err,Err+NY+1);
	Re=umax*ly/niu;
    ostringstream outdata;
    outdata << "tau_f=" << tau_f << ",P=" << P << ".txt";
    ofstream out(outdata.str().c_str());
    out << "Title=\"泊肃叶\"\n"
        << "理论最大速度Umax=" << umax << ",最大相对速度误差maxerr= " << maxErr << ", 松弛时间tau_f= " << tau_f  << ",分子步长dx=" << dx << "\n"
        <<"Re="<< Re <<", 格子重力加速度g="<< g <<", 压力梯度="<<(rho_in-rho_out)/3/lx<<",粘度niu="<<niu<<"\n"
        << "VARIABLES=\"格子纵坐标\",\"x/H=150/40\",\" 理论解\",\"  \"\n" << endl;
	for(i=0;i<=NY;i++){
        out << setw(4) << i << "  " << setprecision(8) << u[test][i][0] << "  " << u_real[i] << endl;
	}
	out.close();
    ostringstream outdata_ex;
    outdata_ex << "tau_f NY and maxErr"<< ".txt";
    ofstream app(outdata_ex.str().c_str(),ofstream::app);
	app <<tau_f<<"  "<<P<<"  "<<maxErr<<"  "<<n<<endl;
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


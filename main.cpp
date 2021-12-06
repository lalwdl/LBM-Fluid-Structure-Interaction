#include "fuction.h"
#include "fuc_lbm.h"

using namespace std;
const int Q=9;          //D2Q9模型
const int NX=800;       //x方向
const int NY=400;       //y方向
const double U=0.1;     //来流速度
const double cir_x=150; //圆心坐标
const double cir_y=200;
const double cir_r=15;  //半径

int e[Q][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};
int k_[9]={0,3,4,1,2,7,8,5,6};
double w[Q]={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};
double rho[NX+1][NY+1],u[NX+1][NY+1][2],u0[NX+1][NY+1][2],f[NX+1][NY+1][Q],u_[2]={0};  //圆的速度
double F[NX+1][NY+1][Q],rate[NX+1][NY+1],B[NX+1][NY+1],om[NX+1][NY+1][Q];
double OM[NX+1][NY+1][2]={0},Ff[2],Tf;
int i,j,k,n;
double aver_rho,Cl,Cd;
double c,Re,dx,dy,Lx,Ly,dt,rho0,rhoo,P0,tau_f,niu,error,hh;

int main()
{
    using namespace std;
    std::cout<<"请指定雷诺数：Re= ";
    std::cin>>Re;
    ostringstream outdata;
    outdata<<"data of Re="<<Re<<".txt";
    ofstream out(outdata.str().c_str());
    out<<"Title=\"other data of cylinder circular flow\"\n"
        <<"Re="<<Re<<", NX= "<<NX<<", NY="<<NY<<"\n"
        <<"cir_x="<<cir_x<<", cir_y="<<cir_y<<", r="<<cir_r<<"\n"
        <<"VARIABLES=\"Time step\",\"Cl\",\"Cd\",\"Tf\",\"aver_rho\"\n"<<endl;

    init();
    for(n=0;;n++)
    {
        evolution();
        out<<setw(6)<<n<<"  "<<setw(10)<<Cl<<"  "<<setw(10)<<Cd<<"  "<<setw(10)<<Tf<<"  "<<setw(5)<<aver_rho<<endl;
        if(n%100==0)
        {
            Error();
            cout<<"The"<<n<<"th computation result:"<<endl;
            cout<<"The u,v of point(NX/2,NY/2) is: "<<setprecision(6)
                <<u[NX/2][NY/2][0]<<", "<<u[NX/2][NY/2][1]<<endl;
            cout<<"The max relative error of uv is: "
                <<setiosflags(ios::scientific)<<error<<endl;
            if(n>=1000)
            {
                if(n%1000==0) output(n);
                if(error<1.0e-7){
                    output(n);
                    break;
                }
                if(error>1){
                    cout<<"结果已经发散，请检查程序"<<endl;
                    break;
                }
            }
            if(n>100000){
                cout<<"100000时间步内未收敛于1e-7,已中断程序"<<endl;
                break;
            }
        }
    }
    out.close();
    return 0;
}


void init()                           //初始化
{
    dx=1.0;
    dy=1.0;
    Lx=dx*double(NX);
    Ly=dy*double(NY);
    dt=dx;
    c=dx/dt;
    rho0=1.0;
    niu=2*U*cir_r/Re;
    tau_f=3.0*niu+0.5;
    //P0=0;//可添加的作用力
    std::cout<<"松弛时间：tau_f= "<<tau_f<<endl;

    for(i=0;i<=NX;i++)
        for(j=0;j<=NY;j++)
        {
            u[i][j][0]=0;
            u[i][j][1]=0;
            rho[i][j]=rho0;
            if((j!=0)&&(j!=NY))
                u[0][j][0]=U;
            for(k=0; k<Q;k++)
                f[i][j][k]=feq(k,rho[i][j],u[i][j]);
            rate[i][j]=solid_rate(i,j,1.0,cir_x,cir_y,cir_r);
            B[i][j]=rate[i][j]*(tau_f-0.5)/(1-rate[i][j]+tau_f-0.5);

        }
}


void evolution()	                //演化：碰撞，迁移，宏观量计算，边界处理
{
    //碰撞
    for(i=0;i<=NX;i++){
        for(j=0;j<=NY;j++){
            for(k=0;k<Q;k++){
                om[i][j][k]=f[i][j][k_[k]]-f[i][j][k]+feq(k,rho[i][j],u_)-feq(k_[k],rho[i][j],u[i][j]);
                //F[i][j][k]=f[i][j][k]+(feq(k,rho[i][j],u[i][j])-f[i][j][k])/tau_f;
                F[i][j][k]=f[i][j][k]+(feq(k,rho[i][j],u[i][j])-f[i][j][k])*(1-B[i][j])/tau_f+B[i][j]*om[i][j][k];//+3*e[k][0]*P0*w[k]
            }
        }
    }

    //迁移
    for(i=1;i<NX;i++)
        for(j=1;j<NY;j++){
            f[i][j][0]=F[i][j][0];
            f[i][j][1]=F[i-1][j][1];
            f[i][j][2]=F[i][j-1][2];
            f[i][j][3]=F[i+1][j][3];
            f[i][j][4]=F[i][j+1][4];
            f[i][j][5]=F[i-1][j-1][5];
            f[i][j][6]=F[i+1][j-1][6];
            f[i][j][7]=F[i+1][j+1][7];
            f[i][j][8]=F[i-1][j+1][8];
        }

    //宏观量计算
    Ff[0]=0;                        //阻力和扭矩累加前置零
    Ff[1]=0;
    Tf=0;
    aver_rho=0;
    for(i=1;i<NX;i++)
        for(j=1;j<NY;j++)
        {
            u0[i][j][0]=u[i][j][0]; //保存误差计算中上一时间步的数据
            u0[i][j][1]=u[i][j][1];

            rho[i][j]=0;            //密度和速度累加前置零
            u[i][j][0]=0;
            u[i][j][1]=0;
            for(k=0;k<Q;k++)        //计算密度和速度
            {
                rho[i][j]+=f[i][j][k];
                u[i][j][0]+=e[k][0]*f[i][j][k];
                u[i][j][1]+=e[k][1]*f[i][j][k];
            }
            u[i][j][0]/=rho[i][j];
            u[i][j][1]/=rho[i][j];
            aver_rho+=rho[i][j];

            if(rate[i][j]!=0){      //计算阻力和扭矩
                OM[i][j][0]=0;
                OM[i][j][1]=0;
                for(k=0;k<Q;k++){
                    OM[i][j][0]+=om[i][j][k]*e[k][0];
                    OM[i][j][1]+=om[i][j][k]*e[k][1];
                }
                Ff[0]+=B[i][j]*OM[i][j][0];
                Ff[1]+=B[i][j]*OM[i][j][1];        //若dx*dx/dt不为1；则需要乘一个额外的系数，以下扭矩计算亦同

                Tf+=(i-cir_x)*B[i][j]*OM[i][j][1]-(j-cir_y)*B[i][j]*OM[i][j][0];
            }
        }
    aver_rho/=((NX-1)*(NY-1));     //平均密度
    Cl=Ff[1]/(aver_rho*U*U*cir_r); //升力系数
    Cd=Ff[0]/(aver_rho*U*U*cir_r); //阻力系数


    //边界处理，均采用非平衡外推格式
    for(j=1;j<NY;j++)	           //左右边界
        for(k=0;k<Q;k++)
        {
            u[NX][j][0]=u[NX-1][j][0];
            u[NX][j][1]=u[NX-1][j][1];
            rho[NX][j]=rho[NX-1][j];
            f[NX][j][k]=feq(k,rho[NX][j],u[NX][j])+f[NX-1][j][k]-feq(k,rho[NX-1][j],u[NX-1][j]);

            u[0][j][0]=U;
            rho[0][j]=rho[1][j];
            f[0][j][k]=feq(k,rho[0][j],u[0][j])+f[1][j][k]-feq(k,rho[1][j],u[1][j]);
        }

    for(i=0;i<=NX;i++)	           //上下边界
        for(k=0;k<Q;k++)
        {
            rho[i][0]=rho[i][1];
            u[i][0][0]=u[i][1][0];
            u[i][0][1]=u[i][1][1];
            f[i][0][k]=feq(k,rho[i][0],u[i][0])+f[i][1][k]-feq(k,rho[i][1],u[i][1]);

            rho[i][NY]=rho[i][NY-1];
            u[i][NY][0]=u[i][NY-1][0];
            u[i][NY][1]=u[i][NY-1][1];
            f[i][NY][k]=feq(k,rho[i][NY],u[i][NY])+f[i][NY-1][k]-feq(k,rho[i][NY-1],u[i][NY-1]);
        }
}

double feq(int k,double rho, double u[2])   //平衡分布函数
{
    double eu,uv,feq;
    eu=(e[k][0]*u[0]+e[k][1]*u[1]);
    uv=(u[0]*u[0]+u[1]*u[1]);
    feq=w[k]*rho*(1.0+3.0*eu+4.5*eu*eu-1.5*uv);
    return feq;
}

void Error()                               //误差计算
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

void output(int m)                          //输出数据
{
    ostringstream name;
    name<<"Re"<<Re<<" cylinder_"<<m<<".dat";
    ofstream out(name.str().c_str());
    out<<"Title=\"circular cylinder flow\"\n"
        <<"VARIABLES=\"X\",\"Y\",\"U\",\"V\"\n"
        <<"ZONE T= \"BOX\", I= "
        <<NX+1<<", J="<<NY+1<<", F=POINT"<<endl;
    for(j=0; j<=NY; j++)
        for(i=0; i<=NX; i++)
        {
            out<<double(i)/Lx*NX/NY<<" "<<double(j)/Ly<<" "
                <<u[i][j][0]<<" "<<u[i][j][1]<<endl;
        }
    out.close();
}

void output_ex1(int m)                      //备用测试输出
{
    ostringstream name;
    name<<"hh"<<m<<".dat";
    ofstream out(name.str().c_str());
    out<<"Title=\"LBM Lid Driven Flow\"\n"
        <<"VARIABLES=\"X\",\"Y\",\"U\",\"V\"\n"
        <<"ZONE T= \"BOX\", I= "
        <<NX+1<<", J="<<NY+1<<", F=POINT"<<endl;

    for(i=0; i<=NX; i++)
        for(j=0; j<=NY; j++)
        {
            out<<double(i)/Lx*NX/NY<<" "<<double(j)/Ly<<" "
                <<rate[i][j]<<" "<<B[i][j]<<endl;
        }
}

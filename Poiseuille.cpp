#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include <algorithm>

void init();	//��ʼ������
double feq(int k,double rho, double u[2]);		//ƽ��̬�ֲ�����
void evolution();			//�ݻ�����
void output(int m);			//���
void output_extra();
void Error();				//������
void output0(int m);

using namespace std;
const int Q=9;          //D2Q9ģ��
const double lx=1;	    //yʵ�ʳ���
const double ly=0.1;	//yʵ�ʳ���
const int NY=20;        //y����������
const int NX=NY*10;     //x����������
double U;				//����Ҷ�������

double e[Q][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};	//��ɢ�ٶ�ģ��(D2Q9ģ��)
double w[Q]={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};		//ƽ��̬�ֲ�������Ȩ��
double rho[NX+1][NY+1],p[NX+1][NY+1],u[NX+1][NY+1][2],u0[NX+1][NY+1][2],f[NX+1][NY+1][Q];
double F[NX+1][NY+1][Q];
int i,j,k,n;
double c,Re,dx,dy,Lx,Ly,dt,rho0,P0,tau_f,niu,error,hh;
double Ch,density,Cden,viscosity,Ct,Cf,f_,g,gradp;
double rho_in,rho_out;

int main()
{
    using namespace std;
    init();				//��ʼ��
    for(n=0;;n++)
    {
        evolution();	//�ֲ��������ݻ�(��ײ��Ǩ�ơ����������ͱ߽紦��)
        if(n%100==0)
        {
            Error();	//��������ʱ�䲽�������ٶȵĵ������
			//����̨��Ļ��ʾ
            cout<<"The"<<n<<"th computation result:"<<endl;
            cout<<"The u,v of point(NX/2,NY/2) is: "<<setprecision(6)
                <<u[NX/2][NY/2][0]<<", "<<u[NX/2][NY/2][1]<<endl;		
            cout<<"The max relative error of uv is: "
                <<setiosflags(ios::scientific)<<error<<endl;
				
			if(error>1) break;
            if(n>=1000)
            {
                if(n%1000==0) output(n);
                if(error<1.0e-8){		//�������
                    output(n);			//����������ٶ�
					output0(n);			//����������������ݣ�����matlab���п��ӻ�(���Բ���)
					output_extra();		//��������ٶȷֲ�����ֵ��;�ȷ��
                    break;
                }
            }
            if(n>100000) break;			//��������
        }
    }
    return 0;
}

void init()
{
    dx=lx/NX; //�ռ䲽��
    dy=ly/NY;
	dt=dx;	//ʱ�䲽��
	c=dx/dt;//�����ٶ�
	Re=0.1;
	tau_f=1; //�������ɳ�ʱ��
	niu=(tau_f-0.5)/3*dt;	//����ճ��
	U=Re*niu/NY;			//�����ٶ�
	rho0=2.7;				//�����ĳ�ʼ���ܶ�
	rho_out=2.6996;			//���ڱ߽���ܶ�
	rho_in=2.7004;			//��ڱ߽���ܶ�  	ѹ���ݶ�ͨ������ڵ��ܶȲ�ʵ��

	//���ӵ�λ������λ��ת��ϵ��
    Ch = ly / NY;           //����ϵ��
    density = 1000;         //ʵ���ܶ�
    Cden = density / rho0;  //�ܶ�ϵ��
    viscosity = 1e-6;       //ʵ��ճ��
    Ct = (tau_f - 0.5) * Ch * Ch / (3 * viscosity);     //ʱ��ϵ��
    Cf = Cden * Ch / (Ct * Ct);                         //�м�ϵ��
    f_ = density * 9.81 / Cf;
    g = f_ / rho0;                                     //�����������ٶ�

    std::cout<<"tau_f= "<<tau_f<<endl;	//����̨��Ļ��ʾ

    for(i=0;i<=NX;i++)   //��ʼ��
        for(j=0;j<=NY;j++)
        {
            u[i][j][0]=0;
            u[i][j][1]=0;	//�����ٶ�
			
            rho[i][j]=rho0;	//�����ܶ�
			p[i][j]=rho[i][j]*c*c/3;	//����ѹ��
        }
}

void evolution()	//����ƽ��̬�ֲ�����
{
    //��ײ
    for(i=0;i<=NX;i++){
        for(j=0;j<=NY;j++){
            for(k=0;k<Q;k++){
                F[i][j][k]=f[i][j][k]+(feq(k,rho[i][j],u[i][j])-f[i][j][k])/tau_f;		
            }
        }
    }
    //Ǩ��
    for(i=0;i<=NX;i++)
        for(j=0;j<=NY;j++)
            f[i][j][0]=F[i][j][0];
    for(j=0;j<=NY;j++){
        for(i=1;i<=NX;i++)
            f[i][j][1]=F[i-1][j][1];//��
        for(i=0;i<NX;i++)
            f[i][j][3]=F[i+1][j][3];//��
    }
    for(i=0;i<=NX;i++){
        for(j=1;j<=NY;j++)
            f[i][j][2]=F[i][j-1][2];//��
        for(j=0;j<NY;j++)
            f[i][j][4]=F[i][j+1][4];//��
    }
    for(j=1;j<=NY;j++){
        for(i=1;i<=NX;i++)
            f[i][j][5]=F[i-1][j-1][5];//����
        for(i=0;i<NX;i++)
            f[i][j][6]=F[i+1][j-1][6];//����
    }
    for(j=0;j<NY;j++){
        for(i=1;i<=NX;i++)
            f[i][j][8]=F[i-1][j+1][8];//����
        for(i=0;i<NX;i++)
            f[i][j][7]=F[i+1][j+1][7];//����
    }
	
	//���������
    for(i=1;i<NX;i++)	
        for(j=1;j<NY;j++)
        {
            u0[i][j][0]=u[i][j][0];			//�ۼ�ǰ����
            u0[i][j][1]=u[i][j][1];
			rho[i][j]=0;
            u[i][j][0]=0;
            u[i][j][1]=0;
            for(k=0;k<Q;k++)
            {
				rho[i][j]+=f[i][j][k];		//�����ܶ�
                u[i][j][0]+=e[k][0]*f[i][j][k];
                u[i][j][1]+=e[k][1]*f[i][j][k];
            }
			p[i][j]=rho[i][j]*c*c/3;		//����ѹ�� 
            u[i][j][0]/=rho[i][j];
            u[i][j][1]/=rho[i][j];			//�����ٶ�
        }

	#if 1//�߽紦��
	for(j=1;j<NY;j++)		//���ұ߽�  ��ƽ�ⷴ����ʽ(Zou-He��ʽ)		
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
        for(j=1; j<NY; j++)   //���ұ߽籸ѡ ��ƽ�����Ƹ�ʽ
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
	for(i=0;i<=NX;i++)	//���±߽�	��ƽ�����Ƹ�ʽ
		for(k=0;k<Q;k++)
		{
			rho[i][0]=rho[i][1];
			f[i][0][k]=feq(k,rho[i][0],u[i][0])+f[i][1][k]-feq(k,rho[i][1],u[i][1]);

			rho[i][NY]=rho[i][NY-1];
			f[i][NY][k]=feq(k,rho[i][NY],u[i][NY])+f[i][NY-1][k]-feq(k,rho[i][NY-1],u[i][NY-1]);
		}
}

double feq(int k,double rho, double u[2])		//ƽ��̬�ֲ�����(������õ���Qian����ĸ�ʽ)
{
    double eu,uv,feq;
    eu=(e[k][0]*u[0]+e[k][1]*u[1]);
    uv=(u[0]*u[0]+u[1]*u[1]);
    feq=w[k]*rho*(1.0+3.0*eu+4.5*eu*eu-1.5*uv);
	
    return feq;
}

void output(int m)  //���
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
                <<u[i][j][0]<<" "<<u[i][j][1]<<" "<<p[i][j]<<endl;	//����������ٶ���ѹ��
        }
}
void output0(int m)  //����ֲ���������(�Ǳ�Ҫ)
{
    ostringstream name;
    name<<"matlab"<<m<<".dat";
	ofstream out(name.str().c_str());
    for(j=0; j<=NY; j++)
        for(i=150; i<170; i++)
        {
            out<<double(i)/NX*NX/NY<<" "<<double(j)/NY<<" "		//����ֲ�����������
                <<u[i][j][0]<<endl;
        }
}
void output_extra()		//��������ٶȷֲ�����ֵ��;�ȷ��(����������ݴ����Ǳ�Ҫ)
{
	double u_real[NY+1],Err[NY+1],maxErr,umax;
	int test=3.0/4*NX;  //ˮƽ��������(�����ˮƽλ��)
	
	umax=(rho_in-rho_out)/3/lx*ly*ly/(8*rho0*niu);	//�����������
	for(i=0;i<=NY;i++){
		u_real[i]=umax*4*(double (i)/NY)*(1-double (i)/NY);	//���۽������ٷֲ�
		Err[i]=sqrt((u_real[i]-u[test][i][0])*(u_real[i]-u[test][i][0])+(0-u[test][i][1])*(0-u[test][i][1]))/umax;		//����ٶ����
	}
	maxErr = *max_element(Err,Err+NY+1);
	Re=umax*ly/niu;
    ostringstream outdata;
    outdata << "tau_f=" << tau_f << ",NY=" << NY << ".txt";
    ofstream out(outdata.str().c_str());
    out << "Title=\"����Ҷ\"\n"
        << "��������ٶ�Umax=" << umax << ",�������ٶ����maxerr= " << maxErr << ", �ɳ�ʱ��tau_f= " << tau_f  << ",���Ӳ���dx=" << dx << "\n"
        <<"Re="<< Re <<", �����������ٶ�g="<< g <<", ѹ���ݶ�="<<(rho_in-rho_out)/3/lx<<",ճ��niu="<<niu<<"\n"
        << "VARIABLES=\"����������\",\"x/H=150/40\",\" ���۽�\",\"  \"\n" << endl;
	for(i=0;i<=NY;i++){
        out << setw(4) << i << "  " << setprecision(8) << u[test][i][0] << "  " << u_real[i] << endl;	//����������ٵ���ֵ������۽�
	}
	out.close();
    ostringstream outdata_ex;
    outdata_ex << "tau_f NY and maxErr"<< ".txt";
    ofstream app(outdata_ex.str().c_str(),ofstream::app);
	app <<tau_f<<"  "<<NY<<"  "<<maxErr<<"  "<<n<<endl;		//������������
	app.close();
}
void Error()  //������
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
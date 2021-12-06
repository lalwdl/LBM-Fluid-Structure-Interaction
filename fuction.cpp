#include "fuction.h"

/*double sol[4][2]={0},node[2][2];//储存线与圆的交点,储存特殊点
int sol_i=0;//i为3或者4的情况
int flag_oncir[4]={0},flag_on_i=0;
int flag_incir[4]={0},flag_in_i=0;
int flag_outcir[4]={0},flag_out_i=0;*/

double sol[4][2],node[2][2];//储存线与圆的交点,储存特殊点
int sol_i=0;//i为3或者4的情况
int flag_oncir[4],flag_on_i;
int flag_incir[4],flag_in_i;
int flag_outcir[4],flag_out_i;

double solid_rate(double sx,double sy,double sl,double cx,double cy,double r)
{
    int ii,jj;
    for(ii=0;ii<4;ii++){
        flag_oncir[ii]=0;flag_incir[ii]=0;flag_outcir[ii]=0;
        for(jj=0;jj<2;jj++)
            sol[ii][jj]=0;
    }

    sol_i=0;
    flag_on_i=0;
    flag_in_i=0;
    flag_out_i=0;

    double rate,node[4][2]={{sx-sl*0.5,sy+sl*0.5},{sx+sl*0.5,sy+sl*0.5},{sx+sl*0.5,sy-sl*0.5},{sx-sl*0.5,sy-sl*0.5}};
    int num_in_cir=0,num_on_cir=0,num_out_cir=0,num_inon_cir;
    int i;
    double s_s=sl*sl;
    double dot_tri[3][2],dot_tri_[3][2],dot_fan[2][2],dot_fan_[2][2];

    for(i=0;i<4;i++){
        if((node[i][0]-cx)*(node[i][0]-cx)+(node[i][1]-cy)*(node[i][1]-cy)==r*r) {num_on_cir+=1;flag_oncir[flag_on_i]=i;flag_on_i++;}
        else if((node[i][0]-cx)*(node[i][0]-cx)+(node[i][1]-cy)*(node[i][1]-cy)<r*r)  {num_in_cir+=1;flag_incir[flag_in_i]=i;flag_in_i++;}
        else if((node[i][0]-cx)*(node[i][0]-cx)+(node[i][1]-cy)*(node[i][1]-cy)>r*r)  {num_out_cir+=1;flag_outcir[flag_on_i]=i;flag_out_i++;}
    }
    num_inon_cir=num_in_cir+num_on_cir;

    line_cir(cx,cy,r,node[0],node[1]);
    line_cir(cx,cy,r,node[2],node[1]);
    line_cir(cx,cy,r,node[3],node[2]);
    line_cir(cx,cy,r,node[3],node[0]);

    switch (num_inon_cir){
        case 0:
            if(sol_i==0) rate=0;
            else if(sol_i==2){
                dot_fan[0][0]=sol[0][0];
                dot_fan[0][1]=sol[0][1];
                dot_fan[1][0]=sol[1][0];
                dot_fan[1][1]=sol[1][1];

                rate=s_fan(r,dot_fan)/s_s;
            }
            break;
            //liang jiao dian
        case 1:
            if(num_in_cir==0) rate=0;
            else{
                dot_tri[0][0]=node[flag_incir[0]][0];
                dot_tri[0][1]=node[flag_incir[0]][1];//if sol i==2
                dot_tri[1][0]=sol[0][0];
                dot_tri[1][1]=sol[0][1];
                dot_tri[2][0]=sol[1][0];
                dot_tri[2][1]=sol[1][1];

                dot_fan[0][0]=sol[0][0];
                dot_fan[0][1]=sol[0][1];
                dot_fan[1][0]=sol[1][0];
                dot_fan[1][1]=sol[1][1];

                rate=(s_tri(dot_tri)+s_fan(r,dot_fan))/s_s;
            }
            break;
        case 2:
            if(num_in_cir==0){
                dot_fan[0][0]=node[flag_oncir[0]][0];
                dot_fan[0][1]=node[flag_oncir[0]][1];
                dot_fan[1][0]=node[flag_oncir[1]][0];
                dot_fan[1][1]=node[flag_oncir[1]][1];

                rate=s_fan(r,dot_fan)/s_s;
            }
            else if(num_in_cir==1){
                dot_tri[0][0]=node[flag_incir[0]][0];
                dot_tri[0][1]=node[flag_incir[0]][1];
                dot_tri[1][0]=sol[0][0];
                dot_tri[1][1]=sol[0][1];
                dot_tri[2][0]=node[flag_oncir[0]][0];
                dot_tri[2][1]=node[flag_oncir[0]][1];

                dot_fan[0][0]=sol[0][0];
                dot_fan[0][1]=sol[0][1];
                dot_fan[1][0]=node[flag_oncir[0]][0];
                dot_fan[1][1]=node[flag_oncir[0]][1];

                rate=(s_tri(dot_tri)+s_fan(r,dot_fan))/s_s;
            }
            else if(num_in_cir==2){

                if(sol_i==4){
                    if(node[flag_outcir[0]][0]==node[flag_outcir[1]][0])
                        qsort(sol,4, sizeof(sol[0]), compare_y);
                    else if(node[flag_outcir[0]][1]==node[flag_outcir[1]][1])
                        qsort(sol,4, sizeof(sol[0]), compare_x);
                    cmp_two(node[flag_outcir[0]],node[flag_outcir[1]]);

                    dot_fan[0][0]=sol[0][0];
                    dot_fan[0][1]=sol[0][1];
                    dot_fan[1][0]=sol[1][0];
                    dot_fan[1][1]=sol[1][1];                    //sol node

                    dot_fan_[0][0]=sol[2][0];
                    dot_fan_[0][1]=sol[2][1];
                    dot_fan_[1][0]=sol[3][0];
                    dot_fan_[1][1]=sol[3][1];

                    dot_tri[0][0]=node[0][0];
                    dot_tri[0][1]=node[0][1];
                    dot_tri[1][0]=sol[0][0];
                    dot_tri[1][1]=sol[0][1];
                    dot_tri[2][0]=sol[1][0];
                    dot_tri[2][1]=sol[1][1];

                    dot_tri_[0][0]=node[1][0];
                    dot_tri_[0][1]=node[1][1];
                    dot_tri_[1][0]=sol[2][0];
                    dot_tri_[1][1]=sol[2][1];
                    dot_tri_[2][0]=sol[3][0];
                    dot_tri_[2][1]=sol[3][1];

                    rate=(s_s-(s_tri(dot_tri)-s_fan(r,dot_fan)+s_tri(dot_tri_)-s_fan(r,dot_fan_)))/s_s;
                }

                else{
                    dot_fan[0][0]=sol[0][0];
                    dot_fan[0][1]=sol[0][1];
                    dot_fan[1][0]=sol[1][0];
                    dot_fan[1][1]=sol[1][1];

                    if((node[flag_incir[0]][0]==sol[0][0])||(node[flag_incir[0]][1]==sol[0][1])){
                        rate=(s_tra(sl,node[flag_incir[0]],sol[0],node[flag_incir[1]],sol[1])+s_fan(r,dot_fan))/s_s;
                    }
                    else{
                        rate=(s_tra(sl,node[flag_incir[0]],sol[1],node[flag_incir[1]],sol[0])+s_fan(r,dot_fan))/s_s;;
                    }
                }
            }
            break;
        case 3:
            if(num_in_cir==1) {
                dot_tri[0][0]=node[flag_incir[0]][0];
                dot_tri[0][1]=node[flag_incir[0]][1];
                dot_tri[1][0]=node[flag_oncir[0]][0];
                dot_tri[1][1]=node[flag_oncir[0]][1];
                dot_tri[2][0]=node[flag_oncir[1]][0];
                dot_tri[2][1]=node[flag_oncir[1]][1];

                dot_fan[0][0]=node[flag_oncir[0]][0];
                dot_fan[0][1]=node[flag_oncir[0]][1];
                dot_fan[1][0]=node[flag_oncir[1]][0];
                dot_fan[1][1]=node[flag_oncir[1]][1];

                rate=(s_tri(dot_tri)+s_fan(r,dot_fan))/s_s;
            }
            else if(num_in_cir==2) {
                dot_tri[0][0]=node[flag_outcir[0]][0];
                dot_tri[0][1]=node[flag_outcir[0]][1];
                dot_tri[1][0]=node[flag_oncir[0]][0];
                dot_tri[1][1]=node[flag_oncir[0]][1];
                dot_tri[2][0]=sol[0][0];
                dot_tri[2][1]=sol[0][1];

                dot_fan[0][0]=node[flag_oncir[0]][0];
                dot_fan[0][1]=node[flag_oncir[0]][1];
                dot_fan[1][0]=sol[0][0];
                dot_fan[1][1]=sol[0][1];

                rate=(s_s-(s_tri(dot_tri)-s_fan(r,dot_fan)))/s_s;
            }
            else if(num_in_cir==3) {
                dot_tri[0][0]=node[flag_outcir[0]][0];
                dot_tri[0][1]=node[flag_outcir[0]][1];
                dot_tri[1][0]=sol[0][0];
                dot_tri[1][1]=sol[0][1];
                dot_tri[2][0]=sol[1][0];
                dot_tri[2][1]=sol[1][1];

                dot_fan[0][0]=sol[0][0];
                dot_fan[0][1]=sol[0][1];
                dot_fan[1][0]=sol[1][0];
                dot_fan[1][1]=sol[1][1];

                rate=(s_s-(s_tri(dot_tri)-s_fan(r,dot_fan)))/s_s;
            }
            break;
        case 4:
            rate=1; break;
    }

    sol_i=0;
    return rate;
}

double s_tri(double dot[][2])
{
    return fabs(dot[0][0]*dot[1][1]+dot[0][1]*dot[2][0]+dot[1][0]*dot[2][1]-dot[0][0]*dot[2][1]-dot[0][1]*dot[1][0]-dot[1][1]*dot[2][0])/2.0;
}

double s_fan(double r,double dot[][2])
{
    double l,angle,stri,sfan,sfan_tri;

    l=pow(((dot[0][0]-dot[1][0])*(dot[0][0]-dot[1][0])+(dot[0][1]-dot[1][1])*(dot[0][1]-dot[1][1])),0.5);
    angle=acos(l/2.0/r);
    angle=(PI/2.0-angle)*2;
    stri=r*r*0.5*sin(angle);
    sfan=angle*0.5*r*r;
    sfan_tri=sfan-stri;

    return sfan_tri;
}

void line_cir(double cx,double cy,double r,double dot1[],double dot2[])//两个点从小到大
{
    double s1,s2;

    if(dot1[0]==dot2[0]){
        if((r*r-(dot1[0]-cx)*(dot1[0]-cx))>0){//相切未考虑。。。
            s1=-sqrt((r*r-(dot1[0]-cx)*(dot1[0]-cx)))+cy;
            s2= sqrt((r*r-(dot1[0]-cx)*(dot1[0]-cx)))+cy;
            if((s1>dot1[1])&&(s1<dot2[1])) {sol[sol_i][0]=dot1[0];sol[sol_i][1]=s1;sol_i++;}
            if((s2>dot1[1])&&(s2<dot2[1])) {sol[sol_i][0]=dot1[0];sol[sol_i][1]=s2;sol_i++;}//一个线两个交点的情况。。。
        }
    }
    else{
        if((r*r-(dot1[1]-cy)*(dot1[1]-cy))>0){
            s1=-sqrt((r*r-(dot1[1]-cy)*(dot1[1]-cy)))+cx;
            s2= sqrt((r*r-(dot1[1]-cy)*(dot1[1]-cy)))+cx;
            if((s1>dot1[0])&&(s1<dot2[0])) {sol[sol_i][0]=s1;sol[sol_i][1]=dot1[1];sol_i++;}
            if((s2>dot1[0])&&(s2<dot2[0])) {sol[sol_i][0]=s2;sol[sol_i][1]=dot1[1];sol_i++;}
        }
    }
}

double s_tra(double sl,double dot1[],double dot2[],double dot3[],double dot4[])
{
    double a,b;

    a=pow(((dot1[0]-dot2[0])*(dot1[0]-dot2[0])+(dot1[1]-dot2[1])*(dot1[1]-dot2[1])),0.5);
    b=pow(((dot3[0]-dot4[0])*(dot3[0]-dot4[0])+(dot3[1]-dot4[1])*(dot3[1]-dot4[1])),0.5);

    return 0.5*(a+b)*sl;
}

int compare_x(const void *p,const void *q)
{
    return ((double *)p)[0] - ((double *)q)[0];
}
int compare_y(const void *p,const void *q)
{
    return ((double *)p)[1] - ((double *)q)[1];
}

void cmp_two(double a[],double b[])
{
    double c;
    if(a[0]==b[0]){
        if(a[1]>b[1])
            c=a[1];a[1]=b[1];b[1]=c;
    }
    if(a[1]==b[1]){
        if(a[0]>b[0])
            c=a[0];a[0]=b[0];b[0]=c;
    }
    node[0][0]=a[0];
    node[0][1]=a[1];
    node[1][0]=b[0];
    node[1][1]=b[1];
}

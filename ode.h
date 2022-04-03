#pragma once
#include<cmath>
#include<vector>
#include<initializer_list>
typedef long double ld;
typedef std::initializer_list<ld> ls;

namespace ode 
{

template<typename function>
//显示欧拉
ld explicit_Euler(function func,ls x,ld y0,long n=1000)//函数,区间,初值,迭代次数
{//error: O(h^2)
ld beg=*x.begin(),ed=*(x.begin()+1);
ld h=(ed-beg)/n;//步长
for(long i=0;i<n;++i)
{
	y0+=h*func(beg,y0);//增值
	beg+=h;//x增加
}
return y0;
}

template<typename function>
//(半)隐式欧拉
ld implicit_Euler(function func,ls x,ld y0,long n=1000)//函数,区间,初值,迭代次数
{
ld beg=*x.begin(),ed=*(x.begin()+1);	
ld h=(ed-beg)/n;//步长
for(long i=0;i<n;++i)
	y0+=h*func(beg+=h,y0);//增值
return y0;
}

template<typename function,typename Diff>
//隐式欧拉方法更擅长求解刚性(stiff)问题,可以减少震荡
//一般过程是每一步用牛顿法求解y_{n+1}
ld implicit_Euler(function func,Diff diff,ls x,ld y0,ld tol=1e-3,long n=1000)//函数,导数,区间,初值,迭代次数
{
ld beg=*x.begin(),ed=*(x.begin()+1);	
ld h=(ed-beg)/n,y1,y2;//步长
for(long i=0;i<n;++i)
{
	beg+=h;
	y1=y0+1,y2=y0;
	while(fabs(y1-y0)>tol)
	{
		y2=y1-(y1-h*func(beg,y1)-y0)/(1-h*diff(beg,y1));
		y1=y2;
	}
	y0=y1;
}
return y0;
}

template<typename function>
//改进欧拉(显示梯形)
ld Heun(function func,ls x,ld y0,long n=1000)//函数,区间,初值,迭代次数
{
ld beg=*x.begin(),ed=*(x.begin()+1);
ld h=(ed-beg)/n,fy;//步长
for(long i=0;i<n;++i)
{
	fy=func(beg,y0);
	y0+=h*(fy+func(beg+h,y0+h*fy))/2;//增值,其中f(t_{n+1},y_{n+1})的y_{n+1}用一阶近似代替
	beg+=h;//x增加
}
return y0;
}

template<typename function,typename Diff>
//2阶泰勒方法,阶数越高越精确,因为需要高阶导数,一般不用的
ld Taylor2(function func,Diff diff,ls x,ld y0,long n=1000)//函数,导数,区间,初值,迭代次数
{
ld beg=*x.begin(),ed=*(x.begin()+1);
ld h=(ed-beg)/n;//步长
for(long i=0;i<n;++i)
{
	y0+=h*func(beg,y0)+h*h*diff(beg,y0)/2;//增值
	beg+=h;//x增加
}
return y0;
}

template<typename function>
//Runge-Kutta
//使用若干个函数点做泰勒展开的近似公式,达到一定的阶数(每一个点的斜率由前一个点的导数值做线性预测)
ld RK4(function func,ls x,ld y0,std::string str="classic",long n=1000)//函数,区间,初值,method,迭代次数
{
ld beg=*x.begin(),ed=*(x.begin()+1);
ld h=(ed-beg)/n,s1,s2,s3,s4;
switch(str[0])
{
case 'c':
case 'C':
	for(long i=0;i<n;++i)
	{//古典RK
	s1=func(beg,y0);
	s2=func(beg+h/2,y0+h/2*s1);
	s3=func(beg+h/2,y0+h/2*s2);
	s4=func(beg+h,y0+h*s3);
	y0+=h*(s1+2*s2+2*s3+s4)/6;
	beg+=h;
	}
	break;
case 'k':
case 'K':
	for(long i=0;i<n;++i)
	{//Kutta方法
	s1=func(beg,y0);
	s2=func(beg+h/3,y0+h/3*s1);
	s3=func(beg+h*2/3,y0-h*s1/3+h*s2);
	s4=func(beg+h,y0+h*(s1-s2+s3));
	y0+=h*(s1+3*s2+3*s3+s4)/8;
	beg+=h;
	}
	break;
case 'g':
case 'G':
	for(long i=0;i<n;++i)
	{//Gill方法
	s1=func(beg,y0);
	s2=func(beg+h/2,y0+h/2*s1);
	s3=func(beg+h/2,y0+h*s1*(sqrt(2)-1)/2.+h*s2*(2-sqrt(2))/2.);
	s4=func(beg+h,y0-h*s2*sqrt(2)/2.+h*s3*(2+sqrt(2))/2.);
	y0+=h*(s1+(2-sqrt(2))*s2+(2+sqrt(2))*s3+s4)/6;
	beg+=h;
	}
	break;
default :
	std::cerr<<"\nRK4: error method.\n";
	exit(0);
}
return y0;
}

template<typename function>
//Runge-Kutta-Fehlberg
//自适应步长
ld RKF45(function func,ls x,ld y0,ld T=1e-4,ld h=0.1)//函数,区间,初值,容差,初始步长
{
ld t=*x.begin(),ed=*(x.begin()+1);
ld s1,s2,s3,s4,s5,s6,z,e,y1,d;
int n=0,i=0;
bool flag=0;
	while(i<100000)//最多循环10万次
	{
	d=ed-t;
	if(fabs(d)<h)//最后一步要刚好覆盖区间
	{
		flag=1;
		h=d;
	}
	s1=func(t,y0);
	s2=func(t+h/4.,y0+h/4.*s1);
	s3=func(t+h*3./8.,y0+h*s1*3./32.+h*s2*9./32.);
	s4=func(t+h*12./13.,y0+h*s1*1932./2197.-h*s2*7200./2197.+h*s3*7296./2197.);
	s5=func(t+h,y0+h*s1*439./216.-h*s2*8.+h*s3*3680./513.-h*s4*845./4104.);
	s6=func(t+h/2,y0-h*s1*8./27.+h*s2*2-h*s3*3544./2565.+h*s4*1859./4104.-h*s5*11./40.);
	//5阶
	z=y0+h*(s1*16./135.+s3*6656./12825.+s4*28561./56430.-s5*9./50.+s6*2./55.);
	//4阶
	y1=y0+h*(s1*25./216.+s3*1408./2565.+s4*2197./4104.-s5/5.);
	if(flag)
	return z;//说明在循环内已完成所有步长
	e=fabs((z-y1)/y0);//相对精度
	if(e<T)//容差内更新步长
	{
		y0=z;
		t+=h;
		n=0;
		++i;
	}
	else if(n++<2)h*=0.8*pow(T/e,0.2);//精度不够,重做该步
	else h/=2.,n=0;//重复失败,减半处理
	}
return y0;
}

template<typename function>
//Adams多步方法选取之前计算的多个点进行拉格朗日插值,用外插值的积分预测下一个点
//同理,隐式方法采用内插值
//Adams4步方法
ld Adams4(function func,ls x,ld y0,long n=1000)
{
	ld beg=*x.begin(),ed=*(x.begin()+1);
	ld h=(ed-beg)/n,s1,s2,s3,s4,s[5];
	long i=0;
	s[0]=y0;
	for(;i<4;++i)
	{//先用RK4得到Adams前四个值
		s1=func(beg,y0);
		s2=func(beg+h/2,y0+h/2*s1);
		s3=func(beg+h/2,y0+h/2*s2);
		s4=func(beg+h,y0+h*s3);
		s[i+1]=s[i]+h*(s1+2*s2+2*s3+s4)/6;
		beg+=h;
	}
	y0=s[4];
	for(;i<n;++i)
	{
		y0+=h*(55*func(beg,s[4])-59*func(beg-h,s[3])+37*func(beg-2*h,s[2])-9*func(beg-3*h,s[1]))/24.;
		s[1]=s[2];
		s[2]=s[3];
		s[3]=s[4];
		s[4]=y0;
		beg+=h;
	}
	return y0;
}

template<typename function>
//改进的Adams算法(内插加外插误差相消)
ld Adams_ipv(function func,ls x,ld y0,long n=1000)
{
	ld beg=*x.begin(),ed=*(x.begin()+1);
	ld h=(ed-beg)/n,s1,s2,s3,s4,s[5];
	long i=0;
	s[0]=y0;
	for(;i<4;++i)
	{//先用RK4得到Adams前四个值
		s1=func(beg,y0);
		s2=func(beg+h/2,y0+h/2*s1);
		s3=func(beg+h/2,y0+h/2*s2);
		s4=func(beg+h,y0+h*s3);
		s[i+1]=s[i]+h*(s1+2*s2+2*s3+s4)/6;
		beg+=h;
	}
	y0=s[4];
	ld p2,p1=0,c2,c1=0,m;
	for(;i<n;++i)
	{
		s1=func(beg,s[4]);
		s2=func(beg-h,s[3]);
		s3=func(beg-2*h,s[2]);
		s4=func(beg-3*h,s[1]);
		p2=y0+h*(55*s1-59*s2+37*s3-9*s4)/24.;//未改进结果
		//m=p2;
		m=p2+(c1-p1)*251./270.;//修正的预估值(误差相消)
		m=func(beg+h,m);
		c2=y0+h*(9*m+19*s1-5*s2+s3)/24.;//修正的矫正值
		//y0=c2;
		y0=c2-(c2-p2)*19./270.;
		//更新参数
		c1=c2;
		p1=p2;
		s[1]=s[2];
		s[2]=s[3];
		s[3]=s[4];
		s[4]=y0;
		beg+=h;
	}
	return y0;
}
	
template<typename function>
//打靶法
std::vector<ld>shoot(function func,ls x,ls yx,ls s={-10,10},ld T=1e-3,long n=1000)
{
ld beg=*x.begin(),ed=*(x.begin()+1);
ld y0=*yx.begin(),yn=*(yx.begin()+1);
ld sp1=*s.begin(),sp2=*(s.begin()+1);
ld h=(ed-beg)/n,z,z1,z2,s1,s2,s3,s4,l1,l2,l3,l4,tem,t,sp0;//步长
std::vector<ld>y(n+1),y1(n+1);
auto pre=[&](ld sp)->ld
{
t=beg;
y[0]=y0;
y1[0]=sp;
for(long i=0;i<n;++i)
{//RK4
	tem=y1[i];
	s1=func(t,y[i],tem);
	l1=tem;
	tem=y1[i]+h/2*l1;
	s2=func(t+h/2,y[i]+h/2*s1,tem);
	l2=tem;
	tem=y1[i]+h/2*l2;
	s3=func(t+h/2,y[i]+h/2*s2,tem);
	l3=tem;
	tem=y1[i]+h*l3;
	s4=func(t+h,y[i]+h*s3,tem);
	l4=tem;
	y1[i+1]=y1[i]+h*(s1+2*s2+2*s3+s4)/6;
	y[i+1]=y[i]+h*(l1+2*l2+2*l3+l4)/6;
	t+=h;//x增加	
}
return y[n]-yn;
};

while(true)//割线法求零点
{
z1=pre(sp1);
z2=pre(sp2);
sp0=sp1-z1*(sp1-sp2)/(z1-z2);
sp2=sp1;
sp1=sp0;
if(fabs(sp1-sp2)<T)
return y;
}
/*
while(true)//改进割线法求零点
//for(int i=0;i<100;++i)
{
z1=pre(sp1);
z2=pre(sp2);
sp0=(sp2*z1-sp1*z2)/(z1-z2);
z=pre(sp0);
if(fabs(z)<T)
return y;
if(z1*z<0)sp2=sp0;
else sp1=sp0;
}
*/
return y;
}
}

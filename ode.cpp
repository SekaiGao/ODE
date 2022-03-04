#pragma GCC optimize(3,"Ofast","inline")
#include<iostream>
#include"matrix.h"
typedef long double ld;
typedef std::initializer_list<ld> ls;
using namespace std;

auto func=[](ld x,ld y)->ld{return x*y+pow(x,3);};//f(x,y)

auto diff=[](ld x,ld y)->ld{return y+3*x*x+x*(x*y+pow(x,3));};//f'(x,y)关于x的全导数


template<typename function>
//显示欧拉
ld explicit_Euler(function func,ls x,ld y0,long n=1000)//函数,区间,初值,迭代次数
{
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
//隐式欧拉
ld implicit_Euler(function func,ls x,ld y0,long n=1000)//函数,区间,初值,迭代次数
{
ld beg=*x.begin(),ed=*(x.begin()+1);	
ld h=(ed-beg)/n;//步长
for(long i=0;i<n;++i)
	y0+=h*func(beg+=h,y0);//增值
return y0;
}

template<typename function>
//改进欧拉(显示梯形)
ld Heun(function func,ls x,ld y0,long n=1000)//函数,区间,初值,迭代次数
{
ld beg=*x.begin(),ed=*(x.begin()+1);
ld h=(ed-beg)/n;//步长
for(long i=0;i<n;++i)
{
	y0+=h*(func(beg,y0)+func(beg+h,y0+h*func(beg,y0)))/2;//增值
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
ld RK4(function func,ls x,ld y0,long n=1000)//函数,区间,初值,迭代次数
{
ld beg=*x.begin(),ed=*(x.begin()+1);
ld h=(ed-beg)/n,s1,s2,s3,s4;
	for(long i=0;i<n;++i)
	{
	s1=func(beg,y0);
	s2=func(beg+h/2,y0+h/2*s1);
	s3=func(beg+h/2,y0+h/2*s2);
	s4=func(beg+h,y0+h*s3);
	y0+=h*(s1+2*s2+2*s3+s4)/6;
	beg+=h;
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
	if(abs(d)<h)//最后一步要刚好覆盖区间
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

int main()
{
	cout<<fixed;
	cout.precision(16);
	ld x=5,y=1;
	ld ty=(2.+y)*exp(x*x/2)-x*x-2.;
	cout<<ty<<endl;

	cout<<explicit_Euler(func,{0,5},1)<<endl;

	cout<<implicit_Euler(func,{0,5},1)<<endl;

	cout<<Heun(func,{0,5},1)<<endl;

	cout<<Taylor2(func,diff,{0,5},1)<<endl;

	cout<<RK4(func,{0,5},1)<<endl;

	cout<<RKF45(func,{0,5},1,1e-16)<<endl;
	return 0;
}
#pragma GCC optimize(3,"Ofast","inline")
#include<iostream>
#include"ode.h"
using namespace std;
auto func=[](ld x,ld y)->ld{return x*y+pow(x,3);};//f(x,y)

auto diff=[](ld x,ld y)->ld{return y+3*x*x+x*(x*y+pow(x,3));};//f'(x,y)关于x的全导数

int main()
{
	cout<<fixed;
	cout.precision(16);
	ld x=5,y=1;
	ld ty=(2.+y)*exp(x*x/2)-x*x-2.;

	cout<<ty<<endl;

	cout<<ode::explicit_Euler(func,{0,5},1)<<endl;

	cout<<ode::implicit_Euler(func,{0,5},1)<<endl;

	cout<<ode::Heun(func,{0,5},1)<<endl;

	cout<<ode::Taylor2(func,diff,{0,5},1)<<endl;

	cout<<ode::RK4(func,{0,5},1)<<endl;

	cout<<ode::RKF45(func,{0,5},1,1e-16)<<endl;

	return 0;
}
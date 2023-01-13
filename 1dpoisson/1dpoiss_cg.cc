//this program for prob uxx=f(x)
// a<x<b u(a)=u1 u(b)=u2
#include<iostream>
#include<cmath>
#include<vector>
#include<numeric>
#include<algorithm>
#include<fstream>
//#include "linsolvers.hh"
#include "cg.hh"
constexpr float val(float t) //for calculatinf f(x)
{ return (2+t)*exp(t);}

template<typename T>
int genvecb(float xo,float xn,float yo,float yn, double h,T& b)
{
  int n=0;
  std::vector<float> x;
  while(xo<(xn+h))
      {
	x.push_back(xo);
	xo+=h;
	n++;		
      }
    x.resize(n);
    b.insert(b.end(),(pow(h,2)*val(*(x.begin()+1)))-yo);
    for_each(x.begin()+2,x.end()-2,[&](auto t)
				   { b.insert(b.end(),pow(h,2)*val(t));});
    b.insert(b.end(),(pow(h,2)*val(x.back()-h))-yn);
    return n-2;
}
void genmatA(const int& dim,std::vector<float>& a,std::vector<int>& ra,std::vector<int>& ca)

{
  ra.assign({0,2});  //ra and ca general for tridiagonal system
  for(int i=0;i<dim-2;i++)
    {
      ra.push_back(3+(ra.back()));
    }
  ra.push_back(2+ra.back());

  for(int i=0;i<dim;i++) //vector a particular for this prob
    {
      if(i==0)
	{
	  ca.assign({i,i+1});
	  a.assign({-2,1});
	}
      else if(i==dim-1)
	{
	  ca.insert(ca.end(),{i-1,i});
	  a.insert(a.end(),{1,-2});
	}
      else
	{
	  ca.insert(ca.end(),{i-1,i,i+1});
	  a.insert(a.end(),{1,-2,1});
	}
      
    }
}
void plot(float& xo,float& xn,float& yo,float &yn,float h,const std::vector<double>& u)
{
  std::ofstream fout;
  fout.open("bvpplot.txt",std::ios::out);
  int n=0;
  std::vector<float> x;
  while(xo<(xn+h))
      {
	x.push_back(xo);
	xo+=h;
	n++;		
      }
    x.resize(n);
    fout<<*(x.begin())<<" "<<yo<<"\n";
    n=0;
    for_each(x.begin()+1,x.end()-1,[u,&n,&fout](auto x1)
				   { fout<<x1<<" "<<u[n]<<"\n";n++;});
    fout<<*(x.end())<<" "<<yn<<"\n";
    fout.close();
				    
  }
int main()
{
  float xo=-1,xn=1,yo=-(1/exp(1)),yn=exp(1),h=0.01;
  std::vector<float> a;std::vector<int> ra,ca;
  std::vector<float> b;
  int dim;
  std::cout<<"b is:";
  dim=genvecb(xo,xn,yo,yn,h,b);
  // for_each(b.begin(),b.end(),[](auto x)
  //			       {std::cout<<x<<" ";});
  //std::cout<<std::endl;
  std::vector<float> u(dim);
  for(int i=0;i<dim;i++)
    {
      u[i]=0;
    }
  std::cout<<"dimension is:"<<dim<<std::endl;
  genmatA(dim,a,ra,ca);
  std::cout<<"matrix generated!";
  std::cout<<std::endl;
  //  gauss_siedel solver;
  //solver.apply(a,ra,ca,b,u0,im);
  cg solver;
  solver.apply(a,b,u,ra,ca,dim);
  std::cout<<"U after gauss siedel/cg:"<<std::endl;
  int k=0;
  for_each(u.begin(),u.end(),[&](auto x)
			     { std::cout<<x<<" ";k++;});
  std::cout<<std::endl<<"k is:"<<k;
  //plot(xo,xn,yo,yn,h,u);
  std::cout<<std::endl;
  return 0;
}

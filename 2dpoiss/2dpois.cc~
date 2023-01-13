 //this program for prob -(uxx+uyy)=f(x)
// given f(x) and boundary values u(0,y),u(1,y),u(x,0),u(x,1) //works for dim 3or more
#include<iostream>
#include<cmath>
#include<vector>
#include<numeric>
#include<algorithm>
#include<fstream>
#include<cmath>
#include "cg.hh"
//#include "linsolvers.hh"

constexpr float valf(float x,float y) //for calculating u(x,0),u(x,1)
{
  return (20*sin(3*M_PI*x)*cos(2*M_PI*y));
}
constexpr float valx(float x,float y) //for calculating u(x,0),u(x,1)
{
  if(y==0)
    return (x*x*x);
  else
    return 1;
}
constexpr float valy(float x,float y) //for calculating u(0,y),u(1,y)
{
  if(x==0)
    return (y*y);
  else
    return 1;
}

template<typename T>
int genvecb(float xo,float xn,float yo,float yn, double h,T& b)
{
  int nx=0,ny=0;
  std::vector<float> x,y;
  while(xo<(xn+h))
    {
      x.push_back(xo);
      // std::cout<<x.back()<<" ";
      xo+=h;
      nx+=1;;		
    }
  if(x.back()>xn)
    {
     x.resize(nx-1);
    }

  //for_each(x.begin(),x.end(),[](auto x1)
  //			     {std::cout<<x1<<" ";});
  nx=x.size();
  // std::cout<<"size is:"<<x.size()<<std::endl;
  while(yo<(yn+h))
    {
      y.push_back((yo));
      yo+=h;
      ny++;		
    }
  if(y.back()>yn)
    {
     y.resize(ny-1);
    }
//for_each(x.begin(),x.end(),[](auto x1)
//  			     {std::cout<<x1<<" ";});
  ny=y.size();
//std::cout<<"xback is:"<<x.back()<<std::endl;
//std::cout<<"yback is:"<<y.back()<<std::endl;
  for(int i=0;i<nx-2;i++)
    {
      std::vector<double> b1;
      double f;
      if(i==0)
	{
	  for(int j=0;j<nx-2;j++)
	    {
	      f=(h*h*valf(*(x.begin()+j+1),*(y.begin()+i+1)));
	      if(j==0)
		b1.push_back(f+valx(*(x.begin()+1),*(y.begin()))+valy(*(x.begin()),*(y.begin()+1)));
	      else if(j==nx-2-1)
		b1.push_back(f+valx(*(x.end()-1),*(y.begin()))+valy(*(x.end()),*(y.begin()+1)));
	      else
		b1.push_back(f+valx(*(x.begin()+j+1),*(y.begin())));
	    }
	  b.push_back(b1);
	  b1.clear();
	}
      else if(i==nx-2-1)
	{
	  for(int j=0;j<nx-2;j++)
	    {
	      f=(h*h*valf(*(x.begin()+j+1),*(y.begin()+i+1)));
	      if(j==0)
		b1.push_back(f+valx(*(x.begin()+1),*(y.end()))+valy(*(x.begin()),*(y.end()-1)));
	      else if(j==nx-2-1)
		b1.push_back(f+valx(*(x.end()-1),*(y.end()))+valy(*(x.end()),*(y.end()-1)));
	      else
		b1.push_back(f+valx(*(x.begin()+j+1),*(y.end())));
	    }
	  b.push_back(b1);
	  b1.clear(); 
	}
      else
	{
	  for(int j=0;j<nx-2;j++)
	    {
	      f=(h*h*valf(*(x.begin()+j+1),*(y.begin()+i+1)));
	      if(j==0)
		b1.push_back(f+valy(*(x.begin()),*(y.begin()+i+1)));
	      else if(j==nx-2-1)
		b1.push_back(f+valy(*(x.end()),*(y.begin()+i+1)));
	      else
		b1.push_back(0);
	    }
	  b.push_back(b1);
	  b1.clear();
	}
    }
  return nx-2;
}
void genmatA(const int& dim,std::vector<std::vector<float>>& a,std::vector<int>& ra,std::vector<int>& ca)

{
  ra.assign({0,2});  //ra and ca general for tridiagonal system
  for(int i=0;i<dim-2;i++)
    {
      ra.push_back(3+(ra.back()));
    }
  ra.push_back(2+ra.back());
  std::vector<float> v1,v2;
  for(int i=0;i<dim;i++)
    {
      if(i==0)
	v1.insert(v1.end(),{4,-1});
      else if(i==dim-1)
	v1.insert(v1.end(),{-1,4});
      else
	v1.insert(v1.end(),{-1,4,-1});
      v2.insert(v2.end(),{-1});
    }

  for(int i=0;i<dim;i++) //vector a particular for this prob
    {
      if(i==0)
	{
	  a.insert(a.end(),{v1,v2});
	  ca.assign({i,i+1});
	  
	}
      else if(i==dim-1)
	{
	  ca.insert(ca.end(),{i-1,i});
	  a.insert(a.end(),{v2,v1});
	}
      else
	{
	  ca.insert(ca.end(),{i-1,i,i+1});
	  a.insert(a.end(),{v2,v1,v2});
	}
      
    }
}

int main()
{
  float xo=0,xn=1,yo=0,yn=1,h=0.2;
  std::vector<std::vector<float>> a;std::vector<int> ra,ca;
  std::vector<std::vector<double>> b,u; std::vector<double> u1;
  int dim;
  dim=genvecb(xo,xn,yo,yn,h,b);
  std::cout<<"dimension is:"<<dim<<std::endl;
  std::cout<<"b is: ";
  for_each(b.begin(),b.end(),[](auto x)
			     {for_each(x.begin(),x.end(),[](auto y)
							 {
							   std::cout<<y<<" ";
							 });
			       std::cout<<std::endl;
			     });
  genmatA(dim,a,ra,ca);
  std::cout<<"matrix A is:";
  for_each(a.begin(),a.end(),[](auto x)
			     { for_each(x.begin(),x.end(),[](auto y)
							  { std::cout<<y<<" ";});
			       std::cout<<std::endl;});
  for_each(ca.begin(),ca.end(),[](auto x)
    { std::cout<<x<<" ";});
  for(int i=0;i<dim;i++)
    {
      for(int j=0;j<dim;j++)
	u1.push_back(0);
      u.push_back(u1); //initial guess
      u1.clear();
    }
  
   cg solver;
  solver.apply(a,b,u,ra,ca,dim);
  std::cout<<std::endl<<"after solving:"<<std::endl;
  for(int i=0;i<dim;i++)
    {
      for(int j=0;j<dim;j++)
	std::cout<<u[i][j]<<" ";
      std::cout<<std::endl;
    }
    std::cout<<std::endl;
  /*  //using gauss_siedel
  std::cout<<std::endl;
  std::vector<float> a1; std::vector<int> ra1,ca1;
  std::vector<double> b1;
  ra1.push_back(0);int s=0;
  int k=0,k1=0,k2=((dim-2)*dim);
  for(int i=0;i<(dim*dim);i++)
    {
      if(i<dim)
	{
	  if(i==0)
	    {s=3;
	      a1.insert(a1.end(),{4,-1,-1});
	      ra1.push_back(s);
	      ca1.insert(ca1.end(),{k,k+1,k+3});
	    }
	  else if(i==dim-1)
	    {
	      s+=3;
	      a1.insert(a1.end(),{-1,4,-1});
	      ra1.push_back(s);
	      ca1.insert(ca1.end(),{k,k+1,k+4});
	    }
	  else
	    {
	      a1.insert(a1.end(),{-1,-1,4,-1});
	      s+=4;
	      ra1.push_back(s);
	      ca1.insert(ca1.end(),{k,k+1,k+2,k+4});
	      k+=1;
	    }
	}
      else if(i>=(dim*dim)-(dim))
	{
	  if(i==(dim*dim)-(dim))
	    {
	      a1.insert(a1.end(),{-1,4,-1});
	      s+=3;
	      ra1.push_back(s);
	      ca1.insert(ca1.end(),{k2,k2+3,k2+4});
	    }
	  else if(i==(dim*dim)-1)
	    {
	      a1.insert(a1.end(),{-1,-1,4});
	      s+=3;
	      ra1.push_back(s);k2+=1;
	      ca1.insert(ca1.end(),{k2,k2+2,k2+3});
	    }
	  else
	    {
	      a1.insert(a1.end(),{-1,-1,4,-1});
	      s+=4;
	      ra1.push_back(s);k2+=1;
	      ca1.insert(ca1.end(),{k2,k2+2,k2+3,k2+4});
	    }
	}
      else
	{
	  if(i%dim==0)
	    {
	      a1.insert(a1.end(),{-1,4,-1,-1});
	      s+=4;
	      ra1.push_back(s);
	      ca1.insert(ca1.end(),{k1,k1+3,k1+4,k1+6});k1+=1;
	    }
	  else if(i%dim==dim-1)
	    {
	      std::cout<<"i is:"<<i<<std::endl;
	      a1.insert(a1.end(),{-1,-1,4,-1});
	      s+=4;
	      ra1.push_back(s);
	      ca1.insert(ca1.end(),{k1,k1+2,k1+3,k1+6});k1=0;
	    }
	  else
	    {
	      a1.insert(a1.end(),{-1,-1,4,-1,-1});
	      s+=5;
	      ra1.push_back(s);
	      ca1.insert(ca1.end(),{k1,k1+2,k1+3,k1+4,k1+6});k1+=1;
	    }
      
	    }
    }
  std::vector<double> u2;
  for(int i=0;i<dim;i++)
    { for(int j=0;j<dim;j++)
	{b1.push_back(b[i][j]);u2.push_back(0);}
	  
    }
  for_each(a1.begin(),a1.end(),[](auto x)
			       { std::cout<<x<<" ";});
  std::cout<<std::endl;
  for_each(ca1.begin(),ca1.end(),[](auto x)
			       { std::cout<<x<<" ";});
  std::cout<<std::endl;
  for_each(ra1.begin(),ra1.end(),[](auto x)
    { std::cout<<x<<" ";});
  std::cout<<std::endl;
  gauss_siedel solver;
  solver.apply(a1,ra1,ca1,b1,u2,(dim*dim)-1);
  for_each(u2.begin(),u2.end(),[](auto x)
    { std::cout<<x<<" ";});
    std::cout<<std::endl;*/
  return 0;
}

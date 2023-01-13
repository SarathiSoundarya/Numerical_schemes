#include<iostream>
#include<cmath>
#include<vector>
#define err 0.0000007

  
template<typename T>
class vect //for 1D
{
public:
  int dim;
  std::vector<T> v;
  vect(int d)
  {
    dim=d;
    for(int i=0;i<dim;i++)
      v.push_back(0);
   }
  vect(const std::vector<T>& arr,int d)
  {
    dim=d;
    for(int i=0;i<dim;i++)
      v.push_back(arr[i]);
  }
 void operator+=(const vect<T>& v1)
  {
    for(int i=0;i<dim;i++)
	 v[i]+=v1.v[i];
  }
  void operator-=(const vect<T>& v1)
  {
    for(int i=0;i<dim;i++)
	 v[i]-=v1.v[i];
   }
  void scalar(double s)
  {
    for(int i=0;i<dim;i++)
      v[i]*=s;
  }
   
  void operator=(vect<T>& v1)
  {
    for(int i=0;i<dim;i++)
      v[i]=v1.v[i];
  }
  long double operator*(const vect<T>& v1)
  {
    long double sum=0.0;
    for(int i=0;i<dim;i++)
      sum+=(v[i]*v1.v[i]);
    return (sum);
  }
  long double norm()
  {
    long double sum=0.0;
    for(int i=0;i<dim;i++)
      sum+=(v[i]*v[i]);
    return sqrt(sum);
  }
  void apadd(const vect<T> &v1,long double alpha,const vect<T> &v2)
  {
    for(int i=0;i<dim;i++)
      v[i]=v1.v[i]+(alpha*v2.v[i]);
  }
  void print(std::vector<T>& ini)
  {
    for(int i=0;i<dim;i++)
      ini[i]=v[i];
  }
};


template<typename T,typename L>
class mat
{
public:
  int dim;
  int check;
  std::vector<T> a;
  std::vector<int> ra;
  std::vector<int> ca;//for tridiagonal
 mat(std::vector<T> a1,std::vector<int> ra1,std::vector<int> ca1,int d)
  {
    dim=d;
    for(int i=0;i<ra1.size();i++)
      ra.push_back(ra1[i]);
    for(int i=0;i<ca1.size();i++)
      {
        ca.push_back(ca1[i]);
	a.push_back(a1[i]);
      }
  }
  void  mulvec(vect<L> &v1)
  {
    
    vect<L> v2(dim);
    int n,l=0,k=0,co=0;
    for(int i=0;i<dim;i++)
      {
	n=ra[l+1]-ra[l];
	l+=1;
	co=0;
	while(co<n)
	  {
	    v2.v[i]+=(a[k]*v1.v[ca[k]]); //multiplication of csr with vec
	    k++;
	    co++;
	  }
      }
    v1=v2;
    
    }
    
};

//FOR 2D
template<typename A>
class vect<std::vector<A>>
{
 
public:
  
  int dim;
  std::vector<std::vector<A>> v;
  vect(int d)
  {
    dim=d;
    std::vector<A> v1; 
    for(int i=0;i<dim;i++)
      {
	for(int j=0;j<dim;j++)
	  v1.push_back(0);
	v.push_back(v1);
      }
  }
    
  vect(std::vector<std::vector<A>> arr,int d)
  {
    dim=d;
    for(int i=0;i<dim;i++)
      v.push_back(arr[i]);
  }
 
  void operator+=(const vect<std::vector<A>>& v1)
  {
    for(int i=0;i<dim;i++)
      {
	for(int j=0;j<dim;j++)
	  v[i][j]+=v1.v[i][j];
      }
  }
  void operator-=(const vect<std::vector<A>>& v1)
  {
    for(int i=0;i<dim;i++)
      {
	for(int j=0;j<dim;j++)
	  v[i][j]-=v1.v[i][j];
      }
  }
  void scalar(double s)
  {
    for(int i=0;i<dim;i++)
      {
	for(int j=0;j<dim;j++)
	  v[i][j]*=s;
      }
  }
  void operator=(const vect<std::vector<A>>& v1)
  {
    for(int i=0;i<dim;i++)
      {
	for(int j=0;j<dim;j++)
	  v[i][j]=v1.v[i][j];
      }
  }
  long double operator*(const vect<std::vector<A>>& v1)
  {
    long double sum=0.0;
    for(int i=0;i<dim;i++)
      {
	for(int j=0;j<dim;j++)
	  sum+=(v[i][j]*v1.v[i][j]);
      }
        
    return (sum);
  }
  long double norm()
  {
    long double sum=0.0;
    for(int i=0;i<dim;i++)
      {
	for(int j=0;j<dim;j++)
	  sum+=(v[i][j]*v[i][j]);
      }
    return sqrt(sum);
  }
  void apadd(const vect<std::vector<A>> &v1,long double alpha,const vect<std::vector<A>> &v2)
  {
    for(int i=0;i<dim;i++)
      {
	for(int j=0;j<dim;j++)
	  v[i][j]=v1.v[i][j]+(alpha*v2.v[i][j]);
      }

  }
  void print(std::vector<std::vector<A>>& ini)
  {
    for(int i=0;i<dim;i++)
      {
	for(int j=0;j<dim;j++)
	  ini[i][j]=v[i][j];
      }
  }
};


template<typename Y,typename X>
class mat<std::vector<Y>,X>
{
public:
  int dim;
  std::vector<std::vector<Y>> a;
  std::vector<int> ra;
  std::vector<int> ca;//for tridiagonal
  mat(std::vector<std::vector<Y>> a1,std::vector<int> ra1,std::vector<int> ca1,int d)
  {
    dim=d;
    for(int i=0;i<ra1.size();i++)
      ra.push_back(ra1[i]);
    for(int i=0;i<ca1.size();i++)
      {
	ca.push_back(ca1[i]);
	a.push_back(a1[i]);
      }
  }
  template<typename W>
  void md(const std::vector<Y>& ak,const W& vk1,W& vk2)
  {
    int n,l=0,k=0,co=0;
    for(int i=0;i<dim;i++)
      {
	n=ra[l+1]-ra[l];
	l+=1;
	co=0;
	while(co<n)
	  {
	    vk2[i]+=(ak[k]*vk1[ca[k]]); //multiplication of csr with vec
	    k++;
	    co++;
	  }
      }
  }
  template<typename W>
  void mo(const std::vector<Y>& ak,const W& vk1,W& vk2)
  {
    for(int i=0;i<dim;i++)
      vk2[i]+=ak[i]*vk1[i];
  }
  void  mulvec(vect<X> &v1)
  {
    vect<X> v2(dim);
    int i=0;
    md(a[0],v1.v[0],v2.v[0]);
    mo(a[1],v1.v[1],v2.v[0]);
    for(int i=1;i<dim-1;i++)
      {
	mo(a[1],v1.v[i-1],v2.v[i]);
	md(a[0],v1.v[i],v2.v[i]);
	mo(a[1],v1.v[i+1],v2.v[i]);
      }
	
    mo(a[1],v1.v[dim-2],v2.v[dim-1]);
    md(a[0],v1.v[dim-1],v2.v[dim-1]);
    /*std::cout<<"after multip:"<<std::endl;
    for(int i=0;i<dim;i++)
      {
	for(int j=0;j<dim;j++)
	  std::cout<<v2.v[i][j]<<" ";
	std::cout<<std::endl;
	}*/
    v1=v2;
       
  }
};
								
										 
									  
  
				       
				   
				     
      
				   




				       
				   
				     
      
				   


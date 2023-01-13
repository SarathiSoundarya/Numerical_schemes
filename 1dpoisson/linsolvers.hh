#include<iostream>
#include<cmath>
#include<vector>
#define err 0.0000007 
template<typename T>
class vect
{
 
public:
  int dim;
  T *ptr;
  long double norm;
  vect(int d)
  {
    dim=d;
    ptr=new T[dim];
    for(int i=0;i<dim;i++)
      ptr[i]=0;
    norm=0;
  }
  vect(std::vector<T> arr,int d)
  {
    dim=d;
    ptr=new T[dim];
    for(int i=0;i<dim;i++)
      ptr[i]=arr[i];
    long double sum=0; //explicit type sum
    for(int i=0;i<dim;i++)
      sum+=(arr[i]*arr[i]);
    norm=sqrt(sum);
  }
  void operator+=(vect<T> v1)
  {
    for(int i=0;i<dim;i++)
      ptr[i]+=v1.ptr[i];
    long double sum=0; //explicit type sum
    for(int i=0;i<dim;i++)
      sum+=(ptr[i]*ptr[i]);
    norm=sqrt(sum);
  }
  void operator-=(vect<T> v1)
  {
    for(int i=0;i<dim;i++)
      ptr[i]=ptr[i]-v1.ptr[i];
    long double sum=0; //explicit type sum
    for(int i=0;i<dim;i++)
      sum+=(ptr[i]*ptr[i]);
    norm=sqrt(sum);
  }
  void scalar(float s)
  {
    for(int i=0;i<dim;i++)
      ptr[i]=ptr[i]*s;
    norm*=s;
  }
  void operator=(vect<T> v1)
  {
    for(int i=0;i<dim;i++)
      ptr[i]=v1.ptr[i];
    norm=v1.norm;
  }
};


template<typename T>
class matrix
{
public:
  int dim;
  std::vector<T> al; //for lower triangle csr
  std::vector<int> rl,cl;
  std::vector<T> au; //for upper triangle csr
  std::vector<int> ru,cu;
  std::vector<T> d; // for diagonal matrix
  matrix(std::vector<T>& a,std::vector<int>& ra,std::vector<int>& ca,int& d1)
  {
    int l=0,k=0,co=0,suml,sumu,r,n;
    dim=d1;
    r=ra.size();
    ru.assign({0});
    rl.assign({0});
    for(int i=0;i<r-1;i++)
      {
	co=0;
	suml=0;sumu=0;
	n=ra[l+1]-ra[l];
	l+=1;
	while(co<n)
	  {
	    if(ca[k]<i)
	      {
		al.push_back(a[k]);
		cl.push_back(ca[k]);
		suml+=1;k++; co++;
	      }
	    else if(ca[k]>i)
	      {
		au.push_back(a[k]);
		cu.push_back(ca[k]);
		sumu+=1; k++;co++;
	      }
	    else
	      {
		d.push_back(a[k]);
		k++;co++;
	      }
	  
	  }
	rl.push_back(suml+rl.back());
	ru.push_back(sumu+ru.back());
      
      }  
  }

  void mulvec(std::vector<T> &v,std::vector<int> &r,std::vector<int> &c,vect<T> v1,int dim) //for multip A*vector
  {
    int n;
    float sum=0;
    int l=0;
    int  k=0,co=0;
    vect<T> vec(dim);
    for(int i=0;i<dim;i++)
      {
	n=r[l+1]-r[l];
	l+=1;
	co=0;sum=0;
	while(co<n)
	  {
	    sum+=(v[k]*v1.ptr[c[k]]); //multiplication of csr with vec
	    k++;
	    co++;
	  }
	vec.ptr[i]=sum;
      }
    v1=vec;
  }
  void addlowdia(std::vector<T>& v,std::vector<int> &ri,std::vector<int> &ci, std::vector<T> &d1,int dim) //adding lower diag & diag
  {
  
    int r=ri.size();
    int k=0,l=2,b=0,n;
    for(int i=0;i<r-1;i++)
      {
	if(i==0)
	  { v.insert(v.begin(),d1[b]);b++; 
	    for(int j=i+1;j<r;j++)
	      ri[j]+=1;
	    ci.insert(ci.begin(),0);
	  }
      
	else
	  {
	    n=ri[l];
	    k=n;
	    v.insert(v.begin()+k,d1[b]);b++;
	    for(int j=i+1;j<r;j++)
	      ri[j]+=1;
	    ci.insert(ci.begin()+k,i);
	    l+=1;
	  }
      }
    /* std::cout<<"inside lower lower siagonal:";
    for_each(v.begin(),v.end(),[](auto x)
      {std::cout<<x;});
      std::cout<<std::endl;
      for_each(ci.begin(),ci.end(),[](auto x)
      {std::cout<<x;}); */
  }

  void lowinv(std::vector<T> &v,std::vector<int> &ri,std::vector<int> &ci ,int dim)//only for lower tridiagonal matrix inverse
  {
    T a1=v[0];
    T b1=v[1];
    int r=ri.size();
    r=r-1;
    long double f=pow(a1,r); //explicit type double
    std::cout<<"f is:"<<f<<std::endl;
    int n=r-1;
    int s=0,k=0,sum=0;
    v.clear();
    ri.clear();
    ci.clear();
    ri.push_back(0);
    for(int i=0;i<r;i++)
      {
	s=0;
	for(int j=0;j<=i;j++)
	  {
	    v.push_back((1/f)*(pow(-1,i+j)*pow(a1,n-i+s)*pow(b1,i-s)));
	    k+=1;
	    s+=1;
	    ci.push_back(j); sum+=1;
	  }
	ri.push_back(sum);
      }
    /*   std::cout<<"inside lowinve:low matric"<<std::endl;
       int j=0;
       for_each(v.begin(),v.end(),[&j](auto x)
				  {std::cout<<x<<" ";
				    j++;});
       std::cout<<"dim is:"<<j<<" columns:"<<std::endl;
       j=0;
       for_each(ci.begin(),ci.end(),[&j](auto x)
				    {std::cout<<x<<" ";
				      j++;});
				      std::cout<<"coldim is:"<<j;*/
  
  }
  
  
};
  

class gauss_siedel
{
public:
  template<typename M,typename X>
  void apply(M A,std::vector<int> ra,std::vector<int> ca,X B,X& U,int dim)
  {
    typedef typename M::value_type real;
    typedef typename X::value_type real1;
    vect<real1> xo(dim);//initial guess
    vect<real1> b(B,dim);
    matrix<real> a(A,ra,ca,dim);
    vect<real1> error(dim);
    a.addlowdia(a.al,a.rl,a.cl,a.d,dim);// L+D
    a.lowinv(a.al,a.rl,a.cl,dim); // (L+D)^-1
    int n=0;
    do
      {
	error=xo;
	a.mulvec(a.au,a.ru,a.cu,xo,dim); //U*x(k)
	xo-=b; //-(b-U*x(k))
	xo.scalar(-1);
	a.mulvec(a.al,a.rl,a.cl,xo,dim);
	//	for(int i=0;i<dim;i++)
	// std::cout<<xo.ptr[i]<<" "<<std::endl;
	error-=xo;
	n++;
      }while(error.norm>err);
    std::cout<<"number of iterations:"<<n<<std::endl;
    for(int j=0;j<dim;j++)
      U.push_back(xo.ptr[j]);
  }
};



 
  
  

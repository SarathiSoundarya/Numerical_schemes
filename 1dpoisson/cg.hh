
#include<iostream>
#include<cmath>
#include<vector>
#define err 0.0000007
#include "vect.hh"

class cg
{
public:
  template<typename T, typename U>
  void apply(std::vector<T> A,std::vector<U> B,std::vector<U>& ini,std::vector<int> ra,std::vector<int> ca,int dim)
  {
    int i=1;
    mat<T,U> am(A,ra,ca,dim);
    vect<U> b(B,dim),xo(ini,dim),p(dim),u(dim),r(dim);
    long double alpha,beta,to,tn;
    u=xo;
    am.mulvec(xo);// Axo
    xo-=b;// Axo-b;
    xo.scalar(-1); //b-Axo 
    r=xo;//residual
    xo=u;//xo
    p=r;
    to=r*r;
    u=p;
    std::cout<<"r is:"<<std::endl;
    std::cout<<"r.norm():"<<r.norm()<<" ";
    do
      {
	am.mulvec(u);//Ap
	alpha=((to)/(p*u));// (rTr)/(pTAp)
	xo.apadd(xo,alpha,p);// xo + (alpha)p
	r.apadd(r,(-1*alpha),u);// r-(alpha)Ap
	tn=r*r;
	std::cout<<"iteration"<<i<<":"<<"norm: "<<r.norm()<<std::endl;
        if(r.norm()<err)
	 break;
	beta=(tn/to); // (r(k+1)Tr(k+1))/(r(k)Tr(k))
	p.apadd(r,beta,p);// r(k+1)+ beta*p
	u=p;
	to=tn;
        i++;
	}while(r.norm()>err);
    xo.print(ini);
  }
};

									
										 
									  
  
				       
				   
				     
      
				   


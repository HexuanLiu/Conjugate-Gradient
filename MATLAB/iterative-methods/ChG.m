function [ x_vec,r_list ] = ChG( A,b,x0,nmax )
%CHG CG, Algorithm 2.2
r=b-A*x0;
x=x0;
p=r;
s=A*p;
alpha=r'*r/(p'*s);

x_vec=[];
r_list=[];

for i=1:nmax
    x=x+alpha*p;
    x_vec=[x_vec,x];
    rp=r;
    r=r-alpha*s;
   
    
    

   % r=r-(r'*rp)/(rp'*rp)*rp;
    
   
    r_list=[r_list,r];
    w=A*r;
    beta=r'*r/(rp'*rp);
    alpha=r'*r/(w'*r-(beta/alpha)*r'*r);
   
  
    
    pp=p;
    p=r+beta*p;
    %p=p-(p'*A*pp)/(pp'*A*pp)*pp;

    
    s=w+beta*s;
end


end


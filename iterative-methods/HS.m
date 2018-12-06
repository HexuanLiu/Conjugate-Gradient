function [ x_vec,r_list ] = HS( A,b,x0,nmax )

x=x0;
r=b-A*x0;
p=r;
beta=r'*r;

x_vec=[];
r_list=[];


for i=1:nmax
    alpha=r'*r/(r'*A*p);
    x=x+alpha*p;
    x_vec=[x_vec,x];
    rp=r;
    r=r-alpha*A*p;
 

%     r=round(r,5);
   % r=r-(r'*rp)/(rp'*rp)*rp;
    r_list=[r_list,r];
    
    beta=r'*r/(rp'*rp);
    pp=p;
    p=r+beta*p;
 
end


end


function [ x_vec,r_list ] = HS_var( A,b,x0,nmax )

x=x0;
r=b-A*x0;
p=r;
x_vec=[];
r_list=[];


for i=1:nmax
    alpha=r'*r/(p'*A*p);
    x=x+alpha*p;
    x_vec=[x_vec,x];
    rp=r;
    r=r-alpha*A*p;
%     r=round(r,5);
% orthogonalize r_k agains r_(k-1)
     r=r-(r'*rp)/(rp'*rp)*rp;
    r_list=[r_list,r];
    
    beta=r'*r/(rp'*rp);
    pp=p;
    p=r+beta*p;
%     % orthogonalize p_k against Ap_(k-1)
%     p=p-(p'*A*rp)/(rp'*A*rp)*rp;
      p=p-(p'*A*pp)/(pp'*A*pp)*pp;
end


end


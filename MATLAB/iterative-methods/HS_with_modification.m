function [ x_vec ] = HS_with_modification(  A,b,x0,nmax )
%HS_with_MODIFICATION replace (3.8) by (3.9), FIG 3.1 

r=b-A*x0;
p=r;
x_vec=[];
x=x0;
for i=1:nmax
    alpha=r'*r/(p'*A*p);
    xp=x;
    x=x+alpha*p;
    x_vec=[x_vec,x];
    rp=r;
    r=r-alpha*A*p;
    beta=r'*r/(rp'*rp);
    p=r+beta/alpha*(x-xp);
end
end


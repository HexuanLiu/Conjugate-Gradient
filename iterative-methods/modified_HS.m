function [ x_vec ] = modified_HS( A,b,x0,nmax )
%MODIFIED_HS, Algorithm 3.1

r=b-A*x0;
x=x0;
p=r;
s=A*p;
x_vec=[];

for i=1:nmax
    alpha=(r'*r)/(p'*s);
    x=x+alpha*p;
    x_vec=[x_vec,x];
    rp=r;
    r=r-alpha*s;
    beta=(r'*r)/(rp'*rp);
    p=r+beta*p;
    s=A*r+beta*s;
end


end


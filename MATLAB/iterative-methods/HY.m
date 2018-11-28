function [ x_vec,r_vec ] = HY( A,b,x0,nmax )
%HY CG: Hegeman and Young (p. 143, applied iterative methods)
r=b-A*x0;
theta=1;
xold=0;
rold=0;
x=x0;
x_vec=[];
r_vec=[];

for i=1:nmax
    mu=r'*r/(r'*A*r);
    rnew=theta*(-mu*A*r+r)+(1-theta)*rold;
    xnew=theta*(mu*r+x)+(1-theta)*xold;
    munew=rnew'*rnew/(rnew'*A*rnew);
    thetanew=(1-munew/mu*(rnew'*rnew)/(r'*r)/theta)^(-1);
    
    mu=munew;
    theta=thetanew;
    
    rold=r;
    xold=x;
    
    r=rnew;
    x=xnew;
    
    r_vec=[r_vec, r];
    x_vec=[x_vec, x];
end


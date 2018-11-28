function [ x_vec, r_vec ] = GV( A,b,x0,nmax )
%GV (pipelined) CG, Algorithm 2.3
r=b-A*x0;
x=x0;
p=r;
s=A*p;
w=A*r;
z=A*w;
alpha=r'*r/(p'*s);
x_vec=[];
r_vec=[];

for i=1:nmax
    x=x+alpha*p;
    x_vec=[x_vec,x];
    
    rp=r;
    r=r-alpha*s;
    % orthogonalize r_k against r_{k-1}
    r=r-(r'*rp)/(rp'*rp)*rp;
    
    r_vec=[r_vec,r];
    w=w-alpha*z;
    
    q=A*w;
    beta=r'*r/(rp'*rp);
    alpha=r'*r/(w'*r-(beta/alpha)*r'*r);
    %pp=p;
    p=r+beta*p;
    %p=p-(p'*A*pp)/(pp'*A*pp)*pp;
    s=w+beta*s;
    z=q+beta*z;
end

end


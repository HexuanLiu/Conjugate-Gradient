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
    if (i==10 | i==20 | i==30 | i==40)
        i;
        r=b-A*x;
    end
    
    r_vec=[r_vec,r];
    w=w-alpha*z;
    
    q=A*w;
    beta=r'*r/(rp'*rp);
    alpha=r'*r/(w'*r-(beta/alpha)*r'*r);
    p=r+beta*p;
    s=w+beta*s;
    z=q+beta*z;
end

end


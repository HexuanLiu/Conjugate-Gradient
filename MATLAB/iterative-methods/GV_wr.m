function [ x_vec, r_vec ] = GV_wr( A,b,x0,nmax )
%GV (pipelined) CG, Algorithm 2.3
r0=b-A*x0;
r=r0;
x=x0;
p=r;
s=A*p;
w=A*r;
z=A*w;
alpha=r'*r/(p'*s);
x_vec=[];
r_vec=[];
count=0;
for i=1:nmax
    x=x+alpha*p;
    x_vec=[x_vec,x];
    rp=r;
    r=r-alpha*s;
    r_vec=[r_vec,r];
    w=w-alpha*z;
  % replace w 
    if (norm(r)/norm(r0)>1e-1)
        count=count+1;
        w=A*r;
    end
    q=A*w;
    beta=r'*r/(rp'*rp);
    alpha=r'*r/(w'*r-(beta/alpha)*r'*r);
    p=r+beta*p;
    s=w+beta*s;
    z=q+beta*z;
    
%     if (norm(r)<1e-6)
%         break
%     end
end
count
end


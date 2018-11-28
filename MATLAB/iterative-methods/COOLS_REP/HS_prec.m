function [ x_vec,r_list, error_list ] = HS_prec( A,Minv, b,x0,nmax, xtrue )

x=x0;
r=b-A*x0;
u=Minv*r;
p=u;
x_vec=[];
r_list=[];
bn=norm(b,2);
error_list=[];

for i=0:nmax
   s=A*p;
   alpha=(r'*u)/(s'*p);
   xnew=x+alpha*p;
   rnew=r-alpha*s;
   unew=Minv*rnew;
   betanew=(rnew'*unew)/(r'*u);
   pnew=unew+betanew*p;
   
   x=xnew;
   r=rnew;
   u=unew;
   beta=betanew;
   p=pnew;
   
   x_vec=[x_vec,x];
   r_list=[r_list,r];
%    relres=norm(b-A*x,2)/bn;
%    relres_list=[relres_list, relres];
    error=sqrt(norm((x-xtrue)'*A*(x-xtrue)));
   error_list=[error_list, error];
end
end


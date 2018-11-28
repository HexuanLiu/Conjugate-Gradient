function [ x_vec,r_list, error_list ] = prec_CGCG( A, Minv, b, x0, nmax, xtrue)
x=x0;
r=b-A*x0;
u=Minv*r;
w=A*u;
alpha=(r'*u)/(w'*u);
beta=0;
gamma=r'*u;
pold=r;
sold=A*pold;

x_vec=[];
r_list=[];
bn=norm(b,2);
error_list=[];
tol=1e-16;


for i=0:nmax
    p=u+beta*pold;
    s=w+beta*sold;
    xnew=x+alpha*p;
    rnew=r-alpha*s;
    unew=Minv*rnew;
    wnew=A*unew;
    gammanew=rnew'*unew;
    delta=wnew'*unew;
    betanew=gammanew/gamma;
    alphanew=(delta/gammanew-betanew/alpha)^(-1);
    
%     relres=norm(b-A*x,2)/bn;
%     relres_list=[relres_list, relres];

    error=sqrt(norm((x-xtrue)'*A*(x-xtrue)));
    error_list=[error_list, error];
    

    pold=p;
    sold=s;
    x=xnew;
    r=rnew;
    u=unew;
    w=wnew;
    gamma=gammanew;
    beta=betanew;
    alpha=alphanew;
    
end
   

end
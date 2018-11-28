function [ x_vec,r_list, error_list ] = prec_pCG( A,Minv, b,x0,nmax, xtrue )
x=x0;
r=b-A*x;

u=Minv*r;
w=A*u;
x_vec=[];
r_list=[];
bn=norm(b,2);
error_list=[];


for i=0:nmax
    gamma=r'*u;
    delta=w'*u;
    m=Minv*w;
    v=A*m;
    if(i>0)
        beta=gamma/gammaold;
        alpha=(delta/gamma-beta/alphaold)^(-1);
    else
        beta=0;
        alpha=gamma/delta;
        zold=v;
        qold=m;
        sold=w;
        pold=u;
    end
    z=v+beta*zold;
    q=m+beta*qold;
    s=w+beta*sold;
    p=u+beta*pold;
    xnew=x+alpha*p;
    rnew=r-alpha*s;
    unew=u-alpha*q;
    wnew=w-alpha*z;
    
    zold=z;
    qold=q;
    sold=s;
    pold=p;
    gammaold=gamma;
    alphaold=alpha;
    x=xnew;
    r=rnew;
    u=unew;
    w=wnew;
    
    x_vec=[x_vec,x];
    r_list=[r_list,r];
%     relres=norm(b-A*x,2)/bn;
%     relres_list=[relres_list, relres];

    error=sqrt(norm((x-xtrue)'*A*(x-xtrue)));
    error_list=[error_list, error];
end
end
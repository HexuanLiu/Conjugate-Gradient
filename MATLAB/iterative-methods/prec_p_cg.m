function [ x_vec,r_vec ] = prec_p_cg( A,Minv,b,x0, nmax)
r=b-A*x0;
u=Minv*r;
w=A*u;
x=x0;

x_vec=[x0];
r_vec=[r0];
for i=0:nmax
    gamma=r'*u;
    delta=w'*u;
    m=Minv*w;
    v=A*m;
    if (i>0)
        beta=gamma/gamma_old;
        alpha=(delta/gamma-beta/alpha_old)^(-1);
    else
        beta=0;
        alpha=gamma/delta;
    end
    z=v+beta*z_old;
    q=m+beta*q_old;
    s=w+beta*s_old;
    p=u*beta*p_old;
    x_new=x+alpha*p;
    r_nw=r-alpha*s;
    u_new=u-alpha*q;
    w_new=w-alpha*z;
end

end
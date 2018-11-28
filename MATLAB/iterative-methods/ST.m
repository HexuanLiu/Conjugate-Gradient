function [ x_vec,r_vec ] = ST( A,b,x0,nmax )
% function of the finite-precision ST (algorithm 2.1)
%   Input: SPD matrix A, right-hand side vector b, initial approximation x0, maximum number of
%   iterations nmax
%   Output: approximate solutions x_vec in each iteration

r=b-A*x0;

x=x0;
xold=x;
rold=r;
eold=0;

x_vec=[];
r_vec=[];

dx=0;
dr=0;

for i=1:nmax
    q=r'*A*r/(r'*r)-eold;
    %%deltax=x-xold;
    %%deltar=r-rold;
    dx_new=1/q*(r+eold*dx);
    xnext=x+dx_new;
    dr_new=1/q*(-A*r+eold*dr);
    rnext=r+dr_new;
    
    dr=dr_new;
    dx=dx_new;
    
    
    xold=x;
    rold=r;
    x=xnext;
    r=rnext;
    
    x_vec=[x_vec,x];
    r_vec=[r_vec,r];
    eold=q*(r'*r)/(rold'*rold);
end


end




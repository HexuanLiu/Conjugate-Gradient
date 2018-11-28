function [ x_vec,r_list ] = HS_A( A,b,x0,nmax )
% function of the finite-precision HS
%   Input: SPD matrix A, right-hand side vector b, initial approximation x0, maximum number of
%   iterations nmax
%   Output: approximate solutions x_vec in each iteration

x=x0;
r=b-A*x0;
p=r;
x_vec=[];
r_list=[];


for i=1:nmax
    alpha=r'*r/(r'*A*p);
    x=x+alpha*p;
    x_vec=[x_vec,x];
    rp=r;
    r=r-alpha*A*p;
    r_list=[r_list,r];
    
    beta=r'*r/(rp'*rp);
    p=r+beta*p;
end


end

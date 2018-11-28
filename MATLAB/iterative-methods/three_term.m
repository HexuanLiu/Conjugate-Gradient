clear all; close all; clc

load('bcsstk03.mat')
sA=Problem.A;
A=full(sA);
N=112;
[V,D]=eig(A);
b=V*ones(112,1);
b=b/norm(b,2);
x=A\b;
x0=zeros(N,1);

iter=1200;
[x_ST,r1]=ST(A,b,x0,iter);
[x_HY,r2]=HY(A,b,x0,iter);
[x_HS,r3]=HS(A,b,x0,iter);
ST_list=[];
HY_list=[];
HS_list=[];

for i=1:iter
    err_ST=sqrt(norm((x_ST(:,i)-x)'*A'*(x_ST(:,i)-x)));
    ST_list=[ST_list,err_ST];
    
    err_HY=sqrt(norm((x_HY(:,i)-x)'*A'*(x_HY(:,i)-x)));
    HY_list=[HY_list,err_HY];
    
    err_HS=sqrt(norm((x_HS(:,i)-x)'*A'*(x_HS(:,i)-x)));
    HS_list=[HS_list,err_HS];

end
semilogy(ST_list)
hold on;
semilogy(HY_list) 
semilogy(HS_list)
legend('ST','HY','HS')
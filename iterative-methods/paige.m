clear all; close all; clc

%% Lanczos
u=A*v;
alpha=v'*v;
w=u-v*alpha;
beta=sqrt(w'*w);
if (beta< 1e-16)
    v=w/beta;
    u=A*v-v*beta;
end


%% example 6.1
U=randn(n,n);
Lam=zeros(n,1);
for i =1:9
    Lam(i)=0.00001*i;
end
Lam(10)=1;
Lam=diag(Lam);

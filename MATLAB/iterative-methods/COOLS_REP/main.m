clear all; close all; clc

load('bcsstk17s.mat')
sA=Problem.A;
A=full(sA);
N=size(A,1);
xtrue=1/sqrt(N)*ones(N,1);
b=A*xtrue;
M=sparse(diag(diag(A)));
Minv=eye(N);
% opts.diagcomp=0.5;
% opts.type = 'nofill';
% opts.michol = 'on';
% L=(ichol(sA)^(-1))';
% Minv=L*L';


x0=zeros(N,1);
nmax=5000;

[x1,r1, err1]=HS_prec(A,Minv,b,x0, nmax,xtrue);
[x2,r2, err2]=prec_CGCG(A,Minv,b,x0, nmax, xtrue);
[x3,r3, err3]=prec_pCG(A,Minv,b,x0, nmax, xtrue);
semilogy(err1)
hold on
semilogy(err2)
semilogy(err3)
legend('HS','CGCG','pCG')
title('bcsstk03')
clear all; 

%%use the model problem
rho=0.8;
n=48;
N=48;
lambda1 = .001; lambdan = 1;
lambda = lambda1*ones(n,1);
for i=2:n, lambda(i) = lambda(1) + ((i-1)/(n-1))*(lambdan-lambda1)*rho^(n-i); end;


Lambda = diag(lambda);
[U,R] = qr(randn(n,n));
A = U*Lambda*U';
b=randn(n,1);
% Make sure A is perfectly symmetric.
for i=1:n-1, for j=i+1:n, A(i,j)=A(j,i); end; end;
x=A\b;
x0=zeros(N,1);

% % n=112
% 
% load('bcsstk03.mat')
% sA=Problem.A;
% 
% A=full(sA);
% N=size(A,2);
% [V,D]=eig(A);
% b=V*ones(N,1);
% b=b/norm(b,2);
% 
% x=A\b;
% x0=zeros(N,1);


%% % form poisson matrix A:
% m=60;
% h=1/(m+1);
% I = speye(m);
% e = ones(m,1);
% T = spdiags([e -4*e e],[-1 0 1],m,m);
% S = spdiags([e e],[-1 1],m,m);
% A = (kron(I,T) + kron(S,I)) / h^2;
% x=A\b;
% x0=zeros(N,1);

%%
iter=150;
% [x_HS,r0]=HS(A,b,x0,iter);
% [x_HSv,r0]=HS_var(A,b,x0,iter);
[x_ChG,r]=ChG(A,b,x0,iter);
[x_GV,r_rec]=GV(A,b,x0,iter);
[x_GVwr,r1]=GV_wr(A,b,x0,iter);

% % s step CG
% %basis_type: string denoting which basis to use. Acceptable values are
% %'monomial', 'newton', or 'chebyshev'. 
% s=2;
% results=cacg(A,b,s,x0,iter,[],'chebyshev');
% x_s=results.xlist';



HS_list=[];
HSV_list=[];
ChG_list=[];
GV_list=[];
GVwr_list=[];
s_list=[];

% [flag,errlist,result] = exact_cg(x, x0, A,b,iter, 1e-16);

for i=1:iter
    

%     err=sqrt(norm((x_HS(:,i)-x)'*A'*(x_HS(:,i)-x)));
%     HS_list=[HS_list,err];
     
%     err1=sqrt(norm((x_HSv(:,i)-x)'*A'*(x_HSv(:,i)-x)));
%     HSV_list=[HSV_list,err1];
     
    err2=sqrt(norm((x_ChG(:,i)-x)'*A'*(x_ChG(:,i)-x)));
    ChG_list=[ChG_list,err2];
    
    err3=sqrt(norm((x_GV(:,i)-x)'*A'*(x_GV(:,i)-x)));
    GV_list=[GV_list,err3];
    
    err4=sqrt(norm((x_GVwr(:,i)-x)'*A'*(x_GVwr(:,i)-x)));
    GVwr_list=[GVwr_list,err4];
    
%     err5=sqrt(norm((x_re(:,i)-x)'*A'*(x_re(:,i)-x)));
%     re_list=[re_list,err5];

%     err6=sqrt(norm((x_s(:,i)-x)'*A'*(x_s(:,i)-x)));
%     s_list=[s_list,err6];
end
% semilogy(HS_list) 

% semilogy(HSV_list) 
semilogy(ChG_list)
hold on
semilogy(GV_list)
semilogy(GVwr_list)
% semilogy(errlist)
% semilogy(re_list)
semilogy(s_list)
legend('CGCG','GV','GVwr')
t = title('bcsstk03', 'Units', 'normalized', 'Position', [0.5, -0.1, 0]); 

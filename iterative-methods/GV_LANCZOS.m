load('bcsstk03.mat')
sA=Problem.A;

A=full(sA);
N=size(A,2);
[V,D]=eig(A);
b=V*ones(N,1);
b=b/norm(b,2);
x=A\b;
x0=zeros(N,1);

% Make sure A is perfectly symmetric.
for i=1:n-1, for j=i+1:n, A(i,j)=A(j,i); end; end;

x0=zeros(n,1);
x=x0;
r=rhs-A*x0;
p=r;
q1=r/norm(r);
Q = zeros(n,J+1);
T = zeros(J+1,J);
Q(:,1) = q1;
q2=A*q1;
alpha1=q2'*q1;
q2=q2-alpha1*q1;
beta1=norm(q2);
alpha_list=zeros(1,J);
beta_list=zeros(1,J);
alpha_list(1)=alpha1;
beta_list(1)=beta1;
a_list=[];
b_list=[]; 
r_list=[r];


clear all; close all; clc


%% set up precision digits. Defaut=10
digits(64);

%% finite precision CG for A
% set up A, b, x_exact
n=24;
l1=0.1;
kappa=1000;
ln=l1*kappa;

lambda=zeros(n,1);
lambda(1)=l1;
lambda(n)=ln;

b=randn(n,1);



%% b with peak in the middle

d=0.2;
m=11;

bn3=[];
for i=1:n
    bnew=zeros(m,1);
    sumb=b(i)^2;
    bnew(6)=sqrt(sumb*(1-d));
    sumb=sumb*d;
    for j=1:4
        bnew(6-j)=sqrt(sumb*(1-d)/2);
        bnew(6+j)=sqrt(sumb*(1-d)/2);
        sumb=sumb*d;
    end
    bnew(1)=sqrt(sumb/2);
    bnew(11)=sqrt(sumb/2);
    bn3=[bn3; bnew];
end

%% check if the sum of sqaures are equal

s1=0;
for i=1:length(b)
    s1=s1+b(i)^2;
end



s3=0;
for i=1:length(bn3)
    s3=s3+bn3(i)^2;
end


for rho=[0.4, 0.6, 0.8, 0.9, 1.0]


    for i=2:n-1
        lambda(i)=l1+((i-1)/(n-1))*(ln-l1)*rho^(n-i);
    end
%     A=diag(lambda);
%     A=vpa(A);
%     b=vpa(b);
%     % exact solution
%     x=vpa(A\b);
%     
%     x0=vpa(zeros(24,1));

%     iter=40;
%     [flag, errlist, xiter]=conjgrad(x, x0, A, b, iter,1e-12);
    


%% exact CG for A

%     iter=40;
%     [flag, errlist_exact, xiter]=exact_cg(x, x0, A, b, iter,1e-12);
       



 %% exact CG for A hat
delta=1e-12;
m=11;

l=zeros(n,m);
for i=1:n
    for j=1:m
        l(i,j)=lambda(i)+(j-(m+1)/2)/(m-1)*delta;
    end
end

nlambda=[];
for i=1:n
    nlambda=[nlambda, l(i,:)];
end
Ahat=diag(nlambda);

% % original setup: b evenly spread out over m
% bn=[];
% for i=1:n
%     bnew=sqrt(b(i)^2/m);
%     bn=[bn; bnew*ones(m,1)];
% end


% iter=50;
% [flag, errlist_beven, xiter]=exact_cg(xhat, x0, Ahat, bn, iter,1e-12);

 


%  % b decrease by 1/d
% bn1=[];
% d=0.2;
% for i=1:n
%     bnew=zeros(m,1);
%     sumb=b(i)^2;
%     for j=1:m
%         if (j<m)
%             bnew(j)=sqrt(sumb*(1-d));
%             sumb=sumb*d;
%         end
%         if (j==m)
%             bnew(j)=sqrt(sumb);
%         end
%     end
%     bn1=[bn1; bnew];
% end


 %% b increase by 1/d
% m=11;
% n=24;
% bn2=[];
% d=0.2;
% for i=1:n
%     bnew=zeros(m,1);
%     sumb=b(i)^2;
%     for j=1:m
%         if (j<m)
%             bnew(m-j+1)=sqrt(sumb*(1-d));
%             sumb=sumb*d;
%         end
%         if (j==m)
%             bnew(m-j+1)=sqrt(sumb);
%         end
%     end
%     bn2=[bn2; bnew];
% end
% 




Ahat=vpa(Ahat);
bn3=vpa(bn3);
xhat=vpa(Ahat\bn3);

    x0=vpa(zeros(m*n,1));

    iter=50;
    [flag, errlist_bd, xiter]=exact_cg(xhat, x0, Ahat, bn3, iter,1e-12);

    semilogy(errlist_bd)
    ylim([1e-12 1])
    hold on
end

legend('0.4', '0.6','0.8','0.9','1.0')
xlabel('Iteration')
ylabel('A norm of error')



function [ xlist, rlist ] = sstep( A, b, x0, nmax, s)
%SSTEP Summary of this function goes here
%   Detailed explanation goes here
r0=b-A*x0;
r=r0;
x=x0;
xlist=[x];
rlist=[norm(r)];

for i=[0:nmax]
    k=i*s;
    T=[r];
    temp=r;
    for j=1:s
        temp=A*temp;
        T=[T temp];
    end
    R=T(:,1:end-1);
    Q=T(:,2:end);
    
    if i==0
        P=R;
    else
        C=-Q'*P;
        B=W\C;
        P=R+P*B;
    end
    W=Q'*P;
  
    g=P'*r;
    a=W\g;
    x=x+P*a;
    xlist=[xlist,x];
    r=b-A*x;
    rlist=[rlist,norm(r)];
end


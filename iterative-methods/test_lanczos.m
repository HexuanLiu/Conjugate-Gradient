clear all; close all; clc

% x = chebfun('x');
% f = abs(x-0.5);
% [p,err] = remez(f,10);
% LW = 'linewidth'; FS = 'fontsize'; fs = 14;
% figure, plot(f,'b',p,'r',LW,1.6)
% title('Function and best polynomial approximation',FS,fs)
z=24;
lamb=zeros(20,z);
for i=1:20
    for j=1:z
        lamb(i,j)=1/(2^i)+j/(z*2^i);
    end
end
A=diag(reshape(lamb, 24*20,1));
x=rand(480,1);
b=A*x;

err_list=[];
for iter=1:480
    x_approx=lanczos(A,b, @(x)1/x, iter);
    err=norm(x-x_approx, inf);
    err_list=[err_list,err];
end
plot((1:480), err_list)
    
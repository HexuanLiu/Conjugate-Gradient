function [x_list, r_list ] = recurrence( A,b,x0,nmax )

x_list=[];
r_list=[];
r0=b-A*x0;
r_list=[r_list,r0];


p0=r0;
a0=r0'*r0/(p0'*A*p0);
r1=r0-a0*A*r0;
r_list=[r_list,r1];
pp=p0;
ap=a0;
x1=x0+a0*p0;
x=x1;
x_list=[x_list,x0,x1];

elist=[norm(r0),norm(r1)];

for j=2:nmax
    r=r_list(:,j);
    rp=r_list(:,j-1);
    b=r'*r/(rp'*rp);
    p=r+b*pp;
    a=r'*r/(p'*A*p);
    
    x=x+a*p;
    x_list=[x_list,x];
    
    rnew=r-a*A*r-a*b/ap*(rp-r);
    
    err=norm(rnew);
    elist=[elist,err];
    
    r_list=[r_list,rnew];
    
    pp=p;
    ap=a;
end
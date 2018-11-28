function [flag, errlist, x] = conjgrad(xexact, x0, A,b,iter, tol)
    
    errlist=[];
    flag=0;
    x=x0;
  
    err1=0;
    
    r=b-A*x;
    p=r;
    for k=1:iter
        
        if (sqrt(norm((xexact-x)'*A'*(xexact-x))))<tol
            flag=1;
            return
        end
        y=A*p;
        s=r'*r;
        a=s./(p'*y);
        x=x+a*p;
        r=r-a*y;
        b=(r'*r)./s;
        p=r+b*p;
        if (k==1)
            err1=sqrt(norm((xexact-x)'*A'*(xexact-x)));
        end
        err=sqrt(norm((xexact-x)'*A'*(xexact-x)))/err1;
        errlist=[errlist,err];
 
    end

 end


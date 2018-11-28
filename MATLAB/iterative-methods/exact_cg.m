function [flag,errlist, x] = exact_cg(xexact, x0, A,b,iter, tol)
    
    errlist=[];
    err1=0;
    m=0;
    flag=0;
    x=x0;
    rlist=[];
    r=b-A*x;
    rlist=[rlist; r];
    m=m+1;

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
        
        % re-orthogonalize
        for kount=1:2
            for j=1:m
                rj=rlist(:,j);
                r=r-((r'*rj)/(rj'*rj))*rj;
            end
        end
        
        
        rlist=[rlist, r];
        m=m+1;
        b=(r'*r)./s;
        
        pp=p;
        p=r+b*p;
     
        
        if (k==1)
            err1=sqrt(norm((xexact-x)'*A'*(xexact-x)));
        end
        err=sqrt(norm((xexact-x)'*A'*(xexact-x)));
        errlist=[errlist,err];
    end

 end


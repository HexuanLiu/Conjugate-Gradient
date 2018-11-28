%  Algorithm 2.  Conjugate Gradient Algorithm (CG).
%                (For Hermitian positive definite matrices)

%  User enters HPD matrix A, b, x0, and itmax.
%  Returns residual norms in array resid and norms of updated
%  residuals in array resest.

Asave = A;
numdec = input('Enter number of decimal digits: ');
digits(numdec);
A = vpa(A); b = vpa(b); x0 = vpa(x0);
A = (A+A')/2;    % Make sure that A is symmetric.
n = length(b);
resid = vpa(zeros(itmax+1,1));
resest = vpa(zeros(itmax+1,1));

xk = x0;
rk = b - A*xk;
%rknrm = norm(rk),
rknrm = sqrt(rk'*rk); drknrm = (rknrm),
resid(1) = rknrm;
resest(1) = rknrm;
pk = rk;
rkrk = rk'*rk;
aks = zeros(itmax,1); bks = zeros(itmax,1);
errlist=[];

for k=1:itmax,
  Apk = A*pk;
  aknum = rkrk;  akden = Apk'*pk;
  ak = aknum/akden;
  aks(k) = ak;
  xk = xk + ak*pk;
  rk = rk - ak*Apk;
%  rknrm = norm(rk);
  rknrm = sqrt(rk'*rk); drknrm = (rknrm);
  resest(k+1) = rknrm;
  bmAxk = b - A*xk; resid(k+1) = sqrt(bmAxk'*bmAxk);
  rkrk = rk'*rk;
  bk = rkrk/aknum;
  bks(k) = bk;
  pk = rk + bk*pk;
  
  err=sqrt(norm((xk-x)'*A'*(xk-x)));
  errlist=[errlist,err];
end;
A = Asave;

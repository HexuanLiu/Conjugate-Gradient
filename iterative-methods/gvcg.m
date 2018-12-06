%  Pipelined CG

%  User enters HPD matrix A, b, x0, and itmax.
%  Returns residual norms in array resid, norms of updated residuals in array resest.
%  If user sets flag=1 and supplies the true solution x_true, will also return 
%  the A-norm of the error at each step in array errA.
%  Looks at recurrence for zk = (-1)^k rk/||rk|| to see how closely it is satisfied
%  and how nearly orthogonal successive vectors are.  This determines size of
%  intervals in G. paper.  Returns arrays Tk and Zk.

n = length(b);
resid = zeros(itmax+1,1); resest = zeros(itmax+1,1);
xkdiff = zeros(itmax,1); e1 = zeros(itmax,1); e1(1) = 1;
if flag==1, errA = zeros(itmax+1,1); errAest = zeros(itmax+1,1); end;
Tk = zeros(itmax+1,itmax); Zk = zeros(n,itmax+1);

%  Initialization
xk = x0;
rk = b - A*xk;
resid(1) = norm(rk); resest(1) = resid(1);
if flag==1, diff = x_true - x0; errA(1) = sqrt(diff'*A*diff); errAest(1) = sqrt(rk'*(A\rk)); end;
Zk(:,1) = rk/resest(1);
pk = rk;
sk = A*pk; 
wk = sk; zk = A*wk; 
aknum = rk'*rk;
akden = sk'*pk; ak = aknum/akden;
Tk(1,1) = 1/ak;
start = 1;

%  Iteration
for k=1:itmax,
  xk = xk + ak*pk;
%      Check to see if xk = x0 + Zk*(Tk\(resid(1)*e1)).
    xkchk = x0 + Zk(:,1:k) * (Tk(1:k,1:k)\(resid(1)*e1(1:k)));
    xkdiff(k) = norm(xk-xkchk)/norm(xk);
  rk = rk - ak*sk;
%  rk = b-A*xk;    %  Try residual replacement!
  wk = wk - ak*zk;
%  if fix(k/10)*10==k, wk = A*rk; end;    %  See what happens if you make this modification!
%  gknorm = norm(wk - A*rk)/norm(A*rk);
%  if gknorm > 1.e-8, k, gknorm, pause, wk = A*rk; end;      %  Try adaptively replacing wk.
  if resest(k)/resest(start) < 1.e-2, start = k; wk = A*rk; k, pause, end;
%  wk = A*rk;
  
  resest(k+1) = norm(rk); resid(k+1) = norm(b-A*xk);
  if flag==1, diff = x_true - xk; errA(k+1) = sqrt(diff'*A*diff); errAest(k+1) = sqrt(rk'*(A\rk)); end;
  Zk(:,k+1) = (-1)^k*rk/resest(k+1);
  Tk(k+1,k) = resest(k+1)/(ak*resest(k));
  if k < itmax, Tk(k,k+1) = Tk(k+1,k); end;
  qk = A*wk;
  bknum = rk'*rk; bk = bknum/aknum; aknum = bknum;
  akm1 = ak;
  ak = aknum/((rk'*wk) - (bk/ak)*aknum);
  if k < itmax, Tk(k+1,k+1) = (1/ak + bk/akm1); end;
  pk = rk + bk*pk;
  sk = wk + bk*sk;
  zk = qk + bk*zk;
end;

Gk = A*Zk(:,1:itmax) - Zk*Tk;
gknorm = 0;
for j=1:itmax,
  colnorm = norm(Gk(:,j));
  if colnorm > gknorm, gknorm = colnorm; end;
end;
maxinprod = 0;
for k=2:itmax
  inprod = abs(Tk(k+1,k)*Zk(:,k-1)'*Zk(:,k));
  if inprod > maxinprod, maxinprod = inprod; end;
end;
normA = norm(full(A));
eps1 = gknorm/normA, eps2 = maxinprod/normA
eps3 = max(xkdiff)

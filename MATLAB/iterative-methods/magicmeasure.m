clear
clf
%  Run the Lanczos algorithm using finite precision arithmetic.
%  Then look for a weight function for which exact Lanczos 
%  generates the same tridiagonal matrix.

%  Use model problem with parameters supplied by user.

n = input('Enter n: ');
rho = input('Enter rho: ');
J = input('Enter number of steps to run: ');

%  Plot eigenvalues.

lambda1 = .001; lambdan = 1;
lambda = lambda1*ones(n,1);
for i=2:n, lambda(i) = lambda(1) + ((i-1)/(n-1))*(lambdan-lambda1)*rho^(n-i); end;
figure(1)
subplot(3,1,1)
plot(lambda,zeros(n,1),'*-'), axis([lambda(1) lambda(n) -1 1]);
title('Eigenvalues')

% finite precision cg

% %  Run J steps of Lanczos algorithm for A = U*diag(lambda)*U'
% %  with random starting vector q1.
% 
Lambda = diag(lambda);
[U,R] = qr(randn(n,n));
A = U*Lambda*U';
rhs=randn(n,1);
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

% HS CG
for j=1:J
    a=(r'*r)/(p'*A*p);
    xnew=x+a*p;
    rnew=r-a*A*p;
    b=(rnew'*rnew)/(r'*r);
    pnew=rnew+b*p;  
    
    x=xnew;
    r=rnew;
    p=pnew;
    
    a_list=[a_list,a];
    b_list=[b_list,b];
    r_list=[r_list,r];
    size(r_list);
end


% 
for k=2:J
   alpha_list(k)=1/a_list(k)+b_list(k-1)/a_list(k-1);
   beta_list(k)=norm(r_list(:,k+1))/(norm(r_list(:,k))*a_list(k));
   %beta_list(k)=sqrt(b_list(k))/a_list(k);
end
for k=1:J
    Q(:,k+1)=(-1)^(k)*r_list(:,k+1)/norm(r_list(:,k+1));
    T(k,k)=alpha_list(k);
    T(k,k+1)=beta_list(k);
    T(k+1,k)=beta_list(k);
end
T=T(1:J+1,1:J);
%
r1=rhs-A*x0;
q1 = r1/norm(r1);
Q1 = zeros(n,J+1); Q1(:,1) = q1;
T1 = zeros(J+1,J);

for j=2:J+1,
  Q1(:,j) = A*Q1(:,j-1);
  if j > 2, Q1(:,j) = Q1(:,j) - betaj*Q1(:,j-2); end;
  alphaj = Q1(:,j)'*Q1(:,j-1);
  T1(j-1,j-1) = alphaj;
  Q1(:,j) = Q1(:,j) - alphaj*Q1(:,j-1);
  betaj = norm(Q1(:,j));
  T1(j,j-1) = betaj;
  if j < J+1, T1(j-1,j) = betaj; end;
  Q1(:,j) = Q1(:,j)/betaj;
end;

%  Plot eigenvalues of T in bins about eigenvalues of A.

eigA = sort(eig(A));
normA = max([abs(eigA(1)); abs(eigA(n))]);

width = 100*eps*normA;                                       % Bin width.
egroup(1,1) = eigA(1);  egroup(1,2) = eigA(1); ngroups = 1;  % Group eigenvalues that are too close.
for i=2:n,
  if eigA(i)-width <= egroup(ngroups,2)+width,
    egroup(ngroups,2) = eigA(i);
  else
    ngroups = ngroups+1;
    egroup(ngroups,1) = eigA(i); egroup(ngroups,2) = eigA(i);
  end;
end;

edges = zeros(2*ngroups+2,1);
edges(1) = -Inf; edges(2*ngroups+2) = Inf;
for i=1:ngroups,
  edges(2*i) = egroup(i,1) - width;
  edges(2*i+1) = egroup(i,2) + width;
end;

eigT = sort(eig(T(1:J,1:J)));
NT = histc(eigT,edges);
subplot(3,1,2)
bar(NT)
title(['Step ',int2str(J)])
pause(.1)

%  Check.
errorth = norm(eye(J,J) - Q(:,1:J)'*Q(:,1:J),'fro');
fprintf('\n')
fprintf('Norm of I - QJtrans*QJ = %12.4d\n',errorth)
normFJ = norm(A*Q(:,1:J) - Q*T,'fro');
fprintf('Norm of FJ = %12.4d\n',normFJ)
%%
%--------------------------------------------------------------------------------------------
% exact Lanczos
%  Now continue 3-term recurrence using multiple precision arithmetic and orthogonalizing
%  future vectors against each other and against the unconverged Ritz vectors from the
%  finite precision computation.

Asave = A; Tsave = T; Qsave = Q; eigAsave = eigA; normAsave = normA;

digits(64);
A = vpa(A); T = vpa(T); Q = vpa(Q);
eigA = sort(eig(A));
normA = abs(eigA(n));  if abs(eigA(1)) > normA, normA = abs(eigA(1)); end;
[S, Theta] = eig(T(1:J,1:J));
for i=1:J,
  if S(J,i) < 0, S(:,i) = -S(:,i); end;   % Make bottom entries of S positive.
end;
Y = Q(:,1:J)*S;

%  Print out lower bound on distance of eigenvalues of any extension of T
%  to eigenvalues of A.

[eigT,indx] = sort(diag(Theta));
Y = Y(:,indx); S = S(:,indx);
lower_bound = vpa(0);

%  Lower bound based on interlacing of roots of orthogonal polynomials.
i = 1;
for l=1:J,
  if i==1 & eigT(l) < eigA(1),
    dist = eigA(1)-eigT(l);
    if dist > lower_bound, lower_bound = dist; end;
  elseif i==n & eigT(l) > eigA(n),
    dist = eigT(l)-eigA(n);
    if dist > lower_bound, lower_bound = dist; end;
  else
    while i < n & eigA(i) <= eigT(l), i=i+1; end;
%      eigA(i-1) <= eigT(l) < eigA(i)
    if l < J & eigT(l+1) < eigA(i),
      dist1 = eigT(l) - eigA(i-1); dist2 = eigA(i) - eigT(l+1);
      dist = dist1; if dist2 < dist, dist = dist2; end;
      if dist > lower_bound, lower_bound = dist; end;
      i = i-1;
    end;
  end;
end;
fprintf('Lower bound on maxdist = %12.4d\n\n',double(lower_bound))
    
%  Try to find a set of nearly orthonormal Ritz vectors that T(J+1,J)*Q(:,J+1) is
%  almost orthogonal to and that T(J+1,J)*Q(:,J) is almost a linear combination of.

%    First replace clustered Ritz vectors by cluster vectors.

%cluster_width = vpa(sqrt(eps)*normA);
cluster_width = vpa(width);
JJ = 0; clusterflag = 0;
for l=1:J,
  if clusterflag==0,
    if l==J | eigT(l+1)-eigT(l) > cluster_width,        % Not part of a cluster
      JJ = JJ+1; YY(:,JJ) = Y(:,l);
    else
      l1 = l; clusterflag = 1;                          % Start of a cluster
    end;
  else
    if l==J | eigT(l+1)-eigT(l) > cluster_width,        % End of a cluster
      JJ = JJ+1; 
      YY(:,JJ) = Y(:,l1:l)*S(J,l1:l)'/sqrt(S(J,l1:l)*S(J,l1:l)');
      clusterflag = 0;
    end;
  end;
end;

%    Now order Ritz vectors according to the size of the projection of Q(:,J+1)
%    in the direction of each.

proj = YY'*Q(:,J+1);
for l=1:JJ,
  proj(l) = proj(l)/sqrt(YY(:,l)'*YY(:,l));
end;
[proj_sort,indx] = sort(abs(proj));
YY = YY(:,indx);

%    Eliminate any with too large a projection onto Q(:,J+1).

m = JJ;
while m >= 1 & proj_sort(m) > sqrt(eps)*normA/T(J+1,J), m=m-1; end;

%    Form an orthonormal basis for the remaining Ritz vectors.
%    If projection of Q(:,J+1) onto this space is small enough, 
%    we are done.  Otherwise, eliminate vectors one at a time
%    until it is.

[W,Sigma,X] = svd(YY(:,1:m));
W = W(:,1:m); Sigma = Sigma(1:m,1:m);
sigma = diag(Sigma);
%dsigma_min = min(double(sigma)), dsigma_max = max(double(sigma)),
WTqJp1 = W'*Q(:,J+1);
err1 = T(J+1,J)*sqrt(WTqJp1'*WTqJp1);
ImWWTqJ = Q(:,J) - W*W'*Q(:,J);
err2 = T(J+1,J)*sqrt(ImWWTqJ'*ImWWTqJ);
%derr1 = double(err1), derr2 = double(err2),
%pause,
while m >= 1 & ( err1 > sqrt(eps)*normA | (err1 > 100*err2 & err1 > sqrt(eps)*normA/100) ),
  ll = m;
  [Wll,Sigmall,Xll] = svd(YY(:,1:m-1));
  Wll = Wll(:,1:m-1); Sigmall = Sigmall(1:m-1,1:m-1);
  sigmall = diag(Sigmall);
  dsigmall_min = min(double(sigmall)); dsigmall_max = max(double(sigmall));
  WllTqJp1 = Wll'*Q(:,J+1);
  err1ll = T(J+1,J)*sqrt(WllTqJp1'*WllTqJp1);
  ImWllWllTqJ = Q(:,J) - Wll*Wll'*Q(:,J);
  err2ll = T(J+1,J)*sqrt(ImWllWllTqJ'*ImWllWllTqJ);
  for lll=m-1:-1:1,
    if lll==1,
      YYY = YY(:,2:m);
    else
      YYY = YY(:,[1:lll-1,lll+1:m]);
    end;
    [Wlll,Sigmalll,Xlll] = svd(YYY);
    Wlll = Wlll(:,1:m-1); Sigmall = Sigmalll(1:m-1,1:m-1);
    sigmalll = diag(Sigmalll);
    dsigmalll_min = min(double(sigmalll)); dsigmalll_max = max(double(sigmalll));
    WlllTqJp1 = Wlll'*Q(:,J+1);
    err1lll = T(J+1,J)*sqrt(WlllTqJp1'*WlllTqJp1);
    ImWlllWlllTqJ = Q(:,J) - Wlll*Wlll'*Q(:,J);
    err2lll = T(J+1,J)*sqrt(ImWlllWlllTqJ'*ImWlllWlllTqJ);
    if err1lll^2 + err2lll < err1ll^2 + err2ll,
      ll = lll; err1ll = err1lll; err2ll = err2lll;
      dsigmall_min = dsigmalll_min; dsigmall_max = dsigmalll_max;
    end;
  end;
%  ll, dsigmall_min, dsigmall_max,
%  derr1ll = double(err1ll), derr2ll = double(err2ll),
%  pause,
%    Eliminate column ll from YY and decrement m by 1.
  if ll==1, 
    YY = YY(:,2:m);
  elseif ll==m,
    YY = YY(:,1:m-1);
  else
    YY = YY(:,[1:ll-1,ll+1:m]);
  end;
  err1 = err1ll; err2 = err2ll;
  m = m-1;
end;


%    Now add Ritz vectors, one at a time, until the two error terms match.

%m = 0; err1 = vpa(0); err2 = T(J+1,J);
%derr1 = double(err1), derr2 = double(err2), pause,
%while err2 > err1 & m < JJ,
%  err1prev = err1; err2prev = err2;
%  m = m+1;
%  [W,Sigma,X] = svd(YY(:,1:m));
%  W = W(:,1:m); Sigma = Sigma(1:m,1:m);
%  sigma = diag(Sigma); 
%  dsigma_min = min(double(sigma)), dsigma_max = max(double(sigma)),
%  WTQJp1 = W'*Q(:,J+1);
%  err1 = T(J+1,J)*sqrt(WTQJp1'*WTQJp1);
%  err1 = vpa(0.1)*err1;
%  ImWWTQJ = Q(:,J) - W*W'*Q(:,J);
%  err2 = T(J+1,J)*sqrt(ImWWTQJ'*ImWWTQJ);
%  derr1 = double(err1), derr2 = double(err2), pause,
%end;

%if m > 0,
%  errprev = err2prev;  if err1prev > errprev, errprev=err1prev; end;
%  err = err2; if err1 > err, err = err1; end;
%  if errprev < err, m=m-1; end;
%end;
fprintf('m = %3i\n',m)
if m > 0,
  [W,Sigma,X] = svd(YY(:,1:m));
  W = W(:,1:m);
end;
 
%  Continue 3-term recurrence with perturbation terms designed to make new vectors
%  orthogonal to each other and to the columns of W.
  
N = n+J-m;
F = vpa(zeros(n,N));
F(:,1:J) = A*Q(:,1:J) - Q*T;   % Rounding errors from fp computation.

%    Redo bins for plotting.

widthsave = width; egroupsave = egroup; ngroupssave = ngroups; edgessave = edges;

width = vpa(width);
egroup(1,1) = eigA(1);  egroup(1,2) = eigA(1); ngroups = 1;  % Group eigenvalues that are too close.
for i=2:n,
  if eigA(i)-width <= egroup(ngroups,2)+width,
    egroup(ngroups,2) = eigA(i);
  else
    ngroups = ngroups+1;
    egroup(ngroups,1) = eigA(i); egroup(ngroups,2) = eigA(i);
  end;
end;

edges = vpa(zeros(2*ngroups+2,1));
edges(1) = -Inf; edges(2*ngroups+2) = Inf;
for i=1:ngroups,
  edges(2*i) = egroup(i,1) - width;
  edges(2*i+1) = egroup(i,2) + width;
end;

%    Orthogonalize Q(:,J+1) against columns of W and renormalize.

vJp1 = A*Q(:,J) - T(J,J)*Q(:,J) - T(J-1,J)*Q(:,J-1);
if m > 0, vJp1 = vJp1 - W*(W'*vJp1); end;
betaJp1 = sqrt(vJp1'*vJp1); 
if m==0, 
  F(:,J) = vpa(zeros(n,1));
else
  F(:,J) = W*W'*(T(J+1,J)*Q(:,J+1)+F(:,J));
end;
Q(:,J+1) = vJp1/betaJp1;
T(J,J+1) = betaJp1; T(J+1,J) = betaJp1;

betaj = T(J+1,J);
for j=J+2:N+1,
  Q(:,j) = A*Q(:,j-1);
  Q(:,j) = Q(:,j) - betaj*Q(:,j-2);
  alphaj = Q(:,j)'*Q(:,j-1);
  T(j-1,j-1) = alphaj;
  Q(:,j) = Q(:,j) - alphaj*Q(:,j-1);
  if j==J+2,
    for l=1:m,
      F(:,j-1) = F(:,j-1) + (Q(:,j-1)'*A*W(:,l) - T(J+1,J)*Q(:,j-2)'*W(:,l))*W(:,l);
    end;
  else
    for l=1:m,
     F(:,j-1) = F(:,j-1) + (Q(:,j-1)'*A*W(:,l))*W(:,l);
    end;
    F(:,j-1) = F(:,j-1) + Q(:,J+1)*(T(J+1,J)*Q(:,j-1)'*Q(:,J));
  end;
  Q(:,j) = Q(:,j) - F(:,j-1);

%    At this point, the new Lanczos vector should be perfectly orthogonal to W
%    and to the other new Lanczos vectors.  If it is not because we used only
%    64 decimal place arithmetic, orthogonalize one more time to be sure.

  if m > 0,
    newpert = W*W'*Q(:,j);
  else
    newpert = zeros(n,1);
  end;
%  checkWorth = double(sqrt(newpert'*newpert)),
  Q(:,j) = Q(:,j) - newpert; F(:,j-1) = F(:,j-1) + newpert;
  newpert = Q(:,J+1:j-1)*Q(:,J+1:j-1)'*Q(:,j);
%  checkQorth = double(sqrt(newpert'*newpert)), pause,
  Q(:,j) = Q(:,j) - newpert; F(:,j-1) = F(:,j-1) + newpert;

  betaj = sqrt(Q(:,j)'*Q(:,j));
  T(j,j-1) = betaj;
  if j < N+1, T(j-1,j) = betaj; end;
%  if j==N+1, betaNp1 = double(betaj), end;
  if j < N+1, Q(:,j) = Q(:,j)/betaj; end;

%    Plot eigenvalues of T in bins about eigenvalues of A.

  eigT = sort(eig(T(1:j-1,1:j-1)));
  subplot(3,1,3)
%  NTnew = histc(eigT,edges);
  NTnew = zeros(length(edges),1);
  k = 1;
  for i=1:length(eigT),
    while k < length(edges) & edges(k) <= eigT(i),
      k = k+1;
    end;
    if edges(k) > eigT(i) & k > 1, k = k-1; end;
    if edges(k) <= eigT(i) & k < length(edges) & eigT(i) < edges(k+1),
      NTnew(k) = NTnew(k)+1;
    end;
  end;
  bar(NTnew),
  title(['Step ',int2str(j-1)])
  pause(1),
end;
title('Eigenvalues of Matrix in Equivalent Exact Arithmetic Computation')

%  Check.

Qorth = [W Q(:,J+1:N)];
errorth = norm(double(eye(n,n) - Qorth'*Qorth),'fro');
fprintf('Norm of I - [W,Q]trans*[W,Q] = %12.4d\n',errorth)
AQmQT = A*Q(:,1:N) - Q(:,1:N)*T(1:N,1:N);
%errF = norm(double(F-AQmQT),'fro'),
normF = norm(double(F),'fro');
fprintf('Norm of perturbations = %12.4d\n',normF)

%  Check distance from eigenvalues of T to those of A.
eigA = sort(eig(A)); eigT = sort(eig(T(1:N,1:N)));
maxdist = 0; imax = 0; iimax = 0;
for i=1:length(eigT),
  disti = abs(eigT(i)-eigA(1)); ii=1;
  for iii=2:n,
    distiii = abs(eigT(i)-eigA(iii));
    if distiii < disti, disti=distiii; ii=iii; end;
  end;
  if disti > maxdist, maxdist = disti; imax=i; iimax=ii; end;
end;
fprintf('Maxdist from eval of T to nearest eval of A = %12.4d\n\n',double(maxdist))

%  This routine loads the bcsstk03 matrix and finds its eigenvalues.
%  It then uses the Remez algorithm to compute a sharp upper bound for 
%  the A-norm of the error in the CG algorithm and the initial vector 
%  b for which this bound will be attained.

%  Load A and find its eigenvalues.

load bcsstk03.mat;
A = Problem.A;
Afull = full(A);
eigA = eig(A);
n = length(eigA);
figure(1)
semilogx(real(eigA),imag(eigA),'xk'), shg, hold on

nhat = n; eigAhat = eigA;

%  Ask user for step numbers k.

kmin = input('Enter kmin: ');
kmax = input('Enter kmax: ');

%  Loop over steps k.

errAk = zeros(kmax-kmin+1,1);
for k=kmin:kmax,
  
%    Use eigenvalues closest to extreme points of Chebyshev polynomial as initial guess.

  if k==1, indices = zeros(k+1,1); indices(1) = 1; indices(k+1) = nhat; end;
  if k > 1,
    indices(1) = 1; indices(k+1) = nhat;
    for kk=k:-1:3,
      indices(kk) = indices(kk-1);
    end;
    kk = 3;
    while kk < k+1 & indices(kk)==kk-1,
      kk = kk+1;
    end;
    indices(2) = kk-1;

%    Keep indices from previous k and add rightmost point that is not already there.
%    kk = k-1;
%    while indices(kk)==nhat-(k-kk),
%      indices(kk+1) = indices(kk);
%      kk = kk-1;
%    end;
%    indices(kk+1) = nhat-(k-kk);
%    chebypts = cos(pi*[k-1:-1:1]'/k);
%    shiftedcheby = (eigAhat(nhat)-eigAhat(1))/2 * chebypts + (eigAhat(nhat)+eigAhat(1))/2;
%    llstart = 2;
%    for l=2:k,
%%        Find eigenvalue of Ahat that is closest to shiftedcheby(l-1).
%      if eigAhat(llstart) >= shiftedcheby(l-1),
%        indices(l) = llstart; llstart = llstart+1;
%      else
%        while eigAhat(llstart) < shiftedcheby(l-1),
%          llstart = llstart+1;
%        end;
%        indices(l) = llstart;
%        if abs(eigAhat(llstart-1)-shiftedcheby(l-1)) < abs(eigAhat(llstart)-shiftedcheby(l-1)),
%          llstart = llstart-1; indices(l) = llstart;
%        end;
%        llstart = llstart+1;
%      end;
%    end;
  end;
%  indices,

%    Evaluate Bk at each eigenvalue of Ahat and find where its absolute value is maximal.
%    If it is at one of the eigenvalues in indices, we are done.  Otherwise, make an
%    exchange and continue.

  flag = 0;
  while flag==0,

  Bkmax = 0; indxmax = 0; signBkx = zeros(nhat,1); Bkxvals = zeros(nhat,1);
  for m=1:nhat,
    x = eigAhat(m);
    Bkx = 0;
    for j=1:k+1,
      termj = (-1)^(j-1);
      for i=1:k+1,
        if i~=j, 
          termj = termj*(eigAhat(indices(i))-x)/(eigAhat(indices(i))-eigAhat(indices(j)));
        end;
      end;
      Bkx = Bkx + termj;
    end;
    signBkx(m) = sign(Bkx);
    Bkxvals(m) = Bkx;
    Bkx = abs(Bkx);
    if Bkx > Bkmax,
      Bkmax = Bkx; indxmax = m;
    end;
  end;
%  figure(2)
%  plot(eigAhat,Bkxvals), shg
%  Bkmax, indxmax, pause,

  l = 1;
  while indices(l) < indxmax & l<k+1,
    l = l+1;
  end;
  if indices(l)==indxmax,   % We have found the right indices.  Evaluate Bkden.
    Bkden = 0;
    for j=1:k+1,
      termj = (-1)^(j-1);
      for i=1:k+1,
        if i~=j, 
          termj = termj*eigAhat(indices(i))/(eigAhat(indices(i))-eigAhat(indices(j)));
        end;
      end;
      Bkden = Bkden + termj;
    end;
    figure(1)
    semilogx(eigAhat(indices),k*ones(k+1,1),'o'), shg, pause(.1)
    flag = 1; errk = 1/Bkden, errAk(k-kmin+1) = errk;
  else
    if l==1, 
      indices(l) = indxmax;
    else
      if signBkx(indices(l-1))==signBkx(indxmax),
        indices(l-1) = indxmax;
      else
        indices(l) = indxmax;
      end; 
    end;
%    indices,
  end;
  end;

end;

%  Form right-hand side vector b.

etilde0sq = zeros(nhat,1);
avg = prod(eigAhat(indices).^(1/(kmax+1)));
for l=1:kmax+1,
  etilde0sq(indices(l)) = 1/eigAhat(indices(l));
  for j=1:kmax+1,
    if j~=l,
      for i=j+1:kmax+1,
        if i~=l,
          etilde0sq(indices(l)) = etilde0sq(indices(l))*...
                                  (eigAhat(indices(i))-eigAhat(indices(j)))/avg;
        end;
      end;
      etilde0sq(indices(l)) = etilde0sq(indices(l));
    end;
  end;
end;
%  Check formula in 1979 paper!
%wsq = zeros(nhat,1);
%for l=1:kmax+1,
%  for j=1:kmax,
%    if j~=l,
%      termj = eigAhat(indices(j));
%      for i=j+1:kmax+1,
%        if i~=l,
%          termj = termj*(eigAhat(indices(i))-eigAhat(indices(j)));
%        end;
%      end;
%      wsq(indices(l)) = wsq(indices(l)) + termj;
%    end;
%  end;
%end;
%[etilde0sq(indices), wsq(indices), wsq(indices)./etilde0sq(indices)]
etilde0 = sqrt(etilde0sq);
[V,Lam] = eig(Afull);
b = V*Lam^(1/2)*etilde0;
b = b/norm(Afull^(-1/2)*b);

figure(2)
semilogy(errAk)


function [x,y] = fastLP(A, Aeq, c, b, beq, K, bsense)
% Set parameters
%
  [m, n] = size(A);
 [me, ~] = size(Aeq);
 
 BigA = [A; Aeq];
 Bigb = [b; beq];
 
 m = m + me;
 
%  br = K * Bigb; % Set the initial remaining inventory
 d = Bigb / n;
%  br = K * Bigb;
 
step = 1 / sqrt(K * n); % Set the step-size
%
 x = zeros(n,1); % Set initial primal solutions
% y = ones(m,1)./(sum(A'))';
% y = (sum(c)/n)*y;
if bsense == 0
  y = ones(m, 1) / exp(1.5); %* mean(abs(c)./sum(abs(A))'); % Set initial dual solutions
else  
    y = zeros(m, 1);
end
%
% Start the outer loop
 for k=1:K
  p=randperm(n); %Randomly permute variable-order in each round

% Start the inner loop
  for i=1:n
   ii=p(i);
   aa=BigA(:,ii);
   % step = 1 / sqrt(((k - 1) * n + i));
   %
   % Set the primal increment
    if c(ii)-aa'*y >= 0
        xk = 1;
    else
        xk = 0;
    end % End if 
   % xk = (sign(c(ii)-aa'*y)+1)/2;
   %
   % Update the dual soluton
   y = y.*exp(step * (xk * aa - d));
   % y = y + step * (xk * aa - d);
   % y(1:m - me) = max(0, y(1:m - me));
   %y=max(0,y+step*(xk*aa-d)); 
   %
   % Update the remaining inventory and primal solution

   x(ii) = x(ii) + xk;
   %
%    br = br - xk * aa;
%    if (i * k) < (n * K 
%        d = br / (n * K - i * k);
%    end % End if
  end
 end
  x=x/K; % Output the primal solution
%
%  This program solves the linear program
%
%      max   c'*x
%      s.t.  Ax <= b, 0 <= x <= 1.
%
%  Input 
%      A: (sparse) inequality constraint matrix.
%      b (>0): inequality right-hand column vector
%      K: boosting (positive) integer
%  
%  Output
%     x      : approximate (fractional) solution betwee 0 and 1,
%     y>=0   : price vector





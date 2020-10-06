function [x,y] = fastLP(A, Aeq, c, b, beq, K, bsense)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m, n] = size(A);
[me, ~] = size(Aeq);

BigA = [A; Aeq];
Bigb = [b; beq];

m = m + me;

%  br = K * Bigb; % Set the initial remaining inventory
d = Bigb / n;
%  br = K * Bigb;

step = 1/ sqrt(K*n); % Set the step-size

x = zeros(n,1);
if bsense == 0
    y = ones(m, 1); % Set initial dual solutions
else
    y = zeros(m, 1);
end

% Start the outer loop
for k = 1:K
    p = randperm(n);
    
    % Start the inner loop
    for i = 1:n
        ii = p(i);
        aa = BigA(:,ii);
        % Set the primal increment
        if c(ii)-aa'*y >= 0
            xk = 1;
        else
            xk = 0;
        end % End if
        
        % Update the dual soluton
        y = y + step * (xk * aa - d);
        y(1:m - me) = max(0, y(1:m - me));
        
        % Update the remaining inventory and primal solution
        x(ii) = x(ii) + xk;
        % br = br - xk * aa;
        % if (i * k) < (n * K
        %     d = br / (n * K - i * k);
        % end % End if
    end
end
x = x / K; % Output the primal solution

end % End function







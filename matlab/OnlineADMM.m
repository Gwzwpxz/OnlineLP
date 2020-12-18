function [x, y] = OnlineADMM(A, b, c, K, CheckInnerFeas)
% This function implements ADMM based online algorithm

% Get size of data array
[m, n] = size(A);

% Initialize auxiliary resources
d = b / n;

% Initialize primal and dual variables for ADMM
p = zeros(m, 1);
q = zeros(m, 1);
lbd = zeros(m, 1);

% Initialize primal variable for ADMM
x = zeros(n, 1);

% Initialize proximal parameter beta
beta = sqrt(K * n);

% Set remaining inventory
br = K * b;

% Start outer iteration
for k = 1:K
    % counter = 0;
    for i = randperm(n)
        % counter = counter + 1;
        % beta = sqrt((k - 1) * n  + counter);
        aa = A(:, i);
        
        % Estimate x based on two dual variables p and q
        % For SCP
        %         if (c(i) <= p' * aa) && (c(i) <= q' * aa)
        %             x(i) = x(i) + 1;
        %         end % End if
        %
        %         % Implement stochastic ADMM
        %         p = max(q - (-d + lbd) / beta, 0);
        %         q = p + lbd / beta;
        %
        %         if (c(i) <= aa'*p - (norm(aa)^2 - aa' * lbd) / beta)
        %             q = q - aa / beta;
        %         end % End if
        %
        %         lbd = lbd + beta * (p - q);
        
        % For Knapsack
        xk = ((c(i) >= p' * aa) && (c(i) >= q' * aa));
        
        % Implement stochastic ADMM
        p = max(q - (d + lbd) / beta, 0);
        
        % Update q by choosing candidate
        q1 = p + lbd / beta;
        q2 = q1 + aa / beta;
        q3 = q1 + (((- aa' * lbd / beta) + c(i) - aa' * p) * aa) / norm(aa)^2;
        
        % Compute minimum
        [~, minidx] = min([ProxObj(aa, c(i), p, beta, -lbd, q1), ...
            ProxObj(aa, c(i), p, beta, -lbd, q2), ...
            ProxObj(aa, c(i), p, beta, -lbd, q3)]);
        
        % Apply admm update for variable q
        if minidx == 1
            q = q1;
        elseif minidx == 2
            q = q2;
        else
            q = q3;
        end % End if
        
        lbd = lbd + beta * (p - q);
        
        % Check feasibility of inner iteration and update
        if (~ CheckInnerFeas) || (min(br - xk * aa) >= 1e-08)
            br = br - xk * aa;
            x(i) = x(i) + xk;
        end % End if
        
    end % End for
end % End for

x = x / K;
y = p;

end % End function
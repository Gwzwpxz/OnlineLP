function [x, y] = OnlineProx(A, b, c, K, CheckInnerFeas, Momentum)
% This function implements proximal point based online algorithm
% Momentum Parameter here is to control path averaging

% Get size of data array
[m, n] = size(A);

% Initialize auxiliary resources
d = b / n;

% Initialize primal and dual variables for proximal point method
p = zeros(m, 1);
p_pre = p;
paverage = zeros(m, 1);

% Initialize primal variable for ADMM
x = zeros(n, 1);

% Get stepsize
gamma = 1;% sqrt(K * n);

% Set remaining inventory
br = K * b;

for k = 1:K
    counter = 0;
    for i = randperm(n)
        counter = counter + 1;
        % gamma = sqrt((k - 1) * n  + counter);
        
        aa = A(:, i);
        xk = (c(i) >= aa' * p);
        
        % Compute proximal candidates
        p1 = p - d / gamma;
        p2 = p1 + aa / gamma;
        p3 = p1 + (((aa' * d / gamma) + c(i) - aa' * p) * aa) / norm(aa)^2;
        
        % Projection onto non-negative orthant
        p1 = max(p1, 0);
        p2 = max(p2, 0);
        p3 = max(p3, 0);
        
        p_pre = p;
        % Compute minimum
        [~, minidx] = min([ProxObj(aa, c(i), p, gamma, d, p1), ...
            ProxObj(aa, c(i), p, gamma, d, p2), ...
            ProxObj(aa, c(i), p, gamma, d, p3)]);
        
        % Apply proximal point update
        if minidx == 1
            p = p1;
        elseif minidx == 2
            p = p2;
        else
            p = p3;
        end % End if
        
        if Momentum < 1
            p = max(p + Momentum * (p - p_pre), 0);
        else
            paverage = 0.1 * paverage + 0.9 * p;
            % p = paverage / ((k - 1) * n + counter);
        end % End if
        
        % p = 0.999 * p + 0.001 * paverage / ((k - 1) * n + counter);
        
        % Check feasibility of inner iteration and update
        if (~ CheckInnerFeas) || (min(br - xk * aa) >= 1e-08)
            br = br - xk * aa;
            x(i) = x(i) + xk;
        end % End if
        
    end % End for
end % End for

x = x / K;
if Momentum < 1
    y = p;
else
    y = paverage;
end % End if 

end % End function

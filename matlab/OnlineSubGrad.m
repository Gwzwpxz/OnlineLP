function [x, y] = OnlineSubGrad(A, b, c, K, CheckInnerFeas, Metric, Momentum)
% This function implements sub-gradient based online algorithm

% Parameter Adaptive is not in use for now
Adaptive = 0;

% Get size of data array
[m, n] = size(A);

% Initialize auxiliary resources
d = b / n;

% Set stepsize
step = 1 / sqrt(K * n);

% Initilize primal and dual solutions
x = zeros(n, 1);

% Initialize the dual variable based one metric used
if Metric == "L2"
    y = zeros(m, 1);
else
    % Initialize parameters for simplex projection
    y = ones(m, 1) / exp(1);
    delta = mean(c) / min(d);
end % End if

% Initialize momentum parameter
y_pre = y;

% Set remaining inventory
br = K * b;

% Start the outer loop
for k = 1:K
    
    % Random permutation
    p = randperm(n);
    
    % Start the inner loop
    for i = 1 : n
        ii = p(i);
        aa = A(:,ii);
        
        % Set primal increment
        infeas = c(ii) - aa' * y;
        xk = (infeas >= 0);
        infeas = max(infeas, 0);
        
        % step = 1 / sqrt((k - 1) * n + i);
        
        % Update the dual soluton
        
        % Use alternative descent stepsizes
        subgrad = (d - xk * aa);
        % step = (d' * y + infeas) / (norm(subgrad)^2 + 1e-08);
%         if ~ (step == 2 / sqrt(n * K))
%             step
%         end % End if
        
        if Metric == "L2"
            y = max(0, y + Momentum * (y - y_pre));
            y_pre = y;
            y = max(0, y - step * subgrad);
        else
            y = max(0, y + Momentum * (y - y_pre));
            y_pre = y;
            y = y.*exp(- step * subgrad);
            
            % Do projection back to the simplex
            if (sum(y) > delta)
                % Projection
                y = y / sum(y) * delta;
            end % End if
                
        end % End if
        
        % Check feasibility of inner iteration and update
        if (~ CheckInnerFeas) || (min(br - xk * aa) >= 0)
            br = br - xk * aa;
            x(ii) = x(ii) + xk;
        end % End if
        
        if Adaptive && (i * k < n * K)
            d = br / (n * K - i * k);
        end % End if
        
    end % End for
end % End for

x = x / K;
end % End function

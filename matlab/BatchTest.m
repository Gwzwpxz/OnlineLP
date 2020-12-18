clear;
clc;

% L2 Approx
ntest = 10;
batchrange = [1, 2, 4, 6, 8, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 768, 1024, 3072, 4096];
ninterval = length(batchrange);
L2LPapprox = zeros(ninterval, 1);
L2LPRegret = zeros(ninterval, 1);

m = 64;
n = 4096 * 3;

% A = randn(m, n) + 1;
% b = n^(1/3) * (rand(m, 1) + 1) / 3;
% c = sum(A, 1)'/m - rand(n, 1) * m;
% A = randn(m, n) + 1;
% b = n * (rand(m, 1) + 1) / 3;
% c = (sum(A, 1)' - rand(n, 1) * m)/m;

A = randi(1000, m, n) / 100;
b = 0.25 * sum(A, 2);
c = sum(A, 1)'/m + rand(n, 1) * 5;

data.A = A;
data.b = b;
data.c = c;

% Do comparison
model.A = sparse(A);
model.rhs = b;
model.sense = '<';
model.modelsense = 'max';
model.ub = ones(n, 1);
model.obj = c;
sol = gurobi(model);
fval = sol.objval;

params.CheckInnerFeas = true;
params.BoostingParam = 10;
params.Xmax = 1;
params.Momentum = 0;
params.SubAlg = "Batch";
params.Metric = "L2";

for i = 1:ninterval
    
    % Generate data
    batchsize = batchrange(i);
    
    params.Batch = batchsize;
    LPObjs = zeros(ntest, 1);
    Regrets = zeros(ntest, 1);
    
    for j = 1:ntest
        
        % Solve problem
        [optx, info] = OLLP(data, params);
        
        LPSolution = optx.Lpx;
        
        % Verify feasibility
        if ~ (sum(A * LPSolution <= b) == m)
            error("LP solution is not feasible");
        end % End if
        
        % Relaxed LP objective
        LPobjectiveValue = info.LPobj;
        
        LPObjs(j) = (LPobjectiveValue / fval * 100);
        Regrets(j) = fval - LPobjectiveValue;
        
    end % End for
    
    % The approximation ratio is
    L2LPapprox(i) = mean(LPObjs);
    L2LPRegret(i) = mean(Regrets);
    
end % End for
clear;
clc;

% KL Approx

ninterval = 10;
ntest = 10;
mrange = floor(logspace(1, 3, ninterval))';
KLLPapprox = zeros(ninterval, 1);
KLLPRegret = zeros(ninterval, 1);
L2LPapprox = zeros(ninterval, 1);
L2LPRegret = zeros(ninterval, 1);

for i = 1:ninterval
    
    % Generate data
    m = mrange(i);
    n = 10000;
    
    % A = randn(m, n) + 1;
    % b = n^(1/3) * (rand(m, 1) + 1) / 3;
    % c = sum(A, 1)'/m - rand(n, 1) * m;
    A = randn(m, n) + 1;
    b = n * (rand(m, 1) + 1) / 3;
    c = (sum(A, 1)' - rand(n, 1) * m)/m;
    
    data.A = A;
    data.b = b;
    data.c = c;
    
    LPObjs = zeros(ntest, 1);
    Regrets = zeros(ntest, 1);
    
    params.CheckInnerFeas = true;
    params.BoostingParam = 1;
    params.Xmax = 1;
    params.Momentum = 0;
    params.Metric = "KL";
    params.SubAlg = "SubGrad";
    
    % Do comparison
    model.A = sparse(A);
    model.rhs = b;
    model.sense = '<';
    model.modelsense = 'max';
    model.ub = ones(n, 1);
    model.obj = c;
    
    sol = gurobi(model);
    fval = sol.objval;
    
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
    KLLPapprox(i) = mean(LPObjs);
    KLLPRegret(i) = mean(Regrets);
    
    params.Metric = "L2";
    
    for j = 1:ntest
        
        % Solve problem
        [optx, info] = OLLP(data, params);
        
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
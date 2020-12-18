% This is a test script to compare performance of different parameter
% settings

% Get parameter struct
params = SetDefaultParam();

% Adjust parameters
params.CheckInnerFeas = true;
params.BoostingParam = 100;
params.Xmax = 1;
params.Momentum = 0;
params.Metric = "L2";
params.SubAlg = "SubGrad";

% Generate data
m = 5;
n = 2000;

A = randi(1000, m, n) / 100;
b = sum(A, 2) * 0.1 / sqrt(n);
c = sum(A, 1)'/m + rand(n, 1) * 500;

data.A = A;
data.b = b;
data.c = c;

% Solve problem
[x, fval, ~, ~, y] = linprog(-c, A, b, [], [], zeros(n, 1), ones(n, 1) * params.Xmax);
y = y.ineqlin;
fval = - fval;

% Subgradient with L2 metric
params.SubAlg = "SubGrad";
params.Metric = "L2";
params.Momentum = 0.1;
[optxSubGradL2, infoSubGradL2] = OLLP(data, params);

% Subgradient with KL divergence
params.SubAlg = "SubGrad";
params.Metric = "KL";
[optxSubGradKL, infoSubGradKL] = OLLP(data, params);

% Proximal point with no Momentum
params.SubAlg = "Proximal";
params.Momentum = 0;
[optxProx, infoProx] = OLLP(data, params);

% Proximal point with Momentum
params.SubAlg = "Proximal";
params.Momentum = 0.1;
[optxProxMomentum, infoProxMomentum] = OLLP(data, params);

% Proximal point with path averaging
params.SubAlg = "Proximal";
params.Momentum = 2;
[optxProxAvg, infoProxAvg] = OLLP(data, params);

% ADMM
params.SubAlg = "ADMM";
[optxADMM, infoADMM] = OLLP(data, params);

disp("Subgradient L2: Approx " + infoSubGradL2.LPobj * 100 / fval + "% Norm: " + norm(y - infoSubGradL2.PriceVec));
disp("Subgradient KL: Approx " + infoSubGradKL.LPobj * 100 / fval + "% Norm: " + norm(y - infoSubGradKL.PriceVec));
disp("Proximal: Approx " + infoProx.LPobj * 100 / fval + "% Norm: " + norm(y - infoProx.PriceVec));
disp("Proximal Momentum: Approx " + infoProxMomentum.LPobj * 100 / fval + "% Norm: " + norm(y - infoProxMomentum.PriceVec));
disp("Proximal Avg: Approx " + infoProxAvg.LPobj * 100 / fval + "% Norm: " + norm(y - infoProxAvg.PriceVec));
disp("ADMM: Approx " + infoADMM.LPobj * 100 / fval + "% Norm: " + norm(y - infoADMM.PriceVec));

normlist = [infoSubGradL2.PriceVec infoSubGradKL.PriceVec, infoProx.PriceVec, infoProxMomentum.PriceVec, infoProxAvg.PriceVec, infoADMM.PriceVec];
normvallist = sqrt(sum((y - normlist).^2, 1));

[~, idx] = min(normvallist);
bestmethod(i) = idx;
PriceVector = normlist(:, idx);



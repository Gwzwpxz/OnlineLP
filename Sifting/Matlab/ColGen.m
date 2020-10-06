%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Sifting Test                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;

file = "scpm1";
xlp = load ("..\presolved_sol\" + file + "_presolved_sol.mat");
xlp = xlp.xlp';

data = load(file + ".mat");
A = data.A;
b = data.b';
Aeq = data.Aeq;
beq = data.beq';
c = - double(data.c');
K = 1;
[m, n] = size(A);
Xmax = 1; 

[x, y] = fastLP(A, Aeq, c, b, beq, K, 0);

% Number of right basic variables found
% sum(x(xlp > 0) > 0) / sum(xlp > 0)

thresh = 1;

% This quantity is number of correct elements selected by online algorithm
nCorrect = sum(xlp(x >= thresh) > 0);

% This quantity is number of all elements selected by online algorithm
nAll = sum(x >= thresh);

% This is percentage of correct elements selected by online algorithm
% This should be large

nCorrectPercentage = nCorrect / sum(xlp > 0);

% This is efficiency for finding correct elements
% This should be large

nCorrecctEfficiency = nCorrect / nAll;


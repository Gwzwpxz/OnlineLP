function [optx, info] = OLLP(data, params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        OnLine Linear Progamming Solution Routine Matlab Interface       %
%                                                                         %
%                 Version 0.2     Nov 1st, 2020                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function serves as a driver of optimization algorithm in OLLP      %
% algorithm with Extended Simple Algorithm from reference                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Li, Xiaocheng , C. Sun , and Y. Ye . "Simple and Fast Algorithm for  %
%    Binary Integer and Online Linear Programming." (2020).               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given input data, OLLP approximates an online binary integer programming
% problem (and alternatively a general integer program) by ESA algorithm
% and returns approximate solutions before and after being rounded
%
% The problem is expected to be represented by
%
%       maximize     c^T * x
%       subject to   Ax <= b
%                     x in {0, 1} (or x in [M])
%
% The input format is specified as follows.
%
% Input
%      data:
%      A struct containing two components X and Y to be tested
%           A             : Inequality Constraint Matrix of size m * n
%           b             : Inequality RHS vector, expected non-negative
%           c               : Linear coefficient vector of size n * 1
%      params:
%      A struct containing parameters for the optimization algorithm
%           CheckInnerFeas: Boolean variable determining whether to check
%                           constraint feasibility in inner iterations
%                           Turned off by default
%           BoostingParam : Boosting parameter in ESA
%                           Set 50 by default (must be >= 1)
%           Xmax          : An integer parameter implying an upperbound
%                           on x in integer programming formulation
%                           Set 1 by default (must be >= 1)
% Usage
%       [Ipx, Lpx, info] = OLLPSolve(data, params);
%
% Output
%       Ipx                 : Approximate solution to IP problem
%       Lpx                 : Approximate solution to LP problem
%
%       info:
%       A struct containing information about optimization
%           OLLPtime        : Total running time of OLLP algorithm
%           LPtime          : Time spent solving relaxed LP
%           Rtime           : Time spent rounding the solution
%           Probintro       : Description of problem data
%           PriceVec        : Price vector of ESA algorithm
%           LPobj           : Objective of relaxed LP
%           IPobj           : Objective of relaxed IP
%
% Function References:
% 1. CheckInput 2. SetDefaultParam 3. OLLPSolve 4. SolRound

% Match the number of arguments
Probintro = CheckInput(data);

if nargin == 2
elseif nargin == 1 && Probintro.isValid
    params = SetDefaultParam();
else
    error("Invalid number of arguments");
end % End if

% Unpack Parameters
CheckInnerFeas = params.CheckInnerFeas;
BoostingParam = params.BoostingParam;
Xmax = params.Xmax;
SubAlg = params.SubAlg;
Metric = params.Metric;
Momentum = params.Momentum;
Batch = params.Batch;

% Unpack data
A = data.A;
b = data.b;
c = data.c;

% Start Algorithm
tic;
[Lpx, PriceVec] = OLLPSolve(data, CheckInnerFeas, BoostingParam, ...
    Xmax, SubAlg, Metric, Momentum, Batch);
LPtime = toc;

% Recovery of integer solution
tic;
Ipx = SolRound(repmat(data.A, 1, Xmax), b, Lpx);

n = length(c);
for i = 1:(Xmax - 1)
    Lpx(1:n) = Lpx(1:n) + Lpx(n * i + 1:n * i + n);  
    Ipx(1:n) = Ipx(1:n) + Ipx(n * i + 1:n * i + n);  
end % End for

Lpx = Lpx(1:n);
Ipx = Ipx(1:n);

Rtime = toc;
OLLPtime = LPtime + Rtime;

% Summarize all information

optx.Ipx = Ipx;
optx.Lpx = Lpx;

info.OLLPtime = OLLPtime;
info.LPtime = LPtime;
info.Rtime = Rtime;
info.Probintro = Probintro;
info.PriceVec = PriceVec;
info.LPobj = c' * Lpx;
info.IPobj = c' * Ipx;

% Begin routine
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
disp("%%%%%%%%%%%%%%%% OLLP Solution Routine %%%%%%%%%%%%%%%%");
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
disp("");

% Print solution details
disp("Parameter Console");

if CheckInnerFeas
    disp("Check inner feasibility: Yes");
else
    disp("Check inner feasibility: No");
end % End if
disp("Range of x: From 0 to " + Xmax);
disp("Boosting Parameter: " + BoostingParam);

disp("Sub-Algorithm: " + SubAlg);
disp("Metric " + Metric);

if Momentum
    disp("Momentum parameter: " + Momentum);
end % End

disp("");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("%%%%%%%%%%%%%  Solution Summary   %%%%%%%%%%%%%%%%");

disp("OLLP Solution Time: " + OLLPtime * 1000 + " ms");
disp("LP  Solution  Time: " + LPtime * 1000 + " ms");
disp("Rounding      Time: " + Rtime * 1000 + " ms");
disp("Approximate LP objective: " + info.LPobj);
disp("Approximate IP objective: " + info.IPobj);

disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
disp("%%%%%%%%%%% OLLP Routine Successfully Ends %%%%%%%%%%%%");
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");

end % End of function


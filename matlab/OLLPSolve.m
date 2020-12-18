function [Lpx, PriceVec] = OLLPSolve(data, CheckInnerFeas, BoostingParam, ...
    Xmax, SubAlg, Metric, Momentum, Batch)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        OnLine Linear Progamming Solution Routine Matlab Interface       %
%                                                                         %
%                 Version 0.2     Nov 1st, 2020                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function serves as the main body of optimization algorithm in OLLP %
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
%           A         : Inequality Constraint Matrix of size m * n
%           b         : Inequality RHS vector, expected non-negative
%           c         : Linear coefficient vector of size n * 1
%
%      CheckInnerFeas : Boolean variable determining whether to check
%                       constraint feasibility in inner iterations
%                       Turned off by default
%      BoostingParam  : Boosting parameter in ESA
%                       Set 50 by default (must be >= 1)
%      Xmax           : An integer parameter implying an upperbound
%                       on x in integer programming formulation
%                       Set 1 by default (must be >= 1)
%      Adaptive(XX)   : Boolean variable specifying whether adaptive
%                       update is to be switched on
%      Type(XX)       : An integer specifying the type of problem to solve
%                       0 : Default knapsack problem
%                       1 : Set covering problem by transformation into MKP
%                       2 : Set covering problem by novel dual
%                           initialization
%                       Following methods are not yet implemented
%                       3 : General LP with indicator penalty
%                       4 : General LP with barrier penalty
%
%                       Above two parameters are now not in use
%
%      SubAlg         : A string specifying sub-algorithm used in online
%                       algorithm. Candidate algorithms include
%                       "SubGrad" : Default sub-gradient method
%                       "Proximal": Proximal Point Method
%                       "ADMM"    : Alternating directionb method of
%                                   multipliers
%      Metric         : A string specifying the Bregman divergence to be
%                       used. Currently supported metrics include
%                       "L2"      : Euclidean norm
%                       "KL"      : KL divergence
%      Momentum       : A real number between 0 and 1 specifying the
%                       momentum parameter in dual descent process. If
%                       set 0, then momentum is shut down
%
% Note that some parameters above are not compatible and we will ignore
% INCOMPATIBLE parameters based on certain priority.
%
% Usage
%       [Lpx, PriceVec] = OLLPSolve(data, BoostingParam, CheckInnerFeas);
%
% Output
%       Ipx                 : Approximate solution to IP problem
%       Lpx                 : Approximate solution to LP problem
%
%       info:
%       A struct containing information about optimization
%           OLLPtime        : Total running time of online algorithm
%           LPtime          : Time spent solving relaxed LP
%           Rtime           : Time spent rounding the solution
%           Probintro       : Description of problem data
%           PriceVec        : Price vector of online algorithm
%
% Function References:
% 1. CheckInput 2. SetDefaultParam 3. OLLPSolve 4. SolRound

% Check parameters
if (BoostingParam <= 0)
    error("Boosting parameter must be a positive integer");
elseif (Xmax < 1)
    error("Maximum of x must be a positive integer");
else
end % End if

% Unpack data
A = repmat(data.A, 1, Xmax);
b = data.b;
c = repmat(data.c, Xmax, 1);

K = BoostingParam;

if SubAlg == "SubGrad"
    [x, y] = OnlineSubGrad(A, b, c, K, CheckInnerFeas, Metric, Momentum);
elseif SubAlg == "Proximal"
    warning("Metric parameter ignored for proximal sub-algorithm");
    [x, y] = OnlineProx(A, b, c, K, CheckInnerFeas, Momentum);
elseif SubAlg == "ADMM"
    warning("Metric and Momentum parameter ignored for ADMM sub-algorithm");
    [x, y] = OnlineADMM(A, b, c, K, CheckInnerFeas);
else
    [x, y] = OnlineBatch(A, b, c, K, CheckInnerFeas, Metric, Momentum, Batch);
end % End if

% Get the relaxed LP solution
Lpx = x;
PriceVec = y;

end % End function


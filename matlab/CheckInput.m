function [Probintro, isValid] = CheckInput(data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        OnLine Linear Progamming Solution Routine Matlab Interface       %
%                                                                         %
%                 Version 0.1   July 28th, 2020                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function serves as the input check of OLLP solution routine        %
% referred from                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Li, Xiaocheng , C. Sun , and Y. Ye . "Simple and Fast Algorithm for  %
%    Binary Integer and Online Linear Programming." (2020).               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given input data, this function checks data feasibility and returns a
% boolean variable representing feasibility of data
% 
% The input format is specified as follows.
%
% Input
%      data: 
%      A struct containing two components X and Y to be tested
%           A             : Inequality Constraint Matrix of size m * n
%           b             : Inequality RHS vector, expected non-negative
%           c             : Linear coefficient vector of size n * 1
% Usage
%      [Probinfo] = OLLPSolve(data);
%
% Output
%       Probinfo          : A struct containing summaries of problem
%                           informationm including number of columns and
%                           rows, number of nonzeros and sparsity
%
% Function References:
% None

% Initialize validity
isValid = true;

% Unpack data
A = data.A;
b = data.b;
c = data.c;

% Unpack size information

disp("Checking size of problem ...");
[m, n] = size(A);
if ~(length(b) == m)
    isValid = false;
    error("Size mismatch: b, expected " + m + ", got " + length(b));
elseif ~(length(c) == n)
    isValid = false;
    error("Size mismatch: c, expected " + n + ", got" + length(c));
elseif min(b) <= 0
    isValid = false;
    error("Resource b must take nonnegative values");
else
    disp("Size Check Passed ......");
end % End for

if (min(c) < 0)
    warning("Linear coefficient c has negative values");
end % End if

% Create problem information struct
Probintro.nnz = 0;
Probintro.sps = 0;
Probintro.m = 0;
Probintro.n = 0;
Probintro.valid = isValid;

% Summarize problem introduction
if isValid
    Probintro.nnz = length(nonzeros(A));
    Probintro.sps = Probintro.nnz / (m * n);
    Probintro.m = m;
    Probintro.n = n;
end % End if

% Print problem introduction
disp("Problem data introduction:                           ");
disp("Number of rows: " + m + " Number of columns: " + n + "  ");
disp("Sparsity of A: " + Probintro.sps);

end % End of function




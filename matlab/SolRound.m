function [Ipx] = SolRound(A, b, Lpx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        OnLine Linear Progamming Solution Routine Matlab Interface       %
%                                                                         %
%                 Version 0.1   July 28th, 2020                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function serves as the final rounding process in OLLP solution     %
% routine with Extended Simple Algorithm from reference                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Li, Xiaocheng , C. Sun , and Y. Ye . "Simple and Fast Algorithm for  %
%    Binary Integer and Online Linear Programming." (2020).               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given input data a relaxed LP solution from OLLP, the function performs
% random rounding to get feasible solution to the original integer
% programming represented by
%
%       maximize     c^T * x
%       subject to   Ax <= b
%                     x in {0, 1} (or x in [M])
% 
% The input format is specified as follows.
%
% Input
%      A              : Inequality constraint matrix
%      b              : Inequality RHS resource vector
%      Lpx            : Soluion found by OLLP
%
% Usage
%       [Ipx] = SolRound(Lpx, Xmax);
%
% Output
%       Ipx                 : Approximate solution to IP problem
%
% Function References:
% None

% Unpack dimension data
n = length(Lpx);

% Initialize x with zeros
Ipx = zeros(n, 1);

% Get permutation
p = randperm(n);

% Start rounding
for i = p
    a = A(:, i);
    xfloor = floor(Lpx(i));
    rdx = xfloor + (rand() <= (Lpx(i) - xfloor));
    consp = b - rdx * a;
    if min(consp) >= 0
        b = consp;
        Ipx(i) = rdx;
    end % End if
end % End for

end % End of function

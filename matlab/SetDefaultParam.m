function [params] = SetDefaultParam()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        OnLine Linear Progamming Solution Routine Matlab Interface       %
%                                                                         %
%                 Version 0.2     Nov 1st, 2020                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function serves as default parameter setting of OLLP solution      %
% routine referred from                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Li, Xiaocheng , C. Sun , and Y. Ye . "Simple and Fast Algorithm for  %
%    Binary Integer and Online Linear Programming." (2020).               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function initialize a parameter struct for later modification or in
% case no parameter is specified
%
% The input format is specified as follows.
%
% Input
%       None
% Usage
%       params = SetDefaultParam();
%
% Output
%       params : A struct containing default parameters for OLLP
%
% Function References:
% None

% Set default parameters
params.CheckInnerFeas = false;
params.BoostingParam = 50;
params.Xmax = 1;
params.Adaptive = 0;
params.SubAlg = "SubGrad";
params.Metric = "L2";
params.Momentum = 0;
params.Batch = 1;
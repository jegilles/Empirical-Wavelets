function [bounds,plane,L]=EWT_GSS_BoundariesDetect(f,type)

% =========================================================================
% function [bounds,plane,L]=EWT_GSS_BoundariesDetect(f,type)
%
% This function performs the Gaussian Scale-Space Boundaries detection.
%
% Inputs:
%   -f: original function/histogram
%   -type: type of wanted threshold detection method (otsu,halfnormal,
%          empiricallaw,mean,kmeans)
%
% Outputs:
%   -bounds: detected boundaries
%   -plane: the built scale-space plane
%   -L: set of minima curve lengths
%
% Author: Jerome Gilles
% Institution: UCLA - Department of Mathematics
% Year: 2013
% Version: 1.0
% =========================================================================

% build the Gaussian Scale-Space representation
plane=EWT_PlanGaussianScaleSpace(f);

% find the meaningful boundaries
[bounds,L]=EWT_MeaningfulScaleSpace(f,plane,type);
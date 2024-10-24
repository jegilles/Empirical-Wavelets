function boundaries = EWT_Angles_Detect(f,params)

%================================================================
% function boundaries = EWT_Angles_Detect(f,params)
%
% This function segments f into a certain amount of supports by 
% using different techniques: 
% - middle point between consecutive local maxima,
% - lowest minima between consecutive local maxima,
%
% Moreover some preprocessing are available in order to remove a
% global trend:
% - by substracting a power law approximation
% - by substracting a polynomial approximation
% - by substracting the average of the opening and closing of ff
%   (the size of the structural element is automatically detected)
% - by substracting the opening of ff
%   (the size of the structural element is automatically detected)
%
% Note: the detected boundaries are given in term of indices
%
% Inputs:
%   -ff: the function to segment
%   -params: structure containing the following parameters:
%       -params.log: 0 or 1 to indicate if we want to work with
%                    the log of the ff
%       -params.curvpreproc: 'none','plaw','poly','morpho','tophat'
%       -params.curvreg: 'none','gaussian','closing'
%       -params.curvlengthFilter: size of the filters used in the above
%                             regularization methods
%       -params.curvsigmaFilter: standard deviation for the above Gaussian
%                            regularization
%       -params.curvmethod: 'locmax','locmaxmin'
%       -params.curvN: maximum number of supports (needed for the
%                  locmax and locmaxmin methods)
%       -params.curvdegree: degree of the polynomial (needed for the
%                       polynomial approximation preprocessing)
%
% Outputs:
%   -boundaries: list of detected boundaries
%
% Author: Jerome Gilles
% Institution: UCLA - Department of Mathematics
% Year: 2013
% Version: 1.0
%===============================================================

params2=params;
params2.globtrend = params.curvpreproc;
params2.reg = params.curvreg;
params2.degree = params.curvdegree;
params2.lengthFilter = params.curvlengthFilter;
params2.sigmaFilter = params.curvsigmaFilter;

% Global trend removal
presig = EWT_RemoveTrend(f,params2);

% Regularization
presig = EWT_SpectrumRegularize(presig,params2);

% We perform the detection itself
switch lower(params.curvmethod)
    case 'locmax'
        boundaries = EWT_AnglesLocalMax(presig,params.curvN);
    case 'locmaxmin'
        %% We extract the lowest local minima between to selected local maxima
        boundaries = EWT_AnglesLocalMaxMin(presig,params.curvN);
    case 'scalespace'
        boundaries=EWT_GSS_BoundariesDetect(presig,params.typeDetect);
end
function [boundaries,presig] = EWT_Boundaries_Detect(f,params)

%==========================================================================
% function [boundaries,presig] = EWT_Boundaries_Detect(f,params)
%
% This function segments f into a certain amount of supports by 
% using different technics: 
% - middle point between consecutive local maxima,
% - lowest minima between consecutive local maxima,
% - fine to coarse histogram segmentation algorithm,
% - from an a priori set of boundaries then adjusted in some neighborhood 
%   (the lowest minima can be computed on the original spectrum or a
%   regularized version)
% - automatic scale space detection of meaningful minima
%
% Some preprocessing are available in order to remove a
% global trend of the spectrum:
% - by substracting a power law approximation
% - by substracting a polynomial approximation
% - by substracting the average of the opening and closing of ff
%   (the size of the structural element is automatically detected)
% - by substracting the opening of ff
%   (the size of the structural element is automatically detected)
%
% Moreover, regularized version of the spectrum can be obtained by the
% following methods:
% - Gaussian filtering (its parameters are filter of width 
%   params.lengthFilter and standard deviation params.sigmaFilter)scalesp
% - Average filtering (its parameters are filter of width 
%   params.lengthFilter)
% - Morphological closing, return the upper envelope (structural element 
%   has width params.lengthFilter)
%
% Note: the detected boundaries are given in term of indices
%
% Inputs:
%   -f: the function to segment
%   -params: structure containing the following parameters:
%       -params.log: 0 or 1 to indicate if we want to work with
%                    the log of the ff
%       -params.preproc: 'none','plaw','poly','morpho,'tophat'
%       -params.reg: 'none','gaussian','closing','average'
%       -params.lengthFilter: size of the filters used in the above
%                             regularization methods
%       -params.sigmaFilter: standard deviation for the above Gaussian
%                            regularization
%       -params.detect: 'locmax','locmaxmin','locmaxminf','adaptive',
%                       'adaptivereg','scalespace'
%       -params.N: maximum number of supports (needed for the
%                  locmax and locmaxmin methods)
%       -params.degree: degree of the polynomial (needed for the
%                       polynomial approximation preprocessing)
%       -params.completion: 0 or 1 to indicate if we try to complete
%                           or not the number of modes if the detection
%                           find a lower number of mode than params.N
%       -params.InitBounds: vector of initial bounds (in index domain)
%                           needed for the adaptive and adaptivereg methods
%       -params.typeDetect: (for scalespace method only) 'otsu',
%                           'halfnormal','empiricallaw','mean','kmeans'
%
% Outputs:
%   -boundaries: list of detected boundaries
%   -presig: preprocessed spectrum
%
% Author: Jerome Gilles + Kathryn Heal
% Institution: UCLA - Department of Mathematics
% Year: 2013
% Version: 2.0
%==========================================================================

% Apply the log if needed
if params.log==1
    f=log(f);
end

% Global trend removal
presig = EWT_RemoveTrend(f,params);

% Regularization
presig = EWT_SpectrumRegularize(presig,params);
%plot(presig);

% Boundaries detection
switch lower(params.detect)
    case 'locmax'
        boundaries = EWT_LocalMax(presig,params.N);
        %% Mid-point between two consecutive local maxima computed on the regularized
        %% spectrum
    case 'locmaxmin'
        %% We extract the lowest local minima between two selected local maxima
        boundaries = EWT_LocalMaxMin(presig,params.N);
    case 'locmaxminf'
        boundaries = EWT_LocalMaxMin(presig,params.N,f);
        %% We extract the lowest local minima on the original spectrum between 
        %% two local maxima selected on the regularized signal
    case 'adaptivereg'
        boundaries = EWT_Adaptive_Bounds_Adapt(presig,params);
        %% epsilon-neighborhood method entirely computed on the regularized signal
    case 'adaptive'
        boundaries = EWT_Adaptive_Bounds_Adapt(presig,params,f);
        %% epsilon-neighborhood method: local maxima computed on the regularized signal
        %% and lowest minima computed on the original spectrum
    case 'scalespace'
        boundaries=EWT_GSS_BoundariesDetect(presig,params.typeDetect);
end


%% If asked, and needed perform the completion of the number of modes
if params.completion==1
    boundaries=EWT_Boundaries_Completion(boundaries,params.N-1);
end

function [bounds,L]=EWT_MeaningfulScaleSpace(f,plane,type)

% =========================================================================
% function [bounds,L]=EWT_MeaningfulScaleSpace(f,plane,type)
%
% This function calls the different methods to detect the meaningful scale
% threshold from the set of initial minima curve lengths and then find the
% meaningful boundaries (it manages the merging in the scale-space plane)
%
% Inputs:
%   f: input function
%   plane: the scale-space representation of f
%   type: the wanted detection method (otsu,halfnormal,empiricallaw,mean,
%         kmeans)
%
% Outputs:
%   bounds: the detected meaningful boundaries
%   L: set of minima curve lengths
%
% Author: Jerome Gilles
% Institution: UCLA - Department of Mathematics
% Year: 2013
% Version: 1.0
% =========================================================================

[L,ind]=EWT_LengthScaleCurve(plane);
%PlotStats(L);

switch lower(type)
    case 'otsu'
        [bounds,th]=EWT_OtsuMethod(L,ind);
    case 'halfnormal'
        [bounds,th]=EWT_HalfNormalLaw(L,ind,size(plane,2));
    case 'empiricallaw'
        [bounds,th]=EWT_EmpiricalLaw(L,ind);
    case 'mean'
        [bounds,th]=EWT_MeanTh(L,ind);
    case 'kmeans'
        [bounds,th]=EWT_kmeansDetect(L,ind);
end
bounds=EWT_RemoveMerge(f,plane,bounds,th);
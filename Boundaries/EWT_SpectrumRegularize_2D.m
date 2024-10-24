function ff = EWT_SpectrumRegularize_2D(f,params )

%================================================================================
% function ff = EWT_SpectrumRegularize_2D(f,params )
%
% This function performs the regularization indicated in the variable params.reg
% on f. The available options are:
%   -'none': return the original f.
%   -'gaussian': apply a gaussian filter on f of width params.lengthFilter
%   and standard deviation params.sigmaFilter.
%   -'average': apply an average filter on f of width params.lengthFilter
%
% Inputs:
%   -f: the function to regularize
%   -params: parameter structure
%
%  Output:
%   -ff: regularized function
%
% Author: Jerome Gilles
% Institution: SDSU - Department of Mathematics & Statistics
% Year: 2022
% Version: 1.0
%================================================================================

switch lower(params.reg)
    case 'gaussian'
        ff = imgaussfilt(f,params.sigmaFilter);
    case 'average'
        filter = fspecial('average',[params.lengthFilter,params.lengthFilter]); 
        ff = conv2(f,filter,'same');
    otherwise
        ff = f;
end


function ff = EWT_SpectrumRegularize( f, params )

%================================================================================
% function ff = EWT_SpectrumRegularize( f, params )
%
% This function performs the regularization indicated in the variable params.reg
% on f. The available options are:
%   -'none': return the original f.
%   -'gaussian': apply a gaussian filter on f of width params.lengthFilter
%   and standard deviation params.sigmaFilter.
%   -'average': apply an average filter on f of width params.lengthFilter
%   -'closing': compute the upper envelope via a morphological closing
%   operator, the structural element has width params.lengthFilter
%
% Inputs:
%   -f: the function to regularize
%   -params: parameter structure
%
%  Output:
%   -ff: regularized function
%
% Author: Kathryn Heal and Jerome Gilles
% Institution: UCLA - Department of Mathematics
% Year: 2013
% Version: 1.0
%================================================================================

switch lower(params.reg)
    case 'none'  
        ff=f;
    case 'gaussian'
        filter = fspecial('gaussian',[1,params.lengthFilter], params.sigmaFilter ); 
        ff = conv(f,filter,'same');
    case 'average'
        filter = fspecial('average',[1,params.lengthFilter] ); 
        ff = conv(f,filter,'same');
    case 'closing'
        ff = EWT_FunctionClosing(f,params.lengthFilter);
end


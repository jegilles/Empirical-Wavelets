function [supports,num_regions] = EWWT_Get_Boundaries(absff,centers,params)
%==========================================================================
% function [supports,num_regions] = EWWT_Get_Boundaries(absff,centers,params)
%
% This function partitions the Fourier domain based on a set of points
% which represent mode locations
%   
% Input:
%   -absff: magnitude spectrum of an image
%   -centers: detected mode locations on magnitude spectrum
%   -params: empirical wavelet parameters, notably a flag on whether image
%            is real or complex
% Output:
%   -supports: label matrix representing detected partitioning of the Fourier 
%       spectrum
%   -num_regions: array holding numbers for each region
%
% Author: Basile Hurat - Jerome Gilles
% Institution: SDSU - Department of Mathematics
% Year: 2022
% Version: 1.0
%==========================================================================

local_max = zeros(size(absff));
for i = 1:size(centers,1)    
    local_max(centers(i,1),centers(i,2)) = 1;
end
% make the maxima to be minima as expected by the watershed function
flipped_absff = imimposemin(-absff,local_max);

% compute the watershed
supports = watershed(flipped_absff,8);


% Split watershed regions into separate supports, pair for symmetry in the real
% case
if params.complex == 0
    for i = 1:length(supports(:))
        supports(end-i+1) = min(supports(i),supports(end-i+1));
    end
    supports = EWWT_removeBounds(supports,1);
else
    supports = EWWT_removeBounds(supports,0);
end
num_regions = unique(supports(:));
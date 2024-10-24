function boundaries = EWT_Adaptive_Bounds_Adapt(f,params,fm)

%==========================================================================
% function boundaries = EWT_Adaptive_Bounds_Adapt(presig,params,fm)
%
% This function adapt an initial set of boundaries to the studied signal.
% First it computes some neighborhood from the initial boundaries then it
% detects the global minima in each neighborhood.
% If the input fm is provided then the local maxima are computed on the
% input signal f and the local minima on fm otherwise everything will be
% computed from f.
%
% Inputs:
%   -f: input signal
%   -params: input parameters
%            -params.InitBounds: vector of initial bounds (in index domain)
%
% Optional input:
%   -fm: function on which the local minima will be computed
%
% Output:
%   -boundaries: the detected set of boundaries
%
% Author: Kathryn Heal and Jerome Gilles
% Institution: UCLA - Department of Mathematics
% Year: 2013
% Version: 1.0
%==========================================================================

if nargin<3
    fm=f;
end

boundaries=zeros(length(params.InitBounds),1);

% Initial boundaries in indices space (not integer) + endpoints
spect_bounds=[0 params.InitBounds length(f)-1];

for i=2:(length(spect_bounds)-1)
    % Detect an asymetric epsilon-neighborhood (EN)
    neighb_low = round(spect_bounds(i)-round(abs(spect_bounds(i)-spect_bounds(i-1))/2));
    neighb_hi = round(spect_bounds(i)+round(abs(spect_bounds(i+1)-spect_bounds(i))/2));

    [~,imini]=min(fm(neighb_low:neighb_hi));
    boundaries(i-1) = imini+neighb_low-2;
end

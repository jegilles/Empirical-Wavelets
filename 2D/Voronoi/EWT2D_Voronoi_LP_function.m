function y=EWT2D_Voronoi_LP_function(x,tau)

%===============================================================================
% function y=EWT2D_Voronoi_LP_function(x,tau)
%
% This function calculate the Littlewood-Paley value based on the provided
% distance.
%   
% Input:
%   -x: distance from the edge 
%   -tau: transition width
%
% Output:
%   -y: value of the filter.
%
% Author: Basile Hurat - Jerome Gilles
% Institution: SDSU - Department of Mathematics
% Year: 2022
% Version: 1.0
%===============================================================================

if x>tau
    y=1;
elseif x<-tau
    y=0;
else
    y=cos(pi*EWT_beta((tau-x)/(2*tau))/2);
end
function ind = EWT2D_get_voronoi_cell_index(j,dist,sizeImg)

% ==============================================================================
% function ind = EWT2D_get_voronoi_cell_index(i,j,dist,maxima,sizeImg)
%
% This function properly build the voronoi cells by guaranteing a symmetrical
% selection (useful to build symmetric filters for a real transform)
%
% Inputs:
%   -j: coordinates of the current position
%   -dist: vector containing the distances between the current position and each
%   maxima
%   -sizeImg: vector containing the size of the processed spectrum
%
% Outputs:
%   -ind: label for the voronoi cell
%
% Author: Jerome Gilles
% Institution: SDSU - Department of Mathematics
% Year: 2022
% Version: 1.0
% ==============================================================================

% find indices of the minima of dist:
cand = find(dist==min(dist));

if length(cand) ==1 % if there is only one closest maxima 
    ind = cand;
else % if there are several maxima equidistants to the current position
    c = floor(length(dist)/2)+1;

    if (cand(1)<c) && (cand(end)<=c)
        ind = cand(1);
    elseif (cand(1)>=c) && (cand(end)>c)
        ind = cand(end);
    else
        if j<(floor(sizeImg(2)/2)+1)
            ind = cand(1);
        else
            ind = cand(end);
        end
    end
end
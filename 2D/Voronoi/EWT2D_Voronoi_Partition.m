function [labelImage, voronoi_cells] = EWT2D_Voronoi_Partition(centroids,sizeImg)

%===============================================================================
% function [labelImage, voronoi_cells, alldist] = EWT2D_Voronoi_Partition(centroids,sizeImg)
%
% This function return the Voronoi partition corresponding to the provided
% centroids. The cells are properly assigned in order to guarantee a symmetric
% partition in the real case. (Note: we used the square of the Euclidean 
% distance to avoid computing radicals)
%   
% Input:
%   -centroids: vector containing the coordinates of each centroid 
%   -sizeImg: spectrum size
%
% Output:
%   -labelImage: image containing the Voronoi partition tagged from 1 to the 
%                number of cells
%   -voronoi_cells: collection of binary images containing the mask of each 
%                   Voronoi cell
%
% Author: Basile Hurat - Jerome Gilles
% Institution: SDSU - Department of Mathematics
% Year: 2022
% Version: 1.0
%===============================================================================

labelImage = zeros(sizeImg);

%find closest centroid for each voxel
for i = 1:size(labelImage,1)
    for j = 1:size(labelImage,2)
         eucDistance = sqrt((centroids(:,1)-i).^2+(centroids(:,2)-j).^2);
         labelImage(i,j) = EWT2D_get_voronoi_cell_index(j,eucDistance,sizeImg);
    end
end

%% Extract each Voronoi cell mask
voronoi_cells = cell(size(centroids,1),1);
for i = 1:size(centroids,1)
    voronoi_cells{i} = zeros(sizeImg);
    voronoi_cells{i}(labelImage == i) = 1; 
end
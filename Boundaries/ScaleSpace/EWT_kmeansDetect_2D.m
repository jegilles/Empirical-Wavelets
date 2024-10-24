function [bounds,th]=EWT_kmeansDetect_2D(L,ind)
% =========================================================================
% function [bounds,th]=EWT_kmeansDetect_2D(L,ind)
%
% This function classifies the set of minima curve lengths stored in L into
% two classes by a kMeans algorithm. Then it returns the boundaries which 
% are supposed to be the meaningful ones. (Modified by Basile Hurat to work
% for matrices)
%
% Inputs:
%   L: set of minima curve lengths
%   ind: original index of each minima
%
% Outputs:
%   bounds: detected bounds
%   th: detected scale threshold
%
% Author: Jerome Gilles (modified for 2D by Basile Hurat)
% Institution: UCLA - Department of Mathematics
% Year: 2013
% Version: 1.0
% =========================================================================

clust=kmeans(L,2,'start','uniform','emptyaction','singleton','replicates',10); %L2+uniform initialization

[~,i]=max(L);
nc=clust(i);

bounds=ind(clust==nc,:);
th=min(L(clust==nc));
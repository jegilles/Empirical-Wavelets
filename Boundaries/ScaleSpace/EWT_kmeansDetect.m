function [bounds,th]=EWT_kmeansDetect(L,ind)

% =========================================================================
% function [bounds,th]=EWT_kmeansDetect(L,ind)
%
% This function classifies the set of minima curve lengths stored in L into
% two classes by a kMeans algorithm. Then it returns the boundaries which 
% are supposed to be the meaningful ones.
%
% Inputs:
%   L: set of minima curve lengths
%   ind: original index of each minima
%
% Outputs:
%   bounds: detected bounds
%   th: detected scale threshold
%
% Author: Jerome Gilles
% Institution: UCLA - Department of Mathematics
% Year: 2013
% Version: 1.0
% =========================================================================

clust=kmeans(L,2,'start','uniform','emptyaction','singleton','replicates',10); %L2+uniform initialization
%clust=kmeans(L,2,'emptyaction','singleton','replicates',10); %L2+random initialization
%clust=kmeans(L,2,'distance','cityblock','start','uniform','emptyaction','singleton','replicates',10); %L1+uniform initialization
%clust=kmeans(L,2,'distance','cityblock','emptyaction','singleton','replicates',10); %L1+random initialization

[~,i]=max(L);
nc=clust(i);

bounds=ind(clust==nc);
th=min(L(clust==nc));
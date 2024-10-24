function symvorcells = EWT2D_merge_symmetric(in)

%===============================================================================
% function realcells = EWT2D_merge_symmetric(in)
%
% This function pairs the Voronoi cell masks in order to have symmetric
% partitions. The mask contining the frequency origin is assigned first.
%   
% Input:
%   -in: collection of Voronoi cell masks
%
% Output:
%   -symvorcells: collection of paired cells
%
% Author: Basile Hurat - Jerome Gilles
% Institution: SDSU - Department of Mathematics
% Year: 2022
% Version: 1.0
%===============================================================================

symvorcells = cell(floor(length(in)/2)+1,1);

i0 = floor(size(in{1},1)/2)+1;
j0 = floor(size(in{1},2)/2)+1;

% find the cell centered at origin
for n=1:length(in)
   if in{n}(i0,j0) == 1
       symvorcells{1} = in{n};
       ind = n;
   end
end
in(ind) = [];

% find the pairs
curin = 1;
for cur=2:length(symvorcells)
    [i0,j0] = find(in{curin});
    syc = size(in{1})-round([sum(i0) sum(j0)]/length(i0))+1;

    for n=curin+1:length(in)
        if in{n}(syc(1),syc(2))==1
           symvorcells{cur} = in{curin}+in{n}; 
           break
        end        
    end
    curin = curin + 1;
end
function partition = EWWT_removeBounds(partition,sym)
%==========================================================================
%function partition = EWWT_removeBounds(partition,sym)
%
% This function takes a partitioning from the watershed function and
% assigns the watershed line to a region, potentially in a symmetric fashion
%   
% Input:
%   -partition: a label image corresponding to a partition
%   -sym: set to 1 if the partition is symmetric
% Output:
%   -partition: resulting partition without watershed line
%
% Author: Basile Hurat - Jerome Gilles
% Institution: SDSU - Department of Mathematics
% Year: 2022
% Version: 1.0
%==========================================================================

for i = 1:length(partition(:))
    if partition(i) == 0
        lab = checkNeighborhood(partition,i);
        partition(i) = lab;
        if sym == 1
            partition(end-i+1) = lab;
        end
    end    
end



function reg = checkNeighborhood(L,i)
[h,~] = size(L);
m = mod(i,h);
n = ceil(i/h);
for i = -1:1
    for j = -1:1
        if i == j
            continue
        end
        if m+i < 1 ||m+i > size(L,1) || n+j < 1 || n+j > size(L,2)
            continue
        end
        if L(m+i,n+j) ~= 0
            reg = L(m+i,n+j);
        end
    end
end
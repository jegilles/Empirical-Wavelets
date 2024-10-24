function [bounds,th]=EWT_MeanTh(L,ind)

% =========================================================================
% function [bounds,th]=EWT_MeanTh(L,ind)
%
% This function classifies the set of minima curve lengths stored in L into
% two classes. The thresold is computed as the mean of L. Then it returns 
% the boundaries which are supposed to be the meaningful ones.
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

th=ceil(mean(L));

%keep only meaningful minima
Lth=L;
for i=1:length(L)
   if L(i)<th
       Lth(i)=0;
   else
       Lth(i)=1;
   end
end
bounds=ind(Lth>0);
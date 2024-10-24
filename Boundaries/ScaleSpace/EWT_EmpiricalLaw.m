function [bounds,th]=EWT_EmpiricalLaw(L,ind)
% =========================================================================
% function [bounds,th]=EWT_EmpiricalLaw(L,ind)
%
% This function classifies the set of minima curve lengths stored in L into
% two classes by considering length which are epsilon-meaningful for an 
% empirical law. Then it returns the boundaries which are supposed to be 
% the meaningful ones.
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

%build normalized histogram
histo=hist(L,max(L));
chisto=cumsum(histo/sum(histo));

th=find(chisto>(1-1/length(L)));
th=th(1);

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

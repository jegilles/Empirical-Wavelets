function [bounds,th]=EWT_HalfNormalLaw(L,ind,Lmax)

% =========================================================================
% function [bounds,th]=EWT_HalfNormalLaw(L,ind)
%
% This function classifies the set of minima curve lengths stored in L into
% two classes by considering length which are epsilon-meaningful for an 
% half-normal law fitted to the data. Then it returns the boundaries which 
% are supposed to be the meaningful ones.
%
% Inputs:
%   L: set of minima curve lengths
%   ind: original index of each minima
%   Lmax: maximum possible length of a minima curve
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

%first estimate sigma
sigma=sqrt(pi/2)*mean(L);

%compute threshold
th=sqrt(2)*sigma*erfinv(erf(Lmax/(sqrt(2)*sigma))-1/length(L));

%keep only meaningful minima
Lth=L;
for i=1:length(L)
   if L(i)<=th
       Lth(i)=0;
   end
end
bounds=ind(Lth~=0);
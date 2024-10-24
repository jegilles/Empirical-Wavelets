function bound = EpsNeighLocalMaxMin(f,fm)

%=======================================================================
% function bound = EpsNeighLocalMaxMin(f,fm)
%
% This function returns the global minima between the two highest
% local maxima (if the function can't find two maxima then it returns
% the global minima over the entire f). If the input fm is provided 
% then the local maxima are computed on f and the local minima on fm 
% otherwise both are computed on f (this is useful if you want to compute 
% the maxima on a regularized version of your signal while detecting the 
% "true" minima).
%
% Input:
%   -f: input function
% Optional input:
%   -fm: function on which the local minima will be computed
%
% Output:
%   -bound: the global minima between the two highest local maxima
%
% Author: Kathryn Heal and Jerome Gilles
% Institution: UCLA - Department of Mathematics
% Year: 2013
% Version: 1.0
%=======================================================================

if nargin>2
    f2=fm;
else
    f2=f;
end

locmax=zeros(size(f));
locmin=max(f)*ones(size(f));
% check if the endpoints are maxima or minima
if (f(1)>f(2))
    locmax(1)=f(1);
elseif (f2(1)<f2(2))
    locmin(1)=f2(1);
end

% We detect the local maxima and minina
for i=2:length(f)-1
    if ((f(i-1)<f(i)) && (f(i)>f(i+1)))
        locmax(i)=f(i);                       
    end
    
    if ((f2(i-1)>f2(i)) && (f2(i)<f2(i+1)))
        locmin(i)=f2(i);
    end
end

% We keep the two highest maxima and their index
[lmax,Imax]=sort(locmax,1,'descend');
if length(lmax)>=2 % If two local maxima exist we keep them
    Imax=sort(Imax(1:2));
else % If there's only one or no local maxima we keep the entire neighborhood
    Imax=[1 length(f)];
end

% We search the global minimum in the previously detected segment
[lmin,ind]=sort(locmin(Imax(1):Imax(2)));
n=1;
tmp=lmin(1);
if n<length(lmin)
       n=2;
       while ((n<=length(lmin)) && (tmp==lmin(n)))
            n=n+1;
       end
end
bound=Imax(1)+ind(ceil(n/2))-2;

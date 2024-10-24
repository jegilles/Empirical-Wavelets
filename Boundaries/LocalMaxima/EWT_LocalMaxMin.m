function bound = EWT_LocalMaxMin(f,N,fm)

%===================================================================
% function bound = EWT_LocalMaxMin(f,N,fm)
%
% This function segments f into a maximum of N supports by detecting
% the lowest local minima between the N largest local maxima. If the
% input fm is provided then the local maxima are computed on f and 
% the local minima on fm otherwise both are computed on f (this is
% useful if you want to compute the maxima on a regularized version
% of your signal while detecting the "true" minima).
%
% Note: the detected boundaries are given in term of indices
%
% Inputs:
%   -f: the function to segment (must a vertical vector)
%   -N: maximal number of bands
%
% Optional input:
%   -fm: function on which the local minima will be computed
%
% Outputs:
%   -bound: list of detected boundaries
%
% Author: Jerome Gilles
% Institution: UCLA - Department of Mathematics
% Year: 2013
% Version: 2.0
%===================================================================

locmax=zeros(size(f));

if nargin>2
    f2=fm;
else
    f2=f;
end

locmin=max(f2)*ones(size(f2));
% We detect the local maxima and minina
for i=2:length(f)-1
    if ((f(i-1)<f(i)) && (f(i)>f(i+1)))
        locmax(i)=f(i);
    end
    
    if ((f2(i-1)>f2(i)) && (f2(i)<f2(i+1)))
        locmin(i)=f2(i);
    end
end

% We keep the N-th highest maxima and their index
if N~=-1
    N=N-1;
    [lmax,Imax]=sort(locmax,1,'descend');
    if length(lmax)>N
        Imax=sort(Imax(1:N));
    else
        Imax=sort(Imax);
        N=length(lmax);
    end

    % We detect the lowest minima between two consecutive maxima
    bound=zeros(1,N);
    for i=1:N
       if i==1
           a=1;
       else
           a=Imax(i-1);
       end
       [lmin,ind]=sort(locmin(a:Imax(i)));
       tmp=lmin(1);
       n=1;
       if n<length(lmin)
           n=2;
           while ((n<=length(lmin)) && (tmp==lmin(n)))
                n=n+1;
           end
       end
       bound(i)=a+ind(ceil(n/2))-2;
    end
else
    k=1;
    for i=1:length(locmin)
       if locmin(i)<max(f2)
          bound(k) = i-1;
          k=k+1;
       end
    end
end

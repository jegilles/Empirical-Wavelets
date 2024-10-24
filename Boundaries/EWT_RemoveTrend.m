function presig = EWT_RemoveTrend(f,params)

%===========================================================================
% function presig = EWT_RemoveTrend(f,params)
%
% This function remove a global trend from f. The available methods,selected 
% via the field params.preproc, are 
%   -'none': do nothing
%   -'plaw': fit a power law to f and then subtract this power law to f
%   -'poly': fit a polynomial interpolation (the degree is specified by
%   params.degree) to f and then subtract this interpolation to f.
%   -'morpho': compute the lower and upper envelopes of f by the
%   morphological opening and closing operators. Then the average of these
%   envelopes is subtracted to f.
%   -'tophat': compute the morphological Top Hat operator of f.
%   -'opening': compute the morphological Opening operator of f.
%
% Inputs:
%   -f: the function to regularize
%   -params: parameter structure
%
%  Output:
%   -presig: regularized function
%
% Author: Kathryn Heal and Jerome Gilles
% Institution: UCLA - Department of Mathematics
% Year: 2013
% Version: 1.0
%===========================================================================

switch lower(params.globtrend)
    case 'none'
        %% No preprocessing
        presig = f;
    case 'plaw'
        %% power law substraction
        [~,~,presig]=EWT_Powerlaw_Estimator(f);
    case 'poly'
        %% Polynomial interpolation substraction
        w=0:length(f)-1;
        w=w';
        [p,~]=polyfit(w,f,params.degree);
        presig=f-polyval(p,w);
    case 'morpho'
        %% (Opening+Closing)/2 substraction
        % We detect first the size of the structural element
        % as the smallest distance between two consecutive maxima +1
        locmax=zeros(size(f));
        % We detect the local maxima
        for i=2:length(f)-1
            if ((f(i-1)<f(i)) && (f(i)>f(i+1)))
                locmax(i)=f(i);
            end
        end
        sizeel=length(f);
        n=1;np=1;
        while (n<length(locmax)+1)
           if (locmax(n)~=0)
                if sizeel>(n-np)
                    sizeel=n-np;
                end
                np=n;
                n=n+1;
           end
           n=n+1;
        end
        presig=f-(EWT_FunctionOpening(f,sizeel+1)+EWT_FunctionClosing(f,sizeel+1))/2;
    case 'tophat'
        %% Opening substraction (TopHat operator)
        % We detect first the size of the structural element
        % as the smallest distance between two consecutive maxima +1
        locmax=zeros(size(f));
        % We detect the local maxima
        for i=2:length(f)-1
            if ((f(i-1)<f(i)) && (f(i)>f(i+1)))
                locmax(i)=f(i);
            end
        end
        sizeel=length(f);
        n=1;np=1;
        while (n<length(locmax)+1)
           if (locmax(n)~=0)
                if sizeel>(n-np)
                    sizeel=n-np;
                end
                np=n;
                n=n+1;
           end
           n=n+1;
        end
        presig=f-EWT_FunctionOpening(f,sizeel+1);
    case 'opening'
        %% Opening operator
        n0=length(f);
        f=[f(end:-1:2);f];
        % We detect first the size of the structural element
        % as the smallest distance between two consecutive maxima +1
        locmax=zeros(size(f));
        % We detect the local maxima
        for i=2:length(f)-1
            if ((f(i-1)<f(i)) && (f(i)>f(i+1)))
                locmax(i)=f(i);
            end
        end
        sizeel=length(f);
        n=1;np=1;
        while (n<length(locmax)+1)
           if (locmax(n)~=0)
                if sizeel>(n-np)
                    sizeel=n-np;
                end
                np=n;
                n=n+1;
           end
           n=n+1;
        end
        presig=EWT_FunctionOpening(f,sizeel+1);
        presig=presig(n0:end);  
end
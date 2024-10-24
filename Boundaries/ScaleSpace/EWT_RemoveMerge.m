function bounds = EWT_RemoveMerge(f,plane,bounds,th)

%==========================================================================
% function bounds = EWT_RemoveMerge(f,plane,bounds,th)
%
% This function manage local minima which merge at some point in the
% scale-space plane according to the following rules:
% - if the mergin occur before the scale th then we keep only one minima
%   (the lowest one) as they are not individually meaningful
% - if the mergin occur after the scale th then we consider that each 
%   initial minima is meaningful and we keep them
%
% -Inputs:
%   -f: the initial histgoram
%   -plane: the scale-space plane of f
%   -bounds: initial detected bounds
%   -th: scale threshold
%
% -Outputs:
%   -bounds: updated bounds
%
% Author: Jerome Gilles
% Institution: UCLA - Department of Mathematics
% Year: 2013
% Version: 1.0
% =========================================================================

if ~isempty(bounds)

    tagplane=zeros(size(plane));
    indrem=zeros(size(bounds));

    % We tag the first curve
    tag=bounds(1);
    stop=0;
    i=tag;
    j=1;
    while stop~=1
        tagplane(i,j)=tag;
        if i>1
           if plane(i-1,j+1)==1
               i=i-1;
               j=j+1;
           elseif plane(i,j+1)==1
               j=j+1;
           elseif plane(i+1,j+1)==1
               i=i+1;
               j=j+1;
           else
               stop=1;
           end
        else
           if plane(i,j+1)==1
               j=j+1;
           elseif plane(i+1,j+1)==1
               i=i+1;
               j=j+1;
           else
               stop=1;
           end 
        end

        if (j>th) || (j==(size(plane,2)-1))
            stop=1;
        end
    end

    % We address the other curves
    for k=2:length(bounds)
       tag=bounds(k);
       i=tag;
       j=1;
       stop=0;
       retag=0;
       while stop~=1
            tagplane(i,j)=tag;
            if i>1
               if plane(i-1,j+1)==1
                   if (tagplane(i-1,j+1)==bounds(k-1)) && (retag==0) %we found a merge
                       if f(bounds(k-1))<f(bounds(k)) %we keep the previous minima
                           indrem(k)=1;
                           stop=1;
                       else %we keep the current minima and need to retag the rest of the curve
                           indrem(k-1)=1;
                           retag=1;
                       end
                   end
                   i=i-1;
                   j=j+1;
               elseif plane(i,j+1)==1
                   if (tagplane(i,j+1)==bounds(k-1)) && (retag==0) %we found a merge
                       if f(bounds(k-1))<f(bounds(k)) %we keep the previous minima
                           indrem(k)=1;
                           stop=1;
                       else %we keep the current minima and need to retag the rest of the curve
                           indrem(k-1)=1;
                           retag=1;
                       end
                   end
                   j=j+1;
               elseif plane(i+1,j+1)==1
                   i=i+1;
                   j=j+1;
               else
                   stop=1;
               end
            else
               if plane(i,j+1)==1
                   if (tagplane(i,j+1)==bounds(k-1)) && (retag==0) %we found a merge
                       if f(bounds(k-1))<f(bounds(k)) %we keep the previous minima
                           indrem(k)=1;
                           stop=1;
                       else %we keep the current minima and need to retag the rest of the curve
                           indrem(k-1)=1;
                           retag=1;
                       end
                   end
                   j=j+1;
               elseif plane(i+1,j+1)==1
                   i=i+1;
                   j=j+1;
               else
                   stop=1;
               end 
            end

            if (j>th) || (j==size(plane,2)-1)
                stop=1;
            end
       end   
    end
    bounds=bounds(indrem==0);
end
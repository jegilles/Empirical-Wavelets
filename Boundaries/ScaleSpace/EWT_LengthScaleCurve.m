function [Length,indices]=EWT_LengthScaleCurve(plane)

%==========================================================================
% function [Length,indices]=EWT_LengthScaleCurve(plane)
%
% This function computes the length of each minima curve in the scale-space 
% representation and returns their initial indices.
%
% Inputs:
%   -plane: scale-space representation
%
% Outputs:
%   -Length: set of minima curve lengths
%   -indices: original indices of each curve
%
% Author: Jerome Gilles
% Institution: UCLA - Department of Mathematics
% Year: 2013
% Version: 1.0
% =========================================================================

Ncurve=0;
for i=1:size(plane,1)
    if plane(i,1)==1
        Ncurve=Ncurve+1;
    end
end
Length=ones(Ncurve,1);
indices=zeros(Ncurve,1);
ic=1;
%tag=plane;

for i=1:size(plane,1)
    if plane(i,1)==1
        indices(ic)=i;
        i0=i;
        j0=2;
        stop=0;
        if i0==1
           while stop==0
               if plane(i0,j0)==1
 %              if (plane(i0,j0)==1) && (tag(i0,j0)==1)
                   Length(ic)=Length(ic)+1;
                   j0=j0+1;
                   if j0>size(plane,2)
                       stop=1;
                   end
               elseif plane(i0+1,j0)==1
%               elseif (plane(i0+1,j0)==1) && (tag(i0+1,j0)==1)
                   Length(ic)=Length(ic)+1;
                   j0=j0+1;
                   i0=i0+1;
                   if (i0==1) || (j0>size(plane,2))
                       stop=1;
                   end
               else
                   stop=1;
               end
            end
            ic=ic+1;

        elseif i0==size(plane,1)
           while stop==0
               if plane(i0,j0)==1
%               if (plane(i0,j0)==1) && (tag(i0,j0)==1)
                   Length(ic)=Length(ic)+1;
                   j0=j0+1;
                   if j0>size(plane,2)
                       stop=1;
                   end
               elseif plane(i0-1,j0)==1
%               elseif (plane(i0-1,j0)==1) && (tag(i0-1,j0)==1)
                   Length(ic)=Length(ic)+1;
                   j0=j0+1;
                   i0=i0-1;
                   if (i0==size(plane,1)) || (j0>size(plane,2))
                       stop=1;
                   end
               else
                   stop=1;
               end
           end
           ic=ic+1;

        else
            while stop==0
               if plane(i0,j0)==1
%               if (plane(i0,j0)==1) && (tag(i0,j0)==1)
                   Length(ic)=Length(ic)+1;
                   j0=j0+1;
                   if j0>size(plane,2)
                       stop=1;
                   end
               elseif plane(i0-1,j0)==1
%               elseif (plane(i0-1,j0)==1) && (tag(i0-1,j0)==1)
                   Length(ic)=Length(ic)+1;
                   j0=j0+1;
                   i0=i0-1;
                   if (i0==1) || (j0>size(plane,2))
                       stop=1;
                   end
               elseif plane(i0+1,j0)==1
%               elseif (plane(i0+1,j0)==1) && (tag(i0+1,j0)==1)
                   Length(ic)=Length(ic)+1;
                   j0=j0+1;
                   i0=i0+1;
                   if (i0==size(plane,1)) || (j0>size(plane,2))
                       stop=1;
                   end
               else
                   stop=1;
               end
            end
            ic=ic+1;
        end
    end
end
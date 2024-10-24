function tag = EWT_LocalMax2D(in)

%===============================================================================
% function tag = EWT_LocalMax2D(in)
%
% This function detects the local maxima in the image (i.e pixels that are
% larger than their surroounding pixels). It manages the images edges and corner 
% properly by considering the appropriate neighboorhood. If a plateau is
% detected, it returns the pixel those coordinates are the mean coordinates of
% the plateau's pixels.
%
% Input:
%   -in: input image
%
% Output:
%   -tag: image of same size as the input those pixels are 1 where a local
%   minima is detected, 0 otherwise.
%
% Author: Jerome Gilles
% Institution: SDSU - Department of Mathematics
% Year: 2022
% Version: 1.0
%===============================================================================

tag=zeros(size(in));

% check center
for i=2:size(in,1)-1
    for j=2:size(in,2)-1
        if (in(i,j)>in(i-1,j-1)) && (in(i,j)>in(i-1,j)) && (in(i,j)>in(i-1,j+1)) ...
                && (in(i,j)>in(i,j-1)) && (in(i,j)>in(i,j+1)) ...
                && (in(i,j)>in(i+1,j-1)) && (in(i,j)>in(i+1,j)) && (in(i,j)>in(i+1,j+1))
            tag(i,j) = 1;
        end
    end
end

% check vertical edges
for i=2:size(in,1)-1
   if (in(i,1)>in(i-1,1)) && (in(i,1)>in(i-1,2)) && (in(i,1)>in(i,2)) ...
           && (in(i,1)>in(i+1,1)) && (in(i,1)>in(i+1,2))
       tag(i,1) = 1;
   end
   
   if (in(i,end)>in(i-1,end)) && (in(i,end)>in(i-1,end-1)) && (in(i,end)>in(i,end-1)) ...
           && (in(i,end)>in(i+1,end)) && (in(i,end)>in(i+1,end-1))
       tag(i,end) = 1;
   end
end

% check horizontal edges
for j=2:size(in,2)-1
   if (in(1,j)>in(1,j-1)) && (in(1,j)>in(2,j-1)) && (in(1,j)>in(2,j)) ...
           && (in(1,j)>in(1,j+1)) && (in(1,j)>in(2,j+1))
       tag(1,j) = 1;
   end
   
   if (in(end,j)>in(end,j-1)) && (in(end,j)>in(end-1,j-1)) && (in(end,j)>in(end-1,j)) ...
           && (in(end,j)>in(end,j+1)) && (in(end,j)>in(end-1,j+1))
       tag(end,j) = 1;
   end
end

% check four corners
if (in(1,1)>in(1,2)) && (in(1,1)>in(2,1)) && (in(1,1)>in(2,2))
    tag(1,1) = 1;
end

if (in(1,end)>in(1,end-1)) && (in(1,end)>in(2,end)) && (in(1,end)>in(2,end-1))
    tag(1,end) = 1;
end

if (in(end,1)>in(end-1,1)) && (in(end,1)>in(end-1,2)) && (in(end,1)>in(end,2))
    tag(end,1) = 1;
end

if (in(end,end)>in(end-1,end)) && (in(end,end)>in(end,end-1)) && (in(end,end)>in(end-1,end-1))
    tag(end,end) = 1;
end

% detect the plateaux
tagg = bwlabel(tag);
tag = zeros(size(in));

for i=1:max(tagg(:))
   [pixx,pixy] = find(tagg==i);
   if length(pixx)>1  %this is a plateau we keep mean center coordinates
        centx = round(mean(pixx));
        centy = round(mean(pixy));
        tag(centx,centy) = 1;
   else % not a plateau
       tag(pixx,pixy) = 1;
   end
end
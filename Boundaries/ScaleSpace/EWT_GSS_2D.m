function [worm_bound,Th,plane]=EWT_GSS_2D(f,params)
%==========================================================================
% function [worm_bound,Th]=EWT_GSS_2D(f,params)
%
% This function builds a scale-space representation of f via successive
% separable convolutions with a discrete Gaussian Kernel. The function 
% automatically estimates how many scales are needed. It then extract the
% meaningful maxima by selected the long-living worms in the scale-space.
%
% Inputs:
%   -f: the input image
%   -params: various parameters tied to the empirical wavelet transform.
%       Notably, the kernel size, number of iterations, scale-space 
%       step-size, and method for classifying meaningful maxima.
%
% Outputs:
% -worm_bound: collection of the coordinates of the significant worms
% -Th: Detected threshold
%
% Author: Jerome Gilles & Jouie Ames & Basile Hurat
% Institution: SDSU - Department of Mathematics
% Year: 2020
% Version: 1.0
% =========================================================================

if nargin==1    %Default params
    params.t = 0.8;
    params.n = 6;
    params.niter = 4;
    params.edges = 0;
end
t=params.t;
n=params.n;

%% Build the kernel
ker=exp(-t)*besseli(-n:n,t,1);
ker=ker'*ker;
ker = ker/sum(ker(:));

Niter=params.niter*ceil(min(size(f,1),size(f,2))/n); %Number of scale iterations

plane = zeros(size(f,1),size(f,2),Niter+1); %initialize the scale-space domain


%% minimum detection on the original image 
plane(:,:,1) = EWT_LocalMax2D(f);

%% We now go through the scale-space
for i=1:Niter
    %% Perform the filtering
    f = im_fourier_extend(f,n);
    f = imfilter(f,ker,'conv');
    f = im_extract(f,n);
    
    plane(:,:,i+1) = EWT_LocalMax2D(f);
end

% to plot the Scale-space
% [rowInd,colInd,zInd]=ind2sub(size(plane),find(plane));
% scatter3(colInd, rowInd, zInd,'.');
% scatter3(colInd, rowInd, zInd,'.','MarkerEdgeColor',[0 0 0]); % black/white

%% Detecting "worms" and their lifetimes
plane_index = bwconncomp(plane); %detect the worm as connected component objects
n_worm = plane_index.NumObjects; %Number of worms
worm_life = zeros(n_worm,1); %vector to store lifetime of each worm
worm_init_loc = zeros(n_worm,2); %vector to store initial worm location
worm_end_loc = zeros(n_worm,2); %vector to store final worm location

for iii = 1:n_worm
    worm = plane_index.PixelIdxList{iii}; %extract worm iii
    [X, Y, Z] = ind2sub(size(plane),worm); %translate into pixel coordinates
    
    % Filter out useless worm roots when several are present
    [X,Y,Z] = select_worm_root(X,Y,Z,floor(size(f)/2)+1);

    %If a worm is not starting at the bottom, find which worm it is closest
    %to and replace that worm's length with the height of this worm.
    if Z(1) > 1 
        
        candidates = find(worm_life == Z(1) - 1);
        if isempty(candidates)
            continue
        end
        cand_dist = pdist2(worm_end_loc(candidates,:),[X(1) Y(1)]);
        ind = candidates(cand_dist == min(cand_dist));
        worm_life(ind(1)) = Z(end);
        worm_end_loc(ind(1),:) = [X(end), Y(end)]; 
    else
        worm_life(iii) = Z(end); %extract the lifetime
        worm_init_loc(iii,:) = [X(1), Y(1)]; %extract initial worm location
        worm_end_loc(iii,:) = [X(end), Y(end)];
    end
end


switch lower(params.typeDetect)
    case 'otsu'
        [worm_bound,Th,~]= EWT_OtsuMethod_2D(worm_life, worm_init_loc);
    case 'kmeans'
        [worm_bound,Th]= EWT_kmeansDetect_2D(worm_life, worm_init_loc);
    case 'halfnormal'
        [worm_bound,Th]= EWT_HalfNormalLaw_2D(worm_life, worm_init_loc,size(plane,3));
    case 'empiricallaw'
        [worm_bound,Th]= EWT_EmpiricalLaw_2D(worm_life, worm_init_loc);
end
end


function [Xf,Yf,Zf] = select_worm_root(X,Y,Z,cc)
% this function selects the correct worm origin (especially when the worm is the
% result of several original worms that merge at some point.

% find indice where Z=1
    rootidx = find(Z==1);

    if length(rootidx)>1 % worm has several roots
        % we select the root the closest to the center

        % compute distances from the center
        dist = (X(rootidx)-cc(2)).^2+(Y(rootidx)-cc(1)).^2;

        % find the index of the smallest distance
        [~,idxmin] = min(dist);

        % keep only the indices to be removes
        rootidx(idxmin)=[];

        % remove all useless roots
        Xf = X; Xf(rootidx) = [];
        Yf = Y; Yf(rootidx) = [];
        Zf = Z; Zf(rootidx) = [];
    else
        Xf = X;
        Yf = Y;
        Zf = Z;
    end
end
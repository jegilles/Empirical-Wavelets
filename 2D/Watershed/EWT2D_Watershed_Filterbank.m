function mfb = EWT2D_Watershed_Filterbank(L,num_regions, extH, extW, tau)
%==========================================================================
%function mfb = EWT2D_Watershed_Filterbank(L,num_regions, extH, extW, tau)
%
% This function constructs 2D empirical wavelet filterbank from an
% arbitrary partitioning of the Fourier spectrum.
%   
% Input:
%   -L: a label image corresponding to the partitioning of the spectrum.
%   Each pixel of spectrum should have an integer corresponding to a region
%   number.
%   -num_regions: array including number of regions
%   -extH: flag to say wether original image's HEIGHT was even and had to
%   be extended to get symmetric spectrum
%   -extW: flag to say wether original image's WIDTH was even and had to
%   be extended to get symmetric spectrum
%   -tau: transition width in (0,2*pi)
% Output:
%   -mfb: constructed filterbank
%
% Author: Basile Hurat - Jerome Gilles
% Institution: SDSU - Department of Mathematics
% Year: 2022
% Version: 1.0
%==========================================================================

%1. Construct distance function from partition
dist=cell(length(num_regions),1); 
for i=1:length(num_regions)
    dist{i}=2*pi*(-bwdist(L==num_regions(i),'quasi-euclidean')...
        + bwdist(1-(L==num_regions(i)),'quasi-euclidean'))/size(L,1); 
end

%2. Construct filters from distance
mfb=cell(size(dist));
for i=1:length(mfb)     
    mfb{i}=zeros(size(dist{i}));
    mfb{i}(dist{i}>tau) = 1;
    mfb{i}(abs(dist{i}) <= tau) = cos(pi*EWT_beta_2D((tau - dist{i}(abs(dist{i}) <= tau))/(2*tau))/2);
end

%3. if image was extended because of even dimensions, readjust
if extH==1
    for i = 1:length(mfb)
        mfb{i} = mfb{i}(1:end-1,:);
    end
end
if extW == 1
    for i = 1:length(mfb)
        mfb{i} = mfb{i}(:,1:end-1);
    end
end
[H,W] = size(mfb{1});

%4. if image was extended because of even dimensions, fix symmetry
if extH == 1
    s=zeros(1,W);
    if mod(W,2) == 0
        for i = 1:length(mfb)
            mfb{i}(1,floor(W/2)+2:end) = mfb{i}(1,floor(W/2)+2:end) + mfb{i}(1,floor(W/2):-1:2);
            mfb{i}(1,2:floor(W/2)) = mfb{i}(1,end:-1:floor(W/2)+2);
            s = s + mfb{i}(1,:).^2;
        end %normalization to preserve the frame condition
        for i = 1:length(mfb)
            for j = 2:floor(W/2)
                if s(1,j) > 0
                    mfb{i}(1,j)=mfb{i}(1,j)./sqrt(s(1,j));
                end
            end
            for j = floor(W/2)+2:size(mfb{i},2)
                if s(1,j) > 0
                    mfb{i}(1,j)=mfb{i}(1,j)./sqrt(s(1,j));
                end
            end
        end
    else
        for i = 1:length(mfb)
            mfb{i}(1,floor(W/2)+2:end) = mfb{i}(1,floor(W/2)+2:end) + mfb{i}(1,floor(W/2):-1:1);
            mfb{i}(1,1:floor(W/2)) = mfb{i}(1,end:-1:floor(W/2)+2);
            s = s + mfb{i}(1,:).^2;
        end
        for i = 1:length(mfb)
           for j = floor(W/2)+2:size(mfb{i},2)
                if s(1,j) > 0
                    mfb{i}(1,j)=mfb{i}(1,j)./sqrt(s(1,j));
                end
           end
           for j = 1:floor(W/2)
                if s(1,j) > 0
                    mfb{i}(1,j)=mfb{i}(1,j)./sqrt(s(1,j));
                end
           end
        end
    end
end

if extW==1
    s=zeros(H,1);
    if mod(H,2) == 0
        for i = 1:length(mfb)
            mfb{i}(floor(H/2)+2:end,1)=mfb{i}(floor(H/2)+2:end,1)+mfb{i}(floor(H/2):-1:2,1);
            mfb{i}(2:floor(H/2),1)=mfb{i}(end:-1:floor(W/2)+2,1);
            s=s+mfb{i}(:,1).^2;
        end
       %normalization to preserve the tight frame condition
       for i=1:length(mfb)
            for j = 2:floor(H/2)
                if s(j,1) > 0
                    mfb{i}(j,1)=mfb{i}(j,1)./sqrt(s(j,1));
                end
            end
            for j = floor(H/2)+2:size(mfb{i},1)
                if s(j,1) > 0
                    mfb{i}(j,1)=mfb{i}(j,1)./sqrt(s(j,1));
                end
            end
       end
    else
        for i = 1:length(mfb)
            mfb{i}(floor(H/2)+2:end,1)=mfb{i}(floor(H/2)+2:end,1)+mfb{i}(floor(H/2):-1:1,1);
            mfb{i}(1:floor(H/2),1)=mfb{i}(end:-1:floor(H/2)+2,1);
            s=s+mfb{i}(:,1).^2;
        end
        for i = 1:length(mfb)
            for j = 1:floor(H/2)
                if s(j,1) > 0
                    mfb{i}(j,1)=mfb{i}(j,1)./sqrt(s(j,1));
                end
            end

            for j =floor(H/2)+2:size(mfb{i},1) 
                if s(j,1) > 0
                    mfb{i}(j,1)=mfb{i}(j,1)./sqrt(s(j,1));
                end
            end
        end
    end
end
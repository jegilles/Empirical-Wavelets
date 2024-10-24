function mfb = EWT2D_Voronoi_Filterbank(sizeImg,vorocells,tau,extH,extW)

%==========================================================================
% function mfb = EWT2D_Voronoi_Filterbank(f,vorocells,tau,extH,extW)
%
% This function builds the filter bank based on the detected Voronoi cells
%   
% Input:
%   -f:  
%   -tau: transition width
%
% Output:
%   -y: value of the filter.
%
% Author: Basile Hurat - Jerome Gilles
% Institution: SDSU - Department of Mathematics
% Year: 2022
% Version: 1.0
%==========================================================================


dist=cell(length(vorocells),1);
for l=1:length(vorocells)
    dist{l}=2*pi*(-bwdist(vorocells{l})+bwdist(1-vorocells{l}))/norm(size(vorocells{l}),2);
end


mfb=cell(size(dist));
for l=1:length(mfb)
    mfb{l}=zeros(sizeImg);
    for i=1:sizeImg(1)
        for j=1:sizeImg(2)
            mfb{l}(i,j)=EWT2D_Voronoi_LP_function(dist{l}(i,j),tau); %build the filters
        end
    end
    
    % manage the expected symmetry if needed
    if extH==1
        mfb{l}(1,:)=0.5*(mfb{l}(1,:)+mfb{l}(1,end:-1:1));
    end
    if extW==1
        mfb{l}(:,1)=0.5*(mfb{l}(:,1)+mfb{l}(end:-1:1,1));
    end
    
    % remove the extensions
    mfb{l} = UnSymmetrize_Fourier(mfb{l},extH,extW);
end
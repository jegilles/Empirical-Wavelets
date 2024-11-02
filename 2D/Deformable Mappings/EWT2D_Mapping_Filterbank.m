function [mfb,phi_Omega,phi] = EWT2D_Mapping_Filterbank(supports,wavelet,mapping,extH,extW,cpx)

%==========================================================================
%[mfb,phi_Omega,phi] = EWT2D_Mapping_Filterbank(L,num_regions, extH, extW, tau)
%
% This function constructs 2D empirical wavelet filterbank from an
% arbitrary partitioning of the Fourier spectrum using mappings.
%   
% Input:
%   -supports: a label image corresponding to the partitioning of the spectrum.
%   Each pixel of spectrum should have an integer corresponding to a region
%   number.
%   -wavelet: structure containing
%        - 'Lambda ': mother wavelet's Fourier support
%        - 'mother_wavelet': image containing the mother wavelet
%        - 'normalization': 'none' or 'det' for the wavelet normalization
%   -mapping: structure containing  
%        - 'demons': structure containing the demons parameters
%              - 'algo': 'Thirion' or 'Vercauteren' for the type of demons
%              - 'level': number of levels in the multiresolution Demons scheme
%              - 'diffus': diffusion regularization parameter in the demons
%              - 'noise': noise parameter in the Vercautern's demons
%              - 'update': 'additive' or 'compositive' for the Vercauteren's demons update
%              - 'field': 'identity' or 'exponential' for the field in the Vercauteren's demons
%        - 'Omega': Fourier supports if no field 'demons'
%        - 'phi': estimated mappings if no field 'demons'
%        - 'phi_Omega': mapped set Lambda if no field 'demons'
%        - 'interp_function': 'imwarp or 'iminterpolate' if no field 'demons'
%        - 'interp_type': 'linear' or 'spline' if no field demons
%   -extH: flag to say wether original image's HEIGHT was even and had to
%   be extended to get symmetric spectrum
%   -extW: flag to say wether original image's WIDTH was even and had to
%   be extended to get symmetric spectrum
%   -cpx: 1 for complex functions
%
% Output:
%   -mfb: constructed filterbank
%   -phi_Omega: images containing the estimated partition
%   -phi: displacement field which corresponds to the estimated mapping
%
% Author: Charles-Gerard Lucas - Jerome Gilles
% Institution: SDSU - Department of Mathematics
% Year: 2024
% Version: 1.0
%==========================================================================

%1. Get mappings
if ~isfield(mapping,'phi'), mapping = EWT2D_demons_mappings(wavelet.Lambda,supports,mapping.demons); end
phi = mapping.phi; Omega = mapping.Omega; phi_Omega = mapping.phi_Omega; interp_function = mapping.interp_function; interp_type = mapping.interp_type;

%2. Computes filters from reference partition
nb_filters = length(phi);
A=cell(1,nb_filters);
switch wavelet.normalization
    case 'none'
        for k=1:nb_filters, A{k} = 1; end
        Sq = 1;
    case 'det'
        for k=1:nb_filters, A{k} = sqrt(abs(jacobian(phi{k}(:,:,1),phi{k}(:,:,2)))); end
        Sq = sqrt(2);
end

wavelet_phi=cell(1,nb_filters); wavelet_phi_sym=cell(1,nb_filters);
mother_wavelet=wavelet.mother_wavelet;
switch interp_function
    case 'imwarp'
        for k=1:nb_filters
            wavelet_phi{k} = imwarp(mother_wavelet,phi{k}); 
            wavelet_phi_sym{k} = imwarp(flip(flip(mother_wavelet,1),2),phi{k}); 
        end
    case 'iminterpolate'
        for k=1:nb_filters
            wavelet_phi{k} = iminterpolate(mother_wavelet,phi{k}(:,:,1),phi{k}(:,:,2),interp_type);
            wavelet_phi_sym{k} = iminterpolate(flip(flip(mother_wavelet,1),2),phi{k}(:,:,1),phi{k}(:,:,2),interp_type);
        end
end

mfb=cell(1,nb_filters);
if cpx == 0
    [h,w] = size(supports); Mi=floor(h/2)+1; Mj=floor(w/2)+1;
    for k=1:nb_filters
        if Omega{k}(Mi,Mj)
            mfb{k}=A{k}.*wavelet_phi{k};
        else
            mfb{k}=A{k}/Sq.*wavelet_phi{k}+flip(flip(A{k},1),2)/Sq.*flip(flip(wavelet_phi_sym{k},1),2);
        end
    end
else
    for k=1:nb_filters
        mfb{k}=A{k}.*wavelet_phi{k};
    end
end   

%3. if image was extended because of even dimensions, readjust
for l=1:nb_filters
    % manage the expected symmetry if needed
    if isreal(ifft(mother_wavelet))
        if extH==1
            mfb{l}(1,:)=0.5*(mfb{l}(1,:)+mfb{l}(1,end:-1:1));
        end
        if extW==1
            mfb{l}(:,1)=0.5*(mfb{l}(:,1)+mfb{l}(end:-1:1,1));
        end
    end
    
    % remove the extensions
    mfb{l} = UnSymmetrize_Fourier(mfb{l},extH,extW);
end
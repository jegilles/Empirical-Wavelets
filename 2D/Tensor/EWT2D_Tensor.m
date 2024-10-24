function [ewtC,mfbR,mfbC,BR,BC] = EWT2D_Tensor(im,params)

%=========================================================================
%
% function [ewtC,mfbR,mfbC,BR,BC] = EWT2D_Tensor(im,params)
%
% Performs the 2D Empirical Wavelet transform based on the tensor product 
% approach.
%
% Inputs:
%   -im: input image
%   -params: structure containing the following parameters:
%       -params.log: 0 or 1 to indicate if we want to work with
%                    the log of the ff
%       -params.preproc: 'none','plaw','poly','morpho,'tophat'
%       -params.reg: 'none','gaussian','closing'
%       -params.lengthFilter: size of the filters used in the above
%                             regularization methods
%       -params.sigmaFilter: standard deviation for the above Gaussian
%                            regularization
%       -params.detect: 'locmax','locmaxmin','locmaxminf','ftc','adaptive',
%                       'adaptivereg','scalespace'
%       -params.N: maximum number of supports (needed for the
%                  locmax and locmaxmin methods)
%       -params.degree: degree of the polynomial (needed for the
%                       polynomial approximation preprocessing)
%       -params.completion: 0 or 1 to indicate if we try to complete
%                           or not the number of modes if the detection
%                           find a lower number of mode than params.N
%       -params.InitBounds: vector of initial bounds (in index domain)
%                           needed for the adaptive and adaptivereg methods
%       -params.typeDetect: (for scalespace method only) 'otsu',
%                           'halfnormal','empiricallaw','mean','kmeans'
%
% Outputs:
%   -ewtC: structure containing all subband images
%   -mfbR: filter bank used for the rows
%   -mfbC: filter bank used for the columns
%   -BR: detected boundaries used for the rows
%   -BC: detected boundaries used for the columns
%
% Author: Jerome Gilles
% Institution: UCLA - Department of Mathematics
% Year: 2013
% Version: 1.0
% ========================================================================

% ================
%    ROWS FIRST
% ================
% We take the 1D FFT for each row
fftim=fft(im');
% We compute the FFT magnitude mean
meanfft=sum(abs(fftim),2)/size(fftim,2);

% We detect the boundaries
boundaries = EWT_Boundaries_Detect(meanfft(1:round(length(meanfft)/2)),params);
BR = boundaries*2*pi/length(meanfft);

% We build the corresponding filter bank
mfbR=EWT_LP_FilterBank(BR,length(meanfft),0);

% We filter the signal to extract each subband
ewtR=cell(length(mfbR),1);
for k=1:length(mfbR)
    filter=repmat(conj(mfbR{k}),1,size(fftim,2));
    ewtR{k}=real(ifft(filter.*fftim))';
end

% ================
%   COLUMNS NEXT
% ================
% We take the 1D FFT for each column
fftim=fft(im);

% We compute the FFT magnitude mean
meanfft=sum(abs(fftim),2)/size(fftim,2);

% We detect the boundaries
boundaries = EWT_Boundaries_Detect(meanfft(1:round(length(meanfft)/2)),params);
BC = boundaries*2*pi/length(meanfft);

% We build the corresponding filter bank
mfbC=EWT_LP_FilterBank(BC,length(meanfft),0);

% We filter the signal to extract each subband
ewtC=cell(length(mfbR),length(mfbC));
for k=1:length(mfbC)
    filter=repmat(conj(mfbC{k}),1,size(fftim,2));
    for s=1:length(mfbR)
        fftim=fft(ewtR{s});
        ewtC{s,k}=real(ifft(filter.*fftim));
    end
end
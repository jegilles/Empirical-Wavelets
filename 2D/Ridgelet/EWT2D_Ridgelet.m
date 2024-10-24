function [ewtLP,mfb,boundaries]=EWT2D_Ridgelet(f,params)
%==========================================================================
% function [ewtLP,mfb,boundaries]=EWT2D_Ridgelet(f,NbScale)
%
% This function performs the 2D Empirical Ridgelet. The Fourier 
% boundaries are detected using the Pseudo-Polar FFT. 
%
% TO RUN THIS FUNCTION YOU NEED TO HAVE THE MATLAB POLARLAB TOOLBOX OF 
% MICHAEL ELAD: http://www.cs.technion.ac.il/~elad/Various/PolarLab.zip
%
% Input:
%   -f: input image
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
% Output:
%   -ewtLP: cell containing each filtered output subband
%   -mfb: cell containing the set of detected filters in the Fourier domain
%   -boundaries: list of the detected boundaries
%
% Author: Jerome Gilles - Giang Tran
% Institution: UCLA - Department of Mathematics
% Year: 2013
% Version: 1.0
%==========================================================================

% Pseudo Polar FFT of f
PseudoFFT=PPFFT(f);

% Compute the mean spectrum with respect to the angle
meanppfft=fftshift(sum(abs(PseudoFFT),2));

% Detect the boundaries
boundaries = EWT_Boundaries_Detect(meanppfft(1:round(length(meanppfft)/2)),params);
boundaries = boundaries*pi/round(length(meanppfft)/2);

% We build the corresponding filter bank
mfb=EWT_LP_FilterBank(boundaries,size(PseudoFFT,1),0);

% We filter each columns of the PseudoPolar FFT to extract each subband
ewtLP=cell(length(mfb),1);
for k=1:length(mfb)
    mfb{k}=repmat(mfb{k},1,size(PseudoFFT,2));
    ewtLP{k}=real(ifft(conj(mfb{k}).*fftshift(PseudoFFT,1)));
end
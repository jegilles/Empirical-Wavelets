function [ewtLP,mfb,boundaries]=EWT2D_LittlewoodPaley(f,params)

%==========================================================================
% function [ewtLP,mfb,boundaries]=EWT2D_LittlewoodPaley(f,params)
%
% This function performs the 2D Littlewood-Paley EWT. The Fourier 
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

W=size(f,2);
H=size(f,1);

% Build the 2D filter bank
mfb=EWT2D_LP_FilterBank(boundaries,W,H);

% We filter the signal to extract each subband
ff=fft2(f);

ewtLP=cell(length(mfb),1);
for k=1:length(mfb)
    ewtLP{k}=real(ifft2(conj(mfb{k}).*ff));
end
function [ewtc,mfb,Bw,Bt] = EWT2D_Curvelet(f,params)

%==========================================================================
% function [ewtc,mfb,Bw,Bt] = EWT2D_Curvelet(f,params)
%
% This function performs the Empirical Curvelet Transform. The Fourier 
% boundaries and angles are detected using the Pseudo-Polar FFT. 
%
% TO RUN THIS FUNCTION YOU NEED TO HAVE THE MATLAB POLARLAB TOOLBOX OF 
% MICHAEL ELAD: http://www.cs.technion.ac.il/~elad/Various/PolarLab.zip
%
% Input:
%   -f: input image (MUST BE SQUARE)
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
%       -params.curvpreproc: 'none','plaw','poly','morpho','tophat'
%       -params.curvmethod: 'locmax','locmaxmin','scalespace'
%       -params.curvreg: 'none','gaussian','closing'
%       -params.curvlengthFilter: size of the filters used in the above
%                             regularization methods
%       -params.curvsigmaFilter: standard deviation for the above Gaussian
%                            regularization
%       -params.curvN: maximum number of angles
%       -params.curvdegree: degree of the polynomial (needed for the
%                       polynomial approximation preprocessing)
%       -params.option: 1 = scales and angles are independants
%                       2 = scales first then angles per scales
%                       3 = angles first then scales per angles
%
% Output:
%   -ewtc: cell containing each filtered output subband (ewtc{1} is the 
%   lowpass subband the next ewtc{s}{t} are the bandpass filtered images, s
%   corresponds to the scales and t to the direction)
%   -mfb: cell containing the set of empirical filters in the Fourier 
%   domain (the indexation is the same as ewtc above)
%   -Bw: list of the detected scale boundaries
%   -Bt: list of the detected angle boundaries
%
% Author: Jerome Gilles
% Institution: SDSU - Department of Mathematics and Statistics
% Version: 1.0 - 2013
% Version: 2.0 - 2015
%==========================================================================

W=size(f,2);
H=size(f,1);

% Pseudo Polar FFT of f
PseudoFFT=PPFFT(f);

switch params.option
    case 1 
        %% Option 1 - scales and angles independent
        % Compute the mean spectrum with respect to the angle to find the scales
        meanppfft=fftshift(sum(abs(PseudoFFT),2));

        % Detect the boundaries
        boundaries = EWT_Boundaries_Detect(meanppfft(1:round(length(meanppfft)/2)),params);
        Bw = boundaries*pi/round(length(meanppfft)/2);
        
        % Compute the mean spectrum with respect to the magnitude frequency to find
        % the angles
        meanppfft=sum(abs(PseudoFFT),1);

        % Detect the boundaries
        boundaries = EWT_Angles_Detect(meanppfft',params);
        Bt = (boundaries-1)*pi/length(meanppfft)-3*pi/4;
                  
        % Build the filter bank
        mfb = EWT2D_Curvelet_FilterBank(Bw,Bt,W,H,params.option);
        % We filter the signal to extract each subband
        ff=fft2(f);

        ewtc=cell(length(mfb),1);
        % We extract the low frequencies first
        ewtc{1}=ifft2(conj(mfb{1}).*ff);
        for s=2:length(mfb)
            ewtc{s}=cell(length(mfb{s}),1);
            for t=1:length(mfb{s})
                ewtc{s}{t}=ifft2(conj(mfb{s}{t}).*ff);
            end
        end
    case 2
        %% Option 2 - scales first and then angles per scale
        % Compute the mean spectrum with respect to the angle to find the scales
        meanppfft=fftshift(sum(abs(PseudoFFT),2));

        % Detect the boundaries
        boundaries = EWT_Boundaries_Detect(meanppfft(1:round(length(meanppfft)/2)),params);
        Bw = boundaries*pi/round(length(meanppfft)/2);
        
        % We compute the angles limits for each scale
        Bt = cell(length(Bw),1);
        for s=1:length(Bw)-1
            % Compute the mean spectrum with respect to the magnitude frequency 
            % annuli to find the angles
            meanppfft=sum(abs(PseudoFFT(boundaries(s):boundaries(s+1),:)),1);
            % Detect the boundaries
            bounds = EWT_Angles_Detect(meanppfft',params);
            Bt{s} = (bounds-1)*pi/length(meanppfft)-3*pi/4;
        end
        meanppfft=sum(abs(PseudoFFT(boundaries(end):end,:)),1);
        % Detect the boundaries
        bounds = EWT_Angles_Detect(meanppfft',params);        
        Bt{length(Bw)} = (bounds-1)*pi/length(meanppfft)-3*pi/4;

        % Build the filter bank
        mfb = EWT2D_Curvelet_FilterBank(Bw,Bt,W,H,params.option);

        % We filter the signal to extract each subband
        ff=fft2(f);

        ewtc=cell(length(mfb),1);
        % We extract the low frequencies first
        ewtc{1}=ifft2(conj(mfb{1}).*ff);
        for s=2:length(mfb)
            ewtc{s}=cell(length(mfb{s}),1);
            for t=1:length(mfb{s})
                ewtc{s}{t}=ifft2(conj(mfb{s}{t}).*ff);
            end
        end
    case 3
        % Compute the mean spectrum with respect to the angle to find the first scale
        meanppfft=fftshift(sum(abs(PseudoFFT),2));
        LL=round(length(meanppfft)/2);

        % Detect the first scale boundary
        boundariesW = EWT_Boundaries_Detect(meanppfft(round(length(meanppfft)/2):end),params);
        Bw1 = boundariesW(1)*pi/LL;
        
        % Compute the mean spectrum with respect to the magnitude frequency to find
        % the angles
        meanppfft=sum(abs(PseudoFFT(floor(size(PseudoFFT,1)/2)+boundariesW(1):end,:)),1);

        % Detect the boundaries
        boundaries = EWT_Angles_Detect(meanppfft',params);
        Bt = (boundaries-1)*pi/length(meanppfft)-3*pi/4;
        
        Bw=cell(length(boundaries)+1,1);
        Bw{1}=Bw1;
        % we detect the scales per angular sector
        for t=1:length(boundaries)-1
            %average spectrum on the given angular range
            meanppfft=sum(abs(PseudoFFT((boundariesW(1)+floor(size(PseudoFFT,1)/2)):end,boundaries(t):boundaries(t+1))),2);
            % Detect the boundaries
            bounds = EWT_Boundaries_Detect(meanppfft,params);
            Bw{t+1} = (boundariesW(1)+bounds)*pi/LL;
        end
        %last one
        %average spectrum on the given angular range
        meanppfft=sum(abs(PseudoFFT((boundariesW(1)+floor(size(PseudoFFT,1)/2)):end,boundaries(end):end)),2)+sum(abs(PseudoFFT((boundariesW(1)+floor(size(PseudoFFT,1)/2)):end,1:boundaries(1))),2);
        
        % Detect the boundaries
        params.reg='closing';
        bounds = EWT_Boundaries_Detect(meanppfft,params);
        Bw{end} = (boundariesW(1)+bounds)*pi/LL;
        
        % Build the filter bank
        mfb = EWT2D_Curvelet_FilterBank(Bw,Bt,W,H,params.option);

        % We filter the signal to extract each subband
        ff=fft2(f);
        
        ewtc=cell(length(mfb),1);
        % We extract the low frequencies first
        ewtc{1}=ifft2(conj(mfb{1}).*ff);
        for s=2:length(mfb)
            ewtc{s}=cell(length(mfb{s}),1);
            for t=1:length(mfb{s})
                ewtc{s}{t}=ifft2(conj(mfb{s}{t}).*ff);
            end
        end
 end

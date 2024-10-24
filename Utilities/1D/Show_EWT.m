function Show_EWT(ewt,f,rec)

% ====================================================================
% function Show_EWT(ewt,f,rec)
%
% This function plots the successive filtered components (low scale 
% first and then wavelets scales). The original and
% reconstructed signals are plotted on a different graph.
% If f and rec are provided, it also plot the original and reconstructed
% signals on a separate figure
%
% Inputs:
%   -ewt: EWT components
%   -f: input signal  (OPTIONNAL)
%   -rec: reconstructed signal (OPTIONNAL)
%
% Author: Jerome Gilles
% Institution: SDSU - Department of Mathematics & Statistics
% Year: 2019
% Version: 2.0
% =====================================================================

%These lines plot the EWT components
figure;
x=0:1/length(ewt{1}):(length(ewt{1})-1)/length(ewt{1});
l=1;
%don't plot more than 6 components per figure for clarity reasons
if length(ewt)>6
    lm=6;
else
    lm=length(ewt);
end

%check if the components are either all real, all complex or a mix
countreal=0;
countcpx=0;
for k=1:length(ewt)
   if isreal(ewt{k})
       countreal=countreal+1;
   else
       countcpx=countcpx+1;
   end
end

if countcpx==0
    realcpx=1; %all real
elseif countreal==0
    realcpx=-1; %all complex
else
    realcpx=0; %mix real/complex
end

for k=1:length(ewt)
    if realcpx==1 %when all components are real
        hold on; 
        subplot(lm,1,l); plot(x,ewt{k}); 
        title(['ewt(',num2str(k),')'])
        if mod(k,6) == 0
            figure;
            l=1;
        else
            l=l+1;
        end
    elseif realcpx==-1 %when all components are complex
        hold on; 
        subplot(lm,2,l); plot(x,real(ewt{k})); 
        title(['real part of ewt(',num2str(k),')'])
        subplot(lm,2,l+1);plot(x,imag(ewt{k})); 
        title(['imaginary part of ewt(',num2str(k),')'])
        if mod(k,6) == 0
            figure;
            l=1;
        else
        l=l+2;
        end
    else %when we have a mix of real and complex components
        hold on;
        if isreal(ewt{k})
            subplot(lm,2,l); plot(x,ewt{k}); 
            title(['ewt(',num2str(k),')'])            
        else
            subplot(lm,2,l); plot(x,real(ewt{k})); 
            title(['real part of ewt(',num2str(k),')'])
            subplot(lm,2,l+1);plot(x,imag(ewt{k})); 
            title(['imaginary part of ewt(',num2str(k),')'])
        end
        if mod(k,6) == 0
            figure;
            l=1;
        else
        l=l+2;
        end
    end
end

%These lines plot f and its reconstruction
if nargin>1
    figure;
    if isreal(f)
        subplot(2,1,1);plot(x,f); title('f')
        subplot(2,1,2);plot(x,rec); title('rec')
    else
        subplot(2,2,1);plot(x,real(f)); title('real part of f')
        subplot(2,2,2);plot(x,imag(f)); title('imaginary part of f')
        subplot(2,2,3);plot(x,real(rec)); title('real part of rec')
        subplot(2,2,4);plot(x,imag(rec)); title('imaginary part of rec')
    end
end

function ShowCurveletFilters(mfb)

%=======================================================================
% function ShowCurveletFilters(mfb)
% 
% This function displays the curvelet filter bank (one figure per scale
% or angular sector + the lowpass frequency) + the sum of the square of 
% each filter.
%
% Inputs:
%   -mfb: the curvelet filter banks to display
%
% Author: Jerome Gilles
% Institution: SDSU - Department of Mathematics and Statistics
% Version: 1.0 - 2013
% Version: 2.0 - 2015
%=======================================================================

mfb{1}=fftshift(mfb{1});
sum=mfb{1}.^2;
for s=2:length(mfb)
   for a=1:length(mfb{s})
     	mfb{s}{a}=fftshift(mfb{s}{a});
      	sum=sum+mfb{s}{a}.^2;
   end
end

Show_EWT2D_Curvelet(mfb)
figure;imshow(sqrt(sum),[]);
norm(sum-1,Inf)

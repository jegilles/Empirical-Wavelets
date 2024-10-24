function dualmfb=EWT_Dual_Filterbank(mfb)

% =========================================================================
% function dualmfb=EWT_Dual_Filterbank(mfb)
%
% This function generate the dual filter bank of the filter bank mfb,
% i.e dualmfb{k}=mfb{k}/sum_l |mfb{l}|^2
%
% Input parameters:
%   -mfb: filter bank
%
% Output:
%   -dualmfb: dual filter bank
%
% Author: Jerome Gilles
% Institution: SDSU - Department of Mathematics & Statistics
% Year: 2020
% Version: 1.0
% =========================================================================

sum=abs(mfb{1}).^2;

for k=2:length(mfb)
    sum=sum+abs(mfb{k}).^2;
end

dualmfb=cell(size(mfb));

for k=1:length(mfb)
    dualmfb{k}=mfb{k}./sum;
end
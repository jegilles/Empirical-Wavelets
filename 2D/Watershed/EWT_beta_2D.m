function bm=EWT_beta_2D(x)

% function used in the construction of Meyer's wavelet, working directly in the
% matrix case
if x<0
    bm=0;
elseif x>1
    bm=1;
else
    bm=x.^4.*(35-84*x+70*x.^2-20*x.^3);
end
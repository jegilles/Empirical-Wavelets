function [Lambda,mother_wavelet] = EWT2D_mother_wavelet(f,name,params)

[h,w] = size(f); if rem(size(f,1),2)==0, h = h+1; end; if rem(size(f,2),2)==0, w = w+1; end
[I,J] = ndgrid(1:h,1:w); Mi=floor(h/2)+1; Mj=floor(w/2)+1;
X = (I-Mi)/Mi; Y = (J-Mj)/Mj;
switch name
    case 'Gabor'
        if ~isfield(params,'sigma'); params.sigma = sqrt(1/(2*2.5^2*pi)); end
        Lambda = double(X.^2+Y.^2<=1/4);
        mother_wavelet = exp(-(X.^2+Y.^2)/(2*params.sigma^2));
    case 'Shannon'
        Lambda = double(max(abs(X),abs(Y))<=1/2);
        mother_wavelet = exp(-1i*pi*(X+3/2)).*exp(-1i*pi*(Y+3/2)).*Lambda;
    case 'Disk'
        if nargin<3; params.tau = .2; end
        tau = params.tau;
        normXi = sqrt(X.^2+Y.^2); Lambda = double(normXi<=1/2);
        mother_wavelet = zeros(h,w);
        mother_wavelet(normXi<1/2-tau) = 1;
        mother_wavelet(abs(normXi-1/2)<=tau) = ...
            cos(pi/2*EWT_beta_2D((tau-1/2+normXi(abs(normXi-1/2)<=tau))/(2*tau)));
    case 'Square'
        if nargin<3; params.tau = .2; end
        tau = params.tau;
        Lambda = double(max(abs(X),abs(Y))<=1/2);
        mother_wavelet_x = zeros(h,w); mother_wavelet_y = zeros(h,w);
        mother_wavelet_x(abs(X)<=1/2-tau) = 1; mother_wavelet_y(abs(Y)<=1/2-tau) = 1;
        mother_wavelet_x(abs(abs(X)-1/2)<tau) = cos(pi/2*EWT_beta_2D((tau-1/2+abs(X(abs(abs(X)-1/2)<tau)))/(2*tau)));
        mother_wavelet_y(abs(abs(Y)-1/2)<tau) =  cos(pi/2*EWT_beta_2D((tau-1/2+abs(Y(abs(abs(Y)-1/2)<tau)))/(2*tau)));
        mother_wavelet = mother_wavelet_x.*mother_wavelet_y;
end

end
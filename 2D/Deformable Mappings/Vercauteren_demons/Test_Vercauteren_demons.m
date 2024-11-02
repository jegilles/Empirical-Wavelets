close all
clear all
clc

load 'Lambda.mat'
load 'Omega.mat' % Ten sets
Omegan = double(Omega{7});

%% grid search
tic;

max_level = floor(log2(min(size(Lambda))));
levels = max_level-1:max_level;
noise = 4; 
diffus = .3:.01:.5; 
s1 = length(diffus); s2 = length(noise); s3 = length(levels);
phi_tmp = cell([s1,s2,s3]); phi_Omega_tmp = cell([s1,s2,s3]);
crit = zeros([s1,s2,s3]);
for l1=1:s1
    for l2=1:s2
        for l3=1:s3
            params.sigma_i = 1;
            params.sigma_fluid = 1;
            params.sigma_diffusion = diffus(l1);
            params.sigma_x = noise(l2);
            params.nlevel = levels(l3);
            params.interp_type = 'linear';
            params.do_display = 0; % faster if 0
            params.niter = 300; % can be a list
            params.stop_criterium  = 1e-3;
            params.update = 'additive';
            params.field = 'identity'; % 'identity' for additive and 'exp' for diffeomorphic
            [phi_tmp{l1,l2,l3},phi_Omega_tmp{l1,l2,l3}] = Vercauteren_demons(Lambda,Omegan,params);
            crit(l1,l2,l3) = norm(phi_Omega_tmp{l1,l2,l3}./256-Omegan,'fro')^2/numel(Omegan);
        end
    end
end
[~,loc] = min(crit(:));
[ld,ln,ll] = ind2sub(size(crit),loc);
phi_Omegan = phi_Omega_tmp{ld,ln,ll};

toc;

%% plot deformation
figure(1); clf;
subplot(121)
imshow(Omegan); colormap turbo;
subplot(122)
imshow(phi_Omegan); colormap turbo;

%% plot MSE
linewidth = 1.5; fontsize = 12;
figure(2); clf;
for l3=1:s3
    subplot(1,s3,l3)
    plot(diffus,crit(:,1,l3),'Linewidth',linewidth); grid on;
    pbaspect([1 1 1])
    xlabel('$\sigma_{\rm diffusion}$','Interpreter','Latex','Fontsize',fontsize)
    ylabel('MSE','Interpreter','Latex','Fontsize',fontsize)
    title(['Number of levels: ',num2str(levels(l3))],'Interpreter','Latex')
    set(gca,'Linewidth',linewidth,'Fontsize',fontsize)
end
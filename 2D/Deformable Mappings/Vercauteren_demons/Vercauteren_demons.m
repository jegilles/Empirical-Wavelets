function [phi,Mp] = Vercauteren_demons(M,F,params)

% Demons Registration
%
% Simple matlab code for 2D image registration using the diffeomorphic log-demons algorithm 
% Code is provided in order to help the understanding of the Demons algorithm
%
% Original Algorithm in:
% [1] - Symmetric Log-Domain Diffeomorphic Registration: A Demons-Based Approach
%       Vercauteren, Pennec, Perchant, Ayache -- MICCAI 2008, 754-761
% [2] - Diffeomorphic demons: Efficient non-parametric image registration,
%       Vercauteren, Pennec, Perchant, Ayache -- NeuroImage 2009, (45)1:61-72
%
% For a more recent work/survey, exploiting global shape characteristics 
% (instead of the conventional local gradient-based approaches), consider citing 
%
% [1] - Spectral Log-Demons: Diffeomorphic Image Registration with Very Large Deformations
%       Lombaert, Grady, Pennec, Ayache, Cheriet -- IJCV 2014, (107)3:254-271

%figure(1); clf; colormap gray; set(gcf,'renderer','painter')

%% Parameters
if ~isfield(params,'niter'), params.niter = 500; end
if ~isfield(params,'sigma_fluid'), params.sigma_fluid     = 1; end % regularize update      field
if ~isfield(params,'sigma_diffusion'), params.sigma_diffusion = .2; end% regularize deformation field
if ~isfield(params,'sigma_i'), params.sigma_i         = 1; end % weight on similarity term
if ~isfield(params,'sigma_x'), params.sigma_x         = 2; end % weight on spatial uncertainties (maximal step)
if ~isfield(params,'nlevel'), params.nlevel         = 8;   end % multiresolution
if ~isfield(params,'do_display'), params.do_display      = 0; end   % display iterations
if ~isfield(params,'interp_type'), params.interp_type     = 'linear'; end   % display iterations
if ~isfield(params,'stop_criterium'), params.stop_criterium  = 0.01; end   % stop criterium
if ~isfield(params,'update'), params.update = 'additive'; end   % stop criterium
if ~isfield(params,'field'), params.field = 'exp'; end   % stop criterium

niter = params.niter;
sigma_fluid  = params.sigma_fluid ;
sigma_diffusion = params.sigma_diffusion;
sigma_i = params.sigma_i;
sigma_x = params.sigma_x;
nlevel = params.nlevel;
do_display  = params.do_display;
interp_type = params.interp_type;
stop_criterium = params.stop_criterium;
update = params.update;
field = params.field;

%% normalize data
F = double(F);
F = 256*(F-min(F(:)))/range(F(:));              % normalize intensities
M = double(M); 
M = 256*(M-min(M(:)))/range(M(:));              % normalize intensities

if nlevel == 1
    
    %% Register
    disp(['Register...']);
    opt = struct('niter',niter, 'sigma_fluid',sigma_fluid, 'sigma_diffusion',sigma_diffusion, 'sigma_i',sigma_i, 'sigma_x',sigma_x, 'do_display',do_display, 'do_plotenergy',1,'interp_type',interp_type);
    [Mp,sx,sy,vx,vy] = register(F,M,opt);

else
    
    %% Multiresolution
    vx = zeros(size(M)); % deformation field
    vy = zeros(size(M));
    for k=nlevel:-1:1
        if params.do_display
            disp(['Register level: ' num2str(k) '...']);
        end

        % downsample
        scale = 2^-(k-1);
        Fl = imresize(F,scale);
        Ml = imresize(M,scale);
        vxl = imresize(vx*scale,scale);
        vyl = imresize(vy*scale,scale);

        if length(params.niter)>1
            niter = params.niter(k);
        else
            niter = params.niter;
        end
        %niter = 2^(3+nlevel-k+1);

        % register
        opt = struct('niter',niter,...
            'sigma_fluid',sigma_fluid,...
            'sigma_diffusion',sigma_diffusion,...
            'sigma_i',sigma_i,...
            'sigma_x',sigma_x,...
            'vx',vxl, 'vy',vyl,...
            'do_display',do_display, ...
            'do_plotenergy',1,...
            'interp_type',interp_type, ...
            'stop_criterium',stop_criterium, ...
            'update',update, ...
            'field',field);
        [Mp,sxl,syl,vxl,vyl] = register(Fl,Ml,opt);

        % upsample
        vx = imresize(vxl/scale,size(M));
        vy = imresize(vyl/scale,size(M));
    end

    switch field
        case 'exp'
            [sx,sy] = expfield(vx,vy,opt.interp_type);  % deformation field
        case 'identity'
            sx = vx; sy = vy;
    end
    phi(:,:,1) = sx; phi(:,:,2) = sy;
    
end

end

%% Register two images
function [Mp,sx,sy,vx,vy] = register(F,M,opt)

    if nargin<3;  opt = struct();  end
    if ~isfield(opt,'sigma_fluid');      opt.sigma_fluid     = 1.0;              end
    if ~isfield(opt,'sigma_diffusion');  opt.sigma_diffusion = 1.0;              end
    if ~isfield(opt,'sigma_i');          opt.sigma_i         = 1.0;              end
    if ~isfield(opt,'sigma_x');          opt.sigma_x         = 1.0;              end
    if ~isfield(opt,'niter');            opt.niter           = 250;              end
    if ~isfield(opt,'vx');               opt.vx              = zeros(size(M));   end
    if ~isfield(opt,'vy');               opt.vy              = zeros(size(M));   end
    if ~isfield(opt,'stop_criterium');   opt.stop_criterium  = 0.01;             end
    if ~isfield(opt,'imagepad');         opt.imagepad        = 1.2;              end
    if ~isfield(opt,'do_display');       opt.do_display      = 0;                end
    if ~isfield(opt,'do_plotenergy');    opt.do_plotenergy   = 1;                end
    if ~isfield(opt,'interp_type');      opt.interp_type     = 'linear';         end
    if ~isfield(opt,'uptdate');          opt.update          = 'additive';         end
    if ~isfield(opt,'field');            opt.update          = 'exp';         end

    %% padded image
    [F,lim] = imagepad(F,opt.imagepad);
    [M,lim] = imagepad(M,opt.imagepad);
    
    %% T is the deformation from M to F
    vx = imagepad(opt.vx,opt.imagepad);
    vy = imagepad(opt.vy,opt.imagepad);
    e  = zeros(1,opt.niter);
    e_min = 1e+100;      % Minimal energy
    
    %% Iterate update field
    for iter=1:opt.niter

        % Find update
        [ux,uy] = findupdate(F,M,vx,vy,opt.sigma_i,opt.sigma_x,opt.interp_type,opt.field);

        % Regularize update
        ux    = imgaussian(ux,opt.sigma_fluid);
        uy    = imgaussian(uy,opt.sigma_fluid);

        % Compute step (e.g., max half a pixel)
        step  = opt.sigma_x;
        
        switch opt.update
            case 'additive'
                % Update velocities (demons) - additive
                vx = vx + step*ux;
                vy = vy + step*uy;
            case 'compositive'
                % Update velocities (demons) - composition
                [vx,vy] = compose(vx,vy,step*ux,step*uy,opt.interp_type);
        end
        
        % Regularize velocities
        vx = imgaussian(vx,opt.sigma_diffusion);
        vy = imgaussian(vy,opt.sigma_diffusion);
        
        % Get Transformation
        switch opt.field
            case 'exp'
                [sx,sy] = expfield(vx,vy,opt.interp_type);  % deformation field
            case 'identity'
                sx = vx; sy = vy;
        end

        % Compute energy
        e(iter) = energy(F,M,sx,sy,opt.sigma_i,opt.sigma_x,opt.interp_type);
        if opt.do_display
            disp(['Iteration: ' num2str(iter) ' - ' 'energy: ' num2str(e(iter))]);
        end
        
        if e(iter)<e_min
            sx_min = sx; sy_min = sy; % update best fields
            vx_min = vx; vy_min = vy; % update best fields
            e_min  = e(iter);
        end
        
        % Stop criterium
        if iter>1 && abs(e(iter) - e(max(1,iter-5))) < e(1)*opt.stop_criterium
            break;
        end
        

        if opt.do_display
            % display deformation
            subplot(2,4,7); showvector(ux,uy,4,3,lim); title('Update');
            %subplot(2,4,8); showgrid  (sx,sy,4,lim); title('Transformation');
            drawnow;
            
            % Display registration
            Mp     = iminterpolate(M,sx,sy,opt.interp_type);
            diff   = (F-Mp).^2;
            showimage(F,'Fixed', M,'Moving', Mp,'Warped', diff,'Diff', 'lim',lim,'nbrows',2); drawnow;
           
            % Plot energy
            risk(iter) = norm(F-Mp,'fro')^2/numel(F);
            if opt.do_plotenergy
                subplot(2,2,3)
                hold on;
                plot(1:iter,risk(1:iter),'r-'); xlim([0 opt.niter]);
                xlabel('Iteration'); ylabel('Energy');
                hold off;
                drawnow
            end
        end

    end
    
    %% Get Best Transformation
    vx = vx_min;  vy = vy_min;
    sx = sx_min;  sy = sy_min;
    
    %% Transform moving image
    Mp = iminterpolate(M,sx,sy,opt.interp_type);
    
    %% Unpad image
    Mp = Mp(lim(1):lim(2),lim(3):lim(4));
    vx = vx(lim(1):lim(2),lim(3):lim(4));
    vy = vy(lim(1):lim(2),lim(3):lim(4));
    sx = sx(lim(1):lim(2),lim(3):lim(4));
    sy = sy(lim(1):lim(2),lim(3):lim(4));

end

%% Find update between two images
function [ux,uy] = findupdate(F,M,vx,vy,sigma_i,sigma_x,interp_type,field)

    % Get Transformation
        switch field
            case 'exp'
                [sx,sy] = expfield(vx,vy,interp_type);  % deformation field
            case 'identity'
                sx = vx; sy = vy;
        end

    % Interpolate updated image
    M_prime = iminterpolate(M,sx,sy,interp_type); % intensities at updated points

    % image difference
    diff = F - M_prime;
    
    % moving image gradient
    [gy,gx] = gradient(M_prime);   % image gradient
    normg2  = gx.^2 + gy.^2;       % squared norm of gradient
    area    = size(M,1)*size(M,2); % area of moving image
    
    % update is Idiff / (||J||^2+(Idiff^2)/sigma_x^2) J, with Idiff = F(x)-M(x+s), and J = Grad(M(x+s));
    scale = diff ./ (normg2 + diff.^2*sigma_i^2/sigma_x^2);
    scale(normg2==0) = 0;
    scale(diff  ==0) = 0;
    ux = gx .* scale;
    uy = gy .* scale;
    
    % Zero non overlapping areas
    ux(F==0)       = 0; uy(F==0)       = 0;
    ux(M_prime==0) = 0; uy(M_prime==0) = 0;

end

%% Apply gaussian filter to image
function I = imgaussian(I,sigma)

    if sigma==0; return; end % no smoothing
    
    % Create Gaussian kernel
    radius = ceil(3*sigma);
    [x,y]  = ndgrid(-radius:radius,-radius:radius); % kernel coordinates
    h      = exp(-(x.^2 + y.^2)/(2*sigma^2));
    h      = h / sum(h(:));
    
    % Filter image
    I = imfilter(I,h);

end

%% Exponentiate vector field
%  Changed: Dec 6th, 2011
%
function [vx,vy] = expfield(vx, vy,interp_type)

    % Find n, scaling parameter
    normv2 = vx.^2 + vy.^2;
    m = sqrt(max(normv2(:)));
    n = ceil(log2(m/0.5)); % n big enough so max(v * 2^-n) < 0.5 pixel)
    n = max(n,0);          % avoid null values
    
    % Scale it (so it's close to 0)
    vx = vx * 2^-n;
    vy = vy * 2^-n;

    % square it n times
    for i=1:n
        [vx,vy] = compose(vx,vy, vx,vy,interp_type);
    end

end

%% Compose two vector fields
function [vx,vy] = compose(ax,ay,bx,by,interp_type)

    [x,y] = ndgrid(0:(size(ax,1)-1), 0:(size(ax,2)-1)); % coordinate image
    x_prime = x + ax; % updated x values
    y_prime = y + ay; % updated y values
    
    % Interpolate vector field b at position brought by vector field a
    bxp = interpn(x,y,bx,x_prime,y_prime,interp_type,0); % interpolated bx values at x+a(x)
    byp = interpn(x,y,by,x_prime,y_prime,interp_type,0); % interpolated by values at y+a(y)

    % Compose
    vx = ax + bxp;
    vy = ay + byp;
    
end

%% Jacobian
function det_J = jacobian(sx,sy)

    % Gradients
    [gx_y,gx_x] = gradient(sx);
    [gy_y,gy_x] = gradient(sy);
    
    % Add identity
    gx_x = gx_x + 1;  % zero displacement should yield a transformation T = Identity (points keep their positions)
    gy_y = gy_y + 1;  % adding identity matrix here
    
    % Determinant
    det_J = gx_x.*gy_y - ...
            gy_x.*gx_y;
end


%% Get energy
function e = energy(F,M,sx,sy,sigma_i,sigma_x,interp_type)

    % Intensity difference
    Mp     = iminterpolate(M,sx,sy,interp_type);
    diff2  = (F-Mp).^2;
    area   = size(M,1)*size(M,2);
    
    % Transformation Gradient
    jac = jacobian(sx,sy);
    
    % Three energy components
    e_sim  = sum(diff2(:)) / area;
    %e_dist = sum((cx(:)-sx(:)).^2 + (cy(:)-sy(:)).^2) / area;
    e_reg = sum(jac(:).^2) / area;
    
    % Total energy
    e      = e_sim + (sigma_i^2/sigma_x^2) * e_reg;

end

%% Pad image
function [I,lim] = imagepad(I,scale)

    if nargin<2; scale = 2; end; % default, pad image twice as big
    
    Ip  = zeros(ceil(size(I)*scale));
    lim = bsxfun(@plus, floor(size(I)*(scale-1)/2), [[1 1];size(I)]); % image limits
    Ip(lim(1):lim(2),lim(3):lim(4)) = I;                              % padded image
    I = Ip;

end

%% Display vector field
%  Changed: Dec 6th, 2011
%
function showvector(ux,uy,downsample,scale,lim)

    if nargin<3; downsample = 1; end;
    
    sizex = size(ux,1);
    sizey = size(uy,2);
    
    ux  = ux(1:downsample:end, 1:downsample:end);
    uy  = uy(1:downsample:end, 1:downsample:end);
    
    if nargin<4; scale = 3;                     end; % Scale vector to show small ones
    if nargin<5; lim   = [0 sizex-1 0 sizey-1]; end; % Display whole image

    [y,x] = ndgrid((0:downsample:(sizex-1))+downsample/2, (0:downsample:(sizey-1))+downsample/2); % coordinate image
    quiver(x,y,ux,uy,scale);                  % show vectors
    daspect([1 1 1]);
    axis([lim(3) lim(4) lim(1) lim(2)]);      % which vector to show
    axis off;
    set(gca,'YDir','reverse');
    
end

%% Display two images
%  Changed: Dec 6th, 2011
%
function showimage(varargin)

    % Check parameters
    nb_args   = size(varargin,2);
    nb_images = nb_args;
    nb_rows   = 1;
    row       = 1;
    crange    = [0 1]*256; % default image intensities
    
    for i=1:nb_args
        if ischar(varargin{i})
            if isequal(varargin{i},'lim')
                lim       = varargin{i+1};
                nb_images = nb_images-2;
            elseif isequal(varargin{i},'nbrows')
                nb_rows   = varargin{i+1};
                nb_images = nb_images-2;
            elseif isequal(varargin{i},'row')
                row       = varargin{i+1};
                if row>nb_rows; nb_rows = row; end
                nb_images = nb_images-2;
            elseif isequal(varargin{i},'caxis')
                crange    = varargin{i+1};
                nb_images = nb_images-2;
            else
                nb_images = nb_images-1;
            end
        end
    end
    
    % Display images
    iter_image = 1;
    for iter_arg=1:nb_args
        if ~ischar(varargin{iter_arg})
            I = varargin{iter_arg};
            subplot(nb_rows,nb_images,(row-1)*nb_images + iter_image);
            imagesc(I,crange);
            daspect([1 1 1]);
            if exist('lim'); axis([lim(3) lim(4) lim(1) lim(2)]); end
            axis off;
            if iter_arg+1<=nb_args && ischar(varargin{iter_arg+1})
                title(varargin{iter_arg+1});
            end
            iter_image = iter_image+1;
        end
        if iter_image>nb_images
            break;
        end
    end
        
end

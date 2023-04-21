function w = get_fov_weights(x, sigma, rad, neigh_rad)


% build a 3D-matrix of patches
patches = foveated_patches(x, rad);

% initialize
[Nr Nc] = size(x);
w = zeros(Nr, Nc, 2*neigh_rad+1, 2*neigh_rad+1);
w_max = realmin * ones(Nr,Nc);

% scroll the search area
for i = -neigh_rad : neigh_rad
    for j = -neigh_rad : neigh_rad
        
        % skip the "central" block (of the search area)
        if i == 0 && j == 0
            continue;
        end
        
        % translate the image
        r_idx = max(1,1-i) : min(Nr,Nr-i);
        c_idx = max(1,1-j) : min(Nc,Nc-j);
        p_shift = patches;
        p_shift(r_idx,c_idx,:) = patches(r_idx+i,c_idx+j,:);
        
        % compute the distances
        D = sum( (p_shift-patches).^2, 3 );
        D(D==0) = Inf;
        
        % compute the weights
        p = exp( -D / sigma^2 );
        
        % update the max
        w_max = max(p, w_max);
        
        % save the weigths
        ip = i + (neigh_rad+1);
        jp = j + (neigh_rad+1);
        w(:,:,ip,jp) = p;
    end
end

% handle the "central" block (of the search area)
ip = 0 + (neigh_rad+1);
jp = 0 + (neigh_rad+1);
w(:,:,ip,jp) = w_max + ~w_max;  % the term "~w_max" is 1 only when w_max == 0 (otherwise is 0).

% reshape in 3D
w = reshape(w, [Nr Nc (2*neigh_rad+1)^2]);





function BIGfov = foveated_patches(z, rad)

% build the kernels
[v_kernels k_kernel] = make_kernels(rad);

% get the size
[Nr Nc] = size(z);
Nk = numel(k_kernel);
Nv = length(v_kernels);

% compute the blurred images
zBlurred = zeros([Nr+2*rad Nc+2*rad Nv]);
for jj = 1:Nv
    v_rad = (length(v_kernels{jj}) - 1) / 2;
    z_ext = padarray(z, [rad rad] + v_rad, 'symmetric');
    zBlurred(:,:,jj) = conv2(z_ext, v_kernels{jj}, 'valid');
end

% get the set of values assumed by the k_kernel (sorted in ascending order)
K = unique(k_kernel);   
K = K(end:-1:1);

% build a 3D matrix with all foveated patches
BIGfov = zeros(Nr, Nc, Nk);
idx = 1;
for c = -rad:rad
    for r = -rad:rad
        
        % select the blurred image
        jj = find(K == k_kernel(r+rad+1,c+rad+1));
        zb = zBlurred( rad+r+(1:Nr), rad+c+(1:Nc), jj);
        
        % translate
        BIGfov(:,:,idx) = zb;
        idx = idx+1;
    end
end





function [v_kernels k_kernel] = make_kernels(rad)

% build the "k" kernel
[u2,u1] = meshgrid( -rad:rad, -rad:rad );
linfdistance = max( abs(u1), abs(u2) );
K_tilde = cumsum( (2*(rad:-1:0)+1).^-2 );
K_tilde = K_tilde( [rad,rad:-1:1] );
k_kernel = interp1(0:rad,K_tilde, linfdistance,'linear', 'extrap');
k_kernel = k_kernel / sum(k_kernel(:));

% get the set of values assumed by the k_kernel (sorted in ascending order)
K = unique(k_kernel);   
K = K(end:-1:1);

% kernel params
p = 1-exp(-2*pi);                           % PSF mass over disk of unit radius
sigma = sqrt(1/(-2*log(1-p))*max(K)./K);    % std.dev. of the blurring kernels (here max(K) = \kappa_0)
widths = 2*ceil(3*sigma)+1;                 % support widths of the blurring kernels

% build the "blurring" kernels (used by the foveation operator)
v_kernels = cell( size(K) );
jj=0;
for kappa = K'
    jj=jj+1;
    v_kernels{jj} = (sqrt(pi*kappa)*2*sigma(jj)) * fspecial('gaussian', widths(jj), sigma(jj));
end
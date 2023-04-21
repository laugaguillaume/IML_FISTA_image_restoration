function [X,param,Ar,Ac,A1,A2,A3] = create_data(choice_image,choice_blur,param)

switch choice_image
    case 'JWST 256 gray'
        I = im2double(imread('SMACS0723_JWST.jpeg')); targetSize = [2048 2048];
        r = centerCropWindow2d(size(I),targetSize);
        X = imcrop(I,r); X = X(1:8:end,1:8:end,:);
        X = rgb2gray(X);
    case 'JWST 256'
        I = im2double(imread('SMACS0723_JWST.jpeg')); targetSize = [2048 2048];
        r = centerCropWindow2d(size(I),targetSize);
        X = imcrop(I,r); X = X(1:8:end,1:8:end,:);
    case 'JWST 512'
        I = im2double(imread('SMACS0723_JWST.jpeg')); targetSize = [2048 2048];
        r = centerCropWindow2d(size(I),targetSize);
        X = imcrop(I,r); X = X(1:4:end,1:4:end,:);
    case 'JWST 1024'
        I = im2double(imread('SMACS0723_JWST.jpeg')); targetSize = [2048 2048];
        r = centerCropWindow2d(size(I),targetSize);
        X = imcrop(I,r); X = X(1:2:end,1:2:end,:);
    case 'JWST 2048'
        I = im2double(imread('SMACS0723_JWST.jpeg')); targetSize = [2048 2048];
        r = centerCropWindow2d(size(I),targetSize);
        X = imcrop(I,r);
    case 'IR Pillars 2048'
        I1 = im2double(imread('NIR-MIR-Pillars.jpg'));
        X = I1;
        I1 = [X ; flip(X,1)];
        x0 = 0 ; y0 = 0; tsize = 2048;
        r = [x0 y0 tsize tsize];
        X = imcrop(I1,r); X = double(X);
        % Credits: SCIENCE: NASA, ESA, CSA, STScI 
        % https://webbtelescope.org/contents/media/images/01GK2KKTR81SGYF24YBGYG7TAP.html
    case 'Jupiter 1024'
        I = im2double(imread('Jupiter_JWST.png')); targetSize = [1024 1024];
        r = centerCropWindow2d(size(I),targetSize);
        X = imcrop(I,r);
    case 'Jupiter 1024 gray'
        I = im2double(imread('Jupiter_JWST.png')); targetSize = [1024 1024];
        r = centerCropWindow2d(size(I),targetSize);
        X = imcrop(I,r); X = rgb2gray(X);
    case 'ImageNet 2' 
        I = im2double(imread('ILSVRC2012_test_00000164.JPEG')); targetSize = [512 512];
        r = centerCropWindow2d(size(I),targetSize);
        X = imcrop(I,r);
end
%rng(1,'twister');
%%
[param.N,param.M,Q] = size(X);
X=double(X);
X = X./max(max(X)); % normalization 


% PSF Blur parameters
if param.N == 4096
    param.dim = [80, 80]; % dimension of the PSF array
    param.s = 14.5455; % standard deviation of the Gaussian blur
elseif param.N == 2048
    param.dim = [20, 20]; % dimension of the PSF array
    param.s = 3.6364;
elseif param.N == 1024
    % For big images
    param.dim = [20, 20]; % dimension of the PSF array
    param.s = 3.6364; % standard deviation of the Gaussian blur
elseif param.N == 512
    param.dim = [10,10]; % dimension of the PSF array
    param.s = 1.8182; % standard deviation of the Gaussian blur
elseif param.N==256 % 51.200
    % For small images
    param.dim = [5,5]; % dimension of the PSF array [5,5]
    param.s = 0.9091; % standard deviation of the Gaussian blur
end
 
param.BC = 'reflexive';% Blur boundary conditions % periodic in HNO = circular in imfilter

% PSF defined functions using HNO Package
% Generate Gaussian PSF parameter to blur image
Acolor = [1 0 0 ; 0 1 0 ; 0 0 1]; % cross channel blurring
switch choice_blur
    case 'Gauss small'
        if param.N == 2048
            param.dim = [20, 20]; % dimension of the PSF array
            param.s = 3.6364; % standard deviation of the Gaussian blur
        elseif param.N == 512
            param.dim = [5,5]; % dimension of the PSF array [5,5]
            param.s = 0.9091; 
        end
        [PSF,center]=psfGauss(param.dim,param.s);
        if Q == 3
            Z1 = imfilter(X(:,:,1),PSF,'symmetric');
            Z2 = imfilter(X(:,:,2),PSF,'symmetric');
            Z3 = imfilter(X(:,:,3),PSF,'symmetric');
            
            Z1 = Acolor(1)*Z1 + Acolor(2)*Z2 + Acolor(3)*Z3;
            Z2 = Acolor(4)*Z1 + Acolor(5)*Z2 + Acolor(6)*Z3;
            Z3 = Acolor(7)*Z1 + Acolor(8)*Z2 + Acolor(9)*Z3;
            
            Z = cat(3,Z1,Z2,Z3);    
        elseif Q == 1
            Z = imfilter(X,PSF,'symmetric');
        end
        Pbig         = padPSF(PSF,[param.N,param.M]);
        [Ar,Ac      ]= kronDecomp(Pbig,center,param.BC);
        Ar(Ar<1e-10) = 0; Ac(Ac<1e-10) = 0;
        Ar = sparse(Ar); Ac = sparse(Ac);
    case 'Gauss big'
        if param.N == 2048
            param.dim = [40, 40]; % dimension of the PSF array
            param.s = 7.2727; % standard deviation of the Gaussian blur
        elseif param.N == 512
            param.dim = [10,10]; % dimension of the PSF array
            param.s = 1.8182; 
        end
        [PSF,center]=psfGauss(param.dim,param.s);
        if Q == 3
            Z1 = imfilter(X(:,:,1),PSF,'symmetric');
            Z2 = imfilter(X(:,:,2),PSF,'symmetric');
            Z3 = imfilter(X(:,:,3),PSF,'symmetric');
            
            Z1 = Acolor(1)*Z1 + Acolor(2)*Z2 + Acolor(3)*Z3;
            Z2 = Acolor(4)*Z1 + Acolor(5)*Z2 + Acolor(6)*Z3;
            Z3 = Acolor(7)*Z1 + Acolor(8)*Z2 + Acolor(9)*Z3;
            
            Z = cat(3,Z1,Z2,Z3);    
        elseif Q == 1
            Z = imfilter(X,PSF,'symmetric');
        end
        Pbig         = padPSF(PSF,[param.N,param.M]);
        [Ar,Ac      ]= kronDecomp(Pbig,center,param.BC);
        Ar(Ar<1e-10) = 0; Ac(Ac<1e-10) = 0;
        Ar = sparse(Ar); Ac = sparse(Ac);
    case 'Inpainting Blur'
        param.dim = [3,3];
        param.s = 1;
        [PSF,center]=psfGauss(param.dim,param.s);
        if Q == 3
            Z1 = imfilter(X(:,:,1),PSF,'symmetric');
            Z2 = imfilter(X(:,:,2),PSF,'symmetric');
            Z3 = imfilter(X(:,:,3),PSF,'symmetric');
            A1 = rand(param.N,param.N)<0.8;
            A2 = rand(param.N,param.N)<0.8;
            A3 = rand(param.N,param.N)<0.8;
            Z1 = A1.*Z1; Z2 = A2.*Z2; Z3 = A3.*Z3;
            Z = cat(3,Z1,Z2,Z3);    
        elseif Q == 1
            A1 = rand(param.N,param.N)<0.8;
            Z = imfilter(X,PSF,'symmetric');
            Z = A1.*Z;
        end
        Pbig         = padPSF(PSF,[param.N,param.M]);
        [Ar,Ac      ]= kronDecomp(Pbig,center,param.BC);
        Ar(Ar<1e-10) = 0; Ac(Ac<1e-10) = 0;
        Ar = sparse(Ar); Ac = sparse(Ac);
    case 'Inpainting 50'
        if Q == 3
            A1 = rand(param.N,param.N)<0.5;
            A2 = rand(param.N,param.N)<0.5;
            A3 = rand(param.N,param.N)<0.5;
            Z1 = A1.*X(:,:,1); Z2 = A2.*X(:,:,2); Z3 = A3.*X(:,:,3);
            A1 = sparse(A1); A2 = sparse(A2); A3= sparse(A3);
            Z = cat(3,Z1,Z2,Z3);    
        elseif Q == 1
            A1 = sparse(rand(param.N,param.N)<0.5);
            Z = A1.*X;
        end
    case 'Inpainting 90'
        if Q == 3
            A1 = rand(param.N,param.N)<0.1;
            A2 = rand(param.N,param.N)<0.1;
            A3 = rand(param.N,param.N)<0.1;
            Z1 = A1.*X(:,:,1); Z2 = A2.*X(:,:,2); Z3 = A3.*X(:,:,3);
            A1 = sparse(A1); A2 = sparse(A2); A3= sparse(A3);
            Z = cat(3,Z1,Z2,Z3);    
        elseif Q == 1
            A1 = sparse(rand(param.N,param.N)<0.2);
            Z = A1.*X;
        end
end
% Config blur for optimization
switch choice_blur
    case 'Gauss small'
        A1 = [];
        A2 = [];
        A3 = [];
    case 'Gauss big'
        A1 = [];
        A2 = [];
        A3 = [];
    case 'Inpainting Blur'
        Ablur = kron(Ar,Ac);
        if Q==3
            A1 = sparse(A1(:).*Ablur);
            A2 = sparse(A2(:).*Ablur);
            A3 = sparse(A3(:).*Ablur);
        elseif Q == 1
            A1 = sparse(A1(:).*Ablur);
            A2 = [];
            A3 = [];
        end
    case 'Inpainting 50'
        Ar = []; Ac = []; 
        if Q == 1
            A2 = [];
            A3 = [];
        end
    case 'Inpainting 90'
        Ar = []; Ac = []; 
        if Q == 1
            A2 = [];
            A3 = [];
        end
end
        
Z = Z + param.std*randn(size(Z));
param.Z = Z;
switch choice_blur
    case 'Inpainting 50'
        param.PSF = [];
        param.center = [];
        param.Acolor = [];
    case 'Inpainting 90'
        param.PSF = [];
        param.center = [];
        param.Acolor = [];
    otherwise
        param.PSF = PSF;
        param.center = center;
        param.Acolor = Acolor;
end






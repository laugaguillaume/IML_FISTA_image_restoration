%% Cope with Octave/matlab differences.
if is_octave
    pkg load fits;
    fitsread = @(filename) read_fits_image(filename);
end

global mu mtf wgt dat;

%% Add source directory to the path.
srcdir = '../src/';
addpath(srcdir);

%% Read the data.
datadir = '../../data/saturn/';
dat = fitsread([datadir, 'saturn.fits']);
psf = fitsread([datadir, 'saturn_psf.fits']);
%colormap('viridis');
plimg(dat, 1); title('Raw data');
plimg(psf, 2); title('PSF');

%% Compute the Modulation Transfer Function (MTF) as the FFT of the PSF after
%% proper zero-padding and centering.
mtf = fft2(fftshift(zeropad(psf, size(dat))));

%% Regularization level.
mu = 0.01;

%% Assume 70% of bad pixels by default.
if ~exist("fraction_of_bad_pixels")
    fraction_of_bad_pixels = 0.7;
end

%% Array of pixelwise weights and corresponding operator.
wgt = ones(size(dat));
if fraction_of_bad_pixels > 0
    wgt(rand(size(wgt)) < fraction_of_bad_pixels) = 0;
end

%% Functions to apply the convolution by the PSF operator and its adjoint.
H = @(x) real(ifft2(mtf.*fft2(x)));
Ht = @(x) real(ifft2(conj(mtf).*fft2(x)));
W = @(x) wgt.*x;
lhs = @(x) Ht(W(H(x))) + mu*apply_DtD(x);

%% RHS "vector" of the normal equations.  *MUST* be recomputed whenever
%% the weights change.
rhs = Ht(W(dat));

%% Show good pixels.
plimg(W(dat), 3); title('Good pixels in data');

%% Iterative deconvolution by the linear conjugate gradients (function name
%% "lhs" or handle @lhs both work).
[x1, status] = optm_conjgrad(lhs, rhs, [], 'maxiter', 50, 'verb', 1);
fprintf('# %s\n\n', optm_reason(status));
plimg(x1, 4); title('Result of deconvolution by linear conjugate gradient');



%% Iterative deconvolution by a quasi-Newton method (without bounds).
[x2, f, g, status] = optm_vmlmb(@fg, zeros(size(dat)), 'fmin', 0, ...
                                'maxiter', 50, 'verb', 1);
plimg(x2, 5); title('Result of deconvolution by variable metric method');
fprintf('# %s\n\n', optm_reason(status));

%% Iterative deconvolution by a quasi-Newton method (with lower bound).
[x3, f, g, status] = optm_vmlmb(@fg, zeros(size(dat)), 'fmin', 0, 'lower', 0, ...
                                'maxiter', 50, 'verb', 1);
fprintf('# %s\n\n', optm_reason(status));
plimg(x3, 6); title('Result of deconvolution by variable metric method and positivity constraints');

%% Iterative deconvolution by a quasi-Newton method (with lower and upper bounds).
[x4, f, g, status] = optm_vmlmb(@fg, zeros(size(dat)), 'fmin', 0, 'lower', 0, 'upper', 1e3, ...
                                'maxiter', 50, 'verb', 1);
fprintf('# %s\n\n', optm_reason(status));
plimg(x4, 7); title('Result of deconvolution by variable metric method and lower and upper bounds');

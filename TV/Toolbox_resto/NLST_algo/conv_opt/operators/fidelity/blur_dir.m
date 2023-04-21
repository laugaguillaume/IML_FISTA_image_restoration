function y = blur_dir(x, psf, valid, border)

% convolution
y = conv_dir(x, psf, border);

% decimation
y = select_dir(y, valid);
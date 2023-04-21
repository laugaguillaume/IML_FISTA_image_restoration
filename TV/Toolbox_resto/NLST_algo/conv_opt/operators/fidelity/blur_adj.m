function y = blur_adj(x, psf, valid, border)

% adjoint selection
yt = select_adj(x, valid);

% adjoint convolution
y = conv_adj(yt, psf, border);
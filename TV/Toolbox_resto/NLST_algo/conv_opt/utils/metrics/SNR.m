function snr = SNR(x_true, x_est)

snr = 20 * log10( norm( x_true(:) ) / norm( x_est(:)-x_true(:) ) );
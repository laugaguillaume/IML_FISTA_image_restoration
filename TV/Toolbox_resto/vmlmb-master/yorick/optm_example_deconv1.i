// optm_example_deconv1.i -
//
// Bound constrained regularized deconvolution of image with VMLMB in Yorick.
//
// In this example, simple operators are built to implement the various linear
// mappings.  All objective functions are quadratic.
//
//-----------------------------------------------------------------------------

// Load the data and the code.
SRCDIR = dirname(current_include()) + "/";
include, SRCDIR+"invpb.i"; // before "saturn_data.i"
include, SRCDIR+"saturn_data.i";
include, SRCDIR+"optm.i";

// Compute the Modulation Transfer Function (MTF) as the FFT of the PSF after
// proper zero-padding and centering and define functions to apply the
// convolution by the PSF operator and its adjoint.
func H__1(x) { return ifft(mtf*fft(x)); }
func Ht_1(x) { return ifft(conj(mtf)*fft(x)); }
func H__2(x) { return FFT(mtf*FFT(x, 0), 2); }
func Ht_2(x) { return FFT(conj(mtf)*FFT(x, 0), 2); }
if (is_func(xfft_new)) {
    // Use FFTW.
    FFT = xfft_new(dims=dims, real=1n, planning=XFFT_ESTIMATE);
    mtf = FFT(fftshift(pad_array(psf, dims)));
    H  = H__2;
    Ht = Ht_2;
} else {
    // Use Yorick's fft.
    mtf = fft(fftshift(pad_array(psf, dims)));
    H  = H__1;
    Ht = Ht_1;
}

// Weighting operator.
func W(x) { return wgt*x; }

// Regularization level.
mu = 0.01;

// Function to apply the LHS "matrix" of the normal equations.
func A(x) { return Ht(W(H(x))) + mu*DtD(x); }

// RHS "vector" of the normal equations.  *MUST* be recomputed whenever
// the weights change.
b = Ht(W(dat));

// Initial solution.
x0 = array(double, dims);

// Iterative deconvolution by the linear conjugate gradients.
local status;
x1 = optm_conjgrad(A, b, x0, status, maxiter=50, verb=1);
if (!batch()) {
    plimg, x1, cmin=0, fig=4,
        title="Result of deconvolution by linear conjugate gradient";
}
printf, "\n";

// Function to compute the objective function and its gradient.
func fg(x, &g)
{
    r = H(x) - dat;
    Wr = W(r);
    muDtDx = mu*DtD(x);
    g = Ht(Wr) + muDtDx;
    g += g; // cheap 2×
    return optm_inner(r, Wr) + optm_inner(x, muDtDx);
}

// Iterative deconvolution by a quasi-Newton method (without bounds).
local f, g, status;
x2 = optm_vmlmb(fg, x0, f, g, status,
                fmin=0, maxiter=50, verb=1);
if (!batch()) {
    plimg, x2, cmin=0, fig=5,
        title=("Result of deconvolution by\n" +
               "variable metric method");
}
printf, "\n";

// Iterative deconvolution by a quasi-Newton method (with lower bound).
x3 = optm_vmlmb(fg, x0, f, g, status, lower=0,
                fmin=0, maxiter=50, verb=1);
if (!batch()) {
    plimg, x3, cmin=0, fig=6,
        title=("Result of deconvolution by\n" +
               "variable metric method and positivity constraints");
}
printf, "\n";

// Iterative deconvolution by a quasi-Newton method (with lower and upper bounds).
x4 = optm_vmlmb(fg, x0, f, g, status, lower=0, upper=1e3,
                fmin=0, maxiter=50, verb=1);
if (!batch()) {
    plimg, x4, cmin=0, fig=7,
        title=("Result of deconvolution by\n" +
               "variable metric method and lower and upper bounds");
}
printf, "\n";

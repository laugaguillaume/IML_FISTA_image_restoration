// optm_example_deconv2.i -
//
// Bound constrained regularized deconvolution of image with VMLMB in Yorick.
//
// In this example, operators are built using the tools provded by "invpb.i".
// All objective functions are quadratic.
//
//-----------------------------------------------------------------------------

// Load the data and the code.
SRCDIR = dirname(current_include()) + "/";
include, SRCDIR+"invpb.i"; // before "saturn_data.i"
include, SRCDIR+"saturn_data.i";
include, SRCDIR+"optm.i";

// Direct model and precision matrix.
H = convolution_operator(fftshift(pad_array(psf, dims)));
W = scaling_operator(wgt);

// Likelihood term.
lkl = quadratic_objective_function(y=dat, H=H, W=W);

// Quadratic regularization.
rgl = quadratic_objective_function(W=DtD);

// Regularization level.
mu = 0.01;

// Composite objective function: objfun = lkl + mu⋅rgl
objfun = composite_objective_function(lkl, mu, rgl);
fg = always_provide_gradient(objfun);

// LHS matrix and RHS vector of the normal equations.
A = lhs_of_normal_equations(H=H, W=W, mu=mu, R=DtD);
b = rhs_of_normal_equations(A, dat);

// Initial solution.
x0 = array(double, dims);

// Iterative deconvolution by the linear conjugate gradients (function name
// "lhs" or handle @lhs both work).
local status;
x1 = optm_conjgrad(A, b, x0, status, maxiter=50, verb=1);
if (!batch()) {
    plimg, x1, cmin=0, fig=4,
        title="Result of deconvolution by linear conjugate gradient";
}
printf, "\n";

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

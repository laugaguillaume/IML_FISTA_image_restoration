// optm_test_deconv.i -
//
// Bound constrained regularized deconvolution of image with VMLMB in Yorick.
//
// In this example, operators are built using the tools provided by "invpb.i".
// Likelihood is quadratic, regularization is quadratic or hyperbolic.
//
//-----------------------------------------------------------------------------

// Load the data and the code.
SRCDIR = dirname(current_include()) + "/";
include, SRCDIR+"invpb.i"; // before "saturn_data.i"
include, SRCDIR+"saturn_data.i";

// Options to use.
if (is_void(method)) {
    method = "optimpack";
}
if (is_void(regul)) {
    regul = "quadratic";
}

if (batch()) {
    argv = get_argv();
    argc = numberof(argv);
    for (argi = 2; argi <= argc; ++argi) {
        if (argv(argi) == "--method" && argi < argc) {
            method = argv(++argi);
            continue;
        }
        if (argv(argi) == "--regul" && argi < argc) {
            regul = argv(++argi);
            continue;
        }
        error, "unknown option or missing value";
    }
}

if (method == "lbfgsb") {
    include, "lbfgsb.i";
} else if (method == "optimpack") {
    include, "opky.i";
} else if (method == "optimpacklegacy") {
    include, "optimpacklegacy.i";
} else if (method == "optimpackrosetta") {
    include, SRCDIR+"optm.i";
} else {
    error, "invalid value for --method";
}

// Direct model and precision matrix.
H = convolution_operator(fftshift(pad_array(psf, dims)));
W = scaling_operator(wgt);

// Likelihood term.
lkl = quadratic_objective_function(y=dat, H=H, W=W);

// Regularization and regularization parameters.
if (regul == "quadratic") {
    mu = 0.01;
    rgl = quadratic_objective_function(W=DtD);
} else if (regul == "hyperbolic") {
    mu = 0.2;
    eps = 6.0;
    rgl = isotropic_edge_preserving(eps=eps);
} else {
    error, "option --regul must be \"quadratic\" or \"hyperbolic\"";
}

// Composite objective function: objfun = lkl + mu⋅rgl
objfun = composite_objective_function(lkl, mu, rgl);
fg = always_provide_gradient(objfun);

// Initial solution.
x0 = array(double, dims);

#if 0
// LHS matrix and RHS vector of the normal equations.
A = lhs_of_normal_equations(H=H, W=W, mu=mu, R=DtD);
b = rhs_of_normal_equations(A, dat);

// Iterative deconvolution by the linear conjugate gradients (function name
// "lhs" or handle @lhs both work).
local status;
x1 = optm_conjgrad(A, b, x0, status, maxiter=200, maxeval=1000, verb=1);
if (!batch()) {
    plimg, x1, cmin=0, fig=4,
        title="Result of deconvolution by linear conjugate gradient";
}
printf, "\n";
#endif

// Convergence setting.
ftol = 0;
gtol = [1e-3, 0];
xtol = 0;
maxiter =  200;
maxeval = 1000;

// Iterative deconvolution by a quasi-Newton method (without bounds).
local f, g, status;
if (method == "lbfgsb") {
    x2 = lbfgsb(fg, x0, f, g, mem=5,
                ftol=ftol, gtol=gtol, xtol=xtol,
                maxiter=maxiter, maxeval=maxeval, verb=1);
} else if (method == "optimpacklegacy") {
    x2 = opl_vmlmb(fg, x0, f, g, mem=5, fmin=0,
                   maxiter=maxiter, maxeval=maxeval, verb=1);
} else if (method == "optimpack") {
    x0_ = x0;
    x2 = opk_minimize(fg, x0_, f, g, mem=5, fmin=0,
                      maxiter=maxiter, maxeval=maxeval, verb=1);
} else if (method == "optimpackrosetta") {
    x2 = optm_vmlmb(fg, x0, f, g, status, mem=5, fmin=0,
                    ftol=ftol, gtol=gtol, xtol=xtol,
                    maxiter=maxiter, maxeval=maxeval, verb=1);
}
if (!batch()) {
    plimg, x2, cmin=0, fig=5,
        title=("Result of deconvolution by\n" +
               "variable metric method");
}
printf, "\n";

// Iterative deconvolution by a quasi-Newton method (with lower bound).
if (method == "lbfgsb") {
    x3 = lbfgsb(fg, x0, f, g, mem=5, lower=0,
                ftol=ftol, gtol=gtol, xtol=xtol,
                maxiter=maxiter, maxeval=maxeval, verb=1);
} else if (method == "optimpacklegacy") {
    x3 = opl_vmlmb(fg, x0, f, g, mem=5, fmin=0, xmin=0,
                   maxiter=maxiter, maxeval=maxeval, verb=1);
} else if (method == "optimpack") {
    x0_ = x0;
    x3 = opk_minimize(fg, x0_, f, g, mem=5, fmin=0, lower=0,
                      maxiter=maxiter, maxeval=maxeval, verb=1);
} else if (method == "optimpackrosetta") {
    x3 = optm_vmlmb(fg, x0, f, g, status, mem=5, fmin=0, lower=0,
                    ftol=ftol, gtol=gtol, xtol=xtol,
                    maxiter=maxiter, maxeval=maxeval, verb=1);
}
if (!batch()) {
    plimg, x3, cmin=0, fig=6,
        title=("Result of deconvolution by\n" +
               "variable metric method and positivity constraints");
}
printf, "\n";

// Iterative deconvolution by a quasi-Newton method (with lower and upper bounds).
if (method == "lbfgsb") {
    x4 = lbfgsb(fg, x0, f, g, mem=5, lower=0, upper=1e3,
                ftol=ftol, gtol=gtol, xtol=xtol,
                maxiter=maxiter, maxeval=maxeval, verb=1);
} else if (method == "optimpacklegacy") {
    x4 = opl_vmlmb(fg, x0, f, g, mem=5, fmin=0, xmin=0, xmax=1e3,
                   maxiter=maxiter, maxeval=maxeval, verb=1);
} else if (method == "optimpack") {
    x0_ = x0;
    x4 = opk_minimize(fg, x0_, f, g, mem=5, fmin=0, lower=0, upper=1e3,
                      maxiter=maxiter, maxeval=maxeval, verb=1);
} else if (method == "optimpackrosetta") {
    x4 = optm_vmlmb(fg, x0, f, g, status, mem=5, fmin=0, lower=0, upper=1e3,
                    ftol=ftol, gtol=gtol, xtol=xtol,
                    maxiter=maxiter, maxeval=maxeval, verb=1);
}
if (!batch()) {
    plimg, x4, cmin=0, fig=7,
        title=("Result of deconvolution by\n" +
               "variable metric method and lower and upper bounds");
}
printf, "\n";

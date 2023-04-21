/*
 * invpb.i --
 *
 * Useful functions for inverse problems in Yorick.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (c) 2015-2021, Éric Thiébaut <eric.thiebaut@univ-lyon1.fr>
 * All rights reserved.
 */

if (is_func(xfft_new) == 3) {
    // Function xfft_new is an autoload object, make sure to include "xfft.i"
    // to have constants defined.
    include, "xfft.i", 3;
    if (is_func(xfft_new) == 3) {
        error, "cannot load XFFT plugin";
    }
}

//-----------------------------------------------------------------------------
// SIGNAL PROCESSING

local fftshift, ifftshift;
/* DOCUMENT fftshift(x);
         or ifftshift(x);

     Circular-shift of a periodic signal `x` along its dimensions.  Along a
     dimension of length `n`, `fftshift(x)` shifts a signal centered at index 1
     to index `n/2 + 1` while `ifftshift(x)` does the opposite.

   SEE ALSO: fft, convolution_operator.
 */
func fftshift(a)
{
  dims = dimsof(a);
  rank = numberof(dims) - 1;
  if (rank < 1) return a;
  return roll(a, dims(2:0)/2);
}

func ifftshift(a)
{
  dims = dimsof(a);
  rank = numberof(dims) - 1;
  if (rank < 1) return a;
  return roll(a, -(dims(2:0)/2));
}

func ifft(x, setup=, cmplx=)
/* DOCUMENT ifft(x, setup=, cmplx=);

     Compute inverse FFT of `x` assuming a real-valued result unless keyword
     `cmplx` is true.

   SEE ALSO: fft, fft_setup, fftshift.
 */
{
    z = fft(unref(x), -1, setup=setup);
    s = 1.0/numberof(z);
    if (cmplx) return s*unref(z);
    return s*double(unref(z));
}

func fftfreq(dim, step)
/* DOCUMENT k = fftfreq(dim);
         or f = fftfreq(dim, step);

     Generate Discrete Fourier Transform (DFT) frequency indices or
     frequencies.

     With a single argument, the function returns a vector of `dim` values (of
     type `long`) set with the frequency indices:

         k = [0, 1, 2, ..., n-1, -n, ..., -2, -1]   if dim = 2*n
         k = [0, 1, 2, ..., n,   -n, ..., -2, -1]   if dim = 2*n + 1

     depending whether `dim` is even or odd.  These rules are compatible to
     what is assumed by fftshift (which to see) in the sense that:

         fftshift(fftfreq(dim)) = [-n, ..., -2, -1, 0, 1, 2, ...]

     With two arguments, `step` is the sample spacing in the direct space and
     the result is a floating point vector with `dim` elements set with the
     frequency bin centers in cycles per unit of the sample spacing (with zero
     at the start).  For instance, if the sample spacing is in seconds, then
     the frequency unit is cycles/second.  This is equivalent to:

         fftfreq(dim)/(dim*step)

   SEE ALSO: fft, fftshift, indgen.
*/
{
    if (dim <= 2) {
        f = (dim == 1 ? [0] : [0,-1]);
    } else {
        n = dim/2;
        f = grow(indgen(0:dim-n-1),indgen(-n:-1));
    }
    return (is_void(step) ? f : (1.0/(double(dim)*step))*f);
}

func goodfftdim(len)
/* DOCUMENT goodfftdim(len);

     Return the smallest integer which is greater or equal `max(len)` and which
     is a multiple of powers of 2, 3 and/or 5.

   SEE ALSO fftfreq, fft.
*/
{
    len = long(max(len));
    best = 2*len;
    for (i5 = 1; i5 < best; i5 *= 5) {
        for (i3 = i5; i3 < best; i3 *= 3) {
            /* last loop (power of 2) is exited as soon as N >= LEN */
            for (i2 = i3; i2 < len; i2 *= 2)
                ; /* empty loop body */
            if (i2 == len) return len;
            best = min(best, i2);
        }
    }
    return best;
}

func fftplot(a, scale=, legend=, hide=, top=, cmin=, cmax=, stairs=,
             type=, width=, color=, smooth=,
             marks=, marker=, mspace=, mphase=)
/* DOCUMENT fftplot, a;

     Plot 1-D or 2-D FFT array `a`, taking care of "rolling" `a` and setting
     correct world boundaries.  Keyword `scale` can be used to indicate the
     "frequel" scale along both axes (`scale` is a scalar) or along each axis
     (`scale` is a 2-element vector: `scale=[xscale,yscale]`); by default,
     `scale=1` for a 1-D argument and `scale=[1,1]` for a 2-D argument.

   KEYWORDS legend, hide, top, cmin, cmax.

   SEE ALSO pli, plg, fftshift. */
{
    dims = dimsof(a);
    rank = numberof(dims) - 1;
    if (rank < 1 || rank > 2 || identof(a) > Y_DOUBLE) {
        error, "expecting 1D or 2D array of reals";
    }
    if (is_void(scale)) {
        scale = 1.0;
    } else if ((is_scalar(scale) || (is_vector(scale) &&
                                     numberof(scale) == rank)) &&
               identof(scale) <= Y_DOUBLE && noneof(scale == 0)) {
        scale = double(scale);
    } else {
        error, "bad value for the `scale` keyword";
    }

    if (rank == 1) {
        scale = scale(1);
        dim1 = dims(2);
        max1 = dim1 + (min1 = -(dim1/2)) - 1;
        y = fftshift(a);
        x = scale*indgen(min1:max1);
        if (stairs) {
            plh, y, x, legend=legend, hide=hide,
                type=type, width=width, color=color,
                marks=marks, marker=marker, mspace=mspace, mphase=mphase;
        } else {
            plg, y, x, legend=legend, hide=hide,
                type=type, width=width, color=color, smooth=smooth,
                marks=marks, marker=marker, mspace=mspace, mphase=mphase;
        }
    } else {
        scale1 = scale(1);
        scale2 = scale(0);
        dim1 = dims(2);
        dim2 = dims(3);
        max1 = dim1 + (min1 = -(dim1/2)) - 1;
        max2 = dim2 + (min2 = -(dim2/2)) - 1;
        pli, fftshift(bytscl(a, top=top, cmin=cmin, cmax=cmax)),
            scale1*(min1 - 0.5), scale2*(min2 - 0.5),
            scale1*(max1 + 0.5), scale2*(max2 + 0.5),
            legend=legend, hide=hide;
    }
}

//-----------------------------------------------------------------------------
// LINEAR OPERATORS

local linear_operators;
/* DOCUMENT Linear Operators.

     Linear operators in "invpb.i" are Yorick function-like objects such that:

        op(x, adj) -> A⋅x     if adj is false;
                      A'⋅x    else.

     where `op` is the operator implementing the linear mapping `A` and `A⋅x`
     denotes applying `A` to `x` while `A'⋅x` denotes applying the adjoint of
     `A` to `x`.

     A missing argument is set to `[]` in Yorick which is considered as false,
     hence `op(x)` also yields `A⋅x`.

   SEE ALSO: DtD, identity_operator, scaling_operator, convolution_operator,
             generalized_matrix_operator.
 */

func identity(x) { return x; }
func identity_operator(nil)
/* DOCUMENT Id = identity_operator();

     Give the identity operator `Id` such that:

         Id(x, adj) -> x    (whatever adj).


     Note that `identity_operator()` always returns the same function.
     So you can use `op == identity_operator()` to check whether `op` is the identity operator.

     The `identity` function just returns its argument:

         identity(x) -> x

   SEE ALSO: scaling_operator, linear_operators.
 */
{
    if (nil != []) error, "no argument should be provided";
    return _apply_identity;
}
// Protect _apply_identity from being redefined.
__apply_identity__ = _apply_identity;
func _apply_identity(arg, adj) { return arg; }
if (is_func(__apply_identity__) == 1) _apply_identity = __apply_identity__;

func scaling_operator(scl)
/* DOCUMENT op = scaling_operator(scl);

     Build a scaling operator such that:

         op(x, adj) -> scl*x          if adj is false;
                       conj(scl)*x    else.

     If `scl == []` or if `scl == 1` everywhere, the identity operator is
     returned.

   SEE ALSO: identity, linear_operators.
 */
{
    if (is_void(scl) || (min(scl) == 1 && max(scl) == 1)) {
        return _apply_identity;
    } else {
        op = h_new(scl=scl);
        h_evaluator, op, (structof(scl) == complex ?
                          "_apply_complex_scaling" : "_apply_scaling");
        return op;
    }
}
func _apply_scaling(op, arg, adj) { return op.scl*arg; }
func _apply_complex_scaling(op, arg, adj) {
    return (adj ? conj(op.scl) : op.scl)*arg; }

// H = convolution_operator(fftshift(zeropad(psf, dimsof(dat))));
func convolution_operator(psf)
/* DOCUMENT op = convolution_operator(psf);

     Build an operator that implements circulant convolution by the point
     spread function PSF:

         op(x, adj) -> ifft(fft(psf)*fft(x))          if adj is false;
                       ifft(fft(conj(psf))*fft(x))    else.

     If the xfft plugin is available, it is used to compute the fast Fourier
     transform (FFT); otherwise, Swarztrauber's FFT is used.

   SEE ALSO: fft, fft_setup, xfft_new, linear_operators.
 */
{
    if (is_func(xfft_new)) {
        fft_ = xfft_new(dims = dimsof(psf),
                        planning = XFFT_ESTIMATE,
                        real = (structof(psf) != complex));
        op = h_new(mtf = fft_(psf), fft = fft_);
        h_evaluator, op, "_apply_xfft_convolution";
    } else {
        plan = fft_setup(dimsof(psf));
        op = h_new(mtf = fft(psf, +1, setup=plan), plan = plan);
        h_evaluator, op, "_apply_fft_convolution";
    }
    return op;
}

// Convolution by Yorick's FFT.
func _apply_fft_convolution(op, arg, adj)
{
    z = (adj ? conj(op.mtf) : op.mtf)*fft(arg, +1, setup=op.plan);
    fft_inplace, z, -1, setup=op.plan;
    if (structof(arg) == complex) {
        return (1.0/numberof(arg))*unref(z);
    } else {
        return (1.0/numberof(arg))*double(unref(z));
    }
}

// Convolution by XFFT.
func _apply_xfft_convolution(op, arg, adj)
{
    fft = op.fft;
    return fft((adj ? conj(op.mtf) : op.mtf)*fft(arg, 0), 2);
}

func generalized_matrix_operator(A) { return closure("mvmult", A); }
/* DOCUMENT op = generalized_matrix_operator(A);

     Build a scaling operator such that:

         op(x, adj) -> mvmult(A, x, adj)

   SEE ALSO: mvmult, linear_operators.
 */

local lhs_of_normal_equations, rhs_of_normal_equations;
/* DOCUMENT A = lhs_of_normal_equations(H=Id, W=Id, mu=0, R=Id);
        and b = rhs_of_normal_equations(A, y);

     Build a linear operator `A` and compute an array `b` corresponding to the
     left-hand-side "matrix" `A` and to the right-hand-side "vector" `b` of the
     so-called "normal equations":

         A⋅x = b

     where `x` are the unknowns and:

         A = H'⋅W⋅H + µ⋅R
         b = H'⋅W⋅y

     Solving these "normal equations" amounts to minimizing the objective
     function:

         f(x) = ‖H⋅x - y‖²_W + µ⋅‖x‖²_R
              = (H⋅x - y)'⋅W⋅(H⋅x - y) + µ⋅x'⋅R⋅x

     The components of the left-hand-side "matrix" `A` are specified by
     keywords.  Unspecifed operators are assumed to be the identity `Id`.  If
     `mu` is unspecifed, `mu = 0` is assumed.  Linear mappings `W` and `R`
     shall be self-adjoint and nonnegative.

   SEE ALSO: mvmult, linear_operators, quadratic_objective_function.
 */
func lhs_of_normal_equations(nil, H=, W=, mu=, R=)
{
    if (nil != []) error, "all arguments given by keywords";
    Id = identity_operator();
    if (H == []) H = Id;
    if (W == []) W = Id;
    if (R == []) R = Id;
    if (mu == []) mu = 0.0;
    lhs = h_new(H=H, W=W, mu=mu, R=R);
    h_evaluator, lhs, "_apply_lhs_of_normal_equations";
    return lhs;
}
func rhs_of_normal_equations(lhs, y) { return lhs.H(lhs.W(y), 1n); }
func _apply_lhs_of_normal_equations(lhs, x, adj)
{
    z = lhs.H(lhs.W(lhs.H(x)), 1n);
    mu = lhs.mu;
    if (mu == 1) {
        z = unref(z) + lhs.R(x);
    } else if (mu != 0) {
        z = unref(z) + mu*lhs.R(x);
    }
    return z;
}

func DtD(x, adj)
/* DOCUMENT DtD(x);
         or DtD(x, adj);

     Apply `D'⋅D` to `x` `D` is a (multi-dimensional) finite difference
     operator and `D'` is its adjoint (transpose).

     Argument `adj` is to specify whether to apply the adjoint of `D'⋅D` but,
     as `D'⋅D` is self-adjoint, `adj` is ignored.

   SEE ALSO: linear_operators.
 */
{
  /* Create an array to store the result.  Manage to use the same type as the
     `dif` operator. */
  if (is_real(x)) {
    type = structof(x);
  } else if (is_integer(x)) {
    type = long;
  } else {
    error, "bad data type";
  }
  dims = dimsof(x);
  r = array(type, dims);
  rank = numberof(dims) - 1; // number of dimensions (up to 10 in Yorick)
  p = 1:-1;
  q = 2:0;
  if (rank >= 1) {
    if (dims(2) >= 2) {
      dx = x(dif,..);
      r(q,..) += dx;
      r(p,..) -= dx;
    }
    if (rank >= 2) {
      if (dims(3) >= 2) {
        dx = x(,dif,..);
        r(,q,..) += dx;
        r(,p,..) -= dx;
      }
      if (rank >= 3) {
        if (dims(4) >= 2) {
          dx = x(,,dif,..);
          r(,,q,..) += dx;
          r(,,p,..) -= dx;
        }
        if (rank >= 4) {
          if (dims(5) >= 2) {
            dx = x(,,,dif,..);
            r(,,,q,..) += dx;
            r(,,,p,..) -= dx;
          }
          if (rank >= 5) {
            if (dims(6) >= 2) {
              dx = x(,,,,dif,..);
              r(,,,,q,..) += dx;
              r(,,,,p,..) -= dx;
            }
            if (rank >= 6) {
              if (dims(7) >= 2) {
                dx = x(,,,,,dif,..);
                r(,,,,,q,..) += dx;
                r(,,,,,p,..) -= dx;
              }
              if (rank >= 7) {
                if (dims(8) >= 2) {
                  dx = x(,,,,,,dif,..);
                  r(,,,,,,q,..) += dx;
                  r(,,,,,,p,..) -= dx;
                }
                if (rank >= 8) {
                  if (dims(9) >= 2) {
                    dx = x(,,,,,,,dif,..);
                    r(,,,,,,,q,..) += dx;
                    r(,,,,,,,p,..) -= dx;
                  }
                  if (rank >= 9) {
                    if (dims(10) >= 2) {
                      dx = x(,,,,,,,,dif,..);
                      r(,,,,,,,,q,..) += dx;
                      r(,,,,,,,,p,..) -= dx;
                    }
                    if (rank >= 10) {
                      if (dims(11) >= 2) {
                        dx = x(,,,,,,,,,dif,..);
                        r(,,,,,,,,,q,..) += dx;
                        r(,,,,,,,,,p,..) -= dx;
                      }
                      if (rank > 10) {
                        error, "too many dimensions";
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return r;
}

//-----------------------------------------------------------------------------
// OBJECTIVE FUNCTIONS

local objective_functions;
/* DOCUMENT Objective Functions

     An objective function in "invpb.i" is a Yorick function-like object, say
     `fn`, such that:

         fx = fn(x);

     yields the value of the objective function at `x`; while:

         local gx;
         fx = fn(x, gx, grad=1);

     yields the value of the objective function at `x` and stores its gradient
     in variable `gx` (when keyword `grad` is set true).

   SEE ALSO: quadratic_objective_function.
 */

func quadratic_objective_function(nil, y=, H=, W=)
/* DOCUMENT fn = quadratic_objective_function(H=, y=, W=);

     Build a function-like object `fn` implementing the quadratic objective
     function:

         f(x) = ‖H⋅x - y‖²_W
              = (H⋅x - y)'⋅W⋅(H⋅x - y)

     Keywords `H` and `W` may be specified as linear operators, the identity
     is assumed if unspecified.

     Keyword `y` may be specified as an array, zero is assumed if unspecified.

     For example, the quadratic objective function `rgl(x) = ‖D⋅x‖²` may be
     built by:

         rgl = quadratic_objective_function(W=DtD);

     where `DtD` implements `D'⋅D`.

   SEE ALSO: DtD, objective_functions, linear_operators.
 */
{
    if (nil != []) error, "all arguments given by keywords";
    Id = identity_operator();
    if (H == Id) H = [];
    if (W == Id) W = [];
    if (is_array(y) && allof(!y)) y = [];
    fn = h_new(y=y, H=H, W=W);
    h_evaluator, fn, "_call_quadratic_objective_function";
    return fn;
}

func _call_quadratic_objective_function(fn, x, &g, grad=)
{
    // Compute residuals: `r = H⋅x - y`.
    local y, H, r;
    eq_nocopy, y, fn.y;
    eq_nocopy, H, fn.H;
    if (H != []) {
        if (y != []) {
            r = H(unref(x)) - y;
        } else {
            r = H(unref(x));
        }
    } else {
        if (y != []) {
            r = unref(x) - y;
        } else {
            eq_nocopy, r, x;
        }
    }

    // Compute weighted residuals `Wr = W⋅r = W⋅(H⋅x - y)`.
    local W, Wr;
    eq_nocopy, W, fn.W;
    if (W != []) {
        Wr = W(r);
    } else {
        eq_nocopy, Wr, r;
    }

    // Compute gradient: `g = 2⋅H'⋅W⋅(H⋅x - y) = 2⋅H'⋅Wr`.
    if (grad) {
        g = Wr + Wr; // temporary
        if (H != []) {
            g = H(unref(g), 1n);
        }
    }

    // Compute objective function `f = (H⋅x - y)'⋅W⋅(H⋅x - y) = r'⋅Wr`.
    return sum(unref(r)*unref(Wr));
}

func always_provide_gradient(fn) {
    return closure("_call_always_provide_gradient", fn); }
func _call_always_provide_gradient(fg, x, &g) { return fg(x, g, grad=1n); }
/* DOCUMENT fg = always_provide_gradient(fn);

     Convert objective function `fn` into a function `fg` that always provides
     the gradient as expected by some optimization methods.  The returned
     object is callable as follows:

         local gx;
         fx = fg(x, gx);

      to store the gradient of the objective function in `gx` and yield the
      value of the objective function.

   SEE ALSO: objective_functions.
 */

func inner_product(x, y)
/* DOCUMENT inner_product(x, y);
     Compute the inner product of `x` and `y` which must be real-valued arrays.
 */
{
    return sum(x*y);
}

func composite_objective_function(..)
/* DOCUMENT fn = composite_objective_function([a1,] f1, [a2,] f2, ...);

     Build a composite objective function which behaves as:

         fn(x) = a1*f1(x) + a2*f2(x) + ...

     where arguments `a1`, `a2`, etc. are optional multipliers (assumed to be 1
     if omitted) while arguments `f1`, `f2`, etc. are objective functions.

   SEE ALSO: objective_functions.
 */
{
    local a, f;
    args = [];
    i = 0;
    while (more_args()) {
        ++i;
        eq_nocopy, a, next_arg();
        if (is_scalar(a)) {
            // A multiplifer and a function shall be specified.
            if (!(is_real(a) || is_integer(a)) || a < 0) {
                error, sprintf("multiplier #%d is not anonnegative number", i);
            }
            eq_nocopy, f, (more_args() ? next_arg() : []);
            if (f == []) {
                error, sprintf("missing function #%d", i);
            }
        } else if (a == []) {
            // Skip empty argument.
            --i;
            continue;
        } else {
            // Only a function has been specified.
            eq_nocopy, f, a;
            a = 1;
        }
        args = _cat(args, a, f);
    }
    fn = h_new(args = args, argc = i);
    h_evaluator, fn, _call_composite_objective_function;
    return fn;
}

func _call_composite_objective_function(fn, x, &g, grad=)
{
    args = fn.args;
    argc = fn.argc;
    f = 0.0;
    cnt = 0;
    while (--argc >= 0) {
        // Pop multiplier and function out of the list.
        local ai, fi;
        eq_nocopy, ai, _car(args); args = _cdr(args);
        eq_nocopy, fi, _car(args); args = _cdr(args);
        if (ai != 0) {
            local tmp;
            f += ai*fi(x, tmp, grad=grad);
            ++cnt;
            if (grad) {
                if (cnt > 1) {
                    // Increment gradient.
                    if (ai != 1) {
                        g = unref(g) + ai*unref(tmp);
                    } else {
                        g = unref(g) + unref(tmp);
                    }
                } else {
                    // Store in gradient.
                    if (ai != 1) {
                        g = ai*unref(tmp);
                    } else {
                        g = unref(tmp);
                    }
                }
            }
        }
    }
    if (cnt == 0 && grad) {
        g = array(floating_point_type(x), dimsof(x));
    }
    return f;
}

//-----------------------------------------------------------------------------
// REGULARIZATIONS

// Edge preserving regularization.
func isotropic_edge_preserving(nil, eps=, mu=)
/* DOCUMENT fn = isotropic_edge_preserving(eps=, mu=1);

     Build a function-like object `fn` implementing the hyperbolic smoothness
     quadratic objective function defined by:

         f(x) = µ Σᵢ fᵢ(x)

     with:

         fᵢ(x) = 2⋅ϵ⋅[sqrt(‖Dᵢ⋅x‖² + ϵ²) - ϵ]

     where `ϵ > 0` is the value of the parameter `eps` and `Dᵢ` is a linear
     operator which yields local finite differences around `i`-th pixel.  These
     functions have the following properties:

         min_x fᵢ(x) = 0       minimum is when Dᵢ⋅x = 0;

         fᵢ(x) ≈ ‖Dᵢ⋅x‖²       when ‖Dᵢ⋅x‖² ≪ ϵ;
         fᵢ(x) ≈ 2⋅ϵ⋅‖Dᵢ⋅x‖    when ‖Dᵢ⋅x‖² ≫ ϵ;

         ∂fᵢ(x)/∂x ≈ 2⋅ϵ⋅Dᵢ'⋅Dᵢ⋅x/sqrt(‖Dᵢ⋅x‖² + ϵ²)
         f(x) = 2⋅ϵ⋅[sqrt(‖Dᵢ⋅x‖² + ϵ²) - ϵ]

     Keywords `eps` and `mu` may be used to specify the respective values of `ϵ`
     and `µ`.

   SEE ALSO: DtD, objective_functions, linear_operators.
 */
{
    if (!is_func(rgl_totvar)) {
        include, "totvar.i", 3;
        if (!is_func(rgl_totvar)) {
            error, "TOTVAR plugin not installed";
        }
    }
    fn = h_new(eps=eps, mu=(is_void(mu) ? 1 : mu));
    h_evaluator, fn, "_call_isotropic_edge_preserving";
    return fn;
}

func _call_isotropic_edge_preserving(fn, x, &g, grad=)
{
    eps = fn.eps;
    phi = 2.0*mu*eps;
    if (grad) {
        f = phi*rgl_totvar(x, g, threshold=eps);
        g = phi*unref(g);
    } else {
        f = phi*rgl_totvar(x, threshold=eps);
    }
    return f;
}

//-----------------------------------------------------------------------------
// UTILITIES

func pad_array(x, dims, val)
/* DOCUMENT res = pad_array(arr, dims);
         or res = pad_array(arr, dims, val);

     Zero-pad array `arr` to dimensions `dims`.  The contents `arr` is
     approximately centered in the result.  If optional argument `val` is
     specified, this value (instead of zero) is used for padding.

   SEE ALSO: crop_array, array.
 */
{
    if (is_void(val)) {
        // Zero-pad keeping the same element type..
        val = structof(x);
    } else if (is_scalar(val)) {
        // Promote type (shall also works with strings, clash if incompatible).
        val += structof(x)();
    } else {
        error, "padding value must be a scalar";
    }
    if (! is_integer(dims) || dims(1) != (rank = numberof(dims) - 1) ||
        (numberof(dims) > 1 && min(dims) <= 0)) {
        error, "bad dimension list";
    }
    xdims = dimsof(x);
    if (rank < 1 || numberof(dims) != numberof(xdims)) {
        error, "bad number of dimensions";
    }
    diff = (dims - xdims)(2:);
    if (anyof(diff < 0)) {
        error, "output dimensions must be larger or equal input dimensions";
    }
    offs = diff/2;
    i = 1 + offs;
    j = xdims(2:) + offs;
    z = array(val, dims);
    if (rank == 1) {
        z(i(1):j(1)) = x;
    } else if (rank == 2) {
        z(i(1):j(1),i(2):j(2)) = x;
    } else if (rank == 3) {
        z(i(1):j(1),i(2):j(2),i(3):j(3)) = x;
    } else if (rank == 3) {
        z(i(1):j(1),i(2):j(2),i(3):j(3)) = x;
    } else if (rank == 4) {
        z(i(1):j(1),i(2):j(2),i(3):j(3),i(4):j(4)) = x;
    } else if (rank == 5) {
        z(i(1):j(1),i(2):j(2),i(3):j(3),i(4):j(4),i(5):j(5)) = x;
    } else if (rank == 6) {
        z(i(1):j(1),i(2):j(2),i(3):j(3),i(4):j(4),i(5):j(5),i(6):j(6)) = x;
    } else if (rank == 7) {
        z(i(1):j(1),i(2):j(2),i(3):j(3),i(4):j(4),i(5):j(5),i(6):j(6),i(7):j(7)) = x;
    } else if (rank == 8) {
        z(i(1):j(1),i(2):j(2),i(3):j(3),i(4):j(4),i(5):j(5),i(6):j(6),i(7):j(7),i(8):j(8)) = x;
    } else if (rank == 9) {
        z(i(1):j(1),i(2):j(2),i(3):j(3),i(4):j(4),i(5):j(5),i(6):j(6),i(7):j(7),i(8):j(8),i(9):j(9)) = x;
    } else if (rank == 10) {
        z(i(1):j(1),i(2):j(2),i(3):j(3),i(4):j(4),i(5):j(5),i(6):j(6),i(7):j(7),i(8):j(8),i(9):j(9),i(10):j(10)) = x;
    } else {
        error, "too many dimensions";
    }
    return z;
}

func crop_array(x, dims)
/* DOCUMENT res = crop_array(arr, dims);

     Crop array `arr` to a smaller dimensions size `dims`.  The central part of
     `arr` is returned.

   SEE ALSO: pad_array, array.
 */
{
    if (! is_integer(dims) || dims(1) != (rank = numberof(dims) - 1) ||
        (numberof(dims) > 1 && min(dims) <= 0)) {
        error, "bad dimension list";
    }
    xdims = dimsof(x);
    if (rank < 1 || numberof(dims) != numberof(xdims)) {
        error, "bad number of dimensions";
    }
    diff = (xdims - dims)(2:);
    if (anyof(diff < 0)) {
        error, "output dimensions must be smaller or equal input dimensions";
    }
    offs = diff/2;
    i = offs + 1;
    j = offs + dims(2:);
    if (rank ==  1) return x(i(1):j(1));
    if (rank ==  2) return x(i(1):j(1),i(2):j(2));
    if (rank ==  3) return x(i(1):j(1),i(2):j(2),i(3):j(3));
    if (rank ==  4) return x(i(1):j(1),i(2):j(2),i(3):j(3),i(4):j(4));
    if (rank ==  5) return x(i(1):j(1),i(2):j(2),i(3):j(3),i(4):j(4),i(5):j(5));
    if (rank ==  6) return x(i(1):j(1),i(2):j(2),i(3):j(3),i(4):j(4),i(5):j(5),i(6):j(6));
    if (rank ==  7) return x(i(1):j(1),i(2):j(2),i(3):j(3),i(4):j(4),i(5):j(5),i(6):j(6),i(7):j(7));
    if (rank ==  8) return x(i(1):j(1),i(2):j(2),i(3):j(3),i(4):j(4),i(5):j(5),i(6):j(6),i(7):j(7),i(8):j(8));
    if (rank ==  9) return x(i(1):j(1),i(2):j(2),i(3):j(3),i(4):j(4),i(5):j(5),i(6):j(6),i(7):j(7),i(8):j(8),i(9):j(9));
    if (rank == 10) return x(i(1):j(1),i(2):j(2),i(3):j(3),i(4):j(4),i(5):j(5),i(6):j(6),i(7):j(7),i(8):j(8),i(9):j(9),i(10):j(10));
    error, "too many dimensions";
}

func check_gradient(f, x0, g, number=, tiny=, dir=)
/* DOCUMENT check_gradient, f, x, g;
         or check_gradient(f, x);

     Compare gradient function with gradient estimated by finite differences.
     `f` is the function, `g` is the gradient and `x` the parameters.  Argument
     `g` can be an array (the computed gradient for `f` at `x0`) or a function
     (such that `g(x)` yields the gradient of `f` at `x`).

     The number of gradient values to check may be specified by keyword
     `number` (the subset of parameters is randomly chosen).  The default is to
     compute the finite difference gradient for all parameters.

     The finite difference step size cand be specified by keyword `tiny` which
     takes the form `tiny=rtol` or `tiny=[atol,rtol]` where `rtol` and `atol`
     are the absolute and relative step sizes (both must be strictly positive
     and `atol` is assumed to be 1e-20 if omitted).  Actual step size is `dx =
     ±(atol + rtol*abs(x))`.  By default, `tiny=[1e-20, 1e-5]`, hence
     `atol=1e-20` and `rtol=1e-5`.

     Keyword `dir` can be used to specify which kind of finite differences to
     use to estimate the gradient: forward (`dir > 0`), backward (`dir < 0`) or
     centered (`dir = 0` or nil).  The default is to use centered finite
     differences which are more precise (of order `rtol^3`) but twice more
     expensive to compute.

     When called as a subroutine, a summary of the gradient differences is
     printed.  Otherwise, the gradient computed by finite differences is
     returned.
*/
{
    local g0, g1;
    subroutine = am_subroutine();
    n = numberof(x0);
    if (is_void(number) || number >= n) {
        list = indgen(n);
        number = n;
    } else if (subroutine) {
        list = long(1 + n*random(number));
    } else {
        error, "all variables must be considered";
    }
    if (is_void(tiny)) {
        atol = 1e-20;
        rtol = 1e-5;
    } else if (numberof(tiny) == 1 && min(tiny) > 0) {
        atol = 1e-20;
        rtol = double(tiny(1));
    } else if (numberof(tiny) == 2 && min(tiny) > 0) {
        atol = double(tiny(1));
        rtol = double(tiny(2));
    } else {
        error, "invalid value(s) for TINY";
    }
    x1 = x0;
    f0 = f(x0);
    if (subroutine) {
        if (is_array(g)) {
            eq_nocopy, g0, g;
        } else {
            /* Assume a function or alike. */
            g0 = g(x0);
        }
        if (! same_dims(dimsof(g0), dimsof(x0))) {
            error, "the gradient and the variables must have the same dimensions";
        }
        g0 = g0(list);
        g1 = array(double, number);
    } else {
        g1 = array(double, dimsof(x0));
    }
    if (dir) {
        /* Forward or backward differences. */
        atol = abs(atol)*sign(dir);
        rtol = abs(rtol)*sign(dir);
        for (i = 1; i <= number; ++i) {
            l = list(i);
            x = x0(l);
            x1(l) = x + (dx = (atol + rtol*abs(x)));
            g1(i) = (f(x1) - f0)/dx;
            x1(l) = x;
        }
    } else {
        /* Centered differences. */
        for (i = 1; i <= number; ++i) {
            l = list(i);
            x = x0(l);
            dx = atol + rtol*abs(x);
            s = dx/2;
            x1(l) = x + s;
            f1 = f(x1);
            x1(l) = x - s;
            f2 = f(x1);
            x1(l) = x;
            g1(i) = (f1 - f2)/dx;
        }
    }
    if (! subroutine) return g1;

    abs_err = abs(g1 - g0);
    max_abs_err = max(abs_err);
    avg_abs_err = avg(abs_err);
    rms_abs_err = sqrt(avg(abs_err*abs_err));
    u = abs(g0) + abs(g1);
    rel_err = 2*(g1 - g0)*abs_err/(u + !u);
    max_rel_err = max(rel_err);
    avg_rel_err = avg(rel_err);
    rms_rel_err = sqrt(avg(rel_err*rel_err));

    write, format="GRADIENT CHECK WITH: tiny=[%.1e,%.1e],  number=%d\n",
        double(atol), double(rtol), number;
    write, format="ABSOLUTE ERROR: max=%.1e,  avg=%.1e,  rms=%.1e\n",
        max_abs_err, avg_abs_err, rms_abs_err;
    write, format="RELATIVE ERROR: max=%.1e,  avg=%.1e,  rms=%.1e\n",
        max_rel_err, avg_rel_err, rms_rel_err;
}

func finite_difference_gradient(f, x, sel=, atol=, rtol=, order=)
{
    if (is_void(atol)) atol = avg(abs(x))*1e-5;
    if (is_void(rtol)) rtol = 1e-5;
    if (is_void(sel)) {
        n = numberof(x);
        sel = indgen(n);
        g = array(double, dimsof(x));
    } else {
        n = numberof(sel);
        g =  array(double, dimsof(sel));
    }
    x += x + 0.0; // make a copy
    if (is_void(order) || order == 2) {
        for (j = 1; j <= n; ++j) {
            i = sel(j);
            x_i = x(i);
            h = atol + rtol*abs(x_i);
            x(i) = x_i - h; fm1 = f(x);
            x(i) = x_i + h; fp1 = f(x);
            x(i) = x_i;
            g(j) = (fp1 - fm1)/(2*h);
        }
    } else if (order == 4) {
        for (j = 1; j <= n; ++j) {
            i = sel(j);
            x_i = x(i);
            h = atol + rtol*abs(x_i);
            x(i) = x_i - 2*h; fm2 = f(x);
            x(i) = x_i -   h; fm1 = f(x);
            x(i) = x_i +   h; fp1 = f(x);
            x(i) = x_i + 2*h; fp2 = f(x);
            x(i) = x_i;
            g(j) = (fm2 - 8*fm1 + 8*fp1 - fp2)/(12*h);
        }
    } else {
        error, "keyword `order` must be 2 or 4";
    }
    return g;
}

func same_dims(adims, bdims)
/* DOCUMENT same_dims(adims, bdims);

     Check whether `adims` and `bdims` are the same dimensions lists.

   SEE ALSO: dimsof.
 */
{
  return (numberof(adims) == numberof(bdims) && allof(adims == bdims));
}

func floating_point_type(arg)
/* DOCUMENT T = floating_point_type(arg);

     Determine floating-point type corresponding to `arg`; that is `float` if
     `arg` is `float` or is a single precision floating-point scalar or array,
     otherwise, `double`.

   SEE ALSO: structof, array.
 */
{
    if (is_array(arg)) {
        arg = structof(arg);
    }
    if (arg == float) {
        return float;
    }
    if (arg == double || arg == long || arg == int ||
        arg == short || arg == char) {
        return double;
    }
    error, "invalid argument type";
}

local sprintf, printf, fprintf;
/* DOCUMENT printf, fmt, a1, a2,  ...;
         or fprintf, io, fmt, a1, a2,  ...;
         or str = sprintf(fmt, a1, a2, ...);

     Use format string `fmt` to print arguments `a1`, `a2`, etc. to standard
     or `io` output stream or to build a string `str`.

   SEE ALSO: write.
 */
func printf(fmt, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)
{
    if (!(is_string(fmt) && is_scalar(fmt))) {
        error, "usage: printf, fmt, a1, a2, ...;";
    }
    if (a1 == []) {
        write, format="%s", fmt;
    } else if (a2 == []) {
        write, format=fmt, a1;
    } else if (a3 == []) {
        write, format=fmt, a1, a2;
    } else if (a4 == []) {
        write, format=fmt, a1, a2, a3;
    } else if (a5 == []) {
        write, format=fmt, a1, a2, a3, a4;
    } else if (a6 == []) {
        write, format=fmt, a1, a2, a3, a4, a5;
    } else if (a7 == []) {
        write, format=fmt, a1, a2, a3, a4, a5, a6;
    } else if (a8 == []) {
        write, format=fmt, a1, a2, a3, a4, a5, a6, a7;
    } else if (a9 == []) {
        write, format=fmt, a1, a2, a3, a4, a5, a6, a7, a8;
    } else if (a10 == []) {
        write, format=fmt, a1, a2, a3, a4, a5, a6, a7, a8, a9;
    } else {
        write, format=fmt, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10;
    }
}

func fprintf(io, fmt, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)
{
    if (!(is_string(fmt) && is_scalar(fmt))) {
        error, "usage: fprintf, io, fmt, a1, a2, ...;";
    }
    if (a1 == []) {
        write, io, format="%s", fmt;
    } else if (a2 == []) {
        write, io, format=fmt, a1;
    } else if (a3 == []) {
        write, io, format=fmt, a1, a2;
    } else if (a4 == []) {
        write, io, format=fmt, a1, a2, a3;
    } else if (a5 == []) {
        write, io, format=fmt, a1, a2, a3, a4;
    } else if (a6 == []) {
        write, io, format=fmt, a1, a2, a3, a4, a5;
    } else if (a7 == []) {
        write, io, format=fmt, a1, a2, a3, a4, a5, a6;
    } else if (a8 == []) {
        write, io, format=fmt, a1, a2, a3, a4, a5, a6, a7;
    } else if (a9 == []) {
        write, io, format=fmt, a1, a2, a3, a4, a5, a6, a7, a8;
    } else if (a10 == []) {
        write, io, format=fmt, a1, a2, a3, a4, a5, a6, a7, a8, a9;
    } else {
        write, io, format=fmt, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10;
    }
}

func sprintf(fmt, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)
{
    if (!(is_string(fmt) && is_scalar(fmt))) {
        error, "usage: sprintf(fmt, a1, a2, ...);";
    }
    if (a1 == []) {
        return fmt;
    } else if (a2 == []) {
        return swrite(format=fmt, a1);
    } else if (a3 == []) {
        return swrite(format=fmt, a1, a2);
    } else if (a4 == []) {
        return swrite(format=fmt, a1, a2, a3);
    } else if (a5 == []) {
        return swrite(format=fmt, a1, a2, a3, a4);
    } else if (a6 == []) {
        return swrite(format=fmt, a1, a2, a3, a4, a5);
    } else if (a7 == []) {
        return swrite(format=fmt, a1, a2, a3, a4, a5, a6);
    } else if (a8 == []) {
        return swrite(format=fmt, a1, a2, a3, a4, a5, a6, a7);
    } else if (a9 == []) {
        return swrite(format=fmt, a1, a2, a3, a4, a5, a6, a7, a8);
    } else if (a10 == []) {
        return swrite(format=fmt, a1, a2, a3, a4, a5, a6, a7, a8, a9);
    } else {
        return swrite(format=fmt, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
    }
}

// A simple function to plot an image.
func plimg(img, fig=, colors=, cmin=, cmax=, title=, xlabel=, ylabel=)
{
    dims = dimsof(img);
    if (numberof(dims) == 3) {
        // Gray-scaled image.
        width = dims(2);
        height = dims(3);
    } else if (numberof(dims) == 4 && numberof(dims) == 3 && structof(img) == char) {
        // RGB image.
        width = dims(3);
        height = dims(4);
    } else {
        error, "invalid image format";
    }
    if (!is_void(fig)) {
        window, fig;
    }
    fma;
    if (!is_void(colors)) {
        cmap, colors;
    }
    pli, img, 0.5, 0.5, width + 0.5, height + 0.5, cmin=cmin, cmax=cmax;
    if (!is_void(title)) pltitle, title;
    if (!is_void(xlabel) || !is_void(ylabel)) {
        xytitle, xlabel, ylabel;
    }
}
errs2caller, plimg;

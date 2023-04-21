// saturn_data.i -
//
// Load Saturn data for tests.
//
// Global variables:
// - `DATADIR` the directory where are stored the data files,
// - `dat` the data,
// - `psf` the PSF,
// - `dims` the dimensions of the problem,
// - `wgt` the weights (0 or 1 everywhere),
// - `fraction_of_bad_pixels`.
//
//-----------------------------------------------------------------------------

if (!batch()) {
    // Needed for plotting.
    require, "invpb.i";
}
// Directory where are stored the data files.
DATADIR = dirname(current_include()) + "/../data/saturn/";

// Read the data.
dat = transpose(fits_read(DATADIR+"saturn.fits"));
psf = transpose(fits_read(DATADIR+"saturn_psf.fits"));

// Size of the problem.
dims = dimsof(dat);

// Assume 70% of bad pixels by default in interactive mode, 0% in batch mode.
if (is_void(fraction_of_bad_pixels)) {
    fraction_of_bad_pixels = (batch() ? 0.0 : 0.7);
}

// Array of pixelwise weights and corresponding operator.
wgt = (dat > min(dat)); // eliminate real bad pixels
if (fraction_of_bad_pixels > 0) {
    // add random bad pixels
    wgt &= (random(dimsof(wgt)) > fraction_of_bad_pixels);
}
wgt = double(wgt); // convert booleans

if (!batch()) {
    plimg, dat, cmin=0, fig=1, title="Raw data";
    plimg, psf, cmin=0, fig=2, title="PSF";
    plimg, (wgt > 0)*dat, cmin=0, fig=3, title="Good pixels in data";
}

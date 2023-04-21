# IML_FISTA
## Matlab code to run IML FISTA
 Run demo_IMLFISTA to compare FISTA and IML FISTA on the same optimization problem \
 Be careful of the image chosen as CPU time increases rapidly with dimension
 ## Important note :
 To function, the code uses a Matlab toolbox to compute proximity operators of TV based norms written by Giovanni Chierchia [prox-repository](http://proximity-operator.net). \
 The information transfers are based on wavelet transforms and thus a function of the [Wavelab850](https://statweb.stanford.edu/~wavelab/) toolbox is included here 
 # Problems possible
 see create_data.m \
 Inpainting with 50 % missing pixels \
 Inpainting with 90 % missing pixels \
 Small gaussian blur \
 High gaussian blur 
 
 # Images possible
 JWST deep field in sizes 256 x 256 x 3 to 2048 x 2048 x 3 + gray version in 256 x 256

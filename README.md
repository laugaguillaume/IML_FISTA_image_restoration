# IML_FISTA
## Matlab code to run IML FISTA
 Run demo_IMLFISTA.m to compare FISTA and IML FISTA on the same optimization problem \
 Be careful of the image chosen as CPU time increases rapidly with dimension
 ## Read before running the code:
 To function, the code uses a Matlab toolbox to define proximity operators of TV based norms written by Giovanni Chierchia [prox-repository](http://proximity-operator.net). \
 The blur operators are constructed using the HNO package which is available here: [HNO](http://www.imm.dtu.dk/~pcha/HNO/). The HNO folder contains every function of this package. \
 The information transfers are based on wavelet transforms and we use the MakeONFilter function of the [Wavelab850](https://statweb.stanford.edu/~wavelab/) toolbox to construct the filters. Only this function is included here: filters available are ones of orthogonal wavelets.
 # Problems possible
 see create_data.m \
 Inpainting with 50 % missing pixels \
 Inpainting with 90 % missing pixels \
 Small gaussian blur \
 High gaussian blur 
 
 # Images possible
 JWST deep field in sizes 256 x 256 x 3 to 2048 x 2048 x 3 + gray version in 256 x 256 \
 Jupiter taken by JWST \
 (links to images in create_data.m)

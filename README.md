# IML_FISTA
## Matlab code to run IML FISTA
 Run demo_IMLFISTA.m to compare FISTA and IML FISTA on the same optimization problem \
 Be careful of the image chosen as CPU time increases rapidly with dimension
 ## Read before running the code:
 To function, the code uses a Matlab toolbox to define proximity operators of TV based norms written by Giovanni Chierchia [prox-repository](http://proximity-operator.net). \
 The blur operators are constructed using the HNO package which is available here: [HNO](http://www.imm.dtu.dk/~pcha/HNO/). The HNO folder contains every function of this package. \
 The information transfers are based on wavelet transforms and we use the MakeONFilter function of the [Wavelab850](https://statweb.stanford.edu/~wavelab/) toolbox to construct the filters. Only this function is included here: filters available are ones of orthogonal wavelets.
 ## Potential problem with TV toolbox
Matlab may refuse to run functions of the TV toolbox and return an error message saying the author of the C++ code is unknown. I obtained the code from the people who wrote the toolbox so I trust them but if you don't you may download the code directly from their website.
 # Available degradations
 see create_data.m \
 Inpainting with 50 % missing pixels \
 Inpainting with 90 % missing pixels \
 Small gaussian blur \
 High gaussian blur 
 
 # Available images
 JWST deep field in sizes 256 x 256 x 3 to 2048 x 2048 x 3 + gray version in 256 x 256. Credits: IMAGE: NASA, ESA, CSA, STScI [Source](https://webbtelescope.org/contents/media/images/2022/035/01G7DCWB7137MYJ05CSH1Q5Z1Z)\
 Jupiter taken by JWST in size 1024 x 1024. Credits: NASA, ESA, CSA, Jupiter ERS Team; image processing by Judy Schmidt. [Source](https://blogs.nasa.gov/webb/2022/08/22/webbs-jupiter-images-showcase-auroras-hazes/)\
 Near IR image of Pillars of Creations in size 2048 x 2048. Credits: SCIENCE: NASA, ESA, CSA, STScI. [Source](https://webbtelescope.org/contents/media/images/01GK2KKTR81SGYF24YBGYG7TAP.html) \

# shellblur-model
In this repository, you can find a XSPEC model to fit the line-profile considering different geometries


First, you have to initialize the model in xspec with the command:
initpackage shellbluri shellblur.dat ./

To use shellblur, remember that it is a convolution model, so the syntaxis is:

mo TBabs*shellblur*vpshock

Do not forget to download the table from:

https://www.dropbox.com/s/ts1jfrg68nx38

Author: Jonathan Quirola
Email: jquirola@astro.puc.cl

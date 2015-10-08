library(Rcpp)
library(RcppArmadillo)
# Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
# Sys.setenv("PKG_LIBS" = "-fopenmp")
sourceCpp("BSFG_functions_c.cpp")
a1(diag(2))

a2(1:2)



# compileAttributes('adsf')
# detach(package:adsf,unload=T)
# library.dynam.unload("adsf", system.file(package = "adsf"))
# install.packages('adsf',repos=NULL,type='source')
# library(adsf)
# a3(diag(2))
# rcpparma_hello_world()
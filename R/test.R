library(Rcpp)
sourceCpp("../IRTpp/src/test.cpp")
cpower(3)
Sys.setenv( "PKG_CXXFLAGS"="-std=c++11 -I/home/liberato/IRTpp/src/SICSRepository/SICS/src")
sourceCpp("../IRTpp/src/rcpp_hello_world.cpp")


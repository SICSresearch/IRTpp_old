library(Rcpp)
sourceCpp("../IRTpp/src/test.cpp")
cpower(3)
Sys.setenv( "PKG_LIBS"="-L/home/liberato/git/IRTpp/src/SICSRepository/SICS/src")
Sys.setenv( "PKG_CXXFLAGS"=" -std=c++11 -I/home/liberato/git/IRTpp/src/SICSRepository/SICS/src -fPIC")
compileAttributes()
sourceCpp("../IRTpp/src/rcpp_hello_world.cpp")


% compile Boonen13 C++ mex code
mex -DGALERKIN_ISO -v CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" Boonen13.mex.cpp -I../../../eigen -output ../../Boonen13_gal_iso
mex -DGALERKIN_CONSTANT -v CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" Boonen13.mex.cpp -I../../../eigen -output ../../Boonen13_gal_const
mex -DCOLLOC_CONSTANT -v CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" Boonen13.mex.cpp -I../../../eigen -output ../../Boonen13_col_const

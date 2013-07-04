% compile Boonen13 C++ mex code
% tic
% mex -DGALERKIN -DISOPARAMETRIC -v CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" Boonen13.mex.cpp -I../../../eigen -output ../../Boonen13_gal_iso
% toc

tic
mex -DGALERKIN -DCONSTANT CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" Boonen13.mex.cpp -I../../../eigen -output ../../Boonen13_gal_const
toc

% tic
% mex -DCOLLOCATIONAL -DCONSTANT -v CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" Boonen13.mex.cpp -I../../../eigen -output ../../Boonen13_col_const
% toc

% compile UnitBoonen13 C++ mex code
% mex -DGALERKIN_CONSTANT -v CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" UnitBoonen13.mex.cpp -I../../../eigen -output ../../UnitBoonen13_gal_const

clear;
clc;

%eigendir = fullfile(mytoolboxdir(), 'eigen');
%nihudir = fullfile(mytoolboxdir(), 'nihu', 'src');
%opts = strcat('-I', eigendir, ' -I', nihudir, ' -std=c++11');

%mex('aca_test.mex.cpp', '-v', opts, '-output', 'aca_test');

mex -v CXXFLAGS="\$CXXFLAGS -std=c++11 -I/D/research/toolbox/eigen -I/D/research/toolbox/nihu/src -Wall" aca_test.mex.cpp -o aca_test

movefile(strcat('*.', mexext()), '../');


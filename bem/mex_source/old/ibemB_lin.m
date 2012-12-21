%ibemB_lin  Generate acoustic BEM system matrices
%   B = ibemB_lin(nodes, elements, g3, g4, dist, k, points)
%   Computes the acoustic BEM system matrix B for the case of constant
%   shape functions.
% Parameters:
%   nodes    : N x 3 matrix of xyz coordinates of model vertices
%   elements : matrix of element node indices. Each row describes one
%              element in the form [3 n1 n2 n3 0] for TRIA and
%              [4 n1 n2 n3 n4] for QUAD elements.
%   g3, g4   : Gaussian quadrature structure for TRIA and QUAD elements
%              obtained from GAUSSQUAD2
%   dist     : Distance limits governing integration density parameter
%   k        : acoustic wave number
%   points   : optional, if not given, then the acoustuc surface matrices
%              are computed. If defined, then the field point matrices
%              are returned. Points is a M x 3 matrix containing the xyz
%              coordinates of the field points
%
% Peter Fiala
% 2011
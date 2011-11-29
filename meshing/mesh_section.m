function [mesh, nodeind, elemind] = mesh_section(mesh, limits, selmode)
%MESH_SECTION   Return a rectangular section of a mesh
%   [SECTION, NODEIND, ELEMIND] = MESH_SECTION(MESH, LIMITS, SELMODE)
%   returns a rectangular section of the fe mesh mesh. The section mesh
%   contains those elements that are located within the LIMITS
% Parameters:
%  MESH    : AcouFEM mesh structure
%  LIMITS  : 2x3 matrix [xmin ymin zmin; xmax ymax zmax] (can contain -Inf
%            and +Inf for no selection)
%  SELMODE : 'all': An element is selected if all of its nodes are within
%                   the limints.
%            'any': An element is selected if any of its nodes are within
%                   the limits.
%            'none': An element is selected if none of its nodes are within
%                   the limits.
%            'nall': An element is selected if not all of its nodes are
%                   within the limits.
%  SECTION : AcouFEM mesh structure containing the desired section. The
%            returned mesh contains all the nodes of the original mesh,
%            only the element structure is changed.
%  NODEIND : The indices of those nodes that are located within the limits.
%  ELEMIND : The indices of those elements that are located within the
%            limits.
%
% Example:
%   mesh = create_sphere_boundary(2, 20);
%   tol = 1e-3;
%   mesh = mesh_section(mesh, [-Inf, -Inf, 1-tol; +Inf, 1.4+tol, 2.5+tol], 'any');
%   plot_mesh(mesh);
%   
% See also: FE_SELECT

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last updated: 02.12.2009.

%% Parameter check
error(nargchk(2, 3, nargin, 'struct'));

d = find(size(limits) == 3);
if d == 1
    limits = limits.';
end
if nargin < 3
    selmode = 'all';
end

%% selection of appropriate nodes and elements
expression = sprintf('x>=%f & x<=%f & y>=%f & y<=%f & z>=%f & z<=%f', limits);
[nodeind elemind] = mesh_select(mesh, expression, 'ind', selmode);
mesh.Elements = mesh.Elements(elemind,:);
end

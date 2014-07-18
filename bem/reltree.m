function [tr, rr, rs, ns] = reltree(tree, xr, fr, xs, fs, ns, symm)
%RELTREE Convert cluster tree to relative tree
%   [tr, rr, rs, ns] = reltree(tree, xr, fr, xs, fs, ns, symm) Converts the
%   cluster tree to relative tree. In the relative tree, the father
%   children distances are stored rather than cluster centers.
% Parameters:
%   tree : Cluster tree obtained by clustertree
%   xr   : Nrx3 xyz receiver coordinates
%   fr   : Nrx1 leaf level cluster indices to each receiver
%   xs   : Nsx3 xyz source coordinates
%   fs   : Nsx1 leaf level cluster indices to each source
%   ns   : Nsx3 xyz source normals
%   symm : (-1, 0, +1) = (positive symm, no symm [default], antisymm)
% Output:
%   tr   : tree structure containing cluster-to-cluster distances
%   rr   : 3xNr receiver node distances from leaf clusters
%   rs   : 3xNs source node distances from leaf clusters
%   ns   : 3xNs source normals (only transposition)

% Peter Fiala
% 2009

%% Parameter check
if nargin < 7
    symm = 0;
end

%% Allocating the new tree structure
c = cell(length(tree),1);
tr = struct('father', c, 'interlist', c, 'r', c);

%% Computing relative tree entries
for l = 1 : length(tree)
    tr(l).father = int32(tree(l).father);
    tr(l).interlist = int32(tree(l).interlist).'; %transp. for faster C routine
    if symm
        tr(l).iminterlist = int32(tree(l).iminterlist).';
    end
    if l > 1
        tr(l).r = (tree(l-1).coord(tree(l).father,:) - tree(l).coord).';
    end
end

%% source and receiver distances
rs = (tree(end).coord(fs,:) - xs).';
rr = (tree(end).coord(fr,:) - xr).';
ns = ns.';
end

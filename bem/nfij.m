function [i, j] = nfij(nearfield, source, receiver)
%NFIJ Generate near field sparse matrix locations of a FMBEM model
%  [i, j] = nfij(nearfield, source, receiver) generates the location
%  indices of the sparse matrices Hnf and Gnf.
% Parameters:
%   nearfield : Nearfield structure of the leaf level of the cluster
%               tree. This entry is usually obtained as tree(end).nearfield
%               where tree is the cluster tree generated by clustertree
%   source    : Source node index matrix of the leaf level, usually
%               obtained as tree(end).nodsou
%   receiver  : Receiver node index matrix of the leaf level, usually
%               obtained as tree(end).nodrec. This parameter is optional,
%               the default value is receiver = source
%
% See also: clustertree

% Peter Fiala
% 2009

%% Parameter check
if nargin == 2
    receiver = source;
end

%% Build sparse matrix
[cls, ~, clr] = find(nearfield);
[ir, ~, jr] = find(receiver(clr,:));
[is, ~, j] = find(source(cls(ir),:));
i = jr(is);
i = i(:);
j = j(:);
end

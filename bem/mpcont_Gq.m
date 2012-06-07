function [Gq tleaf ttrans tsampup tshiftup tsampdn tshiftdn trec] = mpcont_Gq(rr, rs, q, tree, intdata, k, fathersou, fatherrec, symm)
%MPCONT_GQ Multipole contribution of the Gq integral operator
%  Gq = mpcont_Gq(rr, rs, q, rtree, intdata, k, fs, fr, symm)
%   Computes the multipole transfer between sources and receivers.
% Parameters:
%   rr      : 3xN matrix, xyz coordinates of receivers.
%   rs      : 3xM matrix, xyz coordinates of sources.
%   q       : 3xM source amplitudes
%   rtree   : relative cluster tree obtained by reltree
%   intdata : integration parameter structure array built by integpar
%   k       : scalar real wave number
%   fs      : Mx1 father cluster indices of source elements
%   fr      : Nx1 father cluster indices of receiver elements
%   symm    : symmetry parameter (-1, 0 or +1) (antisymm, no symm, symm)
%
% [Gq tleaf ttrans tsampup tshiftup tsampdn tshiftdn trec] = mpcont_Gq(...)
%    also returns computation times spent on each level of the tree.
%  tleaf    : time for computing leaf level far dield signatures
%  ttrans   : (depth+1)x1 times for translation
%  tsampup  : (depth+1)x1 times for interpolation
%  tshiftup : (depth+1)x1 times for shiftup
%  tsampdn  : (depth+1)x1 times for filtering
%  tshiftdn : (depth+1)x1 times for downshift
%  trec     : time for recovering the solution
%  
% See also: mpcont_Hp, clustertree, reltree, integpar

% Peter Fiala
% 2009

%% parameter check
if nargin < 9
    symm = 0;
end

%% timing initialization
depth = length(tree)-1;
if nargout > 1
    ttrans   = zeros(depth+1,1);
    tsampup  = zeros(depth+1,1);
    tshiftup = zeros(depth+1,1);
    tsampdn  = zeros(depth+1,1);
    tshiftdn = zeros(depth+1,1);
end

%% initialization
% trivial solution (important because of initial step of GMRES)
if norm(q) == 0
    Gq = zeros(size(rr,2),1);
    return;
end
% excitation should be strictly complex because of MEX code LeafGq
q = complex(real(q), imag(q));
% Structure to store the near field signatures
Sig = struct('Nq', cell(depth+1,1));

%% Computing far field signature at leaf level
l = depth+1;
tstart = tic;
Fq = leafGq(rs, q, uint32(fathersou)-1, size(tree(l).father,1), intdata(l).S, k);
tleaf = toc(tstart);

%% Translation and upward shift
mindepth = mininterdepth(tree, symm); % highest tranalation level
for l = depth+1 : -1 : mindepth
    % current tree level
    T = tree(l);
    I = intdata(l);
    %% translation
    tstart = tic;
    if ~isempty(T.interlist)
        Sig(l).Nq = trans(Fq, T.interlist-1, I.Dindex-1,...
            I.Pindex-1, I.Perm-1, I.M);
    else
        Sig(l).Nq = zeros(size(Fq));
    end
    if symm
        if ~isempty(T.iminterlist)
            Sig(l).Nq = Sig(l).Nq + symm *...
                trans(Fq(I.Mz,:), T.iminterlist-1, I.imDindex-1,...
                I.imPindex-1, I.Perm-1, I.M);
        end
    end
    ttrans(l) = toc(tstart);
    %% upward shift
    if l > mindepth
        % receiver level
        T2 = tree(l-1);
        I2 = intdata(l-1);
        % interpolate
        tstart = tic;
        Fq2  = sphresamp(Fq, I.L, I2.L, I.Aup);
        tsampup(l) = toc(tstart);
        % shift up
        tstart = tic;
        Fq = shiftup(T.r, Fq2, T.father-1, length(T2.father), I2.S, k);
        tshiftup(l) = toc(tstart);
    end
end

%% Downward cycle
for l = mindepth : depth
    I = intdata(l);     % source level
    T2 = tree(l+1);     % response level
    I2 = intdata(l+1);
    % filter near field signature
    tstart = tic;
    Nq  = sphresamp(Sig(l).Nq, I.L, I2.L, I.Adn);
    tsampdn(l) = toc(tstart);
    % shift down
    tstart = tic;
    Sig(l+1).Nq = Sig(l+1).Nq + ...
        shiftdown(T2.r, Nq, T2.father-1, I2.S, k);
    tshiftdn(l) = toc(tstart);
end

%% response at leaf level
l = depth+1;
I = intdata(l);
tstart = tic;
Gq = recover(rr, Sig(l).Nq, uint32(fatherrec)-1, I.S, I.W, k);
trec = toc(tstart);
end

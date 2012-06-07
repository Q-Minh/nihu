function [Hp tleaf ttrans tsampup tshiftup tsampdn tshiftdn trec] = mpcont_Hp(rr, rs, ns, p, tree, intdata, k, fs, fr, symm)
%MPCONT_HP Multipole contribution of the integral operator Hp
%  Hp = mpcont_Hp(rr, rs, ns, p, rtree, intdata, k, fs, fr, symm)
%   Computes the multipole transfer between sources and receivers.
% Parameters:
%   rr      : 3xN matrix, xyz coordinates of receivers.
%   rs      : 3xM matrix, xyz coordinates of sources.
%   ns      : 3xM matrix, xyz coordinates of source normals
%   p       : Mx1 source amplitudes
%   rtree   : relative cluster tree obtained by reltree
%   intdata : integration parameter structure array built by integpar
%   k       : scalar real wave number
%   fs      : Mx1 father cluster indices of source elements
%   fr      : Nx1 father cluster indices of receiver elements
%   symm    : symmetry parameter (-1, 0 or +1) (antisymm, no symm, symm)
%
% [Hp tleaf ttrans tsampup tshiftup tsampdn tshiftdn trec] = mpcont_Hp(...)
%    also returns computation times spent on each level of the tree.
%  tleaf    : time for computing leaf level far dield signatures
%  ttrans   : (depth+1)x1 times for translation
%  tsampup  : (depth+1)x1 times for interpolation
%  tshiftup : (depth+1)x1 times for shiftup
%  tsampdn  : (depth+1)x1 times for filtering
%  tshiftdn : (depth+1)x1 times for downshift
%  trec     : time for recovering the solution
%  
% See also: mpcont_Gq, clustertree, reltree, integpar

% Peter Fiala
% 2009

%% parameter check
if nargin < 10
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
if norm(p) == 0
    Hp = zeros(size(rr,2),1);
    return;
end
% excitation should be strictly complex because of MEX code LeafHp
p = complex(real(p), imag(p));
% Structure to store the near field signatures
Sig = struct('dNp', cell(depth+1,1));

%% Computing far field signature at leaf level
l = depth+1;
tstart = tic;
dFp = leafHp(rs, ns, p, int32(fs)-1, size(tree(l).father,1), intdata(l).S, k);
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
        Sig(l).dNp = trans(dFp, T.interlist-1, I.Dindex-1,...
            I.Pindex-1, I.Perm-1, I.M);
    else
        Sig(l).dNp = zeros(size(dFp));
    end
    if symm
        if ~isempty(T.iminterlist)
            Sig(l).dNp = Sig(l).dNp + symm *...
                trans(dFp(I.Mz,:), T.iminterlist-1, I.imDindex-1,...
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
        dFp = sphresamp(dFp, I.L, I2.L, I.Aup);
        tsampup(l) = toc(tstart);
        % shift up
        tstart = tic;
        dFp = shiftup(T.r, dFp, T.father-1, length(T2.father), I2.S, k);
        tshiftup(l) = toc(tstart);
    end
end

%% Downward shift
for l = mindepth : depth
    I = intdata(l);     % source level
    T2 = tree(l+1);     % response level
    I2 = intdata(l+1);
    % filter near field signature
    tstart = tic;
    dNp = sphresamp(Sig(l).dNp, I.L, I2.L, I.Adn);
    tsampdn(l) = toc(tstart);
    % shift down
    tstart = tic;
    Sig(l+1).dNp = Sig(l+1).dNp + shiftdown(T2.r, dNp, T2.father-1, I2.S, k);
    tshiftdn(l) = toc(tstart);
end

%% Response at leaf level
l = depth+1;
I = intdata(l);
tstart = tic;
Hp = recover(rr, Sig(l).dNp, int32(fr)-1, I.S, I.W, k);
trec = toc(tstart);
end

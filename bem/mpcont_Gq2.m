function Gq = mpcont_Gq2(nodrec, nodsou, q, tree, intdata, k, fathersou, fatherrec, symm)
%MPCONT_GQ2 Compute multipole contribution of the Gq integral operator
%   Computes the multipole transfer between a set of elements.
%   The transfer is computed for each pair of
%   nodes in the matrix NODES that are in the far field of each other.
%
% See also: MPCONT_HP GEO2GAUSS CLUSTERTREE INTEGPAR

% Peter Fiala
% 2009

%% parameter check
if nargin < 9
    symm = 0;
end

%% initialization
% trivial solution (important because of initial step of GMRES)
if norm(q) == 0
    Gq = zeros(size(nodrec,1),1);
    return;
end
% excitation should be strictly complex because of MEX code LeafGq
q = complex(real(q), imag(q));
% depth of the tree
depth = length(tree)-1;
% Structure to store the near field signatures
Sig = struct('Nq', cell(depth+1,1));

%% Computing far field signature at leaf level
l = depth+1;
Fq = leafGq(nodsou, q, fathersou-1, tree(l).coord, intdata(l).S.', k);

%% Translation and upward shift
if symm
    mindepth = 2;
else
    mindepth = 3;
end
for l = depth+1 : -1 : mindepth
    % translation
    load(sprintf('data/M_%d.mat', l-1), 'M');
    T = tree(l);
    I = intdata(l);
    Sig(l).Nq = trans(Fq, T.interlist-1, T.Dindex-1, M);
    if symm
        if ~isempty(T.iminterlist)
            Sig(l).Nq = Sig(l).Nq + symm * trans(Fq(:,I.perm), T.iminterlist-1, T.imDindex-1, M);
        end
    end
    % upward shift
    if l > mindepth
        T2 = tree(l-1);
        I2 = intdata(l-1);
        Fq2  = sphresamp(Fq.', I, I2).';
        Fq = shiftup(T.coord, Fq2, T.father-1, T2.coord, I2.S.', k);
    end
end

%% Downward cycle
for l = mindepth : depth
    T = tree(l);
    I = intdata(l);
    T2 = tree(l+1);
    I2 = intdata(l+1);
    Nq  = sphresamp(Sig(l).Nq.', I, I2).';
    Sig(l+1).Nq = Sig(l+1).Nq + ...
        shiftdown(T.coord, Nq, T2.father-1, T2.coord, I2.S.', k);
end

%% response at leaf level
l = depth+1;
I = intdata(l);
Gq = recover(tree(l).coord, Sig(l).Nq, fatherrec-1, nodrec, I.S.', I.W, k);

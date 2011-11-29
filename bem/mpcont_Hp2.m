function Hp = mpcont_Hp2(nodrec, nodsou, normsou, p, tree, intdata, k, fathersou, fatherrec, symm)
%MPCONT_HP2 Compute multipole contribution of the Hp integral operator
%   Computes the multipole transfer between a set of elements.
%   The transfer is computed for each pair of
%   nodes in the matrix NODES that are in the far field of each other.
%   Parameters:
%
% See also: MPCONT_GQ GEO2GAUSS CLUSTERTREE INTEGPAR

% Peter Fiala
% 2009

%% parameter check
if nargin < 10
    symm = 0;
end

%% initialization
% trivial solution (important because of initial step of GMRES)
if norm(p) == 0
    Hp = zeros(size(nodrec,1),1);
    return;
end
% excitation should be strictly complex because of MEX code LeafHp
p = complex(real(p), imag(p));
% depth of the tree
depth = length(tree)-1;
% Structure to store the near field signatures
Sig = struct('dNp', cell(depth+1,1));

%% Computing far field signature at leaf level
l = depth+1;
dFp = leafHp(nodsou, normsou, p, fathersou-1, tree(l).coord, intdata(l).S.', k);

%% Translation and upward shift
if symm
    mindepth = 2;
else
    mindepth = 3;
end
for l = depth+1 : -1 : mindepth
    % current tree level
    T = tree(l);
    I = intdata(l);
    % translation
    load(sprintf('data/M_%d.mat', l-1), 'M');
    Sig(l).dNp = trans(dFp, T.interlist-1, T.Dindex-1, M);
    if symm
        if ~isempty(T.iminterlist)
            Sig(l).dNp = Sig(l).dNp + ...
                symm * trans(dFp(:,I.perm), T.iminterlist-1, T.imDindex-1, M);
        end
    end
    % upward shift
    if l > mindepth
        % receiver level
        T2 = tree(l-1);
        I2 = intdata(l-1);
        % interpolate
        dFp2 = sphresamp(dFp.', I, I2).';
        % shift up
        dFp = shiftup(T.coord, dFp2, T.father-1, T2.coord, I2.S.', k);
    end
end

%% Downward cycle
for l = mindepth : depth
    % source level
    T = tree(l);
    I = intdata(l);
    % response level
    T2 = tree(l+1);
    I2 = intdata(l+1);
    % filter near field signature
    dNp = sphresamp(Sig(l).dNp.', I, I2).';
    % shift down
    Sig(l+1).dNp = Sig(l+1).dNp + ...
        shiftdown(T.coord, dNp, T2.father-1, T2.coord, I2.S.', k);
end

%% Response at leaf level
l = depth+1;
I = intdata(l);
Hp = recover(tree(l).coord, Sig(l).dNp, fatherrec-1, nodrec, I.S.', I.W, k);

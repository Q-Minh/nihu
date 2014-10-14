function faces = get_faces(elements)
%GET_FACES Extract all faces of a FE mesh
%   FACES = GET_FACES(ELEMENTS) extracts all the faces of a fe mesh given
%   by the matrix model.ELEMENTS.
%   The structure of the matrix ELEMENTS is
%   [ElemID ElemType MatID PropID NodID1 NODID2 ... NodIDn]
%   where ElemType can be 23, 24, 34, 36 or 38 for standard elements and
%   122, 133 or 134 for infinite elements.
%   The structure of the matrix FACES is
%   [ElemID FaceType NodID1 ... NodIDn]
%   where ElemID refers to the parent element of the face.
%   FaceType can be 12, 23 or 24

%   Copyright 2008-2012 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% zero padding
elements = [elements zeros(size(elements,1), 12-size(elements,2))];

line = elements(elements(:,2) == 12,[1 (5:6)]);
tria = elements(elements(:,2) == 23,[1 (5:7)]);
quad = elements((elements(:,2) == 24) | (elements(:,2) == 122),[1 (5:8)]);
tetra = elements(elements(:,2) == 34,[1 (5:8)]);
penta = elements((elements(:,2) == 36) | (elements(:,2) == 133),[1 (5:10)]);
hexa = elements((elements(:,2) == 38) | (elements(:,2) == 134),[1 (5:12)]);

linefaces = [
    0 1 
    0 2 
    ].';
triafaces = [
    0 1 2 
    0 2 3 
    0 3 1 
    ].';
quadfaces = [
    0 1 2 
    0 2 3 
    0 3 4 
    0 4 1 
    ].';
tetrafaces = [
    0 1 3 2
    0 2 3 4
    0 1 4 3
    0 1 2 4
    ].';
pentatriafaces = [
    0 1 2 3
    0 4 6 5
    ].';
pentaquadfaces = [
 %   0 1 2 5 4
 %   0 2 3 6 5
 %   0 3 1 4 6
    0 1 4 5 2
    0 2 5 6 3
    0 3 6 4 1
    ].';
hexafaces = [
%    0 1 4 3 2
%    0 5 6 7 8
%    0 1 2 6 5
%    0 2 3 7 6
%    0 3 4 8 7
%    0 4 1 5 8
    0 1 2 3 4
    0 5 8 7 6
    0 1 5 6 2
    0 2 6 7 3
    0 3 7 8 4
    0 4 8 5 1

    ].';

linefac = reshape(line(:,linefaces(:)+1).',1+1,2*size(line,1)).';
triafac = reshape(tria(:,triafaces(:)+1).',2+1,3*size(tria,1)).';
quadfac = reshape(quad(:,quadfaces(:)+1).',2+1,4*size(quad,1)).';
tetrafac = reshape(tetra(:,tetrafaces(:)+1).',3+1,4*size(tetra,1)).';
pentatriafac = reshape(penta(:,pentatriafaces(:)+1).',3+1,2*size(penta,1)).';
pentaquadfac = reshape(penta(:,pentaquadfaces(:)+1).',4+1,3*size(penta,1)).';
hexafac = reshape(hexa(:,hexafaces(:)+1).',4+1,6*size(hexa,1)).';
nl = size(linefac,1);
ntr = size(triafac,1);
nq = size(quadfac,1);
nte = size(tetrafac,1);
np3 = size(pentatriafac,1);
np4 = size(pentaquadfac,1);
nh = size(hexafac,1);

faces = zeros(nl+ntr+nq+nte+np3+np4+nh,6);
faces(1:nl,[1 (3)]) = linefac;
faces(nl+(1:ntr),[1 (3:4)]) = triafac;
faces(nl+ntr+(1:nq),[1 (3:4)]) = quadfac;
faces(nl+ntr+nq+(1:nte),[1 (3:5)]) = tetrafac;
faces(nl+ntr+nq+nte+(1:np3),[1 (3:5)]) = pentatriafac;
faces(nl+ntr+nq+nte+np3+(1:np4),[1 (3:6)]) = pentaquadfac;
faces(nl+ntr+nq+nte+np3+np4+(1:nh),[1 (3:6)]) = hexafac;
faces(faces(:,6) ~= 0, 2) = 24;
faces(faces(:,6) == 0, 2) = 23;
faces(faces(:,5) == 0, 2) = 12;
faces(faces(:,4) == 0, 2) = 1;

faces = faces(:,any(faces, 1));

end

function [f ind] = get_free_faces(faces)
%GET_FREE_FACES  Extract free (boundary) faces from FE faces
%   The structure of the FACES matrix is
%   [ElemID FaceType NodID1 ... NodIDn]
%   where ElemID refers to the parent element of the face.

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

f2 = faces(:,3:end);
ff = sort(f2,2);
[fff, ind] = sortrows(ff);
df1 = find(sum(abs(diff(fff,[],1)),2) == 0)+1;
df2 = find(sum(flipud(abs(diff(flipud(fff),[],1))),2) == 0);
ind2 = setdiff(1:size(faces,1), union(df1,df2));
ind = ind(ind2);
f = faces(ind,:);

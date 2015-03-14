function [f, ind] = get_free_faces(faces)
%GET_FREE_FACES  Extract free (boundary) faces from FE faces
%   The structure of the FACES matrix is
%   [ElemID FaceType NodID1 ... NodIDn]
%   where ElemID refers to the parent element of the face.

%   Copyright 2008-2015 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% last modified (acceleration) 2015.03.14

f2 = faces(:,3:end);
ff = sort(f2,2);
[fff, ind] = sortrows(ff);
df1 = find( all(diff(fff,[],1) == 0, 2) ) +1;
df2 = find( all( flipud(diff(flipud(fff),[],1)) == 0, 2) );
ind2 = setdiff(1:size(faces,1), union(df1,df2));
f = faces(ind(ind2),:);

end

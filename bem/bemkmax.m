function k = bemkmax(model, ratio)
%BEMKMAX Maximal wave numbers of a surface mesh
%   k = bemkmax(model, ratio) returns the maximal wave number for each
%   element of the mesh given by model.
% Parameters:
%   model : acoufem mesh structure
%   ratio : number of elements per wavelength (usually set to 6-8)
%   k     : nEx1 vector with wave numbers for each element

% Peter Fiala
% 2012

%%
if nargin < 2
    ratio = 8;
end

%% preproc
elem = drop_IDs(model);

%% process Quad elements
iQuad = find(elem(:,2) == 24);
if ~isempty(iQuad)
    quad = elem(iQuad,5:8);
    d1 = model.Nodes(quad(:,2),2:4) - model.Nodes(quad(:,1),2:4);
    d2 = model.Nodes(quad(:,4),2:4) - model.Nodes(quad(:,1),2:4);
    Avec = cross(d1,d2);
    d = sqrt(sqrt(dot(Avec,Avec,2)));
    k(iQuad) = (2*pi)./(ratio*d);
end

%% process TRIA elements
iTria = find(elem(:,2) == 23);
if ~isempty(iTria)
    tria = elem(iTria,5:7);
    d1 = model.Nodes(tria(:,2),2:4) - model.Nodes(tria(:,1),2:4);
    d2 = model.Nodes(tria(:,3),2:4) - model.Nodes(tria(:,1),2:4);
    Avec = cross(d1,d2);
    d = sqrt(sqrt(dot(Avec,Avec,2)));
    k(iTria) = (2*pi)./(ratio*d);
end
end

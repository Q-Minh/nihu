function k = mesh_kmax(model, ratio)
%MESH_KMAX Maximal wave numbers of a NiHu mesh
%   k = mesh_kmax(model, ratio) returns the maximal wave number for each
%   element of the mesh given by model.
% Parameters:
%   model : acoufem mesh structure
%   ratio : number of elements per wavelength (usually set to 6-8)
%   k     : nEx1 vector with wave numbers for each element

%   Copyright 2009-2010 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last updated: 06.09.2010.

%% preproc
elem = drop_IDs(model);
k = zeros(size(elem,1),1);

%% process Line elements and infinite elements
iLine = find(elem(:,2) == 12 | elem(:,2) == 122);
if ~isempty(iLine)
    line = elem(iLine,5:6);
    dr = model.Nodes(line(:,2),2:4) - model.Nodes(line(:,5),2:4);
    d = sqrt(dot(dr,dr,2));
    k(iLine) = (2*pi)./(ratio*d);
end

%% process TRIA elements
iTria = find(elem(:,2) == 23 | elem(:,2) == 136);
if ~isempty(iTria)
    tria = elem(iTria,5:7);
    d1 = model.Nodes(tria(:,2),2:4) - model.Nodes(tria(:,1),2:4);
    d2 = model.Nodes(tria(:,3),2:4) - model.Nodes(tria(:,1),2:4);
    Avec = cross(d1,d2);
    d = sqrt(sqrt(dot(Avec,Avec,2)));
    k(iTria) = (2*pi)./(ratio*d);
end

%% process Quad elements
iQuad = find(elem(:,2) == 24 | elem(:,2) == 138);
if ~isempty(iQuad)
    quad = elem(iQuad,5:8);
    d1 = model.Nodes(quad(:,2),2:4) - model.Nodes(quad(:,1),2:4);
    d2 = model.Nodes(quad(:,4),2:4) - model.Nodes(quad(:,1),2:4);
    Avec = cross(d1,d2);
    d = sqrt(sqrt(dot(Avec,Avec,2)));
    k(iQuad) = (2*pi)./(ratio*d);
end

%% process Tetra elements
iTetra = find(elem(:,2) == 34);
if ~isempty(iTetra)
    tetra = elem(iTetra,5:8);
    d1 = model.Nodes(tetra(:,2),2:4) - model.Nodes(tetra(:,1),2:4);
    d2 = model.Nodes(tetra(:,3),2:4) - model.Nodes(tetra(:,1),2:4);
    d3 = model.Nodes(tetra(:,4),2:4) - model.Nodes(tetra(:,1),2:4);
    Vvec = dot(cross(d1,d2,2),d3,2);
    d = (Vvec.^2).^(1/6);
    k(iTetra) = (2*pi)./(ratio*d);
end

%% process Penta elements
iPenta = find(elem(:,2) == 36);
if ~isempty(iPenta)
    penta = elem(iPenta,5:10);
    d1 = model.Nodes(penta(:,2),2:4) - model.Nodes(penta(:,1),2:4);
    d2 = model.Nodes(penta(:,3),2:4) - model.Nodes(penta(:,1),2:4);
    d3 = model.Nodes(penta(:,4),2:4) - model.Nodes(penta(:,1),2:4);
    Vvec = dot(cross(d1,d2,2),d3,2);
    d = (Vvec.^2).^(1/6);
    k(iPenta) = (2*pi)./(ratio*d);
end

%% process Hexa elements
iHexa = find(elem(:,2) == 38);
if ~isempty(iHexa)
    hexa = elem(iHexa,5:12);
    d1 = model.Nodes(hexa(:,2),2:4) - model.Nodes(hexa(:,1),2:4);
    d2 = model.Nodes(hexa(:,4),2:4) - model.Nodes(hexa(:,1),2:4);
    d3 = model.Nodes(hexa(:,5),2:4) - model.Nodes(hexa(:,1),2:4);
    Vvec = dot(cross(d1,d2,2),d3,2);
    d = (Vvec.^2).^(1/6);
    k(iHexa) = (2*pi)./(ratio*d);
end

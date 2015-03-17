function v = mesh_volume(model)
% Returns vector v, the approximate volume for each element
% TODO: NOT TESTED YET
%% preproc
elem = drop_IDs(model);
lsetid = elem(:,2);
v = zeros(size(lsetid));

%% process Line elements
iLine = lsetid == ShapeSet.LinearLine.Id;
if any(iLine)
    line = elem(iLine,5:6);
    dr = model.Nodes(line(:,2),2:4) - model.Nodes(line(:,5),2:4);
    d = sqrt(dot(dr,dr,2));
    v(iLine) = d;
end

%% process TRIA elements
iTria = lsetid == ShapeSet.LinearTria.Id;
if any(iTria)
    tria = elem(iTria,5:7);
    d1 = model.Nodes(tria(:,2),2:4) - model.Nodes(tria(:,1),2:4);
    d2 = model.Nodes(tria(:,3),2:4) - model.Nodes(tria(:,1),2:4);
    Avec = cross(d1,d2);
    d = sqrt(dot(Avec,Avec,2))/2;
    v(iTria) = d;
end

%% process Quad elements
iQuad = lsetid == ShapeSet.LinearQuad.Id;
if any(iQuad)
    quad = elem(iQuad,5:8);
    d1 = model.Nodes(quad(:,2),2:4) - model.Nodes(quad(:,1),2:4);
    d2 = model.Nodes(quad(:,4),2:4) - model.Nodes(quad(:,1),2:4);
    Avec = cross(d1,d2);
    d = sqrt(dot(Avec,Avec,2));
    v(iQuad) = d;
end

%% process Tetra elements
iTetra = lsetid == ShapeSet.LinearTetra.Id;
if any(iTetra)
    tetra = elem(iTetra,5:8);
    d1 = model.Nodes(tetra(:,2),2:4) - model.Nodes(tetra(:,1),2:4);
    d2 = model.Nodes(tetra(:,3),2:4) - model.Nodes(tetra(:,1),2:4);
    d3 = model.Nodes(tetra(:,4),2:4) - model.Nodes(tetra(:,1),2:4);
    Vvec = dot(cross(d1,d2,2),d3,2);
    d = sqrt(Vvec.^2)*1/6;
    v(iTetra) = d;
end

%% process Penta elements
iPenta = lsetid == ShapeSet.LinearPenta.Id;
if any(iPenta)
    penta = elem(iPenta,5:10);
    d1 = model.Nodes(penta(:,2),2:4) - model.Nodes(penta(:,1),2:4);
    d2 = model.Nodes(penta(:,3),2:4) - model.Nodes(penta(:,1),2:4);
    d3 = model.Nodes(penta(:,4),2:4) - model.Nodes(penta(:,1),2:4);
    Vvec = dot(cross(d1,d2,2),d3,2);
    d = sqrt(Vvec.^2)*1/2;
    v(iPenta) = d;
end

%% process Hexa elements
iHexa = lsetid == ShapeSet.LinearHexa.Id;
if any(iHexa)
    hexa = elem(iHexa,5:12);
    d1 = model.Nodes(hexa(:,2),2:4) - model.Nodes(hexa(:,1),2:4);
    d2 = model.Nodes(hexa(:,4),2:4) - model.Nodes(hexa(:,1),2:4);
    d3 = model.Nodes(hexa(:,5),2:4) - model.Nodes(hexa(:,1),2:4);
    Vvec = dot(cross(d1,d2,2),d3,2);
    d = sqrt(Vvec.^2);
    v(iHexa) = d;
end

end

function Model = SInitModel(Type, Model)

switch Type
    case '2.5DCartFin'
        Model.Boundary = get_boundary(Model.Domain);
        Model.PlotDomain = extrude(Model.Domain, [0 0 Model.zVector(2)-Model.zVector(1)], length(Model.zVector)-1);
        Model.PlotBoundary = get_boundary(Model.PlotDomain);
    case '2.5DCartInf'
        Model.Boundary = get_boundary(Model.Domain);
        Model.PlotDomain = extrude(Model.Domain, [0 0 Model.zVector(2)-Model.zVector(1)], length(Model.zVector)-1);
        Model.PlotDomain.Nodes(:,4) = Model.PlotDomain.Nodes(:,4) + Model.zVector(1);
        Model.PlotBoundary = extrude(Model.Boundary, [0 0 Model.zVector(2)-Model.zVector(1)], length(Model.zVector)-1);
        Model.PlotBoundary.Nodes(:,4) = Model.PlotBoundary.Nodes(:,4) + Model.zVector(1);
    case '2.5DCylFull'
        coords = Model.Domain.Nodes(:,2:4);
        nonzeroind = find(sum(coords.^2,1) ~= 0);
        Model.Domain.Nodes(:,[2 4]) = coords(:,nonzeroind);
        Model.Domain.Nodes(:,3) = 0;
        Model.Boundary = get_boundary(Model.Domain);
        Model.PlotDomain = revolve(Model.Domain, [0 0 0], [0 0 1], Model.phiVector(2)-Model.phiVector(1), length(Model.phiVector)-1);
        Model.PlotDomain = merge_coincident_nodes(Model.PlotDomain);
        Model.PlotBoundary = get_boundary(Model.PlotDomain);
    case '2.5DCylFin'
        coords = Model.Domain.Nodes(:,2:4);
        nonzeroind = find(sum(coords.^2,1) ~= 0);
        Model.Domain.Nodes(:,[2 4]) = coords(:,nonzeroind);
        Model.Domain.Nodes(:,3) = 0;
        Model.Boundary = get_boundary(Model.Domain);
        Model.PlotDomain = revolve(Model.Domain, [0 0 0], [0 0 1], Model.phiVector(2)-Model.phiVector(1), length(Model.phiVector)-1);
        Model.PlotBoundary = get_boundary(Model.PlotDomain);
    case '2DCyl'
        coords = Model.Domain.Nodes(:,2:4);
        nonzeroind = find(sum(coords.^2,1) ~= 0);
        Model.Domain.Nodes(:,[2 4]) = coords(:,nonzeroind);
        Model.Domain.Nodes(:,3) = 0;
        Model.Boundary = get_boundary(Model.Domain);
        Model.PlotDomain = Model.Domain;
        Model.PlotBoundary = Model.Boundary;
    otherwise %% 3DCart, 2DCart
        Model.Boundary = get_boundary(Model.Domain);
        Model.PlotDomain = Model.Domain;
        Model.PlotBoundary = Model.Boundary;
end
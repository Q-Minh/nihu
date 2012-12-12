function data = SNormalVelocityExcitation(data, selector, expression)

boundary = data.Model.PlotBoundary;

X = boundary.Nodes(:,2);
Y = boundary.Nodes(:,3);
Z = boundary.Nodes(:,4);
[nodeind] = mesh_select(boundary, selector, 'ind');
x = X(nodeind);
y = Y(nodeind);
z = Z(nodeind);
r = sqrt(x.^2+y.^2);
phi = atan2(y, x);
vhat = zeros(size(X));
vhat(nodeind) = eval(expression);

E = data.Excitation;
if isfield(E, 'NormalVelocity')
    nV = length(E.NormalVelocity);
else
    E.NormalVelocity = [];
    nV = 0;
end
E.NormalVelocity(nV+1).ElemSelector = selector;
E.NormalVelocity(nV+1).vExpression = expression;
E.NormalVelocity(nV+1).vhat = vhat;

data.Excitation = E;

%% Acoustic response in a car compartment
% This tutorial shows how an internal sound radiation problem is solved in
% a modal space. The example considered is sound radiation into a car
% compartment by vibrating walls.

%% Geometry and mesh
% The car mesh is loaded and plotted.
load data/carinterior mesh

figure('Name', 'Car compartment interior mesh');
plot_mesh(mesh);

%%
% The finite element mass and stiffness matrices are computed by the
% toolbox function <matlab:doc('model_mk') model_mk>. 
[M, K] = model_mk(mesh);

figure('Name', 'Structure of mass and stiffness matrices');
subplot(1,2,1);
spy(M);
title('Structure of mass matrix');
subplot(1,2,2);
spy(K);
title('Structure of stiffness matrix');

%% Acoustic modes
% The acoustic mode shapes with rigid wall boundary condition are
% determined by the toolbox function <matlab:doc('fe_modes') fe_modes>.
numModes = 100;                       % max. number of modes to determine
[Phi, Om] = fe_modes(M, K, numModes);

figure('Name', 'Mode frequencies');
plot(Om/2/pi);
title('Frequencies of the first acoustic modes');
xlabel('Mode number');
ylabel('Frequency [Hz]');

figure('Name', 'Mode shape');
n = 100;
plot_mesh(mesh, 'node', mesh.Nodes(:,1), Phi(:,n));
title(sprintf('Acoustic pressure mode #%d at %.2f Hz', n, Om(n)/2/pi));
view(3);

%% Excitation
boundary = get_boundary(mesh);
nNodes = size(mesh.Nodes,1);
q = zeros(nNodes,1);
botsel = sprintf('z <= %g', min(mesh.Nodes(:,4))+2e-2);
[nbot, ebot] = mesh_select(boundary, botsel, 'ind');
q(nbot) = 1;

figure('Name', 'Excitation');
plot_mesh(boundary, 'node', mesh.Nodes(:,1), q);
view([-135 -30]);

boundary.Elements = boundary.Elements(ebot,:);
A = model_a(boundary);

%% Modal solution
freq = 300;
om = 2*pi*freq;
alp = (Om.^2 - om^2) .\ (Phi.'*A*q);
p = Phi * alp;

plot_mesh(mesh, 'node', p);

%% Postprocessing
field = field_bbx(boundingbox(mesh), mesh.Materials(1,4)/freq/10);

Np = fe_interp(mesh, field.Nodes(:,2:4));

figure('Name', 'Response')
plot_mesh(mesh), alpha 0;
plot_mesh(field, 'node', Np * p);
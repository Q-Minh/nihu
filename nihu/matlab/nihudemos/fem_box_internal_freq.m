%% Acoustic internal radiation problem with the finite element method
% This demo solves an acoustic internal radiation problem with the finite
% element method.

%% Geometrical parameters and material properties
% Geometrical parameters and material properties of the problem are defined
% as follows. The element length of the finite element mesh is determined
% using the rule of thumb $L_e = \lambda / 8$, where $\lambda$ denotes the
% minimal wave length.
L = [5 4 3];            % dimensions of the box [m]
f_max = 300;            % maximal frequency
rho = 1.125;            % mass density of air [kg/m3]
c_sound = 340;          % speed of sound
Le = c_sound/f_max/8;   % element length

%% Create mesh
% The finite element mesh is created using the toolbox function
% <matlab:doc('create_brick') create_brick>. The material properties are
% incorporated into the model structure by hand.
box = create_brick(L, ceil(L/Le));
box.Materials(1,3:4) = [rho c_sound];

%% Compute finite element system matrices
% The finite element mass and stiffness matrices (${\bf M}, {\bf K}$) are
% computed using the toolbox function <matlab:doc('model_mk') model_mk>.
% The excitation surface matrix $\bf A$ is computed using the toolbox function
% <matlab:doc('model_a') model_a>.
[vI, vJ, M, K, vDOF] = model_mk(box, 'ind');
bou = get_boundary(box);
[sI, sJ, A, sDOF] = model_a(bou, 'ind');
sI = reindex(sI, sDOF, vDOF);
sJ = reindex(sJ, sDOF, vDOF);

%% create sparse matrices
nNodes = size(box.Nodes,1);
A = sparse(sI, sJ, A, nNodes, nNodes);
M = sparse(vI, vJ, M, nNodes, nNodes);
K = sparse(vI, vJ, K, nNodes, nNodes);

%% Excitation
% The excitation is unit normal velocity on a single wall of the box and
% zero normal velocity on other walls.
% The nodes of the selected wall are selected using the toolbox function
% <matlab:doc('mesh_select') mesh_select>.
v = zeros(nNodes,1);
top_ind = mesh_select(box, sprintf('z > %g-1e-3', L(3)), 'ind');
v(top_ind) = 1;

figure;
title('Velocity excitation');
plot_mesh(box, v);

%%
% The algebraic system is solved on a given frequency by a direct solver
f = f_max * 1;     % frequency
om = 2*pi*f;        % angular frequency

% direct solution
p = (K - om^2 * M) \ (1i*om*A*v);

%%
% The response is plotted as follows:
figure;
title('Pressure response');
plot_mesh(box, abs(p));

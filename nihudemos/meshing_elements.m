%% Basic mesh elements
% This tutorial describes how to create simple meshes using the toolbox
% |NIHU|.

%% Create a straight line
% Simple straight lines can be created using the toolbox function
% <matlab:doc('create_line') create_line>. There are three different ways of meshing straight lines:
%
% * Defining its length $L$ and the number of elements $n_E$
% * Defining its starting and ending positions ${\bf x}_1$ and ${\bf x}_2$,
% as well as the number of elements $n_E$
% * Defining the internal nodes ${\bf x}_i$ of the straight line.
%
% In the first example, the first method is applied.
%
L = 5;              % length of the line
Le = 0.3;           % approximate length of a line element
N = ceil(L/Le);     % number of elements
lin = create_line(L, N);

figure;
plot_mesh(lin);

%%
% The second method can be applied as follows:
R = [                               % start and end locations
    1 0 0
    2 3 1
    ];
Le = 0.3;                           % approx. length of a line element
N = ceil(norm(diff(R,[],1))/Le);    % number of elements
lin = create_line(R, N);

figure;
plot_mesh(lin); view(3);

%%
% In the third method, the line's nodes ${\bf x}_i$ are directly defined.
% In this case, the line is aligned along the $x$ axis.
Cx = (-1:.1:+2)';          % internal nodes
lin = create_line(Cx);

figure;
plot_mesh(lin);


%% Create slab
% Rectangular slabs can be created using the toolbox function
% <matlab:doc('create_slab') create_slab>. There are three ways of defining slabs:
%
% * Define its side lengths $L_x$ and $L_y$, and the number of elements
% $N_x$ and $N_y$
% * Define its corner nodes ${\bf x}_i$, $i = 1 \dots 4$ and the number of
% elements along the two dimensions.
% * Define the internal nodes $x_i$, $y_j$ directly,
% so that the slab's nodes are elements of the Descartes product $x_i \times y_j$.
%
% In the first case, the slab is located at the origin, its sides are
% aligned along the $x$ and $y$ axes.
L = [4 3];                  % side lengths
Le = .3;                    % approximate side length of an element
N = ceil(L/Le);             % number of elements along the two sides
slab = create_slab(L, N);

figure;
plot_mesh(slab);

%%
% The parameters |L| and |N| can be defined as scalars, so that the same
% meshing parameters are used in the two directions.
slab = create_slab(4, 7);

figure;
plot_mesh(slab);

%%
% In the second method, the four corner nodes are defined in the 4x3 matrix
% |C|.
C = [                       % corner coordinates
    0  0 0
    1  0 0
    .7  1 1
    .2 .7 0
    ];
N = [10 7];                 % division parameters
slab = create_slab(C, N);
figure;
plot_mesh(slab); view(3);

%%
% With the third method, non-uniformly meshed slabs can be easily defined
% as follows:
Cx = logspace(1, 2, 30).';    % internal nodes along the x axis
Cy = linspace(10, 30, 10).';  % internal nodes along the y axis
slab = create_slab(Cx, Cy);

figure;
plot_mesh(slab);

%%
% If only |Cx| is defined, the same nodes are applied for both dimensions:
Cx = logspace(1, 2, 20);    % internal nodes along the x and y axes
slab = create_slab(Cx);

figure;
plot_mesh(slab);

%% Create circle
% Circle quadrants can be meshed by using the function
% <matlab:doc('create_circle_quadrant') create_circle_quadrant>. The necessary parameters are the circle's
% radius $R$ and a division parameter $n_R$. The division parameter defines
% the number of elements along the radius of the circle quadrant,
% the same as the number of elements along the perimeter. The circle
% quadrant is located in the $xy$ plane, the center is at the origin, and
% the sides are aligned along the $x$ and $y$ axes.
R = 2;      % radius
nR = 10;    % division parameter
circlequad = create_circle_quadrant(R, nR);

figure;
plot_mesh(circlequad);

%%
% Circles can be meshed by using the function <matlab:doc('create_circle') create_circle>. The
% necessary parameters are the circle's radius $R$ and a division
% parameter $n_R$. The circle is located in the $xy$ plane, the center being
% at the origin.  The number of elements along the perimeter is
% $4 n_R$.
R = 2;      % radius
nR = 10;    % division parameter
circle = create_circle(R, nR);

figure;
plot_mesh(circle);

%% Create brick volume and brick surface
% Bricks volume meshes can be created using the toolbox function
% <matlab:doc('create_brick') create_brick>.
% There are three ways of defining brick volumes:
%
% * Define its side lengths $[L_x, L_y, L_z]$ and division paramters $[N_x N_y N_z]$
% * Define its corner nodes ${\bf x}_i$, $i = 1 \dots 8$,
% * Define its internal nodes $x_i$, $y_j$, $z_k$ along the three dimensions
%
% The first method is demonstrated first, where the brick is located in the
% origin, with sides parallel to the coordinate axes:
L = [5 4 3];        % Brick dimensions [Lx Ly Lz]
Le = .3;            % element length
N = ceil(L / Le);   % Division parameter [Nx Ny Nz]
brick = create_brick(L, N);

figure;
plot_mesh(brick); view(3);

%%
% The second method makes it possible to define skewed bricks, by defining
% the corner nodes
C = [           % 8 corner nodes of the brick
    0   0  0
    1   0  0
    1.5 1  0
    -.1 .9 .1
    0   0  1
    1   0  1.2
    1.5 1  1
    -.1 .7 1
    ];
brick = create_brick(C, [10, 10, 10]);

figure;
plot_mesh(brick); view(3);

%%
% With the third method, non-uniform sampling can easily be defined:
Cx = (.1:.05:1).^2; % nodes along the x axis
Cy = (.1:.05:1).^3; % nodes along the y axis
Cz = (.1:.05:1).^4; % nodes along the z axis
brick = create_brick(Cx, Cy, Cz);

figure;
plot_mesh(brick); view(3);

%%
% Brick boundaries can be meshed using the toolbox function
% <matlab:doc('create_brick_boundary') create_brick_boundary>, by defining the brick dimensions and division
% parameters.
L = [5 4 3];        % Brick dimensions [Lx Ly Lz]
Le = .3;            % element length
N = ceil(L / Le);   % Division parameter [Nx Ny Nz]
brick = create_brick_boundary(L, N);

figure;
plot_mesh(brick); view(3);

%% Create sphere and sphere surface
% Spherical volumes can be created using the toolbox function
% <matlab:doc('create_sphere') create_sphere>, by defining the sphere radius $R$ and a division
% parameter $n_R$. The division parameter defines the number of elements
% along the sphere radius. The number of elements along the perimeter is
% $4 n_R$.
R = 2;          % sphere radius
nR = 10;        % number of elements along the radius
sphere = create_sphere(R, nR);

figure;
plot_mesh(sphere); view(3);

%%
% Sphere boundaries can be created using the toolbox function
% <matlab:doc('create_sphere_boundary') create_sphere_boundary>, by defining the sphere radius $R$ and division
% parameter $n_R$.
R = 2;
nR = 10;
sphere = create_sphere_boundary(R, nR);

figure;
plot_mesh(sphere); view(3);

%% Create Cat's eye
% Cat's eye surface models can be created using the function
% <matlab:doc('create_catseye') create_catseye>, by defining the radius $R$ and division parameter $n_R$.
R = 2;
nR = 10;
cat = create_catseye(R, nR);

figure;
plot_mesh(cat); view([110 35]);
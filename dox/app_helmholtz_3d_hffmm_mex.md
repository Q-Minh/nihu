3D Helmholtz FMBEM with Mex interface {#app_helmholtz_3d_hffmm_mex}
=====================================

\page app_helmholtz_3d_hffmm_mex 3D Helmholtz FMBEM with Mex interface

[TOC]

Description {#app_helmholtz_3d_hffmm_mex_desc}
===========

This application implements a fast multipole solver for the Helmholtz equation in 3D with the following features:

- Collocational formulation
- Homogeneous meshes with quadrilateral elements
- High frequency multilevel fast multipole method
- Burton-Miller technique for mitigating fictitious eigenfrequencies
- Applicable for Neumann and Dirichlet problems
- MEX interface for interacting with Matlab - Compatible with Matlab's built-in iterative 
solvers
 
The application is used for assembling the single and double layer potential matrices (SLP and DLP matrices) and then evaluate fast matrix-vector products with these matrices.
The same matrix-vector products can be called from Matlab's built-in iterative solvers, such as BiCGSTAB or GMRES.

Usage {#app_helmholtz_3d_hffmm_mex_usage}
=====

This MEX application can be executed in various steps.
Each step is represented by a command string that is passed to the MEX function as its first parameter.
The corresponding MEX function has an inner state, therefore initialization and cleanup are necessary.
The typical sequence of the steps and the corresponding commands is as follows.

1. `init` -- Initializes the data structure and the storage object.
2. `set` -- Sets up the parameters for the computation, such as the wave number, or the multipole expansion accuracy.
3. `mesh` -- Defines the surface and field meshes. The field mesh is not necessarily given, in this case the surface problem is solved with singular integral evaluation and the Burton-Miller technique. If the field mesh is also given, the standard collocational method is applied for computing the radiated field.
4. `tree` -- Build the cluster tree and initialize the interaction lists of the clusters.
5. `matrix` -- Assemble the FMM matrices for subsequent fast matrix-vector products. This step includes the precomputation of the discretized operators.
6. `mvp_dlp` and `mvp_slp` -- Evaluate fast matrix-vector product by the SLP or the DLP matrix. These steps are called in each iteration of an iterative solver. The same steps are called once for computing the radiated field.
7. `cleanup` -- Free the memory allocated for the internal storage and reset to the initial state.

Arriving at a solution of a complete Helmholtz problem may require solving the surface problem at first, and then computing the radiated field.
In this case, the whole sequence of steps is repeated twice, first the appilaction is used in surface mode, and then in field mode.
The example in the following section also demonstrates this case.

Example {#app_helmholtz_3d_hffmm_mex_app}
=======

This example demonstrates the usage of the application through a transparent problem.
The field of a point source is computed by prescribing the normal derivative field of the source on a surface mesh as a Neumann boundary condition.
The solution is performed in two steps: first the field on the surface is attained by solving the BEM system of equations using the fast multipole method, and then, the field in the field points is evaluated by means of fast matrix-vector products.
The numerical solution of the transparent problem is finally compared to the analytical solution.

As the first step, the surface and field meshes are created and the wavenumber is set.
The following code snippet creates the sphere surface and a fan-shaped field point mesh next to the sphere.
For the sake of convenience the MEX function to call for each operation is stored in the variable `mex_fun`.

```Matlab
% Name of the MEX function to call
mex_fun = @helmholtz_3d_hffmm_mex;
% Sphere surface with radius R = 1, divided into 30 elements along the radius.
surface = create_sphere_boundary(1, 30);
% Field point mesh: fan shaped, a line revolved around the z axis
field = revolve_mesh( ...
    create_line([1.125, 0, 0; 5.625, 0, 0], 180), pi/200, 100, [0 0 1]);
% Wave number, allow minimum 10 elems / wavelength
k = min(mesh_kmax(surface, 10));
```

Then, the MEX function is initialized and the parameters are set.

```Matlab
mex_fun('init');
mex_fun('set', ...
    'accuracy', 3.0, ...            % Expansion length parameter
    'far_field_order', 5, ...       % Far field quadrature order
    'wave_number', k);              % Wave number
```

The surface mesh is tranferred to the MEX function, and the cluster tree is built in the next steps.
The function `extract_core_mesh` is applied for extracting the mesh in the format accepted by the MEX function.
Each cluster is divided into children clusters in this case as long as the diameter of the cluster is greater than the reciprocal of the wave number, i.e., \f$ D > 1/k \f$.
When the tree is built, the command `print_tree` is utilized to display statistics of the tree.

```Matlab
[s_nodes, s_elems] = extract_core_mesh(surface, 'surface');
mex_fun('mesh', s_nodes, s_elems);
mex_fun('tree', 'divide_diameter', 1/k);
mex_fun('print_tree');
```

The final step of preparation is the assembly of system matrices.
This step also contains the precomputation of some operators that will accelerate the evaluation of matrix--vector products.
The command `matrix` assembles the single layer potential (SLP) and double layer potential (DLP) matrices.
Assembly times are displayed by invoking the command `print_times`.
If parallelization is enabled (as is by default), parallel evaluation is already exploited in the matrix assembly procedure.

```Matlab
mex_fun('matrix');                 % Assemble FMM matrices
mex_fun('print_times');            % Display assembly times
```

The surface problem is solved in the following sequence of steps.
Namely, the solution of the matrix equation 

\f$
	\displaystyle \mathbf{M}_{\mathrm{s}} \mathbf{p}_{\mathrm{s}} = \mathbf{L}_{\mathrm{s}} \mathbf{q}_{\mathrm{s}}
\f$

where the matrices \f$ \mathbf{M}_{\mathrm{s}} \f$ and \f$ \mathbf{L}_{\mathrm{s}} \f$ are the surface DLP and SLP matrices, respectively.
The complex vectors \f$ \mathbf{p}_{\mathrm{s}} \f$ and \f$ \mathbf{q}_{\mathrm{s}} \f$ hold the pressure and its normal derivative at the collocational points.
In this case the Neumann problem is solved, i.e., \f$ \mathbf{q}_{\mathrm{s}} \f$ is known, and \f$ \mathbf{p}_{\mathrm{s}} \f$ is sought.
As the first step of the solution process, the right-hand-side is evaluated by means of a fast matrix--vector product by the SLP matrix.
\f$ \mathbf{q}_{\mathrm{s}} \f$ is defined by the field of a virtual monopole source located inside the volume bounded by the spherical surface.

The function `centnorm` is utilized for extracting element centers (locations of collocational points) and normals at the center.
Then, the function `incident` is used for evaluating the field of the virtual monopole source.
Once the vector \f$ \mathbf{q}_{\mathrm{s}} \f$ is assembled, the command `mvp_slp` is passed to the MEX function to get the right-hand-side vector.

```Matlab
x_src = [.2 .3 .6];
[centers, normals] = centnorm(surface);
[ps_ana, qs] = incident('point', x_src, centers, normals, k);
rhs = mex_fun('mvp_slp', qs);
```

To solve the surface problem, Matlab's built-in GMRES iterative solver is utilized.
In the solution process, matrix--vector product by the DLP matrix are evaluated.
Therefore, the MEX function with the command `mvp_dlp` is passed to the `gmres` method.
The helper function `ensure_comples` is needed to ensure that input vector is always interpreted as a complex valued vector.
The solution vector \f$ \mathbf{p}_{\mathrm{s}} \f$ (or `ps_bem`) in the code is attained in the following code lines.

```Matlab
Afun = @(x)mex_fun('mvp_dlp', ensure_complex(x));
[ps_bem, flag, relres, iter, resvec] = gmres(Afun, rhs, [], 1e-8, 1000);
```

Note that the Dirichlet problem can be solved in a similar manner, with interchanging the role of the SLP and DLP matrices.

In order to check that the FMBEM solution is correct, the numerical and analytical surface results are compared.
The relative error of the numerical solution \f$ \varepsilon \f$ is evaluated as

\f$
	\displaystyle \varepsilon = \frac{||\mathbf{p}_{\mathrm{s}}^{\mathrm{(BEM)}} - \mathbf{p}_{\mathrm{s}}^{\mathrm{(ana)}}||}{|| \mathbf{p}_{\mathrm{s}}^{\mathrm{(ana)}}||}
\f$

```Matlab
err = norm(ps_ana-ps_bem)/norm(ps_ana);
fprintf('Relative error on surface: %.2f %%\n', err*100);
```

Once the surface computation is ready, 

% NiHu Toolbox BEM functions
%
% BEM related meshing functions
%   burtmillgen          - Generate modified Burton and Miller points
%   chiefgen             - Generate random CHIEF points to a bem mesh
%   bemkmax              - Maximal wave numbers of a surface mesh
%   get_spheric_angles   - Outward spheric angles of a surface mesh
%
% Conventional Helmholtz BEM
%   bemHG                - Build matrices of an acoustic bem model
%   bemHG_const          - Generate acoustic BEM system matrices
%   bemHG_lin            - Generate acoustic BEM system matrices
%   bemhg_lin_m          - Generate acoustic BEM system matrices
%
% Fast multipole HelmholtzBEM
%   bemHG_const_sp       - Generate sparse acoustic BEM system matrices
%   clustertree          - Build cluster tree of a fast multipole bem model
%   fm_postproc          - bem post processing with fast multipole algorithm
%   gmres_iter_dirichlet - One GMRES iteration of a Neumann problem
%   gmres_iter_neumann   - One GMRES iteration of a Neumann problem
%   integpar             - Compute integration parameters of a FM tree
%   leafGq               - leafGq
%   leafHp               - leafHp
%   mpcont_Gq            - Multipole contribution of the Gq integral operator
%   mpcont_Hp            - Multipole contribution of the Hp integral operator
%   nfij                 - Generate near field sparse matrix locations of a FMBEM model
%   recover              - RECOVER
%   reltree              - Convert cluster tree to relative tree
%   shiftdown            - SHIFTDOWN
%   shiftup              - SHIFTUP
%   spherequad           - Gaussian quadrature over the unit sphere
%   sphresamp            - Resampling a function over the unit sphere
%   spht                 - Spherical harmonic transform
%   ispht                - Inverse spherical harmonic transform
%   trans                - TRANS
%   translation          - Compute translation between multipoles
%   plot_clusters        - display one level of the cluster tree
%   print_tree_info      - 
%   rad_neumann          - Compute radiation problem with Neumann BC
%
% Iterative solver
%   mygmres              - GMRES iterative solver for FMBEM
%   my_gmres             - GMRES algorithm with right preconditioner
%   my_gmres_lu          - allocating space for the algorithm
%   fgmres               - 
%   excdecomp            - Decomposition of right hand side matrix for iterative solver
%   saileft              - Sparse approximate inverse
%   sairight             - Sparse approximate right inverse
%   saistruct            - Determine sparsity structure of SAI preconditioner
%
% Other functions
%   bempow               - 
%   incident             - Incident pressure and velocity wave field
%   nihu_manager         - 

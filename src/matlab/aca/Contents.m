% MATLAB
%
% Tree construction
%   build_block_tree        - build block structure from a cluster tree
%   build_cluster_tree      - build cluster tree from a set of coordinates
%   build_dual_block_tree   - build block structure from two cluster trees
%   sort_cluster_tree       - cluster-contiguous node numbering
% ACA
%   lowrank_approx          - low rank approximation of a matrix
%   lowrank_approx_block    - low rank approximation of a matrix block
%   max_aca_compression     - maximal compression of the aca algorithm
% bbFMM
%   chebinterp              - Chebyshev anterpolation matrix
%   chebpoly                - values of Chebyshev polynomials
%   chebroots               - roots of Chebyshev polynomials
%   lintrans                - linear transformation from multi-interval to other
%   bb_tree_cheb_nodes      - Chebyshev nodes of a cluster tree
%   bb_M2L                  - Multipole to Local sparse matrix
% helper functions
%   is_admissible_bb        - 
%   is_admissible_dist      - 
%   helmholtz_kernel        - 
%   helmholtz_kernel_sp     - 
%   laplace_kernel          - 
%   laplace_kernel_sp       - 
% postprocessing
%   display_block_structure - visualise the block tree
% unit tests
%   build_tree              - build cluster tree for bbFMM method
%   leaf_contribution       - 
%   bbFMM_test              - 
%   bbmultilevelFMM_test    - 
%   test_aca_dual_1D        - test 1D ACA with different meshes
%   test_aca_dual_2D        - test 2D ACA with different meshes
%   test_aca_dual_3D        - test 3D ACA with different meshes
%   aca_tester              - general tester of the ACA algorithm
%   test_bb_matrices        - 
%   chebinterp_test         - 1D
%   plot_cluster            - 

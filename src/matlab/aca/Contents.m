% Fast methods based on Clustering
%
% Tree construction
%   build_cluster_tree        - build cluster tree from a set of coordinates
%   sort_cluster_tree         - cluster-contiguous node numbering
%   build_block_tree          - build block structure from a cluster tree
%   build_dual_block_tree     - build block structure from two cluster trees
%   build_regular_clustertree - Build cluster tree
% ACA
%   lowrank_approx            - low rank approximation of a matrix
%   lowrank_approx_block      - low rank approximation of a matrix block
%   max_aca_compression       - maximal compression of the aca algorithm
% bbFMM
%   chebinterp                - Chebyshev anterpolation matrix
%   chebpoly                  - values of Chebyshev polynomials
%   chebroots                 - roots of Chebyshev polynomials
%   lintrans                  - linear transformation from multi-interval to other
%   bb_tree_cheb_nodes        - Chebyshev nodes of a cluster tree
%   bb_M2L                    - Multipole to Local sparse matrix
%   bb_M2M                    - Multipole to Multipole sparse matrix
%   bb_P2M                    - Point to Multipole sparse matrix
%   bb_P2P                    - Point to Point sparse matrix
%   bb_matvec                 - Matrix-vector multiplication with bb FMM
%   bb_P2P_regular            - compute P2P sparse matrix of a bb FMM
%   bb_P2M_regular            - Point to Multipole sparse matrix
%   bb_M2M_regular            - Multipole to Multipole sparse matrix
%   bb_M2L_regular            - compute M2L transfer matrices
%   bb_L2L_regular            - Local to Local sparse matrix
%   compress_M2L              - Compress M2L matrices with ACA
%   bb_far_transfer           - 
%   bb_matvec_regular         - 
%   nfij                      - Generate near field sparse matrix locations of a FMBEM model
%   trans2idx                 - assign index to a translation vector
%   transgen                  - Generate unique translation vectors
% helper functions
%   is_admissible_bb          - 
%   is_admissible_dist        - 
%   helmholtz_kernel          - 
%   helmholtz_kernel_sp       - 
%   laplace_kernel            - 
%   laplace_kernel_sp         - 
% postprocessing
%   plot_cluster              - 
%   display_block_structure   - visualise the block tree
% unit tests
%   build_tree                - build cluster tree for bbFMM method
%   leaf_contribution         - 
%   bbFMM_test                - 
%   bbmultilevelFMM_test      - 
%   test_aca_dual_1D          - test 1D ACA with different meshes
%   test_aca_dual_2D          - test 2D ACA with different meshes
%   test_aca_dual_3D          - test 3D ACA with different meshes
%   aca_tester                - general tester of the ACA algorithm
%   test_bb_matrices          - 
%   test_bb_solution          - 
%   chebinterp_test           - 1D
%   test_regular_cluster_tree - 

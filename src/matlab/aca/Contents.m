% ACA
%
% Files
%   build_cluster_tree      - build cluster tree from a set of coordinates
%   build_block_tree        - build block structure from a cluster tree
%   build_dual_block_tree   - build block structure from two cluster trees
%   sort_cluster_tree       - cluster-contiguous node numbering
% ACA
%   lowrank_approx          - low rank approximation of a matrix
%   lowrank_approx_block    - low rank approximation of a matrix block
%   max_aca_compression     - maximal compression of the aca algorithm
% Post processing
%   display_block_structure - visualise the block tree
%   plot_cluster            - 
% Helper functions
%   is_admissible_bb        - 
%   is_admissible_dist      - 
%   helmholtz_matrix        - 
%   helmholtz_matrix_sp     - 
%   laplace_matrix          - 
%   laplace_matrix_sp       - 
% Unit tests
%   aca_tester              - general tester of the ACA algorithm
%   test_dual_1D            - test 1D ACA with different meshes
%   test_dual_2D            - test 2D ACA with different meshes
%   test_dual_3D            - test 3D ACA with different meshes

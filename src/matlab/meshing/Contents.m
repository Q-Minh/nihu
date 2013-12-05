% MESHING
%
% Files
%   adjacency                  - Compute adjacency matrix of a fe mesh
%   boundingbox                - Compute bounding box of a NiHu mesh
%   centnorm                   - Element centers and normals
%
%   create_brick               - Create a brick volume mesh (NiHu / meshing)
%   create_brick_args          - Process brick mesh arguments. (NiHu / meshing)
%   create_brick_base          - Create a basic brick mesh. (NiHu / meshing)
%   create_brick_boundary      - Create a brick surface mesh (NiHu / meshing)
%   create_brick_boundary_base - Create a basic brick mesh. (NiHu / meshing)
%   create_catseye             - Create a Cat's eye surface mesh (NiHu / meshing)
%   create_circle              - Create a circle surface mesh (NiHu / meshing)
%   create_circle_boundary     - Create a circle boundary mesh (NiHu / meshing)
%   create_circle_quadrant     - Create quadrant of a circle surface mesh (NiHu / meshing)
%   create_empty_mesh          - Create an empty NiHu mesh (NiHu / meshing)
%   create_line                - Create a line mesh (NiHu / meshing)
%   create_line_base           - Create basic line mesh. (NiHu / meshing)
%   create_slab                - Create slab mesh (NiHu / meshing)
%   create_slab_args           - Process slab mesh arguments. (NiHu / meshing)
%   create_slab_base           - Create basic slab mesh. (NiHu / meshing)
%   create_slab_boundary       - Create a slab boundary mesh (NiHu / meshing)
%   create_slab_boundary_base  - Create a basic slab mesh. (NiHu / meshing)
%   create_sphere              - Create a sphere volume mesh (NiHu / meshing)
%   create_sphere_boundary     - Create a sphere surface mesh (NiHu / meshing)
%
%   default_mat_prop           - Default Materials and Properties matrices (NiHu / meshing)
%   drop_IDs                   - Get rid of material, property, element and node IDs
%   drop_mesh_IDs              - Get rid of material, property, element and node IDs
%   drop_unused_nodes          - Exclude unused nodes from a mesh
%
%   extrude_mesh               - Extrude 1D and 2D mesh along a given direction (NiHu / meshing)
%   scale_mesh                 - Scale mesh (NiHu / meshing)
%   reflect_mesh               - Reflect mesh to a symmetry plane (NiHu / meshing)
%   repeat_mesh                - Repeat (Copy) mesh along a given vector (NiHu / meshing)
%   revolve_mesh               - Revolve 1D and 2D mesh around a given vector (NiHu / meshing)
%   rotate_mesh                - Rotate mesh around a given vector (NiHu / meshing)
%   translate_mesh             - Translate mesh along a given vector (NiHu / meshing)
%   join_meshes                - Create a single mesh from several submeshes
%   refine_mesh                - refine mesh by dividing its elements into subelements
%
%   field_bbx                  - Cross-shaped field point mesh based on bounding box data
%   fill_polygon               - Create 2D TRI-mesh from polygon mesh
%   flip_elements              - Flip elements of a NiHu mesh
%   gaussquad                  - Gaussian quadrature integration.
%   gaussquad1                 - 1D Gaussian quadrature integration
%   gaussquad2                 - 2D Gaussian quadrature integration
%   gaussquad3                 - 3D Gaussian quadrature integration
%   geo2gauss                  - Gaussian quadrature over a mesh
%   get_boundary               - Extract domain boundary
%   get_faces                  - Extract all faces of a FE mesh
%   get_free_faces             - Extract free (boundary) faces from FE faces
%   merge_coincident_nodes     - Merge coincident nodes in a mesh
%
%   isinside                   - Check whether points are inside a mesh volume
%   mesh_edge_size             - Estimate edge size of a mesh (NiHu / meshing)
%   mesh_kmax                  - Maximal wave numbers of a NiHu mesh (NiHu / meshing)
%   mesh_section               - Return a rectangular section of a mesh (NiHu / meshing)
%   mesh_select                - Element and node selection
%   mesh_volume                - Returns vector v, the volume for each element
%   meshpath                   - Create line mesh from a path of lines and Bezier Curves
%   node_normals               - calculate outward normal vectors in nodes
%   node_partitions            - 
%   node_volume                - Calculates corresponding volume for nodes
%   plot_elem_normals          - Plot element normals of a mesh
%   plot_mesh                  - Plot colored NiHu mesh
%   plot_node_numbers          - Plot node numbers of a NiHu mesh
%   quad2tria                  - Replace Quad elements with TRIA elements
%   quadratic                  - Solve quadratic polynomial equation ax^2 + bx + c = 0
%   quadratise                 - replace linear elements with quadratic elements
%   read_epspath               - Read curves from .eps file
%   rotation_matrix            - Matrix of rotation around a vector
%   shapefun                   - Shape functions over standard elements
%   split_independent_meshes   - Split mesh to independent submeshes
%   surface2wireframe          - 
%   vert2gauss                 - Gaussian quadrature from vertices and elements

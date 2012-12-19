% Toolbox NiHu Meshing
% Version 1.0 (R2006a) 08-Nov-2009
%
% Creating simple meshes
%   create_empty_mesh        - Create an empty NiHu mesh (NiHu / meshing)
%   create_line              - Create a line mesh (NiHu / meshing)
%   create_slab              - Create slab mesh (NiHu / meshing)
%   create_slab_boundary     - Create a slab boundary mesh (NiHu / meshing)
%   create_brick_boundary    - Create a brick surface mesh (NiHu / meshing)
%   create_sphere_boundary   - Create a sphere surface mesh (NiHu / meshing)
%   create_brick             - Create a brick volume mesh (NiHu / meshing)
%   create_circle_quadrant   - Create quadrant of a circle surface mesh (NiHu / meshing)
%   create_circle            - Create a circle surface mesh (NiHu / meshing)
%   create_sphere            - Create a sphere volume mesh (NiHu / meshing)
%   create_catseye           - Create a Cat's eye surface mesh (NiHu / meshing)
%
% Internal mesh creating functions
%   default_mat_prop         - Default Materials and Properties matrices (NiHu / meshing)
%   create_line_base         - Create basic line mesh. (NiHu / meshing)
%   create_slab_args         - Process slab mesh arguments. (NiHu / meshing)
%   create_slab_base         - Create basic slab mesh. (NiHu / meshing)
%   create_slab_boundary_base - Create a basic slab mesh. (NiHu / meshing)
%   create_brick_args        - Process brick mesh arguments. (NiHu / meshing)
%   create_brick_base        - Create a basic brick mesh. (NiHu / meshing)
%   create_brick_boundary_base - Create a basic brick mesh. (NiHu / meshing)
%
% Basic mesh transforms
%   translate_mesh           - Translate mesh along a given vector (NiHu / meshing)
%   scale_mesh               - Scale mesh (NiHu / meshing)
%   reflect_mesh             - Reflect mesh to a symmetry plane (NiHu / meshing)
%   rotate_mesh              - Rotate mesh around a given vector (NiHu / meshing)
%   repeat_mesh              - Repeat (Copy) mesh along a given vector (NiHu / meshing)
%   extrude_mesh             - Extrude 1D and 2D mesh along a given direction (NiHu / meshing)
%   revolve_mesh             - Revolve 1D and 2D mesh around a given vector (NiHu / meshing)
%
% Mesh manipulation
%   mesh_select              - Element and node selection
%   mesh_section             - Return a rectangular section of a mesh
%   drop_IDs                 - Get rid of material, property, element and node IDs
%   drop_mesh_IDs            - Get rid of material, property, element and node IDs
%   drop_unused_nodes        - Exclude unused nodes from a mesh
%   merge_coincident_nodes   - Merge coincident nodes in a mesh
%   join_meshes              - Create a single mesh from several submeshes
%   split_independent_meshes - Split mesh to independent submeshes
%   flip_elements            - Flip elements of a NiHu mesh
%   get_boundary             - Extract domain boundary
%   get_faces                - Extract all faces of a FE mesh
%   get_free_faces           - Extract free (boundary) faces from FE faces
%   quad2tria                - Replace Quad elements with TRIA elements
%   surface2wireframe        - 
%
% Numerical integration over meshes
%   gaussquad1               - 1D Gaussian quadrature integration
%   gaussquad2               - 2D Gaussian quadrature integration
%   gaussquad3               - 3D Gaussian quadrature integration
%   shapefun                 - Shape functions over standard elements
%   vert2gauss               - Gaussian quadrature from vertices and elements
%   geo2gauss                - Gaussian quadrature over a mesh
%   centnorm                 - Element centers and normals
%
% Mesh import and export
%   read_epspath             - Read curves from .eps file
%   meshpath                 - Create line mesh from a path of lines and Bezier Curves
%
% Other
%   boundingbox              - Compute bounding box of a NiHu mesh
%   isinside                 - Check whether points are inside a mesh volume
%   adjacency                - Compute adjacency matrix of a fe mesh
%   field_bbx                - Cross-shaped field point mesh based on bounding box data
%   fill_polygon             - Create 2D TRI-mesh from polygon mesh
%   gaussquad                - Gaussian quadrature integration.
%   invmap                   - Inverse isoparametric mapping
%   mesh_kmax                - Maximal wave numbers of a NiHu mesh
%   plot_node_numbers        - Plot node numbers of a NiHu mesh
%   plot_mesh                - Plot colored NiHu mesh
%   node_normals             - calculate outward normal vectors in nodes
%   plot_elem_normals        - Plot element normals of a mesh
%   mesh_volume              - Returns vector v, the volume for each element
%   node_volume              - Calculates corresponding volume for nodes
%   rotation_matrix          - Matrix of rotation around a vector

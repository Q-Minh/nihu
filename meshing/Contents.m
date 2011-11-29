% Toolbox NiHu Meshing
% Version 1.0 (R2006a) 08-Nov-2009
%
% Creating simple meshes
%   create_line              - Create a line mesh
%   create_slab              - Create slab mesh
%   create_circle_quadrant   - Create quadrant of a circle surface mesh
%   create_circle            - Create a circle surface mesh
%   create_brick_boundary    - Create a brick surface mesh
%   create_sphere_boundary   - Create a sphere surface mesh
%   create_catseye           - Create a Cat's eye surface mesh
%   create_brick             - Create a brick volume mesh
%   create_sphere            - Create a sphere volume mesh
%
% Basic mesh transforms
%   translate_mesh           - Translate mesh along a given vector
%   scale_mesh               - Scale mesh
%   reflect_mesh             - Reflect mesh to a symmetry plane
%   rotate_mesh              - Rotate mesh around a given vector
%   repeat_mesh              - Repeat (Copy) mesh along a given vector
%   extrude_mesh             - Extrude 1D and 2D mesh along a given direction
%   revolve_mesh             - Revolve 1D and 2D mesh around a given vector
%
% Mesh manipulation
%   mesh_select              - Element and node selection
%   mesh_section             - Return a rectangular section of a mesh
%   drop_IDs                 - Get rid of material, property, element and node IDs
%   drop_unused_nodes        - Exclude unused nodes from a mesh
%   merge_coincident_nodes   - Merge coincident nodes in a mesh
%   join_meshes              - Create a single mesh from several submeshes
%   split_independent_meshes - Split mesh to independent submeshes
%   flip_elements            - Flip elements of a fe model
%   get_boundary             - Extract domain boundary
%   get_faces                - Extract all faces of a FE mesh
%   get_free_faces           - Extract free (boundary) faces from FE faces
%   quad2tria                - Replace Quad elements with TRIA elements
%   surface2wireframe        - Extract surface TRIA and QUAD elements
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
%   import_bulk_mesh         - Import NiHu mesh from bulk file (.bdf)
%   export_bulk_mesh         - Export NiHu mesh into bulk file
%
% Other
%   boundingbox              - Compute bounding box of a NiHu mesh
%   isinside                 - Check whether points are inside a mesh

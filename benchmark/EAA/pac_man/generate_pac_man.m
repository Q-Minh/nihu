% generate_pac_man
% Generate input for Pac-Man computations

clear;

%% Create a Pac-Man mesh and field point meshes
Le = 3e-3;
[pac, field, dir_field, dir_ref_field] = create_pac_man(Le);

%% export surface mesh
surf_off_name = sprintf('mesh/pac_man_surf_%03dmm.off', Le*1000);
export_off_mesh(pac, surf_off_name);

%% export field point mesh
if (0)
    export_off_mesh(field, sprintf('mesh/pac_man_field_%03dmm_quad.off', Le*1000));
    
    line_field = field;
    line_field.Elements(:,2) = ShapeSet.LinearLine.Id;
    line_field.Elements(:,5:6) = field.Elements(:,[5 7]);
    line_field.Elements(:,7:end) = [];
    
    field_off_name = sprintf('mesh/pac_man_field_%03dmm.off', Le*1000);
    export_off_mesh(line_field, field_off_name);
end

%% export directivity field point mesh
dir_field_off_name = sprintf('mesh/pac_man_dir_field_%03dmm.off', Le*1000);
export_off_mesh(dir_field, dir_field_off_name);

%% export reference directivity field point mesh
dir_ref_field_off_name = sprintf('mesh/pac_man_ref_field_%03dmm.off', Le*1000);
export_off_mesh(dir_ref_field, dir_ref_field_off_name);

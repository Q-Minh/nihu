head = import_bulk_mesh('mesh\outer_low.nas');
k = min(mesh_kmax(head));
field = create_sphere_boundary(1, 7);
field = translate_mesh(field, mean(boundingbox(head)));

figure;
plot_mesh(field);   % plot the sphere
alpha .5;           % make it transparent
plot_mesh(head);    % plot the head

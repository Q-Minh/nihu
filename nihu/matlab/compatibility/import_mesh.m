function mesh = import_mesh(filename)
%IMPORT_MESH Import a mesh into NiHu
%   MESH = IMPORT_MESH(FILENAME) imports the mesh in file FILENAME into
%   NiHu.
%
% See also: import_off_mesh import_gmsh_mesh import_bulk_mesh

% known mesh format detectors and importers
formats = {
	@detect_off_mesh, @import_off_mesh
	@detect_bulk_mesh, @import_bulk_mesh
	@detect_gmsh_mesh, @import_gmsh_mesh
};

for i = 1 : size(formats,1)
	detector = formats{i,1};
	importer = formats{i,2};
	if detector(filename)
		mesh = importer(filename);
		return;
	end
end

end

function degen_elems_ind = find_degenerate_elements (model)
faces = get_faces(model.Elements);
degen_elem = false(size(model.Elements,1),1);
for i = 1:size(faces,1)
    if size(unique(faces(i,3:6)),2) == 2
        degen_elem(faces(i,1)) = true;
    end
end
degen_elems_ind = find(degen_elem);
end


function export_unv_data(data, fname)
fid = fopen(fname, 'w');
for i = 1 : length(data)
    data(i).write_to_file(fid);
end
fclose(fid);
end % of function

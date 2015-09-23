function res = import_unv_data(fname)
%import_unv_data import field data from universal file
%   res = import_unv_data(fname) imports all field data from the file
%   fname. res is an array of structures, each containing all parameters
%   of a field data.

% open file and read the whole into memory
fid = fopen(fname, 'rt');
data = textscan(fid, '%s', 'delimiter', '\n');
fclose(fid);
data = data{1};

% find universal blocks
blocksep = find(strcmp('-1', data));
blockind = [blocksep(1:2:end-1)+1 blocksep(2:2:end)-1];
blockid = str2double(data(blockind(:,1)));

% filter 2414 blocks
datablock = find(blockid == 2414);

% traverese blocks and process them one-by-one
for i = 1 : length(datablock)
    input = data(blockind(datablock(i),1)+1 : blockind(datablock(i),2));
    res(i) = unv_data();
    res(i).read_from_cellstr(input);
end

end % of function IMPORT_UNV_DATA

function block = read_unv_2414_block(input)
end % of function READ_UNV_2414_BLOCK

function value = get_value_from_cell(id, cell)
value = cell{id == cell2mat(cell(:,1)),2};
end % of function GET_VALUE_FROM_CELL

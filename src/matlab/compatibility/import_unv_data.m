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
blockind = [blocksep(1:2:end)+1 blocksep(2:2:end)-1];
blockid = str2double(data(blockind(:,1)));

% filter 2414 blocks
datablock = find(blockid == 2414);

% traverese blocks and process them one-by-one
for i = 1 : length(datablock)
    input = data(blockind(datablock(i),1)+1 : blockind(datablock(i),2));
    res(i) = read_unv_2414_block(input);
end

end % of function IMPORT_UNV_DATA

function block = read_unv_2414_block(input)
%READ_UNV_2414_BLOCK

block.dataset_label = sscanf(input{1}, '%u');
block.dataset_name = input{2};

dataset_locations = {
    1, 'Data at nodes'
    2, 'Data on elements'
    3, 'Data at nodes on elements'
    5, 'Data at points'
    };

rec3 = sscanf(input{3}, '%u');
block.dataset_location = get_value_from_cell(rec3(1), dataset_locations);

model_types = {
    0, 'Unknown'
    1, 'Structural'
    2, 'Heat transfer'
    3, 'Fluid flow'
    };

analysis_types = {
    0,   'Unknown'
    1,   'Static'
    2,   'Normal mode'
    3,   'Complex eigenvalue first order'
    4,   'Transient'
    5,   'Frequency response'
    6,   'Buckling'
    7,   'Complex eigenvalue second order'
    9,   'Static non-linear'
    10,   'Craig-Bampton constraint modes'
    11,   'Equivalent attachment modes'
    12,   'Effective mass modes'
    13,   'Effective mass matrix'
    14,   'Effective mass matrix'
    };

data_characteristics = {
    0,   'Unknown'
    1,   'Scalar'
    2,   '3 DOF global translation vector'
    3,   '6 DOF global translation & rotation vector'
    4,   'Symmetric global tensor'
    6,   'Stress resultants'
    };

result_types = {
    8, 'Displacement'
    11, 'Velocity'
    117, 'Pressure'
    303, 'Sound Intensity'
    };

data_types = {
    1, 'Integer'
    2, 'Single precision floating point'
    4, 'Double precision floating point'
    5, 'Single precision complex'
    6, 'Double precision complex'
    };

block.id = input(4:8);
rec9 = sscanf(input{9}, '%u', 6);
block.model_type = get_value_from_cell(rec9(1), model_types);
block.analysis_type = get_value_from_cell(rec9(2), analysis_types);
block.data_characteristic = get_value_from_cell(rec9(3), data_characteristics);
block.result_type = get_value_from_cell(rec9(4), result_types);
block.data_type = get_value_from_cell(rec9(5), data_types);
iscplx = ~isempty(strfind(block.data_type, 'complex'));

block.nValDC = rec9(6);

% rec10 = [sscanf(input{10}, '%u', 8); sscanf(input{11}, '%u', 2)];
rec12 = [sscanf(input{12}, '%g', 6); sscanf(input{13}, '%g', 6)];
block.frequency = rec12(2);

input = input(14:end);

stream = strcat(input, {' '});
stream = horzcat(stream{:});
n = block.nValDC;
if iscplx
    n = n * 2;
end
format = ['%u' repmat('%g', 1, n)];
numbers = sscanf(stream, format, [n+1 Inf]);
block.node = numbers(1,:);
block.data = numbers(2:end,:);
if iscplx
    block.data = complex(block.data(1:2:end,:), block.data(2:2:end,:));
end

end % of function READ_UNV_2414_BLOCK

function value = get_value_from_cell(id, cell)
value = cell{id == cell2mat(cell(:,1)),2};
end % of function GET_VALUE_FROM_CELL

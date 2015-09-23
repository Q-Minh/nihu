function export_unv_data(dof, data, freq, datatype, fname)
fid = fopen(fname, 'w');
fprintf(fid, '%6d\n', -1);              % block start
fprintf(fid, '%6d\n', 2414);            % data block
fprintf(fid, '%10d\n', 1);              % dataset label
fprintf(fid, 'Exported by NiHu\n');     % dataset name
fprintf(fid, '%10d\n', 1);              % data at nodes
% textual id lines 1-5
fprintf(fid, '%s\nNONE\nNONE\nNONE\nNONE\n', datestr(now()));

model_type = 1; % structural
analysis_type = 2; % normal mode
data_characteristics = 2; % 3dof global translation vector
result_type = 11; % velocity
data_type = 5; % single precision complex
nvaldc = 3; % number of values per data component
fprintf(fid, '%10d%10d%10d%10d%10d%10d\n', model_type, analysis_type,...
    data_characteristics, result_type, data_type, nvaldc);
% integer analysis type data all zeros
idata = zeros(10,1);
fprintf(fid, '%10d%10d%10d%10d%10d%10d%10d%10d\n%10d%10d\n', idata);
% real analysis type data
rdata = zeros(12,1);
rdata(2) = freq;
fprintf(fid, [repmat('%13.5e', 1, 6) '\n' repmat('%13.5e', 1, 6) '\n'], rdata);
for i = 1 : length(dof)
    fprintf(fid, '%10d\n', dof(i));
    fprintf(fid, [repmat('%13.5e', 1, 6), '\n'], [real(data(i,:)); imag(data(i,:))]);
end
fprintf(fid, '%6d\n', -1);
fclose(fid);
end % of function

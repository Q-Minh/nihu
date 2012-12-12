function AcouFEMLog(string)

fid = fopen('AcouFEM.Log.m', 'a');
fseek(fid, 0, 'eof');
fprintf(fid, '%s\n', string);
fclose(fid);
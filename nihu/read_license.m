function Include = read_license

fid = fopen('license.txt');
Include = textscan(fid, '%s', 'Delimiter', '\n', 'WhiteSpace', '');
Include = Include{1};
fclose(fid);

Include = strcat({'// '}, Include);

end

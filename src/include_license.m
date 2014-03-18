function include_license(Insert, root)
if nargin == 1
    root = pwd();
end

files = dir(fullfile(root, '*.*pp'));

for i = 1 : length(files)
    fprintf(1, 'processing file %s\n', fullfile(root, files(i).name));
    fid = fopen(fullfile(root, files(i).name), 'rt');
    content = textscan(fid, '%s', 'delimiter', '\n', 'WhiteSpace', '');
    content = content{1};
    fclose(fid);
    
    if (strcmp(Insert{1}, content{1}))
%         content = regexprep(content, ...
%             'Copyright \(C\) 2012-2013',...
%             'Copyright \(C\) 2012-2014');
        continue;
    else
        content = [
            Insert
            {''}
            content
            ];
    end
    
    fid = fopen(fullfile(root, files(i).name), 'wt');
    for l = 1 : length(content);
        fprintf(fid, '%s\n', content{l});
    end
    fclose(fid);
end

exclude = {'.', '..', 'ThirdParty'};

dirs = dir(fullfile(root));
dirs = dirs(cell2mat({dirs.isdir}));
for i = 1 : length(dirs)
    if ismember(dirs(i).name, exclude)
        continue;
    end
    
    include_license(Insert, fullfile(root, dirs(i).name));
end

end

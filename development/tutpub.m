function tutpub(directory)
%TUTPUB  Publish all m files in current directory

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

%% Parameter check
if nargin == 0
    directory = '.';
end

%% Recursively enter each subdirectory
d = dir(directory);
nD = length(d);
for iD = 1 : nD
    if d(iD).isdir && ~strcmp(d(iD).name, '.') &&...
            ~strcmp(d(iD).name, '..')
        tutpub(fullfile(directory, d(iD).name));
    end
end

%% Publish all m files
d = dir(fullfile(directory, '*.m'));
nD = length(d);
for iD = 1 : nD
    curdir = pwd;
    cd(directory);
    publish(d(iD).name);
    close all;
    cd(curdir);
end

end

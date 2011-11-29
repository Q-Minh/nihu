function tuttoc(varargin)
%tuttoc  Tutorial table of content file generation
%   tuttoc(directory) Creates a html table of contents for a collection of
%   tutorials published by Matlab. The function searches in the
%   subdirectories of the given root folder for files of type
%   /html/*.htm*
%   and generates the table of contents using the html tag <title>. The
%   result file is called 'index.html' and is placed in the current
%   directory.

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

%% Parameter check
switch nargin
    case 1
        directory = '.';
        title = varargin{1};
        level = 0;
        fid = [];
    case 2
        title = varargin{1};
        directory = varargin{2};
        level = 0;
        fid = [];
    case 3
        directory = varargin{1};
        level = varargin{2};
        fid = varargin{3};
    otherwise
        error('nihu');
end

%% Initialization at level 0
if level == 0
    % open html file for writing
    fid = fopen('index.html', 'w');
    % html header
    fprintf(fid, '<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n');
    fprintf(fid, '<head>\n');
    fprintf(fid, '<meta http-equiv="Content-Type" content="text/html; charset=utf-8">\n');
    fprintf(fid, '<title>Tutorials table of contents</title>\n');
    fprintf(fid, '<style type="text/css">\n');
    fprintf(fid, 'body {background-color: white; margin:10px;}\n');
    fprintf(fid, 'h1 {color: #990000; font_size:x-large;}\n');
    fprintf(fid, 'h2 {color: #990000; font_size:medium;}\n');
    fprintf(fid, 'p,h1,h2,div.content,div {max-width: 600px;}\n');
    fprintf(fid, 'p.footer {text-align: right; font-size: xx-small; font-weight: lighter; font-style: italic; color: gray; }\n');
    fprintf(fid, '</style>\n');
    fprintf(fid, '</head>\n');
    % html body initialization
    fprintf(fid, '<body>\n');
    fprintf(fid, '<div class="content">\n');
    fprintf(fid, '<h1>%s tutorials</h1>\n', title);
    fprintf(fid, '<h2>Contents</h2>\n');
end
fprintf(fid, repmat('\t', 1, level));
fprintf(fid, '<ul>\n');

%% Recursively enter each subdirectory different from 'html'
d = dir(directory);
nD = length(d);
for iD = 1 : nD
    if d(iD).isdir && ~strcmp(d(iD).name, '.') &&...
            ~strcmp(d(iD).name, '..') && ~strcmp(d(iD).name, 'html')
        fi = fopen([directory filesep d(iD).name filesep 'description.txt']);
        if fi ~= -1
            dirinfo = fgetl(fi);
            fclose(fi);
            fprintf(fid, repmat('\t', 1, level+1));
            fprintf(fid, '<li>%s (%s)</li>\n', dirinfo, d(iD).name);
            tuttoc([directory filesep d(iD).name], level+1, fid);
        end
    end
end

%% Build list of 'html/*.htm*' files
d = dir([directory filesep 'html' filesep '*.htm*']);
nD = length(d);
for iD = 1 : nD
    fullname = fullfile(directory, 'html', d(iD).name);
    % open html file and get <title>...</title> information
    f = fopen(fullname, 'r');
    fstring = fscanf(f, '%c', Inf);
    fclose(f);
    [start,stop] = regexp(lower(fstring), '<title>.*?</title>');
    titlestr = fstring(start+7:stop-8);
    %
    fprintf(fid, repmat('\t', 1, level+1));
    fprintf(fid, '<li><a href="%s">%s (%s)</a></li>\n',...
        fullname, titlestr, d(iD).name);
end

%% Close the list and the html file
fprintf(fid, repmat('\t', 1, level));
fprintf(fid, '</ul>\n');
if level == 0
    fprintf(fid, '<p class="footer">\n');
    fprintf(fid, 'Generated %s by <tt>NIHU tuttoc</tt>\n', date);
    fprintf(fid, '</p>\n');
    fprintf(fid, '</div>\n');
    fprintf(fid, '</body>\n');
    fclose(fid);
end

end

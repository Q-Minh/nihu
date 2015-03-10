function [bb, path] = read_epspath(filename)
%READ_EPSPATH Read curves from .eps file
%   [BB, PATH] = READ_EPSPATH(FILENAME) reads straight lines and
%   Bezier-curves from a post script file FILENAME. The curves are
%   returned in a matrix PATH, the BoundingBox data is returned in the
%   4-element vector BB.
%
% See also: MESHPATH

%   Copyright 2008-2012 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified: 2012.12.19.

% Default parameters
scale = [1 1];
translate = [0 0];

% Read eps file
fid = fopen(filename, 'rt');
feps = textscan(fid, '%s', 'delimiter', '\n');
fclose(fid);
feps = feps{1};

stop = find(strcmp('%%EndComments', feps));
Comments = feps(1:stop-1);
Comments = strcat(Comments, {' '});

start = find(strcmp('%%BeginProlog', feps));
if ~isempty(start)
    stop = find(strcmp('%%EndProlog', feps));
    Prolog = feps(start+1:stop-1);
else
    Prolog = [];
end

start = find(strcmp('%%Page: 1 1', feps));
Page = feps(start+1:end);
Page = strcat(Page, {' '});

float = '[-+]?[0-9]*\.?[0-9]+';
floats = ['(' float '\s+)+'];

% Read Bounding Box from Comments
I = skipempty(regexp(Comments, ['%%BoundingBox:\s*(?<bb>' floats ')'], 'names'));
bb = reshape(sscanf(I.bb, '%g'),2,2)';

I = skipempty(regexp(Comments, ['%%HiResBoundingBox:\s*(?<bb>' floats ')'], 'names'));
if ~isempty(I)
    bb = reshape(sscanf(I.bb, '%g'),2,2)';
end

%
if ~isempty(Prolog)
    I = skipempty(regexp(Prolog, '/(?<token>\S+)\s*{\s*(?<command>.*)\s*}\s*bind\s*def', 'names'));
    for l = 1 : length(I)
        % we replace ' l ' to ' lineto '
        Page = strrep(Page, [' ' I(l).token ' '], [' ' I(l).command ' ']);
        % we replace 'l ' to 'lineto ' at the beginnig of lines
        st = find(strncmp(Page, [I(l).token ' '], length(I(l).token)+1));
        for i = 1 : length(st)
            Page{st(i)} = [I(l).command ' ' Page{st(i)}(length(I(l).token)+1:end)];
        end
    end
end

%
Commands = regexp(cell2mat(Page(:)'), ['(?<lw>' floats ')(?<command>scale|translate|moveto|lineto|curveto|6 array astore concat)'], 'names');
nCom = length(Commands);

iEntry = 0;     % number of curves read from the eps file
pos = [0 0];    % starting cursor position
path = zeros(nCom, 9);   % allocate space for 100 curves
for iCom = 1 : nCom
    C = Commands(iCom);
    coord = sscanf(C.lw, '%g');
    switch C.command
        case 'scale'
            scale = coord(end-1:end);
        case 'translate'
            translate = coord(end-1:end);
        case 'moveto'
            pos = coord(end-1:end)';
        case 'lineto'
            iEntry = iEntry+1;
            to = coord(end-1:end)';
            path(iEntry,1:5) = [1 pos to];
            pos = to(end-1:end);
        case 'curveto'
            iEntry = iEntry+1;
            to = coord(end-5:end)';
            path(iEntry,1:9) = [2 pos to];
            pos = to(end-1:end);
        case '6 array astore concat'
            scale = coord([1,4]);
            translate = coord(5:6);
    end
end

% delete zeros lines from path
path = path(1:iEntry,:);

% Scale and translate coordinates
path(:,2:2:end) = path(:,2:2:end)*scale(1);
path(:,3:2:end) = path(:,3:2:end)*scale(2);
path(:,2:2:end) = path(:,2:2:end)+translate(1);
path(:,3:2:end) = path(:,3:2:end)+translate(2);

end

function Items = skipempty(Items)
k = 0;
for i = 1 : length(Items)
    if ~isempty(Items{i})
        k = k+1;
        Items{k} = Items{i};
    end
end
Items = cell2mat(Items(1:k));

end
